import numpy as np,os,glob,sys
#import losoto; from losoto.h5parm import h5parm
from pyrap.tables import table
#PT = '/opt/cep/parseltongue/parseltongue-2.3/bin/ParselTongue'
PT = '/home/njj/root_overflow/parseltongue-2.1/bin/ParselTongue'
DIFMAP = '/pkg/uvf_difmap/difmap'
#DIFMAP = '/home/jackson/uvf_difmap/difmap'
try:
    from AIPSData import AIPSUVData
    from AIPSTask import AIPSTask,AIPSList,AIPSMessageLog
    import Wizardry
    from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
    HAVE_AIPS = True
except:
    HAVE_AIPS = False

def find_difmap_chan (infile,indisk=1,inseq=1):
    try:
        AIPSUVData('DIFMAP','UVDATA',indisk,inseq).zap()
        print ('Removed DIFMAP.UVDATA.%d disk %d'%(inseq,indisk))
    except:
        pass
    fitld=AIPSTask('fitld')
    fitld.datain=infile if '/' in infile else './'+infile
    fitld.outname = 'DIFMAP'
    fitld.outclass = 'UVDATA'
    fitld.outdisk = indisk
    fitld.outseq = inseq
    fitld.go()
    data = WizAIPSUVData('DIFMAP','UVDATA',indisk,inseq)
    for i in data:
        v = i.visibility
        na = np.asarray(np.ravel(v[:,:,0,0]),dtype='bool')
        try:
            a = na | a
        except:
            a = np.copy(na)
    return a

# given an array of OK channels, write difmap select line
def chan2write (output, a):
    nums=sorted(set(a))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    seq = np.ravel(list(zip(edges,edges)))+1 # difmap first chan=1 not 0
    fo=open('dchan_%s'%(output),'w')
    fo.write('select I,')
    for i in range(len(seq)):
        fo.write('%d'%seq[i])
        if i<len(seq)-1:
            fo.write(',')
    fo.write('\n')
    fo.close()

def check_ms (filename,datacolumn="CORRECTED_DATA",flagname="FLAG"):
    # open ms
    ms=table(filename,ack=False)
    # read data, flags
    d=ms.getcol(datacolumn)
    f=ms.getcol(flagname)
    d[f==True]=0.
    # see presence of nonzero vis per chan
    goodchans=[]
    for i in range(d.shape[1]):
        datasum=np.sum(np.abs(d[:,i],dtype=np.float64))
        if datasum!=0.:
            goodchans.append(i)
    return goodchans

def corplt2array():
    f = open('./CORPLT')
    fc = f.readlines()
    f.close()
    ut=[]
    stn = np.array([],dtype='S')
    for l in fc:
        ls = l.split()
        if not len(ls):
            continue
        if ls[0]=='Stn':
            stn = np.append(stn,ls[3])
        if len(ls)==4 and ls[0] in ['Amp','Phs']:
            ut.append(float(ls[1]))  # faster than numpy append
    
    ut = np.unique(np.asarray(ut,dtype='f'))
    amp = np.ones((len(stn),len(ut)))*np.nan
    phs = np.ones((len(stn),len(ut)))*np.nan
    amperr = np.ones((len(stn),len(ut)))*np.nan
    phserr = np.ones((len(stn),len(ut)))*np.nan
    stn_idx = 0
    for l in fc:
        if l[:3] in ['Amp','Phs']:
            t_ut,t_val,t_err = np.asarray(l.split()[1:],dtype='float')
            ut_idx = np.argwhere(ut==t_ut)[0][0]
            if l[:3] == 'Amp':
                amp[stn_idx,ut_idx] = t_val
                amperr[stn_idx,ut_idx] = t_err
            else:
                phs[stn_idx,ut_idx] = t_val
                phserr[stn_idx,ut_idx] = t_err
        if l[:3] == 'Stn':
            stn_idx += 1
    return amp,amperr,phs,phserr,ut,stn

def clean_selfcal_loop (fs,weight):
    fs.write('uvw %s\n'%weight)
    fs.write('if(peak(flux,max)/imstat(rms) > clean_sigma)\n')
    fs.write('	repeat;\\\n')
    fs.write('		peakwin 1.5; clean; selfcal\n')
    fs.write('	until(peak(flux,max)/imstat(rms) < clean_sigma)\n')
    fs.write('end if\n')

def get_isfits (filename):
    os.system('od -c -N 6 %s>isfits_temp'%filename)
    if open('isfits_temp').readline()=='0000000   S   I   M   P   L   E\n':
        return True
    return False
    
def dif_script (infile,pol='XX',aipsno=340,clean_sigma=6,map_size=512,\
                pixel_size=100,obs_length=900,datacolumn='CORRECTED_DATA'):
    isfits = get_isfits(infile)
    if isfits:
        fitsfile = infile
    else:
        if infile[-3:] in ['.MS','.ms']:
            fitsfile = infile[:-3]+'.fits' 
        else:
            fitsfile = infile+'.fits'
        os.system('ms2uvfits in=%s out=%s writesyscal=F'%(infile,fitsfile))
    if not os.path.isfile('dchan_'+fitsfile):
        if isfits:
            if HAVE_AIPS:
                a = find_difmap_chan (fitsfile)
            else:
                pass
                # convert to MS here and use check_ms
        else:
            a = check_ms(infile,datacolumn=datacolumn)
        chan2write('dchan_'+fitsfile,a)
    # Open script and declare variables
    fs = open('dif_script','w')
    fs.write('float clean_sigma; clean_sigma = %f\n'%clean_sigma)
    fs.write('float map_size; map_size = %f\n'%map_size)
    fs.write('float pixel_size; pixel_size = %f\n'%pixel_size)
    fs.write('float obs_length; obs_length = %f\n'%obs_length)
    fs.write('float sn_p;\nfloat sn_a\n')
    # Do the imaging with channels selected by find_difmap_chan
    fs.write('observe %s\n'%fitsfile)
    if not os.path.isfile ('dchan_'+fitsfile):
        os.system('%s find_difmap_chan.py %d %s'%(PT,aipsno,fitsfile))
    fsel = open ('dchan_'+fitsfile)
    fs.write('%s'%fsel.readline().replace('I',pol))
    fsel.close()
    fs.write('mapsize map_size,pixel_size\n')
    fs.write('startmod "",1\n')
    fs.write('peakwin 1.5\nselfcal false,false,0\n')
    clean_selfcal_loop (fs,'2,0')
    clean_selfcal_loop (fs,'2,-1')
    clean_selfcal_loop (fs,'0,-1')
    fs.write('gscale\n')
    fs.write('sn_p = peak(flux,max)/imstat(rms)\n')
    fs.write('repeat;\\\n')
    fs.write('    selfcal true,true,obs_length\n')
    fs.write('    sn_a = peak(flux,max)/imstat(rms)\n')
    fs.write('    if(peak(flux,max)/imstat(rms) > clean_sigma)\n')
    fs.write('        if(sn_a <= sn_p)\n')
    fs.write('            peakwin 1.5; clean; selfcal\n')
    fs.write('        obs_length=obs_length/2\n')
    fs.write('            sn_a = sn_p\n')
    fs.write('        end if\n')
    fs.write('        if(sn_a > sn_p)\n')
    fs.write('            clrmod false,true\n')
    fs.write('        obs_length=obs_length/2\n')
    fs.write('            sn_a = sn_p\n')
    fs.write('        end if\n')
    fs.write('    else\n')
    fs.write('        obs_length=obs_length/2\n')
    fs.write('    end if\n')
    fs.write('until(obs_length < 2)\n')
    fs.write('selfcal true,true,0\n')
    fs.write('delwin\nclean 1000,0.01\n')
    fs.write('device %s.ps/vps\ncorplot\n'%fitsfile.replace('.fits','_auto%s.p'%pol))
    fs.write('device /NULL\nmapl\ncmul=3*imstat(rms)\n')
    fs.write('device %s.ps/vps\nmapl clean,false\n' % \
             fitsfile.replace('.fits','_auto%s'%pol))
    fs.write('save %s\nquit\n' % fitsfile.replace('.fits','_auto%s'%pol))
    fs.close()

    os.system('%s <dif_script'%DIFMAP)

# Specified file may be fits or MS
# If MS, will be inspected for good channels and converted to FITS for mapping
#    (any .MS or .ms extension converted to .fits, otherwise .fits added)
#    Corrected_data column assumed unless otherwise specified
# If FITS and we have AIPS/Parseltongue, then inspected for good channels
#    in AIPS and mapped
# If FITS and we do not have AIPS/Parseltongue, converted to ms for inspection
#    of good channels and the original FITS file mapped
    
def loop3_difmap (infile,datacolumn='CORRECTED_DATA'):
    os.system('rm CORPLT')
    infile = dif_script(infile,'XX')
    ampXX,amperrXX,phsXX,phserrXX,utXX,stnXX = corplt2array()
    os.system('mv CORPLT CORPLT_XX')
    infile = dif_script(infile,'YY')
    ampYY,amperrYY,phsYY,phserrYY,utYY,stnYY = corplt2array()
    os.system('mv CORPLT CORPLT_YY')
    #   nb assumes XX and YY solutions have same number of telescopes and UTs -
    #   combine properly if this is not the case
    amp = np.rollaxis(np.dstack((ampXX,ampYY)),2,0)
    amperr = np.rollaxis(np.dstack((amperrXX,amperrYY)),2,0)
    phs = np.rollaxis(np.dstack((phsXX,phsYY)),2,0)
    phserr = np.rollaxis(np.dstack((phserrXX,phserrYY)),2,0)
    ut,stn = utXX,stnXX
#    h5parmfile = infile.replace('.fits','_auto.h5')
#    data = h5parm(h5parmfile, readonly = False)
#    outSolset = data.makeSolset('sol000')
#    outSolset.makeSoltab(soltype='amplitude',soltabName='amplitude000',axesNames=['pol','ant','time'],axesVals=[['XX','YY'],stn,ut],vals=amp,weights=np.ones_like(amp))
#    outSolset.makeSoltab(soltype='phase',soltabName='phase000',axesNames=['pol','ant','time'],axesVals=[['XX','YY'],stn,ut],vals=phs,weights=np.ones_like(phs))
#    data.close()
    
