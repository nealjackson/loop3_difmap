import numpy as np,os,glob,sys
import losoto; from losoto.h5parm import h5parm
PT = '/opt/cep/parseltongue/parseltongue-2.3/bin/ParselTongue'
DIFMAP = '/home/jackson/uvf_difmap/difmap'

def array2h5 (h5file,amp,amperr,phs,phserr,ut,stn):
    data = h5parm (h5file,readonly=False)
    OutSolset = data.getSolset(solset_out)
    OutSolset.makeSoltab(soltype='amplitude', soltabName='amplitude', \
                   axesNames=['ant', 'time'], \
                   axesVals=[station_names, freqs_ampl, ['XX','YY']], \
                   vals=amps_array[:,-1,:,:], weights=weights_amp[:,-1,:,:])


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

def dif_script (fitsfile,pol='XX',aipsno=340,clean_sigma=6,map_size=512,pixel_size=100,\
                obs_length=900):
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

os.system('rm CORPLT')
dif_script(sys.argv[1],'XX')
ampXX,amperrXX,phsXX,phserrXX,utXX,stnXX = corplt2array()
os.system('mv CORPLT CORPLT_XX')
dif_script(sys.argv[1],'YY')
ampYY,amperrYY,phsYY,phserrYY,utYY,stnYY = corplt2array()
os.system('mv CORPLT CORPLT_YY')
#   nb assumes XX and YY solutions have same number of telescopes and UTs -
#   combine properly if this is not the case
amp = np.rollaxis(np.dstack((ampXX,ampYY)),2,0)
amperr = np.rollaxis(np.dstack((amperrXX,amperrYY)),2,0)
phs = np.rollaxis(np.dstack((phsXX,phsYY)),2,0)
phserr = np.rollaxis(np.dstack((phserrXX,phserrYY)),2,0)
ut,stn = utXX,stnXX
h5parmfile = sys.argv[1].replace('.fits','_auto.h5')
data = h5parm(h5parmfile, readonly = False)
outSolset = data.makeSolset('sol000')
outSolset.makeSoltab(soltype='amplitude',soltabName='amplitude000',axesNames=['pol','ant','time'],axesVals=[['XX','YY'],stn,ut],vals=amp,weights=np.ones_like(amp))
outSolset.makeSoltab(soltype='phase',soltabName='phase000',axesNames=['pol','ant','time'],axesVals=[['XX','YY'],stn,ut],vals=phs,weights=np.ones_like(phs))
data.close()
