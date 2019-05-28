from AIPS import AIPS,AIPSDisk
from AIPSData import AIPSUVData
from AIPSTask import AIPSTask,AIPSList,AIPSMessageLog
import Wizardry
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import numpy as np,sys

def find_difmap_chan (inna, incl, indisk=1, inseq=1):
    INCR=1
    if incl=='':
        output = inna.split('/')[-1] if '/' in inna else inna
        try:
            AIPSUVData('DIFMAP','UVDATA',indisk,inseq).zap()
            print 'Removed DIFMAP.UVDATA.%d disk %d'%(inseq,indisk)
        except:
            pass
        fitld=AIPSTask('fitld')
        fitld.datain=inna if '/' in inna else './'+inna
        fitld.outname = 'DIFMAP'
        fitld.outclass = 'UVDATA'
        fitld.outdisk = indisk
        fitld.outseq = inseq
        fitld.go()
        inna,incl = 'DIFMAP','UVDATA'
    else:
        output = '%s_%s'%(inna,incl)
    data = WizAIPSUVData(inna,incl,indisk,inseq)
    for i in data:
        v = i.visibility
        na = np.asarray(np.ravel(v[:,:,0,0]),dtype='bool')
        try:
            a = na | a
        except:
            a = np.copy(na)
    bstart, bend = [],[]
    for i in range(len(a)):
        if (i==0 and a[i]) or (i!=0 and a[i] and not a[i-1]):
            bstart.append(i+1)
        if (i!=0 and not a[i] and a[i-1]):
            bend.append(i)
        if (a[i] and i==len(a)-1):
            bend.append(len(a))
    fo=open('dchan_%s'%(output),'w')
    fo.write('select I,')
    for i in range(len(bstart)):
        fo.write('%d,%d'%(bstart[i],bend[i]))
        if i!=len(bstart)-1:
            fo.write(',')
    fo.write('\n')
    fo.close()
        

if len(sys.argv)==3:
    AIPS.userno = int(sys.argv[1])
    find_difmap_chan (sys.argv[2],'')
elif len(sys.argv)==4:
    AIPS.userno = int(sys.argv[1])
    find_difmap_chan (sys.argv[2],sys.argv[3])
else:
    print 'Usage: parseltongue find_difmap_chan.py [userno] [filename]'
    exit(0)

