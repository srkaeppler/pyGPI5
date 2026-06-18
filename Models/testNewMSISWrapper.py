import numpy
import pymsis 
import datetime

"""
06/17/2026
new wrapper using pymsis: https://github.com/SWxTREC/pymsis
This is pip installable and the packages are consistent with iricore!

Make the wrapper effectively the same as what I did in MSIS.py 
and deprecate that code.


"""
class MSIS:


    def __init__(self):
        
        return


    def MSIS2(self, tUnix,glat,glong,Heibeg,Heiend,step, CGSorSI = 'SI'):
            # standard interface - hate it or love it...
            #tUnix,Lat,Lon,Heibeg,Heiend,step
            # datetime.datetime.fromtimestamp(timestamp, datetime.UTC)
            # dt1 = datetime.datetime.utcfromtimestamp(tUnix) # deprecated
            dt1 = datetime.datetime.fromtimestamp(tUnix,datetime.UTC) 
            altkm = numpy.arange(Heibeg,Heiend,step)
            
            # pymsis.calculate(dates, lons, lats, alts, f107s=None, f107as=None, aps=None, *, options=None, version=2.1, **kwargs)
            # this is only in SI coming out, so I will need to do the conversion if want it in CGS
            # might have to deprecate that feature
            # version will be MSIS00 to be consistent with past results, should add some keyword stuff
            
            f107,f107a, ap = pymsis.utils.get_f107_ap(dt1)
            # print('f107',f107)
            # print('f107a',f107a)
            # print('ap',ap)
            d = pymsis.calculate(dt1,glong,glat, altkm,version=0,f107s=f107, f107as=f107a,
                                 aps=numpy.array([ap]))
            
            # output part
            # [0- Total mass density (kg/m3),
            # 1-N2 # density (m-3),
            # 2-O2 # density (m-3),
            # 3-O # density (m-3),
            # 4-He # density (m-3),
            # 5-H # density (m-3),
            # 6-Ar # density (m-3),
            # 7-N # density (m-3),
            # 8-Anomalous oxygen # density (m-3),
            # 9-NO # density (m-3),
            # 10-Temperature (K)]

            MSISout=dict()
            MSISout['Altitude']=altkm
            MSISout['Tn']= numpy.ravel(d[0,0,0,:,-1])
            MSISout['MassDensity']= numpy.ravel(d[0,0,0:,0])
            MSISout['nO'] = numpy.ravel(d[0,0,0,:,3])
            MSISout['nN2'] = numpy.ravel(d[0,0,0,:,1])
            MSISout['nO2'] = numpy.ravel(d[0,0,0,:,2])
            MSISout['nN'] = numpy.ravel(d[0,0,0,:,7])
            # print('dshape', d.shape, MSISout['nN'].shape)
            # print('nN at last altitude', d[0,0,0,-1,7])
        
            # see if this breaks anything...
            #MSISout['F107']=f107
            #MSISout['F107a']=f107a
            #MSISout['AP']=ap
            # MSISout['AverageMass'] = MSISout['MassDensity']/
            # MSISout['Nm'] = numpy.zeros(altkm.shape)
            # keep these definitions in
            #MSISout['MassDensity'][Iht] = d[5]
            #     MSISout['AverageMass'][Iht] = d[5]/(d[0]+d[1]+d[2]+d[4]+d[6]+d[7]+d[8])
            #MSISout['Nm'][Iht] = (d[0]+d[1]+d[2]+d[4]+d[6]+d[7]+d[8])


            # for Iht in range(altkm.size):
            #     (d,t)=self.run_MSIS(int(year),int(doy),hrUT,altkm[Iht],glat,glong,ap,\
            #                         f107a,f107,-1, CGSorSI=CGSorSI)
            #     MSISout['MassDensity'][Iht] = d[5]
            #     MSISout['AverageMass'][Iht] = d[5]/(d[0]+d[1]+d[2]+d[4]+d[6]+d[7]+d[8])
            #     MSISout['Texo']=t[0]
            #     MSISout['Tn'][Iht]=t[1]
            #     MSISout['Nm'][Iht] = (d[0]+d[1]+d[2]+d[4]+d[6]+d[7]+d[8])
            #     MSISout['nO'][Iht] = d[1]
            #     MSISout['nN2'][Iht] = d[2]
            #     MSISout['nO2'][Iht] = d[3]
            #     MSISout['nN'][Iht] = d[7]



            return MSISout

if __name__ == '__main__':
    # test case
    # '/Users/srkaeppler/research/data/pygpi5_dev/pyGPI5/Models/nrlmsise00/libnrlmsise-00.so'
    msis = MSIS()
    year=2025
    doy=172
    hrUT=29000/3600.
    g_lat=60
    g_long=45#-70
    CGSorSI = 'SI'
    # outDict = msis.MSIS(doy,hrUT,g_lat,g_long,year, CGSorSI=CGSorSI)
    # print(outDict)

    print("testing MSIS 2 \n \n \n \n \n")
    dt1970 = datetime.datetime(1970,1,1,0,0,0)
    dt1 = datetime.datetime(2025,6,1,13,41,56)
    print(dt1)
    tunix1 = (dt1-dt1970).total_seconds()
    outDict2 = msis.MSIS2(tunix1,g_lat,g_long,50.,100.,1.,CGSorSI='SI')

    print('testnewMSISwrapper nO', outDict2['Altitude'][-30],outDict2['nO'][-30]/1e11)
    print('testnewMSISwrapper Tn', outDict2['Altitude'][-30],outDict2['Tn'][-30])

    print('testnewMSISwrapper nO', outDict2['Altitude'][-20],outDict2['nO'][-20]/1e11)
    print('testnewMSISwrapper Tn', outDict2['Altitude'][-20],outDict2['Tn'][-20])
    
    print('testnewMSISwrapper nO', outDict2['Altitude'][-1],outDict2['nO'][-1]/1e11)
    print('testnewMSISwrapper Tn', outDict2['Altitude'][-1],outDict2['Tn'][-1])

