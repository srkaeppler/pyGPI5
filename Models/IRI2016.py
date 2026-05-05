import numpy
import sys
import datetime
import os
import iricore

"""
Updated on 04/22/2026
my old wrapper which is now labeled as IRI2016.bak will not work anymore because f2py got sunsetted.
I would have to write the wrapper in Meson which I don't want to do.
Found an alternative python package which seems to do most of the heavy lifting 
and will hopefully play well in python3.xx for a while to come
https://github.com/MIST-Experiment/iricore
https://iricore.readthedocs.io/en/latest/#

Objective is to  keep the wrapper interface effectively the same, but now just use this instead

"""

class IRI2016:

    def __init__(self):
        # iricore.update() # this command will just make sure the internal apf107.dat and ig_rz.dat are updated

        # 04/22/2026 - simply setting all of the options as 1 was a bad idea
        # more careful consideration is needed about what are good defaults.
        self.JF = numpy.ones(50) # all of the options
        #self.JF[23] = 1
        self.JMAG = 0 # geographic = 0, geomagnetic =1
        return

    def IRI2016(self,tUnix,glat,glon,AltitudeMin,AltitudeMax,deltaAltitude):
        """

        # OUTF(1,*)  ELECTRON DENSITY/M-3
        # C               OUTF(2,*)  NEUTRAL TEMPERATURE/K
        # C               OUTF(3,*)  ION TEMPERATURE/K
        # C               OUTF(4,*)  ELECTRON TEMPERATURE/K
        # C               OUTF(5,*)  O+ ION DENSITY/% or /M-3 if jf(22)=f
        # C               OUTF(6,*)  H+ ION DENSITY/% or /M-3 if jf(22)=f
        # C               OUTF(7,*)  HE+ ION DENSITY/% or /M-3 if jf(22)=f
        # C               OUTF(8,*)  O2+ ION DENSITY/% or /M-3 if jf(22)=f
        # C               OUTF(9,*)  NO+ ION DENSITY/% or /M-3 if jf(22)=f
        # C                 AND, IF JF(6)=.FALSE.:
        # C               OUTF(10,*)  CLUSTER IONS DEN/% or /M-3 if jf(22)=f
        # C               OUTF(11,*)  N+ ION DENSITY/% or /M-3 if jf(22)=f
        """
        # set the altitude array
        zaltkm = numpy.array([AltitudeMin, AltitudeMax, deltaAltitude])

        # set the date time
        dt = datetime.datetime.fromtimestamp(tUnix)#.astimezone(datetime.timezone.utc)
        print(dt)
        


        # pull the harding move, move to the directory run and move back
        # os.chdir(self.iridir)
        # outf,oarr = self.iri2016.iri_sub(self.JF,self.JMAG,glat,glon\
        #                             ,year,doy,utHrs,AltitudeMin,\
        #                             AltitudeMax,deltaAltitude)
        # call IRI
        jf = iricore.get_jf()
        jf[:] = True
        iriout = iricore.iri(dt,zaltkm,glat,glon,version=16,jf=jf)

        print(iriout.edens[0:-1])
        outDict = dict()
        outDict['Ne'] = iriout.edens[0:-1]
        outDict['Te'] = iriout.etemp[0:-1]
        outDict['Ti'] = iriout.itemp[0:-1]
        outDict['O+'] = iriout.o[0:-1]
        outDict['O2+'] = iriout.o2[0:-1]
        outDict['NO+'] = iriout.no[0:-1]
        outDict['N+'] = iriout.n[0:-1]
        outDict['Altitude'] = iriout.height[0:-1]

        return outDict
    
if __name__ == "__main__":
    t1970 = datetime.datetime(1970,1,1,00,00,00)
    # IYYYY = 2017
    # MMDD = -01
    # DHour = 10.25
    Lat = 45.
    Lon = 228.
    Heibeg = 60.
    Heiend = 500.
    step = 0.5
    t1 = datetime.datetime(2018,1,1,1,0,0)
    iri = IRI2016()
    import pylab as plt
    plt.figure()
    for i in range(24):
        tUnix = (t1-t1970).total_seconds()+i*3600.

        # tUnix,glat,glon,AltitudeMin,AltitudeMax,deltaAltitude
        testDict = iri.IRI2016(tUnix,Lat,Lon,Heibeg,Heiend,step)

        plt.semilogx(testDict['Ne'], testDict['Altitude'])
    plt.xlim([1e8,1e12])
    plt.show()
    print(testDict['Ne'])

    # tmpNO = testDict['O+']
    # print(tmpNO)
    # tmpNO[tmpNO < 0.] = 0.
    # print(tmpNO)
