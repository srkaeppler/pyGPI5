import numpy
import sys
import datetime
import os

class IRI2016:

    def __init__(self):

        self.iridir = '/Users/srkaeppler/research/data/AFOSR_Eregion_Conductivity/Models/iri2016'

#        self.iridir = '/Users/srkaeppler/Dropbox/research/data/NSF_Dregion_ParticlePrecipitation/Models/iri2016'

        sys.path.append(self.iridir)
        self.cwd = os.getcwd()
        import iri2016
        self.iri2016 = iri2016
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
        zaltkm = numpy.arange(AltitudeMin, AltitudeMax, deltaAltitude)

        # set the time variables
        year = int(tUnix/(24.*3600.*365)+1970.)
        doy = -1*int((tUnix/(24.*3600))%365)
        utHrs = ((tUnix/3600.)%24)+25. # see IRIsub.for instructions


        # pull the harding move, move to the directory run and move back
        os.chdir(self.iridir)
        outf,oarr = self.iri2016.iri_sub(self.JF,self.JMAG,glat,glon\
                                    ,year,doy,utHrs,AltitudeMin,\
                                    AltitudeMax,deltaAltitude)
        os.chdir(self.cwd)
        outDict = dict()
        outDict['Ne'] = outf[0,0:zaltkm.shape[0]]
        outDict['Te'] = outf[3,0:zaltkm.shape[0]]
        outDict['Ti'] = outf[2,0:zaltkm.shape[0]]
        outDict['O+'] = outf[4,0:zaltkm.shape[0]]
        outDict['O2+'] = outf[7,0:zaltkm.shape[0]]
        outDict['NO+'] = outf[8,0:zaltkm.shape[0]]
        outDict['N+'] = outf[10,0:zaltkm.shape[0]]
        outDict['Altitude'] = zaltkm

        return outDict

if __name__ == "__main__":
    t1970 = datetime.datetime(1970,1,1,00,00,00)
    # IYYYY = 2017
    # MMDD = -01
    # DHour = 10.25
    Lat = 45.
    Lon = 228.
    Heibeg = 65.
    Heiend = 500.
    step = 0.5
    t1 = datetime.datetime(2010,1,1,1,0,0)
    iri = IRI2016()
    import pylab as plt
    plt.figure()
    for i in range(24):
        tUnix = (t1-t1970).total_seconds()+i*3600.

        # tUnix,glat,glon,AltitudeMin,AltitudeMax,deltaAltitude
        testDict = iri.IRI2016(tUnix,Lat,Lon,Heibeg,Heiend,step)

        plt.semilogx(testDict['Ne'], testDict['Altitude'])
    plt.show()
    print(testDict)

    tmpNO = testDict['O+']
    print(tmpNO)
    tmpNO[tmpNO < 0.] = 0.
    print(tmpNO)
