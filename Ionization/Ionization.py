 # coding: utf-8
import numpy
import sys

# sys.path.append('../Models/')
# import MSIS
from scipy.constants import G as G
from scipy.constants import Boltzmann as kb
import datetime
import matplotlib.pyplot as plt
import scipy.integrate

class Ionization:
    """

    """
    def __init__(self):
        """
        Implementing Fang et al., 2010
        """
        self.MSISDir = '/Users/srkaeppler/research/data/pyGPI5/Models/'
        sys.path.append(self.MSISDir)
        import MSIS
        self.Pij = numpy.zeros([7,4])
        self.Pij = numpy.array([[1.24616e0,1.45903e0, -2.42269e-1, 5.95459e-2], \
                                [2.23976e0, -4.22918e-7, 1.36458e-2, 2.53332e-3], \
                                [1.41754e0, 1.44597e-1, 1.70433e-2, 6.39717e-4],\
                                [2.48775e-1,-1.50890e-1, 6.30894e-9, 1.23707e-3], \
                                [-4.65119e-1, -1.05081e-1, -8.95701e-2, 1.22450e-2], \
                                [3.86019e-1, 1.75430e-3, -7.42960e-4, 4.60881e-4], \
                                [-6.45454e-1, 8.49555e-4, -4.28581e-2, -2.99302e-3], \
                                [9.48930e-1, 1.97385e-1, -2.50660e-3, -2.06938e-3]])
        self.PijMax = numpy.zeros([7,4])
        self.PijMax = numpy.array([[3.49979e-1,-6.18200e-2, -4.08124e-2, 1.65414e-2], \
                                [5.85425e-1, -5.00793e-2, 5.69309e-2, -4.02491e-3], \
                                [1.69692e-1, -2.58981e-2, 1.96822e-2, 1.20505e-3],\
                                [-1.2227e-1,-1.15532e-2, 5.37951e-6, 1.20189e-3], \
                                [1.57018e0, 2.87896e-1, -4.14857e-1, 5.18158e-2], \
                                [8.83195e-1, 4.31402e-2, -8.33599e-2, 1.02515e-2], \
                                [1.90953e0, -4.74704e-2, -1.80200e-1, 2.46652e-2], \
                                [-1.29566e0, -2.10952e-1, 2.73106e-1, -2.92752e-2]])
        self.msis = MSIS.MSIS()
        return

    def GravitationalAcceleration(self, z):
        """
        Input Altitude in meters!
        G = m^3*kg^-1*s^-2
        units of meters
        """
        Me = 5.97237e24 # kg wikipedia source
        Re = 6.371e6 # m

        # probably good enough for altitude
        g = G*Me/(Re+z)**2

        return g


    def CalculateCi(self, E0):
        """
        Fang equation 5
        """
        Ci = numpy.zeros(8)
        x = numpy.log(E0)
        for i in range(8):
            y = self.Pij[i,0]+self.Pij[i,1]*x+self.Pij[i,2]*x**2+self.Pij[i,3]*x**3
            Ci[i] = numpy.exp(y)
        return Ci

    def CalculateCiMaxwellian(self, E0):
        """
        Fang equation 5
        """
        Ci = numpy.zeros(8)
        x = numpy.log(E0)
        for i in range(8):
            y = self.PijMax[i,0]+self.PijMax[i,1]*x+self.PijMax[i,2]*x**2+self.PijMax[i,3]*x**3
            Ci[i] = numpy.exp(y)
        return Ci

    def FangModel(self, E0,Q0, altkm, Hz,MassDensity):
        """
        Input:  E = MonoEnergetic Energy (eV)
                altkm = Altitude (km)
                Hz (m)
                Mass Density (kg/m^3)
        """
        # do all of the conversions here
        zcm = altkm*(1000.*100.) # cm
        EkeV = E0/1000. # KeV
        rhoz = MassDensity*(1e3/1e6) # g/cm^3
        Hzcm = Hz*100.
        deltaE = 0.035 #keV
        # Q0 = Q0/1000.
        Q0 = Q0*6.242e11/1000. # conversion to keV/cm2/s
        # print 'Hzcm', Hzcm
        # print 'EkeV', EkeV

        y = (2./EkeV)*(rhoz*Hzcm/6e-6)**0.7 # equation 1
        Ci = self.CalculateCi(EkeV)
        f = Ci[0]*(y**Ci[1])*numpy.exp(-Ci[2]*(y**Ci[3]))+\
            Ci[4]*(y**Ci[5])*numpy.exp(-Ci[6]*(y**Ci[7])) # equation 4

        # print f
        qz = f*Q0/deltaE/Hzcm # equation 3 which is wrong, look in fang 2008, eq2
        return qz, y, f

    def FangModelMaxwellian(self, Q0,E0, altkm, Hz,MassDensity):
        """
        Input:  E = MonoEnergetic Energy (eV)
                altkm = Altitude (km)
                Hz (m)
                Mass Density (kg/m^3)
            Fang 2008 paper
        """

        #print 'In Fang Model Maxwellian Q0 E0', Q0,E0
        # do all of the conversions here
        zcm = altkm*(1000.*100.) # cm
        EkeV = E0/1000. # KeV
        rhoz = MassDensity*(1e3/1e6) # g/cm^3
        Hzcm = Hz*100.
        deltaE = 0.035 #keV
        # Q0 = Q0/1000.
        Q0 = Q0*6.242e11/1000. # conversion to keV/cm2/s
        # print 'Hzcm', Hzcm
        # print 'EkeV', EkeV

        y = (1./EkeV)*(rhoz*Hzcm/4e-6)**0.606 # equation 4
        Ci = self.CalculateCiMaxwellian(EkeV)
        f = Ci[0]*(y**Ci[1])*numpy.exp(-Ci[2]*(y**Ci[3]))+\
            Ci[4]*(y**Ci[5])*numpy.exp(-Ci[6]*(y**Ci[7])) # equation 4

        # print f
        qz = (Q0*f)/(2.*deltaE*Hzcm) # equation 3 which is wrong, look in fang 2008, eq2
        return qz, y, f


    def FangModelMatrix(self, Q0,E0, altkm, Hz,MassDensity):
        """
        Input:  E = MonoEnergetic Energy (eV)
                altkm = Altitude (km)
                Hz (m)
                Mass Density (kg/m^3)
        """
        # do all of the conversions here
        zcm = altkm*(1000.*100.) # cm
        EkeV = E0/1000. # KeV
        rhoz = MassDensity*(1e3/1e6) # g/cm^3
        Hzcm = Hz*100.
        deltaE = 0.035 #keV
        # Q0 = Q0/1000.
        Q0 = Q0*6.242e11/1000. # conversion to keV/cm2/s
        # print 'Hzcm', Hzcm
        # print 'EkeV', EkeV

        y = (2./EkeV)*(rhoz*Hzcm/6e-6)**0.7 # equation 1
        Ci = self.CalculateCi(EkeV)
        f = Ci[0]*(y**Ci[1])*numpy.exp(-Ci[2]*(y**Ci[3]))+\
            Ci[4]*(y**Ci[5])*numpy.exp(-Ci[6]*(y**Ci[7])) # equation 4

        # print f
        qz = f*Q0/deltaE/Hzcm # equation 3 which is wrong, look in fang 2008, eq2
        return qz, y, f




    def RunMSISFang(self, tUnix, glat,glon,altkm = numpy.arange(80,150,1)):
        """
        This should only have to be run once to calculate what I need
        Input: unix time, longitude and latitude, altitude grid
        Do everything in SI
        """
        year = int(tUnix/(24.*3600.*365)+1970.)
        doy = int((tUnix/(24.*3600))%365)
        utHrs = (tUnix/3600.)%24
        CGSorSI = 'SI'
        outDictSI = self.msis.MSIS(doy,utHrs,glat,glon,year,altkm=altkm, CGSorSI=CGSorSI)
        Tn = outDictSI['Tn'] # in Kelvin
        MassDensity = outDictSI['MassDensity'] # kg/m^3
        AverageMass = outDictSI['AverageMass'] # kg
        # print 'Mass Density', MassDensity
        # print 'AverageMass', AverageMass/1.6726219e-27
        # have to have a what to do if not sent in
        gz = self.GravitationalAcceleration(altkm/1000.)
        Hz = kb*Tn/(gz*AverageMass)


        return MassDensity, Hz

    def MaxwellianFlux(self, E, Q0, E0):
        """
        Input:  E energy array (eV)
                Q0 ergs/cm^2/s
                E0 eV

        Output: ergs/cm^2/s
        """

        NumFlux = Q0*E*numpy.exp(-E/E0)/2./E0**3
        EnergyFlux = E*E*NumFlux

        return NumFlux,EnergyFlux

    def MakeFangModelMatrix(self, EeV,tUnix, glat,glon,altkm = numpy.arange(80,150,1)):
        """
        This should only have to be run once to calculate what I need
        Input: unix time, longitude and latitude, altitude grid
        Do everything in SI
        """
        year = int(tUnix/(24.*3600.*365)+1970.)
        doy = int((tUnix/(24.*3600))%365)
        utHrs = (tUnix/3600.)%24
        CGSorSI = 'SI'
        outDictSI = self.msis.MSIS(doy,utHrs,glat,glon,year,altkm=altkm, CGSorSI=CGSorSI)
        Tn = outDictSI['Tn'] # in Kelvin
        MassDensity = outDictSI['MassDensity'] # kg/m^3
        AverageMass = outDictSI['AverageMass'] # kg
        # print 'Mass Density', MassDensity
        # print 'AverageMass', AverageMass/1.6726219e-27
        # have to have a what to do if not sent in
        gz = self.GravitationalAcceleration(altkm/1000.)
        Hz = kb*Tn/(gz*AverageMass)

        # do all of the conversions here
        zcm = altkm*(1000.*100.) # cm
        EkeV = EeV/1000. # KeV
        rhoz = MassDensity*(1e3/1e6) # g/cm^3
        Hzcm = Hz*100.
        deltaE = 0.035 #keV
        # Q0 = Q0/1000.
        # Q0 = Q0*6.242e11/1000. # conversion to keV/cm2/s
        # print 'Hzcm', Hzcm
        # print 'EkeV', EkeV

        # A is a matrix that is  N altitude elements by M energy elements
        A = numpy.zeros([altkm.shape[0],EkeV.shape[0]-1,])


        for iEnergy in range(EkeV.shape[0]-1):
            y = (2./EkeV[iEnergy])*(rhoz*Hzcm/6e-6)**0.7 # equation 1
            Ci = self.CalculateCi(EkeV[iEnergy])
            f = Ci[0]*(y**Ci[1])*numpy.exp(-Ci[2]*(y**Ci[3]))+\
                Ci[4]*(y**Ci[5])*numpy.exp(-Ci[6]*(y**Ci[7])) # equation 4
            dE = EkeV[iEnergy+1]-EkeV[iEnergy]
            tmpqz = (f*EkeV[iEnergy]*dE)/deltaE/Hzcm
            A[:,iEnergy] = tmpqz

        # print f
        # qz = f*Q0/deltaE/Hzcm # equation 3 which is wrong, look in fang 2008, eq2
        return A

    def Ionization(self,E,EnergyFlux,AltMin,AltMax,AltStep,tUnix,glat,glon, IonizationType='Fang'):
        # E and EnergyFlux == same shape
        """
        Input: E energy array (eV)
                EnergyFlux (erg/cm^2/s or mW/m^2 - same units)
                AltMin (km)
                AltMax (km)
                AltStep (km)
                msisIn dictionary containing what need to run MSIS
        """

        # run msis
        altkm = numpy.arange(AltMin,AltMax,AltStep)
        qZE = numpy.zeros([altkm.shape[0],E.shape[0]])
        qZ = numpy.zeros(altkm.shape[0])
        qZsimp = numpy.zeros(altkm.shape[0])
        if IonizationType == 'Fang':
            MassDensity,Hz = self.RunMSISFang(tUnix,glat,glon,altkm=altkm)
            for iE in range(len(E)):
                tmpqz,y,f = self.FangModel(E[iE],EnergyFlux[iE],altkm,Hz,MassDensity)
                qZE[:,iE] = tmpqz


        # this is slightly different
        # in this case would put in Q0 and E0 for each energy
        if IonizationType == 'Maxwellian':
            MassDensity,Hz = self.RunMSISFang(tUnix,glat,glon,altkm=altkm)
            # FangModelMaxwellian(self, Q0,E0, altkm, Hz,MassDensity):
            for iE in range(len(E)):
                tmpqz,y,f = self.FangModelMaxwellian(EnergyFlux[iE],E[iE],altkm,Hz,MassDensity)
                qZE[:,iE] = tmpqz



        # sum over energy to get the final altitude profile
        # need to multiply by dE
        for ih in range(altkm.shape[0]):
            qZsimp[ih] = scipy.integrate.simps(qZE[ih,:], E/1000.)
        # this needs to be normalized by a dimensionless quantity or I am not doing something
        # ahead of time
        qZ = numpy.sum(qZE, axis=1)

        return qZ, qZE.T,qZsimp

    def IonizationFang2008Maxwellian(self,Q0,E0,AltMin,AltMax,AltStep,tUnix,glat,glon):
        # E and EnergyFlux == same shape
        """
        Input: E0 - characteristic energy in eV
                Q0 (erg/cm^2/s or mW/m^2 - same units)
                AltMin (km)
                AltMax (km)
                AltStep (km)
                msisIn dictionary containing what need to run MSIS
        """

        # def FangModelMaxwellian(self, E0,Q0, altkm, Hz,MassDensity):
        #     """
        #     Input:  E = MonoEnergetic Energy (eV)
        #             altkm = Altitude (km)
        #             Hz (m)
        #             Mass Density (kg/m^3)
        #         Fang 2008 paper
        # run msis
        # print 'IonizationFang Function Q0, E0', Q0, E0
        altkm = numpy.arange(AltMin,AltMax,AltStep)
        MassDensity,Hz = self.RunMSISFang(tUnix,glat,glon,altkm=altkm)
        qZ,y,f = self.FangModelMaxwellian(Q0,E0,altkm,Hz,MassDensity)



        return qZ







if __name__ == '__main__':
    iz = Ionization()
    altkm = numpy.arange(50,400,1)
    Q0 = 1.0 #erg/cm^2/s
    E0 = 1*1e3 # 10 keV
    t1970 = datetime.datetime(1970,1,1,0,0,0)
    t2010 = datetime.datetime(2010,1,2,10,0,0)
    tunix = (t2010-t1970).total_seconds()
    # glat = 45
    # glon = 45
    glat = 45
    glon = 228
    MassDensity,Hz = iz.RunMSISFang(tunix,glat,glon,altkm=altkm)
    qz,y,f = iz.FangModel(E0,Q0,altkm,Hz,MassDensity)
    plt.figure(100)
    plt.semilogx(qz,altkm)
    plt.xlim([1e1,1e5])
    plt.ylim([50,400])
    plt.grid()
    plt.title('Ionization')

    Earr = numpy.logspace(2,5,31)
    plt.figure(1001)
    for iE in Earr:
        qz,y,f = iz.FangModel(iE,Q0,altkm,Hz,MassDensity)
        plt.semilogx(qz,altkm)
        plt.xlim([1e1,1e5])
        plt.ylim([50,400])
        plt.grid()

    plt.figure(101)
    plt.semilogy(f,y)
    plt.xlim([0,1])
    plt.ylim([0.1,10.])
    plt.title('f vs y')

    EeV = numpy.logspace(2,5,num=49)
    Earr = numpy.logspace(2,5,31)
    plt.figure(2001)
    for iE in Earr:
        Q0 = 1.0
        E0 = iE#5.0*1e3
        NumFlux,QeV = iz.MaxwellianFlux(EeV,Q0,E0)
        msisDict = dict()
        msisDict['tUnix'] = tunix
        msisDict['glat'] = glat
        msisDict['glon'] = glon
    # E,EnergyFlux,AltMin,AltMax,AltStep,msisIn, IonizationType='Fang'):
        qZ, qZE,qZSimp = iz.Ionization(EeV,QeV,50,400,1,tunix,glat,glon)
        plt.semilogx(qZSimp, altkm)
        plt.xlim([1e-2,1e5])
        plt.ylim([50,400])
        plt.grid()


    plt.figure(105)
    plt.plot(EeV/1000., NumFlux*(6.242e11*1000))
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1e1,1e7])

    altkm = numpy.arange(50,400,1)
    plt.figure(103)
    for i in range(qZE.shape[0]):
        plt.semilogx(qZE[i,:], altkm)
    plt.semilogx(qZ,altkm, 'k-', lw=2)
    plt.semilogx(qZSimp,altkm,'r-',lw=2)
    plt.xlim([1e2,1e5])
    plt.savefig('LogSpace_201points.png')
    # plt.show()

    #49 points
    NumPts = numpy.array([49,101,201,501])
    color = ['red','green','blue','purple']
    kk = 0
    plt.figure(104)
    for ii in NumPts:
        # print ii
        EeV = numpy.logspace(2,5,num=ii)
        Q0 = 1.0
        E0 = 5.0*1e3
        NumFlux,QeV = iz.MaxwellianFlux(EeV,Q0,E0)
        qZ, qZE,qZSimp = iz.Ionization(EeV,QeV,50,400,1,tunix,glat,glon)
        A = iz.MakeFangModelMatrix(EeV,tunix,glat,glon,altkm=altkm)
        qZMat = numpy.dot(A,NumFlux[0:-1]*6.242e11*1000.*10.) # this factor of 10 is adhoc...

        # print 'qZMat', qZMat.shape
        # print 'qZMat',qZMat
        # print 'qZSimp', qZSimp
        # plt.semilogx(qZ,altkm, '-', lw=2, color=color[kk])
        plt.semilogx(qZSimp,altkm,'--',lw=2,color=color[kk])
        plt.semilogx(qZMat,altkm, '+',color=color[kk] )
        plt.xlim([1e3,1e6])
        kk=kk+1

    # MakeFangModelMatrix(self, EeV,tUnix, glat,glon,altkm = numpy.arange(80,150,1)):
    EeV = numpy.logspace(2,6,num=51)
    A = iz.MakeFangModelMatrix(EeV,tunix,glat,glon,altkm=altkm)
    # print A, A.shape
    Q0 = 2.0
    E0 = 5.0*1e3
    NumFlux,QeV = iz.MaxwellianFlux(EeV,Q0,E0)
    # print 'NumFlux', NumFlux
    # print 'QeV', QeV
    # print scipy.integrate.simps(QeV,EeV)
    # print scipy.integrate.simps(NumFlux*EeV,EeV)
    dE = numpy.diff(EeV)
    # print dE.shape, EeV.shape
    # print numpy.sum(NumFlux[0:-1]*EeV[0:-1]*dE)
    plt.show()
