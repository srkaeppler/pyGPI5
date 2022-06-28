"""
E and D-region chemistry

11-28-2017:
A bunch of these reaction rates come from ionocoeffs_spec.m which was written
by Lehtinen and given to me by Marshall.
I am only doing this in python for convenince, they are the original authors

I will retain (e.g., copy) comments from the other code and also the header
what may be differen is the input/output I choose to use.  I could Research
if there are better reaction rates to use.

I am also making a pretty strategic decision to name the functions after
the variable names found in the GPI paper equations 1-4
"""

import numpy
import scipy.integrate
import copy
import sys
sys.path.append('../Models')
import IRI2016
iri2016 = IRI2016.IRI2016()
import MSIS
msis = MSIS.MSIS()

import matplotlib.pyplot as plt

class Chemistry:

    def __init__(self, SteadyStateTime = 1.e4,ISRIntegrationTime=1.e4):
        self.NeIn = None
        self.Sin = None
        self.Tn = None
        self.Nm = None
        self.DregionChem = None
        self.altkm = None
        self.MSISDict = None
        self.y0 = None
        self.MSISatRef = None # to do cosmic ray calculation
        self.SteadyState = SteadyStateTime
        self.nSteps = 10000
        self.ISRIntTime = ISRIntegrationTime#60. # seconds
        return



    def Recombination_Te(self,Te, CO2 = None, CNO=None,CO=None):
        """

        """

        alpha = numpy.zeros(Te.shape[0])
        # from Schunk and Nagy
        # need to check units
        # units output in cm^3*s^-1
        # Te in Kelvin

        #NO
        a_NO = (4.e-7)*((300./Te)**(0.70))

        #O2+
        a_O2 = (2.4e-7)*((300./Te)**(0.5))

        #O
        a_O = (3.7e-12)*((250./Te)**0.7)


        CNO[CNO<0] = 1e-8
        CO[CO < 0] = 1e-8
        CO2[CO2 < 0] = 1e-8

        if len(CO2) > 0:
            assert (CO2.shape[0] == Te.shape[0], "CO Shape is not the same as Te")
            alpha = (CO2/100.)*a_O2 + alpha
        if len(CNO) > 0:
            assert (CNO.shape[0] == Te.shape[0], "CNO Shape is not the same as Te")
            alpha = (CNO/100.)*a_NO + alpha
        if len(CO) > 0:
            assert (CO.shape[0] == Te.shape[0], "CO Shape is not the same as Te")
            alpha = (CO/100.)*a_O + alpha

        return alpha

    def Recombination_Altitude(self, altkm):
        # units cm^3*s^-1
        alpha = (2.5e-6)*numpy.exp(-altkm/51.2)

        return alpha

    def Calculate_alphaD(self,Nm):
        """
        effective coefficient of dissociative recommbination
        Kelley calls this: dissociative recombination
        % alfad -- eff. coeff. of dissociative recombination
        % --------------------------------------------------
        % [Kelley p. 14]
            alphaD = 4e-7 cm^3s^-1
        [GPI]
        % alfad=1e-7 -- 3e-7
        % [PI]
        %   alfad=3e-7
        % [RI]
        %   alfad = alfa(NO+) N(NO+)/Npos + alfa(O2+) N(O2+)/Npos
        % where
        %   alfa(NO+)=4.1e-7*(300/Tn^.5/Te^.5)
        %   alfa(O2+)=2.1e-7*(300^.7/Tn^.1/Te^.6)
        % We take temperatures (nighttime, [PI])
        %   Te=Tn=200 K
        % The relative concentration of positive ions is
        %   N(NO+)/Npos=0.84 at 100km (IRI-95)
        % (dominates) at night at 45 deg latitude
        % so approximately alfad=5.6e-7
        """
        alphaD=6e-7*numpy.ones(Nm.shape[0])
        return alphaD

    def Calculate_alphaDC(self,Nm):
        """
        effective coefficient of recombation of electrons with positive
        cluster ions
        Kelley calls this: Cluster-ion recombination

         alfadc -- eff. coeff. of recombination of Ne with Nclus
        % -----------------------------------------------------------------
        [Kelley p. 14], alphad = 5e-5 cm^3s^-1
        % [GPI, PI] alfadc=1e-5
        % [RI] alfadc=3e-6*(Tn/Te)^(0--.08)
        """
        alphaDC=1e-5*numpy.ones(Nm.shape[0])
        return alphaDC

    def Calculate_alphaI(self,Nm):
        """
        Effective coefficent of ion-ion recombination
        for all kinds of positive ions with negative ions

        Kelley: ion-ion mutual neutralization
        Kelley Value: 10e-6 cm^3 s^-1

        % alfai -- eff. coeff. of mutual neutralization
        % -------------------------------------------------------
        % [GPI, PI, RI] alfai=1e-7
        % [MH] alfai ~ Nm at h<50 km
        alfai_opt=getvaluefromdict(options,'alfai','MH');
        switch alfai_opt
            case 'MH'
                alfai=1e-7*ones(size(Nm))+1e-24*Nm;
            case 'GPI'
                alfai=1e-7*ones(size(Nm));
            otherwise
                error(['unknown option (alfai) = ' alfai_opt]);
        end
        """
        print Nm.shape[0]
        alphaI = 1e-7*numpy.ones(Nm.shape) + Nm*1e-24
        return alphaI

    def Calculate_Beta(self,Nm):
        """
        Kelley: effective electron attachment rate
        # % beta -- eff. electron attachment rate
        # % ------------------------------------------
        # % Assume N(O2)=0.2 Nm; N(N2)=0.8 Nm, T=Tn=Te=200 K
        # % Simplify:
        # % [GPI,PI,RI] beta=1e-31*N(N2)*N(O2)+k(O2)*(N(O2))^2;
        # % [GPI,PI] k(O2)=1.4e-29*(300/T)*exp(-600/T) == 1e-30
        # %   => beta=5.6e-32*Nm^2
        # % [RI] k(O2)=(see function kO2RI below)=1.5e-30
        # %   => beta=7.6e-32*Nm^2
        """
        beta=5.6e-32*Nm**2;
        return beta

    def Calculate_B(self,Nm):
        """
        % Bcoef -- effective rate of conversion of Npos into Nclus
        % --------------------------------------------------------
        % [GPI,PI,RI] Bcoef=1e-31*Nm^2
        """
        Bcoef=1e-31*Nm**2;
        return Bcoef

    def Calculate_Gamma(self,Nm,Tn = 0, option='GPI', daytime=False):
        """
        % gamma -- eff. electron detachment rate
        % -------------------------------------------

        """
        if option == 'Temp':
            """
            % The new rate [A] of reaction
            %  O2- + O2 -> e + 2 O2
            % and photodetachment [Gur] 0.44 s^{-1}
            """
            gamma=8.61e-10*numpy.exp(-6030./Tn)*Nm
            if daytime:
                gamma = gamma+0.44

        elif option == 'Kozlov':
            """
            % [Koz] Rate of
            %  O2- + O2 -> e + 2 O2
            % [Koz] has photodetachment=0.33 s^{-1}
            """
            gamma=2.7e-10*Nm*numpy.exp(-5590./Tn)*numpy.sqrt(Tn/300)
            if daytime:
                gamma = gamma+0.33

        elif option == 'GPI':
            gamma = 3e-17*Nm
        else:
            gamma = 3e-17*Nm

        return gamma

    def Calculate_GammaX(self,Nm, daytime=False):
        """
        Special stuff Lehtinen put in - this is his IP.
        % gammaX -- detachment rate from the slow negative ions (NO3-)
        """
        if daytime:
            gammaX = 0.002*numpy.ones(Nm.shape[0])
        else:
            gammaX = numpy.zeros(Nm.shape[0])
        return gammaX

    def Calculate_Xbar(self,Nm):
        """
        % Xbar -- the rate of conversion of Nneg (mostly O2-) into NX (N03-)
        % ------------------------------------------------------------------
        % See [M]
        % Xbar=1e-30*N(0_2)*N(M)+3e-10*N(O3), we neglect the ozone.
        """

        Xbar = 0.2e-30*Nm**2
        return Xbar


    def Calculate_Dregion_ReactionRates(self,Nm, Tn=None, \
                                        options = {'GammaType':'GPI', 'daytime':False}, \
                                        iriDict = None):

        # should probably assert that Nm is a 1-D array
        outDict = dict()
        outDict['alphaD'] = self.Calculate_alphaD(Nm)
        outDict['alphaDC'] = self.Calculate_alphaDC(Nm)
        outDict['alphaI'] = self.Calculate_alphaI(Nm)
        outDict['beta'] = self.Calculate_Beta(Nm)
        outDict['B'] = self.Calculate_B(Nm)

        if iriDict is not None:
            outDict['alphaD'] = self.Recombination_Te(iriDict['Te'], \
                                                    CO2 = iriDict['O2+'],\
                                                    CNO = iriDict['NO+'],\
                                                    CO = iriDict['O+'])
            print 'alphaD,', outDict['alphaD']



        if len(Tn) > 1 & (options['GammaType'] != 'Temp'):
            outDict['gamma'] = self.Calculate_Gamma(Nm,Tn=Tn,\
                                                    option=options['GammaType'], \
                                                    daytime=options['daytime'])
        else:
            outDict['gamma'] = self.Calculate_Gamma(Nm)

        outDict['gammaX'] = self.Calculate_GammaX(Nm, daytime=options['daytime'])
        outDict['Xbar'] = self.Calculate_Xbar(Nm)

        return outDict

    def Dregion_Chemistry_5species(self,t,N,params):
        """
        Differential equation solver for 5 species
        This is a direct translation from the code


        % The units in the source papers are m^{-3} [RI] and cm^{-3} [all other]
        %
        % The equations to be solved are (the prime denotes time derivative)
        %
        %  Ne'    = S + gamma*Nneg - beta*Ne - alfad*Ne*Npos - alfadc*Ne*Nclus + gammaX*NX;
        %  Nneg'  = beta*Ne - gamma*Nneg - alfai*Nneg*(Npos+Nclus) - Xbar*Nneg;
        %  Nclus' = - alfadc*Ne*Nclus + Bcoef*Npos - alfai*(Nneg+NX)*Nclus;
        %  NX'    = Xbar*Nneg - gammaX*NX - alfai*NX*(Npos+Nclus);
        %  Npos'  = S - Bcoef*Npos - alfad*Ne*Npos - alfai*(Nneg+NX)*Npos;
        %
        % where
        %
        %  S     - external source of ionization
        %  Ne    - electron density
        %  Npos  - positive ion density
        %  Nneg  - negative ions, from which electrons are easily detached (O2-, ...)
        %  Nclus - positive cluster ion density (hydrated Npos, H+(H2O)n)
        %  NX    - negative ions, from which electrons are not detached (NO3-, ...)
        %          (new in 5-species model)
        %
        % For neutrality, it is necessary that Ne+Nneg+NX=Npos+Nclus, so that the
        % last equation is dependent.

        """
        S=params[0];
        alfad=params[1];
        alfadc=params[2];
        alfai=params[3];
        gamma=params[4];
        beta=params[5];
        Bcoef=params[6];
        gammaX=params[7];
        Xbar=params[8];
        Np=numpy.zeros(4);
        Ne=N[0]; Nneg=N[1]; Nclus=N[2]; NX=N[3];
        Npos=Ne+Nneg+NX-Nclus;
        Np[0] = S + gamma*Nneg - beta*Ne - alfad*Ne*Npos - alfadc*Ne*Nclus + gammaX*NX;
        Np[1] = beta*Ne - gamma*Nneg - alfai*Nneg*(Npos+Nclus) - Xbar*Nneg;
        Np[2] = - alfadc*Ne*Nclus + Bcoef*Npos - alfai*(Nneg+NX)*Nclus;
        Np[3] = Xbar*Nneg - gammaX*NX - alfai*NX*(Npos+Nclus);
        return Np



    def Dregion_Chemistry_4species(self,t,N,params):
        """
        Differential equation solver for 4 species
        This is a direct translation from the GPI papers
        equations 1-4


        % The units in the source papers are m^{-3} [RI] and cm^{-3} [all other]
        %
        % The equations to be solved are (the prime denotes time derivative)
        %
        %  Ne'    = S + gamma*Nneg - beta*Ne - alfad*Ne*Npos - alfadc*Ne*Nclus + gammaX*NX;
        %  Nneg'  = beta*Ne - gamma*Nneg - alfai*Nneg*(Npos+Nclus) - Xbar*Nneg;
        %  Nclus' = - alfadc*Ne*Nclus + Bcoef*Npos - alfai*(Nneg+NX)*Nclus;
        %  NX'    = Xbar*Nneg - gammaX*NX - alfai*NX*(Npos+Nclus);
        %  Npos'  = S - Bcoef*Npos - alfad*Ne*Npos - alfai*(Nneg+NX)*Npos;
        %
        % where
        %
        %  S     - external source of ionization
        %  Ne    - electron density
        %  Npos  - positive ion density
        %  Nneg  - negative ions, from which electrons are easily detached (O2-, ...)
        %  NX    - negative ions, from which electrons are not detached (NO3-, ...)
        %          (new in 5-species model)
        %
        % For neutrality, it is necessary that Ne+Nneg+NX=Npos+Nclus, so that the
        % last equation is dependent.

        """
        S=params[0];
        alfad=params[1];
        alfadc=params[2];
        alfai=params[3];
        gamma=params[4];
        beta=params[5];
        Bcoef=params[6];
        gammaX=params[7];
        Xbar=params[8];
        Np=numpy.zeros(3);
        Ne=N[0]; Nneg=N[1]; Nclus=N[2];
        Npos=Ne+Nneg-Nclus # equation 3
        Np[0] = S + gamma*Nneg - beta*Ne - alfad*Ne*Npos - alfadc*Ne*Nclus
        Np[1] = beta*Ne - gamma*Nneg - alfai*Nneg*(Npos+Nclus)
        Np[2] = - alfadc*Ne*Nclus + Bcoef*Npos - alfai*(Nneg)*Nclus;
        return Np




    def Integrate_ODE(self,NeIn,Sin,ChemistryDict,IntType='5species'):
        """
        Similar to testodeintegrate.py in the matlab folder
        This step could be parallelized
        should have a statement which checks the parameters.
        """

        if IntType == '4species':
            ode15s = scipy.integrate.ode(self.Dregion_Chemistry_4species)
            ode15s.set_integrator('vode', method='bdf', order=15, nsteps=self.nSteps)
            y0 = numpy.zeros(3)
        elif IntType == '5species':
            # configuring
            ode15s = scipy.integrate.ode(self.Dregion_Chemistry_5species)
            ode15s.set_integrator('vode', method='bdf', order=15, nsteps=self.nSteps)
            y0 = numpy.zeros(4)

        args = numpy.zeros(9)
        N = ChemistryDict['B'].shape[0]
        outDict = dict()
        # Ne+Nneg+NX=Npos+Nclus
        outDict['Ne'] = numpy.zeros(N)
        outDict['NnegIon'] = numpy.zeros(N)
        outDict['NX'] = numpy.zeros(N)
        outDict['NposIon'] = numpy.zeros(N)
        outDict['NposCluster'] = numpy.zeros(N)
        print 'NeIn.shape', NeIn['Ne'].shape, N
        for iz in range(N):

            args[0] = Sin[iz];
            args[1] = ChemistryDict['alphaD'][iz];
            args[2] = ChemistryDict['alphaDC'][iz];
            args[3] = ChemistryDict['alphaI'][iz]
            args[4] = ChemistryDict['gamma'][iz]
            args[5] = ChemistryDict['beta'][iz]
            args[6] = ChemistryDict['B'][iz]
            args[7] = ChemistryDict['gammaX'][iz]
            args[8] = ChemistryDict['Xbar'][iz]

            y0[0] = NeIn['Ne'][iz]
            y0[1] = NeIn['NnegIon'][iz]
            y0[2] = NeIn['NposCluster'][iz]
            if IntType == '5species':
                y0[3] = NeIn['NX'][iz]

            ode15s.set_initial_value(y0,0.).set_f_params(args)
            results = ode15s.integrate(self.ISRIntTime)
            # print 'iz', iz, results
            outDict['Ne'][iz] = results[0]
            outDict['NnegIon'][iz] = results[1]
            outDict['NposCluster'][iz] = results[2]
            if IntType == '5species':
                # Ne=N[0]; Nneg=N[1]; Nclus=N[2]; NX=N[3];
                # Npos=Ne+Nneg+NX-Nclus;
                outDict['NX'][iz] = results[3]
                outDict['NposIon'][iz] = (results[0]+results[1]+results[3])-results[2]
            if IntType == '4species':
                # Ne=N[0]; Nneg=N[1]; Nclus=N[2];
                # Npos=Ne+Nneg-Nclus # equation 3
                outDict['NposIon'][iz] = (results[0]+results[1])-results[2]

        return outDict

    def ODE(self,Sin,iz,ChemistryDict, IntType='5species'):

        if IntType == '4species':
            ode15s = scipy.integrate.ode(self.Dregion_Chemistry_4species)
            ode15s.set_integrator('vode', method='bdf', order=15, nsteps=self.nSteps)
            y0 = numpy.zeros(3)
        elif IntType == '5species':
            # configuring
            ode15s = scipy.integrate.ode(self.Dregion_Chemistry_5species)
            ode15s.set_integrator('vode', method='bdf', order=15, nsteps=self.nSteps)
            y0 = numpy.zeros(4)

        args = numpy.zeros(9)
        args[0] = Sin;
        args[1] = ChemistryDict['alphaD'][iz];
        args[2] = ChemistryDict['alphaDC'][iz];
        args[3] = ChemistryDict['alphaI'][iz]
        args[4] = ChemistryDict['gamma'][iz]
        args[5] = ChemistryDict['beta'][iz]
        args[6] = ChemistryDict['B'][iz]
        args[7] = ChemistryDict['gammaX'][iz]
        args[8] = ChemistryDict['Xbar'][iz]


        ode15s.set_initial_value(y0,0.).set_f_params(args)
        results = ode15s.integrate(self.SteadyState)

        return results

    def Binary_Search(self, NeTarget,ChemistryDict, iz):


        # print 'NeTarget, NeTest', NeTarget, NeTest

        # find zero crossing
        S = 1.0
        # Sarr = numpy.zeros(10)
        # FS = numpy.zeros(10)
        kk = 0
        while True:

            # ODE(Sin,iz,ChemistryDict
            NeTest = self.ODE(S,iz, ChemistryDict)[0]

            Sign = NeTest-NeTarget
            Sprev = copy.copy(S)
            if NeTest > NeTarget:
                S = S/2.

            if NeTest < NeTarget:
                S = 2.*S

            if kk==0:
                PrevSign = copy.copy(Sign)

            if numpy.sign(Sign) != numpy.sign(PrevSign):
                # print 'stop'
                break


            kk=+1

            if kk==100:
                print 'could not find bounds'
                break
            PrevSign = copy.copy(Sign)
            Neprev = copy.copy(NeTest)
            # print kk, S, NeTest, NeTarget, Sign, S, Sprev

        # just to check
        NeS = self.ODE(S,iz, ChemistryDict)[0]
        NeSprev = self.ODE(Sprev,iz, ChemistryDict)[0]
        # print 'S', S,NeS, NeTest
        # print 'Sprev',Sprev, NeSprev, Neprev
        # print 'NeTarget', NeTarget

        # now at this point I can do the whole dividing up and finding
        Shigh = numpy.max([S,Sprev])
        Slow = numpy.min([S,Sprev])
        Nehigh = numpy.max([NeS, NeSprev])
        Nelow = numpy.min([NeS, NeSprev])

        # print 'Nehigh', Nehigh, Nelow

        # # go from high to low
        kk = 0
        while True:

            dS = (Shigh-Slow)/2.

            S1 = Shigh-dS
            Ne = self.ODE(S1,iz, ChemistryDict)[0]
            dNe = numpy.abs(Ne-NeTarget)
            dNe2 = (Shigh-Slow)/Shigh
            # print 'S1', S1
            # print 'Ne', Ne
            # print 'Nelow, Nehigh', Nelow, Nehigh
            # print 'dNe',i, dNe, dNe2, Ne, NeTarget
            if dNe2 < 1e-6:
                # print 'stopping'
                break

            if Ne < NeTarget:
                Slow = copy.copy(S1)
                Nelow = copy.copy(Ne)

            if Ne > NeTarget:
                Shigh = copy.copy(S1)
                Nehigh = copy.copy(Ne)

            kk=+1
            if kk == 100:
                print 'too many iterations in binary search'
                break

        # output phase
        NeOut = Ne
        Sout = S1
        # calculate other parameters
        y0 = self.ODE(S1,iz, ChemistryDict)
        return Sout, NeOut, y0

    def Calculate_Background_Ionization(self, altkm,IRIin, ChemistryDict):

        # The altitudes MUST MUST MUST be the same as the IRIGird as the Chemistry

        Sout = numpy.zeros(altkm.shape[0])
        yInitial = numpy.zeros([altkm.shape[0],5])


        for iz in range(altkm.shape[0]):
            # print 'iz, IRIiz', iz, IRIin[iz]
            if IRIin[iz] < 0:

                continue
            else:
                tmpSout, tmpNeOut, y0 = self.Binary_Search(IRIin[iz],ChemistryDict, iz)
                Sout[iz] = tmpSout
                #ZZZ need to generalize


        # extrapolate down to lower altitudes like what is done
        izMin = numpy.where(IRIin > 0)[0]
        q0 = numpy.where(IRIin[izMin] == numpy.min(IRIin[izMin]))
        i0 = izMin[q0][0]
        ScaleHeight = 2. # km
        Sout[0:i0] = Sout[i0]*numpy.exp((altkm[0:i0]-altkm[i0])/ScaleHeight)
        print 'izMin', izMin
        print 'q0',q0
        print 'i0',i0
        # % Extend the source to low altitudes
        # ii=find(S0<=0);
        # i0=min(find(S0>0));
        # HS=2; % from figures in [R]
        # S0(ii)=S0(i0)*exp((z(ii)-z(i0))/HS);
        #

        # % The cosmic-ray source
        # hcr=getvaluefromdict(options,'hcr',15);
    	# % -- for mid-latitudes
        # Scrmax=getvaluefromdict(options,'Scrmax',10);
    	# % -- From NASA-CP-2090 report [R], in 1/cm^3/s;
        # tmp=getNm(z)/getNm(hcr);
        # Scr=tmp.*exp(-tmp+1)*Scrmax;
        # S0=S0+Scr;
        ztmp = self.MSISDict['Nm']/self.MSISatRef['Nm']
        Scr = 10.*ztmp*numpy.exp(-ztmp+1)
        Sout = Sout+Scr
        for iz in range(altkm.shape[0]):
            y0 = self.ODE(Sout[iz],iz, ChemistryDict)
            yInitial[iz,0:4] = y0
            yInitial[iz,-1] = (y0[0]+y0[1]+y0[3])-y0[2]
            print 'iz y0,', iz, y0

        NeIn = dict()
        NeIn['Ne'] = yInitial[:,0]
        NeIn['NnegIon'] = yInitial[:,1]
        NeIn['NposCluster'] = yInitial[:,2]
        NeIn['NX'] = yInitial[:,3]
        NeIn['NposIon'] = yInitial[:,4]


        return Sout, NeIn


    def Calculate_Ionization_From_Ne(self, altkm,NeIn, ChemistryDict):
        """
        This subroutine is designed to determine the ionization profile
        given an input/observed electron density.
        This is part 1 in a two step inversion process
        """
        iriAltGrid = self.altkm
        # to do an irregularly spaced grid, I need to find the correct index to feed into the binary search
        # I can find that by looking at which grid point is nearest and subbing that into iz

        # The altitudes MUST MUST MUST be the same as the IRIGird as the Chemistry

        Sout = numpy.zeros(altkm.shape[0])
        NeOut = numpy.zeros(altkm.shape[0])
        yInitial = numpy.zeros([altkm.shape[0],5])


        for iz in range(altkm.shape[0]):

            if NeIn[iz] < 0:

                continue
            else:
                AltDiff = numpy.abs(iriAltGrid-altkm[iz])
                indx = numpy.where(AltDiff == numpy.min(AltDiff))[0][0]
                tmpSout, tmpNeOut, y0 = self.Binary_Search(NeIn[iz],ChemistryDict, indx)
                Sout[iz] = tmpSout
                NeOut[iz] = y0[0]
                print 'Ne2QZ iz,indx, IRIiz', iz,indx,altkm[iz],iriAltGrid[indx], NeIn[iz],y0[0], tmpSout
        return Sout, NeOut



    def Set_Inital_Ionization(self,tUnix,glat,glon,AltitudeMin,AltitudeMax,deltaAltitude):
        """
        This function will be run outside of the main run routine only because
        I want to give control to the user about how often to update the initial
        ionization profile
        """

       # now run IRI to get the profile in
        iriDict = iri2016.IRI2016(tUnix,glat,glon,AltitudeMin,AltitudeMax,deltaAltitude)
        self.NeIn = iriDict['Ne']/1e6 # needs to be in cm^-3
        self.altkm = iriDict['Altitude']


        # now run MSIS to calculate Dregion Chemistry stuff on Altitude grid
        # set the time variables
        year = int(tUnix/(24.*3600.*365)+1970.)
        doy = -1*int((tUnix/(24.*3600))%365)
        utHrs = (tUnix/3600.)%24

        self.MSISDict = msis.MSIS(doy,utHrs,glat,glon,year,altkm=self.altkm, CGSorSI = 'CGS')
        self.MSISatRef = msis.MSIS(doy,utHrs,glat,glon,year,altkm=numpy.array([15.]), CGSorSI = 'CGS')
        # [ZZZ]needs to be user specified
        options = dict()
        options['GammaType'] = 'Temp'
        options['daytime'] = False
        # self.DregionChem = self.Calculate_Dregion_ReactionRates(self.MSISDict['Nm'], \
        #                                                     Tn = self.MSISDict['Tn'], \
        #                                                     options=options)
        self.DregionChem = self.Calculate_Dregion_ReactionRates(self.MSISDict['Nm'], \
                                                            Tn = self.MSISDict['Tn'], \
                                                             options=options, \
                                                             iriDict = iriDict)

        self.Sin, self.y0 =self.Calculate_Background_Ionization(self.altkm,self.NeIn,self.DregionChem)
        return

    def __call__(self,qz,altkm,TypeName='Dregion'):
        """
        Run routine requires the altitude grid, and whether this is
        E or D-region ionization, which must be specified

        This program takes the height ionization profile and produces the electron density
        Assumptions:
            1. For each interval of time, the chemistry coefficients are fixed and
            only the ionization is changing.

        """
        iriAltGrid = self.altkm
        qin = numpy.zeros(qz.shape[0])
        print 'qz.shape, altkm.shape', qz.shape, altkm.shape
        print 'self.y0, self.altkm', self.y0['Ne'].shape, self.altkm.shape
        y0 = dict()
        for ikeys in self.y0.keys():
            y0[ikeys] = numpy.zeros(qz.shape[0])

        DregionChemDict = dict()
        for ikeys in self.DregionChem.keys():
            DregionChemDict[ikeys] = numpy.zeros(qz.shape[0])
        # y0['NnegIon'] = numpy.zeros(qz.shape[0])
        # y0['NposCluster'] = numpy.zeros(qz.shape[0])
        # y0['NX'] = numpy.zeros(qz.shape[0])
        # y0['NposIon'] = numpy.zeros(qz.shape[0])
        if TypeName == 'Dregion':

            # check if the initial ionization was calculated
            # again put everything onto the same altitude grid and run the model.
            if self.Sin is not None:
                if qz.shape[0] == altkm.shape[0]:

                    for iz in range(qz.shape[0]):
                        AltDiff = numpy.abs(iriAltGrid-altkm[iz])
                        indx = numpy.where(AltDiff == numpy.min(AltDiff))[0][0]
                        qin[iz] = qz[iz]+self.Sin[indx]

                        for ikeys in self.y0.keys():
                            y0[ikeys][iz] = self.y0[ikeys][indx]

                        for ikeys in self.DregionChem.keys():
                            DregionChemDict[ikeys][iz] = self.DregionChem[ikeys][indx]
                        # y0['NnegIon'][iz] = self.y0['NnegIon'][indx]
                        # y0['NposCluster'][iz] = self.y0['NposCluster'][indx]
                        # y0['NX'][iz] = self.y0['NX'][indx]
                        # y0['NposIon'][iz] = self.y0['NposIon'][indx]

                        print altkm[iz],iriAltGrid[indx],qz[iz],self.Sin[indx]

                    print 'qin.shape,', qin.shape, y0['Ne'].shape, self.DregionChem['B'].shape

                    results = self.Integrate_ODE(y0,qin,DregionChemDict,IntType='5species')
                    print results['Ne'].shape
                    print '###############################'
                    #results = self.ODE(self.y0,qin,self.DregionChem,IntType='5species')
                else:
                    raise ValueError("Ionization and altitude sizes do not agree")

            else:
                raise ValueError ('Initial Ionization was not set. Run: Set_Inital_Ionization')

        if TypeName == 'Special':
            if self.Sin is not None:
                if qz.shape[0] == altkm.shape[0]:
                    # results = self.Integrate_ODE(self.y0,qin,self.DregionChem,IntType='5species')
                    results = self.Integrate_ODE(self.y0,qz,self.DregionChem,IntType='5species')
                else:
                    raise ValueError("Ionization and altitude sizes do not agree")

            else:
                raise ValueError ('Initial Ionization was not set. Run: Set_Inital_Ionization')




        return results

if __name__ == "__main__":
    # do validation of previous results as building up the class.
    from scipy.io import loadmat
    dataIn = loadmat('./Matlab/Everything.mat')
    # need to make sure Nm is cm^-3
    chem = Chemistry()

    # run the D-region Chemistry and check against
    options = dict()
    options['GammaType'] = 'Temp'
    options['daytime'] = False
    Nm = dataIn['Nm']
    Tn = dataIn['Tn']
    DregionChem = chem.Calculate_Dregion_ReactionRates(Nm, \
                                                        Tn = Tn,\
                                                        options=options)

    # check values:
    print 'alphaD', numpy.nanmax(DregionChem['alphaD'] - dataIn['alfad'])
    print 'alphaDC', numpy.nanmax(DregionChem['alphaDC'] - dataIn['alfadc'])
    print 'alphaI', numpy.nanmax(DregionChem['alphaI'] - dataIn['alfai'])
    print 'beta', numpy.nanmax(DregionChem['beta']-dataIn['beta'])
    print 'B', numpy.nanmax(DregionChem['B']-dataIn['Bcoef'])
    print 'gamma', numpy.nanmax(DregionChem['gamma']-dataIn['gamma'])
    print 'gammaX', numpy.nanmax(DregionChem['gammaX']-dataIn['gammaX'])
    print 'Xbar', numpy.nanmax(DregionChem['Xbar'] - dataIn['Xbar'])


    # now want to validate previous results, using S0 from Everything.mat
    # this is the same as testodeintegrate in the other
    S0 = dataIn['S'][:,0]#dataIn['S0']
    # specnames={'Ne','Nneg','Nclus','NX','Npos'};
    NeDictIn = {}

    NeDictIn['Ne'] = dataIn['Nspec0'][:,0]
    NeDictIn['NnegIon'] = dataIn['Nspec0'][:,1]
    NeDictIn['NposCluster'] = dataIn['Nspec0'][:,2]
    NeDictIn['NX'] = dataIn['Nspec0'][:,3]
    NeDictIn['NposIon'] = dataIn['Nspec0'][:,4]

    NeOut = chem.Integrate_ODE(NeDictIn,S0,DregionChem)


    # plot the results and check them
    plt.figure(101)
    plt.semilogx(dataIn['Nspec'][:,-1,0], dataIn['z'], 'b-', label='Ne')
    plt.semilogx(NeOut['Ne'], dataIn['z'], 'r-')
    #
    # plt.figure(2)
    # plt.semilogx(dataIn['Nspec'][:,-1,1], dataIn['z'], 'b-', label='NnegIon')
    # plt.semilogx(NeOut['NnegIon'], dataIn['z'], 'r-')
    #
    # plt.figure(3)
    # plt.semilogx(dataIn['Nspec'][:,-1,2], dataIn['z'], 'b-', label='NposCluster')
    # plt.semilogx(NeOut['NposCluster'], dataIn['z'], 'r-')
    #
    # plt.figure(4)
    # plt.semilogx(dataIn['Nspec'][:,-1,3], dataIn['z'], 'b-', label='NX')
    # plt.semilogx(NeOut['NX'], dataIn['z'], 'r-')
    #
    # plt.figure(5)
    # plt.semilogx(dataIn['Nspec'][:,-1,4], dataIn['z'], 'b-', label='NposIon')
    # plt.semilogx(NeOut['NposIon'], dataIn['z'], 'r-')
    # plt.show()

    """
    Test the finding the background distribution
    This validates the binary search routine
    """
    """
    NeIn = dataIn['Nspec0'][:,0]
    altkm = dataIn['z']

    # below 80 km set to -1, because that is where IRI is defined down to at night
    NeIn[0:50] = -1

    TestS = chem.Calculate_Background_Ionization(altkm,NeIn,DregionChem)
    plt.figure(11)
    plt.plot(dataIn['S0'], altkm, 'r-')
    plt.plot(TestS, altkm, 'b-')
    #plt.show()
    """
    """
    Testing initialize ionosphere and msis
    This will allow someone to update the ionosphere and MSIS as needed
    """

    import datetime
    AltMinKm = 50.
    AltMaxKm = 150.
    AltStepKm = 1.0

    t1970 = datetime.datetime(1970,1,1,0,0,0)
    t2010 = datetime.datetime(2010,6,2,01,0,0)
    tUnix = (t2010-t1970).total_seconds()

    glat = 45
    glon = 228

    chem.Set_Inital_Ionization(tUnix,glat,glon,AltMinKm,AltMaxKm,AltStepKm)

    print chem.Sin

    plt.figure(1)
    plt.semilogx(chem.NeIn, chem.altkm)

    plt.figure(2)
    plt.semilogx(chem.Sin, chem.altkm)
    print 'Sin', chem.Sin


    plt.figure(3)
    for ii in chem.y0.keys():
        plt.semilogx(chem.y0[ii], chem.altkm, label='%s'%ii)
    plt.legend()

    # Npos0=Ne0+Nneg0-Nclus0+NX0;
    # plt.figure(31)
    qPos = chem.y0['NposCluster']+chem.y0['NposIon']
    qNeg = chem.y0['Ne']+chem.y0['NnegIon']+chem.y0['NX']
    # plt.semilogx(qPos, chem.altkm, 'b-')
    plt.semilogx(qNeg, chem.altkm, 'k-')
    """
    Test the run function
    """

   # instantiate a glow instance
    sys.path.append('../Ionization')
    import Ionization
    iz = Ionization.Ionization()
    Q0 = 1.0
    E0 = 10.0*1e3
    EeV = numpy.logspace(2,6,num=201)
    NumFlux,QeV = iz.MaxwellianFlux(EeV,Q0,E0)
    qZ, qZE,qZEsimps = iz.Ionization(EeV,QeV,AltMinKm,AltMaxKm,AltStepKm,tUnix,glat,glon)

    plt.figure(4)
    plt.semilogx(qZ, chem.altkm)

    d = chem(qZ,chem.altkm)

    plt.figure(5)
    for ii in d.keys():
        plt.semilogx(d[ii], chem.altkm, label='%s'%ii)
    qNeg = d['Ne']+d['NnegIon']+d['NX']
    # plt.semilogx(qPos, chem.altkm, 'b-')
    plt.semilogx(qNeg, chem.altkm, 'k-')
    plt.legend()
    # plt.figure(2)
    # plt.semilogx(chem.MSISDict['Nm'], chem.altkm)
    # plt.semilogx(Nm,dataIn['z'])
    #
    # plt.figure(3)
    # plt.semilogx(chem.MSISDict['Tn'], chem.altkm)
    # plt.semilogx(Tn, dataIn['z'])
    #
    # print 'gammaX', DregionChem['gammaX'],chem.DregionChem['gammaX']
    # kk = 4
    # for ikey in chem.DregionChem.keys():
    #     plt.figure(kk)
    #     plt.semilogx(chem.DregionChem[ikey], chem.altkm, 'rx')
    #     plt.semilogx(DregionChem[ikey],dataIn['z'])
    #     plt.title(ikey)
    #     print ikey
    #     print chem.DregionChem[ikey].shape, DregionChem[ikey].shape
    #     kk+=1





    plt.show()
