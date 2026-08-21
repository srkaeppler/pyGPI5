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

class Chemistry:

    def __init__(self):
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
        if len(CO2) > 0:
            assert (CO2.shape[0] == Te.shape[0], "CO Shape is not the same as Te")
            alpha = CO2*a_O2 + alpha
        if len(CNO) > 0:
            assert (CNO.shape[0] == Te.shape[0], "CNO Shape is not the same as Te")
            alpha = CNO*a_NO + alpha
        if len(CO) > 0:
            assert (CO.shape[0] == Te.shape[0], "CO Shape is not the same as Te")
            alpha = CO*a_O + alpha

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
        alphaI = 1e-7*numpy.ones(Nm.shape[0])
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

    def Calculate_Gamma(self,Nm,Tn = 0, option='GPI'):

        if option == 'GPI':
            gamma = 3e-17*Nm

        return gamma

    def Dregion_Chemistry(self,Nm, Tn=None):

        # should probably assert that Nm is a 1-D array
        outDict = dict()
        outDict['alphaD'] = self.Calculate_alphaD(Nm)
        outDict['alphaDC'] = self.Calculate_alphaDC(Nm)
        outDict['alphaI'] = self.Calculate_alphaI(Nm)
        outDict['beta'] = self.Calculate_Beta(Nm)
        outDict['B'] = self.Calculate_B(Nm)

        if Tn:
            outDict['gamma'] = self.Calculate_Gamma(Nm,Tn=Tn,option='foo')
        else:
            outDict['gamma'] = self.Calculate_Gamma(Nm)

        return outDict

if __name__ == "__main__":
    # do validation of previous results as building up the class.
    from scipy.io import loadmat
    dataIn = loadmat('./Matlab/Everything.mat')
    # need to make sure Nm is cm^-3
    chem = Chemistry()

    # run the D-region Chemistry and check against
    DregionChem = chem.Dregion_Chemistry(dataIn['Nm'])

    # check values:
    print 'alphaD', numpy.nanmax(DregionChem['alphaD'] - dataIn['alfad'])
    print 'alphaDC', numpy.nanmax(DregionChem['alphaDC'] - dataIn['alfadc'])
    print 'alphaI', numpy.nanmax(DregionChem['alphaI'] - dataIn['alfai'])
    print 'beta', numpy.nanmax(DregionChem['beta']-dataIn['beta'])
    print 'B', numpy.nanmax(DregionChem['B']-dataIn['Bcoef'])
    print 'gamma', numpy.nanmax(DregionChem['gamma']-dataIn['gamma'])
