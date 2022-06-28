#! /usr/bin/env python

"""
xxxxx

Originally written  by:
~M. Nicolls
last revised: xx/xx/2007
with same functionality

Updated into the structure I wanted by SRK


MSIS Documentation:
# Switches: to turn on and off particular variations use these switches.
# 0 is off, 1 is on, and 2 is main effects off but cross terms on.
# Standard values are 0 for switch 0 and 1 for switches 1 to 23. The
# array "switches" needs to be set accordingly by the calling program.
# The arrays sw and swc are set internally.
#   switches[i]:
#    i - explanation
#   -----------------
#   0 - output in centimeters instead of meters
#   1 - F10.7 effect on mean
#   2 - time independent
#   3 - symmetrical annual
#   4 - symmetrical semiannual
#   5 - asymmetrical annual
#   6 - asymmetrical semiannual
#   7 - diurnal
#   8 - semidiurnal
#   9 - daily ap [when this is set to -1 (!) the pointer
#               ap_a in struct nrlmsise_input must
#               point to a struct ap_array]
#   10 - all UT/long effects
#   11 - longitudinal
#   12 - UT and mixed UT/long
#   13 - mixed AP/UT/LONG
#   14 - terdiurnal
#   15 - departures from diffusive equilibrium
#   16 - all TINF var
#   17 - all TLB var
#   18 - all TN1 var
#   19 - all S var
#   20 - all TN2 var
#   21 - all NLB var
#   22 - all TN3 var
#   23 - turbo scale height var


# Array containing the following magnetic values:
    # 0 : daily AP
    # 1 : 3 hr AP index for current time
    # 2 : 3 hr AP index for 3 hrs before current time
    # 3 : 3 hr AP index for 6 hrs before current time
    # 4 : 3 hr AP index for 9 hrs before current time
    # 5 : Average of eight 3 hr AP indicies from 12 to 33 hrs prior to current time
    # 6 : Average of eight 3 hr AP indicies from 36 to 57 hrs prior to current time

# OUTPUT VARIABLES:
#   d[0] - HE NUMBER DENSITY(CM-3)
#   d[1] - O NUMBER DENSITY(CM-3)
#   d[2] - N2 NUMBER DENSITY(CM-3)
#   d[3] - O2 NUMBER DENSITY(CM-3)
#   d[4] - AR NUMBER DENSITY(CM-3)
#   d[5] - TOTAL MASS DENSITY(GM/CM3) [includes d[8] in td7d]
#   d[6] - H NUMBER DENSITY(CM-3)
#   d[7] - N NUMBER DENSITY(CM-3)
#   d[8] - Anomalous oxygen NUMBER DENSITY(CM-3)
#   t[0] - EXOSPHERIC TEMPERATURE
#   t[1] - TEMPERATURE AT ALT*
#
#   O, H, and N are set to zero below 72.5 km
#
#   t[0], Exospheric temperature, is set to global average for
#   altitudes below 120 km. The 120 km gradient is left at global
#   average value for altitudes below 72 km.
#
#   d[5], TOTAL MASS DENSITY, is NOT the same for subroutines GTD7 and GTD7D
#
#   SUBROUTINE GTD7 -- d[5] is the sum of the mass densities of the
#   species labeled by indices 0-4 and 6-7 in output variable d.
#   This includes He, O, N2, O2, Ar, H, and N but does NOT include
#   anomalous oxygen (species index 8).
#
#   SUBROUTINE GTD7D -- d[5] is the "effective total mass density
#   for drag" and is the sum of the mass densities of all species
#   in this model, INCLUDING anomalous oxygen.

"""

import os, ctypes
# import numpy, numpy.fftpack, numpy.interpolate, numpy.optimize
import numpy
import datetime


# this is just how this needs to be defined
class MSIS_OUTPUT(ctypes.Structure):
        _fields_ = [("d", ctypes.ARRAY(ctypes.c_double,9)),
                    ("t", ctypes.ARRAY(ctypes.c_double,2))]

class MSIS_APARRAY(ctypes.Structure):
        _fields_ = [("a",ctypes.ARRAY(ctypes.c_double,7))]

class MSIS_FLAGS(ctypes.Structure):
        _fields_ = [("switches",ctypes.ARRAY(ctypes.c_int,24)),
                    ("sw",ctypes.ARRAY(ctypes.c_double,24)),
                    ("swc",ctypes.ARRAY(ctypes.c_double,24))]

class MSIS_INPUT(ctypes.Structure):
            _fields_ = [("year", ctypes.c_int), # year, doesnt matter
                    ("doy", ctypes.c_int), # day of year
                    ("sec", ctypes.c_double), # seconds in day (UT)
                    ("alt", ctypes.c_double), # altitude (km)
                    ("g_lat", ctypes.c_double), # geodetic latitude
                    ("g_long", ctypes.c_double), # geodetic longitude
                    ("lst", ctypes.c_double), # local apparent solar time (hours)
                    ("f107A", ctypes.c_double), # 81 day average of F10.7 flux (centered on doy)
                    ("f107", ctypes.c_double), # daily F10.7 flux for previous day
                    ("ap", ctypes.c_double), # magnetic index(daily)
                    ("ap_array", ctypes.POINTER(MSIS_APARRAY)) ]

class MSIS:


    def __init__(self):
        self.geophys_dir = '/Users/srkaeppler/research/data/AFOSR_Eregion_Conductivity/Models/AP_KP'

        self.inLib = '/Users/srkaeppler/research/data/AFOSR_Eregion_Conductivity/Models/nrlmsise00/libnrlmsise-00.so'
        if os.path.isfile(self.inLib):
            self.ctype_msis = ctypes.cdll.LoadLibrary(self.inLib)
        else:
            raise ValueError('NRL msis shared object file does not exist or bad location \n')

        return


    def read_geophys(self,year,doy,curtime):

        # Returns:	F107 - for previous day
        #			F107a - 81 day average
        #			AP - Array containing the following magnetic values:
        #				*   0 : daily AP
        #				*   1 : 3 hr AP index for current time
        #				*   2 : 3 hr AP index for 3 hrs before current time
        #				*   3 : 3 hr AP index for 6 hrs before current time
        #				*   4 : 3 hr AP index for 9 hrs before current time
        #				*   5 : Average of eight 3 hr AP indicies from 12 to 33 hrs
        #				*           prior to current time
        #				*   6 : Average of eight 3 hr AP indicies from 36 to 57 hrs
        #				*           prior to current time
        #

        if os.path.exists(os.path.join(self.geophys_dir,str(year)))==False:
            year=year-1
        # if os.path.exists(os.path.join(self.geophys_dir,str(year)))==False:
        #     raise IOError, print('Geophys param directory %s/%s does not exist.' % (geophys_dir, str(year)))

        year=str(year)

        f = open(os.path.join(self.geophys_dir,year))
        lines = f.readlines()
        if os.path.exists(os.path.join(self.geophys_dir,str(int(year)-1)))==True:
            f = open(os.path.join(self.geophys_dir,str(int(year)-1)))
            lines_py = f.readlines()
        else:
            lines_py=[]
        if os.path.exists(os.path.join(self.geophys_dir,str(int(year)+1)))==True:
            f = open(os.path.join(self.geophys_dir,str(int(year)+1)))
            lines_ny = f.readlines()
        else:
            lines_ny=[]

        if len(lines)<doy:
            doy=len(lines)

        print(len(lines))
        # F!07d - previous day
        if doy==1:
            try:
                F107D = float(lines_py[-1][65:70])
            except: # if can't get prev day, just get current day
                F107D = float(lines[doy-1][65:70])
        else:
            F107D = float(lines[doy-2][65:70])

        # AP
        AP = numpy.zeros((7))
        AP[0] = float(lines[doy-1][55:58]) # daily, for current day
        TINDEX = int(curtime/3.0)
        tlines=lines; tind=doy-1
        for i in range(4):
            AP[i+1] = float(tlines[tind][31+TINDEX*3:34+TINDEX*3]) # for current time
            if TINDEX>0:
                TINDEX=TINDEX-1
            else:
                if tind==0:
                    tlines=lines_py
                    tind=len(tlines)-1
                else:
                    tind -= 1
                TINDEX=7
        for i in range(8):
            AP[5] += float(tlines[tind][31+TINDEX*3:34+TINDEX*3]) # for current time
            if TINDEX>0:
                TINDEX=TINDEX-1
            else:
                if tind==0:
                    tlines=lines_py
                    tind=len(tlines)-1
                else:
                    tind=tind-1
                TINDEX=7
        AP[5] /= 8.0
        for i in range(8):
            AP[6] += float(tlines[tind][31+TINDEX*3:34+TINDEX*3]) # for current time
            if TINDEX>0:
                TINDEX=TINDEX-1
            else:
                if tind==0:
                    tlines=lines_py
                    tind=len(tlines)-1
                else:
                    tind=tind-1
                TINDEX=7
        AP[6] /= 8.0

        F107A = 0
        imin=doy-1-40
        imax=doy+40

        if imax>len(lines) and len(lines_ny)==0:
            imin=imin-(imax-len(lines))+len(lines_ny)
            imax=imax-(imax-len(lines))+len(lines_ny)

        lines2=[]
        for aa in range(imin,imax):
            try:
                if aa<0:
                	lines2.append(lines_py[len(lines_py)+aa])
                elif aa>=len(lines):
                	lines2.append(lines_ny[aa-len(lines)])
                else:
                	lines2.append(lines[aa])
            except:
                ''

        F107A=0.0
        for aa in range(len(lines2)):
            try:
                F107A=F107A+float(lines2[aa][65:70])
            except:
                try:
                    F107A=F107A+float(lines2[aa+1][65:70])
                except:
                    F107A=F107A+float(lines2[aa-1][65:70])
        F107A=F107A/len(lines2)

        return F107D, F107A, AP

    def run_MSIS(self,yr,doy,hrUT,altkm,\
                glat,glong,ap_array,f107a,\
                f107,ap=-1, CGSorSI = 'SI'):

        # input structure
        InStruct=MSIS_INPUT()
        InStruct.year=yr
        InStruct.doy=doy
        InStruct.sec=hrUT*3600 # convert from decimal hours into seconds
        InStruct.alt=altkm
        InStruct.g_lat=glat
        InStruct.g_long=glong
        InStruct.lst=hrUT+glong/15.0
        InStruct.f107A=f107a
        InStruct.f107=f107
        InStruct.ap=ap
        InStruct.ap_array=ctypes.POINTER(MSIS_APARRAY)() # null pointer

        # ap_array
        tmp=MSIS_APARRAY()
        for i in range(7):
            tmp.a[i]=ap_array[i]
        InStruct.ap_array.contents = tmp

        # flags structure
        flags=MSIS_FLAGS()
        sw= (ctypes.c_int * 24)()
        sw[0:]=[0,1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] # this is default
        if CGSorSI == 'CGS':
            sw[0] = 0
        if CGSorSI == 'SI':
            sw[0] = 1
        flags.switches=sw

        # output structure
        output=MSIS_OUTPUT()

        # call MSIS
        self.ctype_msis.gtd7(ctypes.byref(InStruct),ctypes.byref(flags),ctypes.byref(output))



        return output.d,output.t

    def MSIS(self, doy,hrUT,glat,glong,year,altkm=numpy.arange(50.0,1000.0,30.0), CGSorSI = 'SI'):

        (f107, f107a, ap)=self.read_geophys(int(year),int(doy),hrUT)

        MSISout={}
        MSISout['Altitude']=altkm
        MSISout['Tn']=numpy.zeros(altkm.shape)
        MSISout['MassDensity']=numpy.zeros(altkm.shape)
        MSISout['AverageMass'] = numpy.zeros(altkm.shape)
        MSISout['Nm'] = numpy.zeros(altkm.shape)
        MSISout['F107']=f107
        MSISout['F107a']=f107a
        MSISout['AP']=ap
        #d[0] - HE NUMBER DENSITY(CM-3)
        #   d[1] - O NUMBER DENSITY(CM-3)
        #   d[2] - N2 NUMBER DENSITY(CM-3)
        #   d[3] - O2 NUMBER DENSITY(CM-3)
        # d[3] - O2 NUMBER DENSITY(CM-3)
        # #   d[4] - AR NUMBER DENSITY(CM-3)
        # #   d[5] - TOTAL MASS DENSITY(GM/CM3) [includes d[8] in td7d]
        # #   d[6] - H NUMBER DENSITY(CM-3)
        # #   d[7] - N NUMBER DENSITY(CM-3)
        MSISout['nO'] = numpy.zeros(altkm.shape)
        MSISout['nN2'] = numpy.zeros(altkm.shape)
        MSISout['nO2'] = numpy.zeros(altkm.shape)
        MSISout['nN'] = numpy.zeros(altkm.shape)


        for Iht in range(altkm.size):
            (d,t)=self.run_MSIS(int(year),int(doy),hrUT,altkm[Iht],glat,glong,ap,\
                                f107a,f107,-1, CGSorSI=CGSorSI)
            MSISout['MassDensity'][Iht] = d[5]
            MSISout['AverageMass'][Iht] = d[5]/(d[0]+d[1]+d[2]+d[4]+d[6]+d[7]+d[8])
            MSISout['Texo']=t[0]
            MSISout['Tn'][Iht]=t[1]
            MSISout['Nm'][Iht] = (d[0]+d[1]+d[2]+d[4]+d[6]+d[7]+d[8])
            MSISout['nO'][Iht] = d[1]
            MSISout['nN2'][Iht] = d[2]
            MSISout['nO2'][Iht] = d[3]
            MSISout['nN'][Iht] = d[7]



        return MSISout

    def MSIS2(self, tUnix,glat,glong,Heibeg,Heiend,step, CGSorSI = 'SI'):
        # standard interface - hate it or love it...
        #tUnix,Lat,Lon,Heibeg,Heiend,step

        dt1 = datetime.datetime.utcfromtimestamp(tUnix)
        altkm = numpy.arange(Heibeg,Heiend,step)
        year = dt1.year
        doy = dt1.timetuple().tm_yday
        hrUT = float(dt1.hour)+float(dt1.minute/60.)+float(dt1.second/3600.)

        (f107, f107a, ap)=self.read_geophys(int(year),int(doy),hrUT)

        MSISout={}
        MSISout['Altitude']=altkm
        MSISout['Tn']=numpy.zeros(altkm.shape)
        MSISout['MassDensity']=numpy.zeros(altkm.shape)
        MSISout['AverageMass'] = numpy.zeros(altkm.shape)
        MSISout['Nm'] = numpy.zeros(altkm.shape)
        MSISout['F107']=f107
        MSISout['F107a']=f107a
        MSISout['AP']=ap
        #d[0] - HE NUMBER DENSITY(CM-3)
        #   d[1] - O NUMBER DENSITY(CM-3)
        #   d[2] - N2 NUMBER DENSITY(CM-3)
        #   d[3] - O2 NUMBER DENSITY(CM-3)
        MSISout['nO'] = numpy.zeros(altkm.shape)
        MSISout['nN2'] = numpy.zeros(altkm.shape)
        MSISout['nO2'] = numpy.zeros(altkm.shape)
        MSISout['nN'] = numpy.zeros(altkm.shape)

        for Iht in range(altkm.size):
            (d,t)=self.run_MSIS(int(year),int(doy),hrUT,altkm[Iht],glat,glong,ap,\
                                f107a,f107,-1, CGSorSI=CGSorSI)
            MSISout['MassDensity'][Iht] = d[5]
            MSISout['AverageMass'][Iht] = d[5]/(d[0]+d[1]+d[2]+d[4]+d[6]+d[7]+d[8])
            MSISout['Texo']=t[0]
            MSISout['Tn'][Iht]=t[1]
            MSISout['Nm'][Iht] = (d[0]+d[1]+d[2]+d[4]+d[6]+d[7]+d[8])
            MSISout['nO'][Iht] = d[1]
            MSISout['nN2'][Iht] = d[2]
            MSISout['nO2'][Iht] = d[3]
            MSISout['nN'][Iht] = d[7]



        return MSISout

if __name__ == '__main__':
    # test case
    msis = MSIS()
    year=2010
    doy=172
    hrUT=29000/3600.
    g_lat=60
    g_long=-70
    CGSorSI = 'SI'
    outDict = msis.MSIS(doy,hrUT,g_lat,g_long,year, CGSorSI=CGSorSI)
    print(outDict)

    dt1970 = datetime.datetime(1970,1,1,0,0,0)
    dt1 = datetime.datetime(2010,6,1,13,41,56)
    tunix1 = (dt1-dt1970).total_seconds()
    outDict2 = msis.MSIS2(tunix1,g_lat,g_long,50.,100.,1.,CGSorSI='SI')
    print(outDict2)
