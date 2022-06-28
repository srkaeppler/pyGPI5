#! /usr/bin/env python

"""
xxxxx

~M. Nicolls
last revised: xx/xx/2007

"""

import os, ctypes 
import scipy, scipy.fftpack, scipy.interpolate, scipy.optimize

from constants import *

class MSIS_FLAGS(ctypes.Structure):
    _fields_ = [("switches",ctypes.ARRAY(ctypes.c_int,24)),
                ("sw",ctypes.ARRAY(ctypes.c_double,24)),
                ("swc",ctypes.ARRAY(ctypes.c_double,24))]
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

class MSIS_APARRAY(ctypes.Structure):
    _fields_ = [("a",ctypes.ARRAY(ctypes.c_double,7))]
# Array containing the following magnetic values:
    # 0 : daily AP
    # 1 : 3 hr AP index for current time
    # 2 : 3 hr AP index for 3 hrs before current time
    # 3 : 3 hr AP index for 6 hrs before current time
    # 4 : 3 hr AP index for 9 hrs before current time
    # 5 : Average of eight 3 hr AP indicies from 12 to 33 hrs prior to current time
    # 6 : Average of eight 3 hr AP indicies from 36 to 57 hrs prior to current time 

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
                ("ap_array", ctypes.POINTER(MSIS_APARRAY)) ] # see above

class MSIS_OUTPUT(ctypes.Structure):
    _fields_ = [("d", ctypes.ARRAY(ctypes.c_double,9)),
                ("t", ctypes.ARRAY(ctypes.c_double,2))] 
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

def compute_mfrac(z,Tinf,T120,z50=185,a=0.01,Hinf=20):
    
    H=Hinf/Tinf*(Tinf-(Tinf-T120)*scipy.exp(-a*(z-120.0)))
    
    qop=2.0/(1+scipy.sqrt(1+8.0*scipy.exp(-(z-z50)/H)))

    return qop

def compute_collfreq(d,Tn, Te=1000.0, mj=30.0):
    
#   nu_in = 0.0
#   nu_in = nu_in + d[1]*scipy.sqrt(0.79/16) # O
#   nu_in = nu_in + d[2]*scipy.sqrt(1.76/28) # N2
#   nu_in = nu_in + d[3]*scipy.sqrt(1.59/32) # O2
#   nu_in = nu_in*2.6e-9

    Ti=Tn
    Tr = (Tn+Ti)/2.0
    nu_in = 0.0
    if mj==16.0: # O
        nu_in = nu_in + d[1]*3.67e-11*scipy.sqrt(Tr)*(1.0-0.064*scipy.log10(Tr))**2.0
    else:
        nu_in = nu_in + 1.0e6*d[1]*9.14e-15/scipy.sqrt(mj*(mj+16.0)) 
    if mj==28.0: # N2
        nu_in = nu_in + d[2]*5.14e-11*scipy.sqrt(Tr)*(1.0-0.069*scipy.log10(Tr))**2.0
    else:
        nu_in = nu_in + 1.0e6*d[2]*1.80e-14/scipy.sqrt(mj*(mj+28.0)) 
    if mj==32.0 and Tn>800.0: # O2
        nu_in = nu_in + d[2]*2.59e-11*scipy.sqrt(Tr)*(1.0-0.073*scipy.log10(Tr))**2.0
    else:
        nu_in = nu_in + 1.0e6*d[3]*1.83e-14/scipy.sqrt(mj*(mj+32.0)) 
    
    nu_en = 0.0
    nu_en = nu_en + d[1]*8.2e-10*scipy.sqrt(Te) # O
    nu_en = nu_en + d[2]*2.33e-11*(1-1.2e-4*Te)*Te # N2
    nu_en = nu_en + d[3]*1.8e-10*(1+3.6e-2*scipy.sqrt(Te))*scipy.sqrt(Te)
    
    return nu_in, nu_en

def call_MSIS(ct_msis,doy,hrUT,glat,glong,year,altkm,ap,f107a,f107,z50,mass=[30]):
    
        
    (d,t,nui,nue)=run_MSIS(ct_msis,int(year),int(doy),hrUT,altkm,glat,glong,ap,f107a,f107,-1,mass)
    (x,t120,xx,xxx)=run_MSIS(ct_msis,int(year),int(doy),hrUT,120.0,glat,glong,ap,f107a,f107,-1,mass)
        
    HEdens=d[0]; Odens=d[1]; N2dens=d[2]; O2dens=d[3]
    ARdens=d[4]; MassDens=d[5]; Hdens=d[6]; Ndens=d[7]; AnomOdens=d[8]
    Texo=t[0]; Tn=t[1]
    
    qOp=compute_mfrac(altkm,Texo,t120[1],z50=z50)
        
    return HEdens,Odens,N2dens,O2dens,ARdens,MassDens,Hdens,Ndens,AnomOdens,Texo,Tn,nui,nue,qOp

def iterate_MSIS(ct_msis,doy,hrUT,glat,glong,year,geophys_path,altkm=scipy.arange(50.0,1000.0,30.0),z50=185.0,mass=[30]):
    
    (f107, f107a, ap)=read_geophys(int(year),int(doy),hrUT,geophys_path)    
    
    MSISout={}
    MSISout['ht']=altkm
    MSISout['nu_en']=scipy.zeros((len(mass),altkm.shape[0])) 
    MSISout['nu_in']=scipy.zeros((len(mass),altkm.shape[0]))   
    MSISout['Tn']=scipy.zeros(altkm.shape)
    MSISout['F107']=f107
    MSISout['F107a']=f107a
    MSISout['AP']=ap
    
    for Iht in range(altkm.size):
        (d,t,nui,nue)=run_MSIS(ct_msis,int(year),int(doy),hrUT,altkm[Iht],glat,glong,ap,f107a,f107,-1,mass)
        
        MSISout['Texo']=t[0]
        MSISout['Tn'][Iht]=t[1]
        MSISout['nu_en'][:,Iht]=nue
        MSISout['nu_in'][:,Iht]=nui
        
    # get O+ frac
    T120=scipy.interpolate.interp1d(MSISout['ht'],MSISout['Tn'],bounds_error=0)(120.0)
    MSISout['qOp']=compute_mfrac(MSISout['ht'],MSISout['Texo'],T120,z50=z50)

    return MSISout
    
def read_geophys(year,doy,curtime,geophys_dir):

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
    
    if os.path.exists(os.path.join(geophys_dir,str(year)))==False:
        year=year-1
    if os.path.exists(os.path.join(geophys_dir,str(year)))==False:
        raise IOError, 'Geophys param directory %s/%s does not exist.' % (geophys_dir, str(year))

    year=str(year)

    f = open(os.path.join(geophys_dir,year))
    lines = f.readlines()
    if os.path.exists(os.path.join(geophys_dir,str(int(year)-1)))==True:
        f = open(os.path.join(geophys_dir,str(int(year)-1)))
        lines_py = f.readlines()
    else:
        lines_py=[]
    if os.path.exists(os.path.join(geophys_dir,str(int(year)+1)))==True:
        f = open(os.path.join(geophys_dir,str(int(year)+1)))
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
    AP = scipy.zeros((7))
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

def run_MSIS(ct_msis,yr,doy,hrUT,altkm,glat,glong,ap_array,f107a,f107,ap=-1,mass=[30]):
    
    # input structure
    input=MSIS_INPUT()
    input.year=yr
    input.doy=doy
    input.sec=hrUT*3600
    input.alt=altkm
    input.g_lat=glat
    input.g_long=glong
    input.lst=hrUT+glong/15.0
    input.f107A=f107a
    input.f107=f107
    input.ap=ap
    input.ap_array=ctypes.POINTER(MSIS_APARRAY)() # null pointer

    # ap_array
    tmp=MSIS_APARRAY()
    for i in range(7):
        tmp.a[i]=ap_array[i]
    input.ap_array.contents = tmp
        
    # flags structure
    flags=MSIS_FLAGS()
    sw= (ctypes.c_int * 24)()
    sw[0:]=[0,1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    flags.switches=sw
    
    # output structure
    output=MSIS_OUTPUT()
    
    # call MSIS
    ct_msis.gtd7(ctypes.byref(input),ctypes.byref(flags),ctypes.byref(output))
    
    # get collision frequencies
    nu_in=scipy.zeros(len(mass))
    nu_en=scipy.zeros(len(mass))
    tn=output.t[1]*1.0
    for ia in range(len(mass)):
        (nu_in[ia],nu_en[ia])=compute_collfreq(output.d,tn,mj=mass[ia])
                
    return output.d,output.t,nu_in,nu_en
    
    
def test_MSIS():
    pathath='/Users/mnicolls/Documents/Work/ISfit/AMISR_fitter_py/lib/nrlmsise00/nrlmsis00_c_version/libnrlmsise-00.dylib'
    ct_msis=ctypes.CDLL(pathath) # spectra library

    # input structure
    input=MSIS_INPUT()
    input.year=0
    input.doy=172
    input.sec=29000
    input.alt=100
    input.g_lat=60
    input.g_long=-70
    input.lst=16.0
    input.f107A=150
    input.f107=150
    input.ap=-1.0
    input.ap_array=ctypes.POINTER(MSIS_APARRAY)()
    
    # ap_array
    tmp=MSIS_APARRAY()
    for i in range(7):
        tmp.a[i]=100.0
    input.ap_array.contents = tmp
    
    # flags structure
    flags=MSIS_FLAGS()
    sw= (ctypes.c_int * 24)()
    sw[0:]=[0,1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    flags.switches=sw
    
    # output structure
    output=MSIS_OUTPUT()
    
    # call MSIS
    ct_msis.gtd7(ctypes.byref(input),ctypes.byref(flags),ctypes.byref(output))
    
    return output
    
        
