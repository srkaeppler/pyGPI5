


import os
import numpy
import sys
import matplotlib.pyplot as plt

# sys.path.append('./iri2016/')
#import iri16duly
# iri16duly.read_ig_rz()
# iri16duly.readapf107()
import iri2016

 #       SUBROUTINE IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
 # 183      &    HEIBEG,HEIEND,HEISTP,OUTF,OARR)
 # 184 C-----------------------------------------------------------------
 # 185 C
 # 186 C INPUT:  JF(1:50)      true/false switches for several options
 # 187 C         JMAG          =0 geographic   = 1 geomagnetic coordinates
 # 188 C         ALATI,ALONG   LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
 # 189 C         IYYYY         Year as YYYY, e.g. 1985
 # 190 C         MMDD (-DDD)   DATE (OR DAY OF YEAR AS A NEGATIVE NUMBER)
 # 191 C         DHOUR         LOCAL TIME (OR UNIVERSAL TIME + 25) IN DECIMAL
 # 192 C                          HOURS
 # 193 C         HEIBEG,       HEIGHT RANGE IN KM; maximal 100 heights, i.e.
 # 194 C          HEIEND,HEISTP        int((heiend-heibeg)/heistp)+1.le.100

# JF = numpy.array([  1,1,1,0,0,0,1,1,1,1,
#       1,1,1,1,1,1,1,1,1,1,
#       0,1,0,1,1,1,1,0,0,0,
#       1,1,0,1,0,1,1,1,0,0,
#       0,0,0,0,0,0,0,0,0,0,
#       ])
JF = numpy.ones(50)
JF[3] = 0.
JF[4] = 0
JF[5] = 0
JF[11]=1 # will pipe out print statements to a text file called messages.txt
JMAG = 0 # =0 geographic   = 1 geomagnetic coordinates
IYYYY = 2017
MMDD = -01
DHour = 10.25
Lat = 45.
Lon = 228.
Heibeg = 80.
Heiend = 500.
step = 0.5

Altitude = numpy.arange(Heibeg,Heiend,step)

outf,oarr = iri2016.iri_sub(JF,JMAG,Lat,Lon,IYYYY,MMDD,DHour,Heibeg,Heiend,step)
# print 'outf', outf
# print 'oarr', oarr

Ne = outf[0,0:Altitude.shape[0]]
plt.figure()
plt.semilogx(Ne,Altitude)
plt.show()
# # OUTF(1,*)  ELECTRON DENSITY/M-3
# # C               OUTF(2,*)  NEUTRAL TEMPERATURE/K
# # C               OUTF(3,*)  ION TEMPERATURE/K
# # C               OUTF(4,*)  ELECTRON TEMPERATURE/K
# Ne = outf[0,:]

# plt.figure()
# plt.semilogx(Ne)
