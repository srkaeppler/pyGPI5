


import os
import numpy

inSo = '/Users/srkaeppler/research/data/NSF_Dregion_ParticlePrecipitation/Models/iri2016/iri2016.so'
iriso = cdll.LoadLibrary(inSo)

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

JF = numpy.array([  1,1,1,0,0,0,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,
      0,1,0,1,1,1,1,0,0,0,
      1,1,0,1,0,1,1,1,0,0,
      0,0,0,0,0,0,0,0,0,0,
      ]).ctypes.data_as(POINTER(c_float))
JF[11]=0 # will pipe out print statements to a text file called messages.txt
JMAG = c_float(0) # =0 geographic   = 1 geomagnetic coordinates
IYYYY = c_float(2010.)
MMDD = c_float(0101.)
DHour = c_float(21.50)
Lat = c_float(45.)
Lon = c_float(45.)
Heibeg = c_float(80.)
Heiend = c_float(150.)
step = c_float(1.0)
OARR = numpy.zeros(100).ctypes.data_as(POINTER(c_float))
OUTF = numpy.zeros([20,1000]).ctypes.data_as(POINTER(c_float))
# thoGrid.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

iriso.iri_sub_.argtype[POINTER(c_float), POINTER(c_float),\
                       POINTER(c_float), POINTER(c_float),\
                       POINTER(c_float), POINTER(c_float),\
                       POINTER(c_float), POINTER(c_float),\
                       POINTER(c_float), POINTER(c_float),\
                       POINTER(c_float), POINTER(c_float)]

# iriso.iri_sub_(byref(JF),byref(JMAG),byref(Lat),byref(Lon),byref(IYYYY),\
#                 byref(MMDD),byref(DHour),byref(Heibeg),byref(Heiend),byref(step),\
#                 byref(OUTF),byref(OARR))

# IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
# # 183      &    HEIBEG,HEIEND,HEISTP,OUTF,OARR)
