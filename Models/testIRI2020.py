import numpy
import datetime
import iri2020
import pylab as plt

# initial test
dt = datetime.datetime(2021, 3,1,21,0) # some random time
glat = 45.
glon = 60.
altmin = 100.
altmax = 1000.
altinc = 1.

# couple of important notes:
# 1. Michael sets the JF matrix in fortran in the iri_driver.f90 code.  
# cannot directly 

# SIMOUT = ["ne", "Tn", "Ti", "Te", "nO+", "nH+", "nHe+", "nO2+", "nNO+", "nCI", "nN+"]
# outDict['Ne'] = outf[0,0:zaltkm.shape[0]]
#         outDict['Te'] = outf[3,0:zaltkm.shape[0]]
#         outDict['Ti'] = outf[2,0:zaltkm.shape[0]]
#         outDict['O+'] = outf[4,0:zaltkm.shape[0]]
#         outDict['O2+'] = outf[7,0:zaltkm.shape[0]]
#         outDict['NO+'] = outf[8,0:zaltkm.shape[0]]
#         outDict['N+'] = outf[10,0:zaltkm.shape[0]]
#         outDict['Altitude'] = zaltkm
iriIn = iri2020.IRI(dt, [altmin,altmax,altinc], glat,glon)

iriDict = dict()
iriDict['Ne'] = iriIn.ne.data
iriDict['Ti'] = iriIn.Ti.data
iriDict['Te'] = iriIn.Te.data
totalmass = iriIn['nO+'].data + iriIn['nH+'].data + iriIn['nHe+'].data +\
             iriIn['nO2+'].data + iriIn['nNO+'].data + iriIn['nN+'].data

# these are technically the concentrations.  I set it up to calculate those instead
iriDict['O+'] = iriIn['nO+'].data/totalmass
iriDict['O2+'] = iriIn['nO2+'].data/totalmass
iriDict['NO+'] = iriIn['nNO+'].data/totalmass
iriDict['N+'] = iriIn['nN+'].data/totalmass
iriDict['Altitude'] = iriIn['alt_km'].data
print(iriIn)

# run a couple of test plots.
plt.figure(101)
plt.plot(iriDict['O+'], iriDict['Altitude'], 'r-')
plt.plot(iriDict['O2+'], iriDict['Altitude'], 'b-')
plt.plot(iriDict['NO+'], iriDict['Altitude'], 'k-')

plt.figure(102)
plt.semilogx(iriDict['ne'], iriDict['Altitude'])
# testDict['CO2']

plt.show()

# print(iriDict['O+'])