FFILE	= irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for
OBJS	= irisub.o irifun.o iritec.o iridreg.o igrf.o cira.o iriflip.o
OBJS1	= $(FFILE).o $(OBJS)
EXE  	= iritest
CINCLUDE =
BIN      =
LIBS     =
GLIB     =
BLIB     =
CC       = gcc -ffast-math
FC       = gfortran
SO 		 = iri2016
PYF 	 = iri2016.pyf
F2PY	= "/Users/srkaeppler/Env3/conductivity/bin/f2py"
#FFLAGS   =-Wall -ansi -D_GNU_SOURCE -D_DARWIN
# I am not sure I know what all this crap does...
#FFLAGS   = -fno-automatic -finit-local-zero -O2 -ffast-math -fPIC
FFLAGS 	  = -w -O3 -fbacktrace -fno-automatic -fPIC -lstdc++ -std=legacy -ffast-math
#-funroll-all-loops -fexpensive-optimizations

$(FFILE):  $(OBJS1)
	@echo "Creating ...."
	$(FC) $(FFLAGS) $(OBJS1) $(BLIB) $(LIBS) -o $(FFILE)

test:
	$(FC) $(FFLAGS) $(BLIB) $(LIBS) -o $(EXE) iritest.for $(FFILE)

mkso:
	$(FC) -c $(FFILE)
	$(FC) $(FFLAGS) -shared -fPIC -fcheck=all  $(OBJS) -o $(SO)

# this follows closely to what Duly and Sabestian did
# not sure what the f77 = "-w" flag does
mkf2py:
	mktemp -d
	$(F2PY) --overwrite-signature irisub.for -m $(SO) -h $(PYF)
	$(FC) $(FFLAGS) -c $(FFILE)
	$(F2PY) --f77flags="-w"  -c $(PYF) $(FFILE)
echotest:
	@echo $(F2PY)

clean:
	@echo "Cleaning ... "
	rm -f *.o
	rm -f *.so
