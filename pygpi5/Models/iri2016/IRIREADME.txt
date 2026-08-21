What I did to get IRI to work
SRK
11/13/2017

I wanted to make a couple of notes since getting IRI 2016 to work took some effort.

1. I tried to wrap using ctypes but that didn't work.
I was getting issues about trying to call a function pointer and I also think
making sure that the right type was going into the function was hard.

2. I decided to try to do f2py following Mike's way of doing it.
If you look in irisub.for, you will see this:
C*********************************************************************
CCC SRK Edits start:
Cf2py intent(in) JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,HEIBEG,HEIEND,HEISTP
Cf2py intent(out) OUTF,OARR
        call read_ig_rz
        call readapf107
CCC SRK Edits end
C*********************************************************************

These are the only edits I made to the code.
Also note you need to call read_ig_rz and readapf107 or else it won't work.

2. When I ran the compilation commands, I did have one issues that required me to use sudo.
However, I found that if I used Tim Duly's little perl fix I didn't have to do sudo again.
It is in the iriflip.for model and here is the fix:
#remove_comments:
#	perl -pi -e 's/![\.|-].*$//g' iriflip.for


3. I relied heavily on Duly's makefile.  However, in the end, when I basically wasn't
getting a F-region peak, that had to do with how JF was defined.
