#===============================================================
# Loscar Makefile
#
# LOSCAR Model: Long-term Ocean-atmosphere-Sediment 
#                 CArbon cycle Reservoir Model
# 
# *** LOSCAR comes with ABSOLUTELY NO WARRANTY ***		 
# *** Use at your own risk. DO NOT DISTRIBUTE  ***
#
# When using LOSCAR, cite as:
# 
# Zeebe, R. E., Zachos, J. C., and Dickens, G. R. Carbon dioxide
# forcing alone insufficient to explain Paleocene-Eocene Thermal
# Maximum warming. Nature Geoscience, doi:10.1038/NGEO578 (2009)
# 
# Richard E. Zeebe
# School of Ocean and Earth 
# Science and Technology 
# Department of Oceanography 
# University of Hawaii at Manoa
# 1000 Pope Road, MSB 504
# Honolulu, HI 96822, USA
# email: loscar.model@gmail.com
#
# updates:
# 
#  11/29/11 gcc 4.6.1. omit -ansi (takes constants as long 
#           double): comparison fails (epslvr) + different results
#           for old results use: -ansi -msse -mfpmath=sse
#  10/17/11 added Paleo option. use: make loscar PALEO=1
#           Default: Modern.
#  04/05/11 added clean and .h dependencies
#  02/13/11 new file
#
#===============================================================

# name of C complier
CC = gcc

# compile + all warnings (-ansi -pedantic -Wextra) -O3? -msse -mfpmath=sse
#CFLAGS = -c -Wall -pedantic -Wextra
CFLAGS = -c

# optimize: -O2 -O3
#CFLAGS = -c -O2 

# 32 bit option required on 64 bit machine?
#CBIT = -m32
CBIT = 

# libraries
LIBS = -lm 

# Paleo version? Default: Modern
ifeq ($(PALEO),1)
MPSET = -D FTYS
else
MPSET =
endif

# OS: uname
UNAME := $(shell uname)

# name of executable (default)
EXENAME = loscar.x	

# name of executable (Windows)
WINDIR := $(shell set WINDIR)

ifneq ($(WINDIR),)
EXENAME = loscar.exe	
endif


# Compile source files
csys.o: csys.c defs.h utils.h
	$(CC) $(CFLAGS) $(CBIT) $(MPSET) csys.c
emiss.o: emiss.c defs.h utils.h common.h
	$(CC) $(CFLAGS) $(CBIT) $(MPSET) emiss.c
initfree.o: initfree.c defs.h utils.h common.h
	$(CC) $(CFLAGS) $(CBIT) $(MPSET) initfree.c
initstart.o: initstart.c defs.h utils.h common.h
	$(CC) $(CFLAGS) $(CBIT) $(MPSET) initstart.c
loscarderivs.o: loscarderivs.c defs.h utils.h common.h
	$(CC) $(CFLAGS) $(CBIT) $(MPSET) loscarderivs.c
matrix.o: matrix.c defs.h utils.h
	$(CC) $(CFLAGS) $(CBIT) $(MPSET) matrix.c
readparms.o: readparms.c defs.h utils.h common.h
	$(CC) $(CFLAGS) $(CBIT) $(MPSET) readparms.c
solver.o: solver.c defs.h utils.h
	$(CC) $(CFLAGS) $(CBIT) $(MPSET) solver.c
thmfun.o: thmfun.c defs.h utils.h
	$(CC) $(CFLAGS) $(CBIT) $(MPSET) thmfun.c
utils.o: utils.c utils.h
	$(CC) $(CFLAGS) $(CBIT) $(MPSET) utils.c
writedat.o: writedat.c defs.h utils.h
	$(CC) $(CFLAGS) $(CBIT) $(MPSET) writedat.c
loscar.o: csys.o emiss.o initfree.o initstart.o matrix.o loscarderivs.o readparms.o solver.o thmfun.o utils.o writedat.o loscar.c defs.h utils.h common.h
	$(CC) $(CFLAGS) $(CBIT) $(MPSET) loscar.c
# Link object files
loscar: loscar.o
	$(CC) $(CBIT) -o $(EXENAME) csys.o emiss.o initfree.o initstart.o matrix.o loscarderivs.o readparms.o solver.o thmfun.o utils.o writedat.o loscar.o $(LIBS)

info: 
	@echo $(UNAME)
	@echo $(WINDIR)
	@echo $(EXENAME)

# clean
clean:
	\rm *.o *~ loscar.x loscar.exe


