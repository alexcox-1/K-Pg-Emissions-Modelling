#===============================================================#
#
# README.TXT
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
#===============================================================#

(1) Installation

Loscar should run under various OS including Linux, Windows, 
and Mac OS (see below for systems tested so far). In order 
to run Loscar, a few applications are required on all 
platforms including tar, make, and a C compiler. Under 
Windows, gzip (or similar) may also be needed to unpack the 
compressed file (Loscar-x.x.x.tar.gz, where x.x.x is the current
version number, e.g. 2.0.x). GNU applications, for instance, 
are freely available including the GNU C compiler gcc.

(1a) Linux (Unix)

Given that a C compiler is installed, installing and running 
Loscar should be fairly simple. Otherwise, most Linux 
distributions provide convenient options to install required 
packages (for instance, under Ubuntu installing gcc is as 
simple as: sudo apt-get install gcc).

When the C compiler and tools such as 'make' are available,
proceed as follows. Save the compressed file (Loscar-x.x.x.tar.gz) 
to the directory of your choice, open a terminal, cd to that 
directory and run (replace x.x.x by actual version number): 

***for the Modern version:

tar -xvf Loscar-x.x.x.tar.gz
cd Loscar-x.x.x
make loscar
./loscar.x preind.inp

Step 3 should produce a number of object files (.o) and the 
executable loscar.x. Step 4 runs the program using the input 
control file preind.inp. If you would like to compile and 
run the Paleo version:

***for the Paleo version:

tar -xvf Loscar-x.x.x.tar.gz
cd Loscar-x.x.x
make loscar PALEO=1
./loscar.x prepetm.inp


*NOTE* If you switch compiling between the Modern and Paleo 
version, first remove object files (*.o) or use the 
cleanup script. Otherwise make will only re-link, but
not re-compile all source files. 

The command line output of the run can be captured in a log 
file by e.g. 

./loscar.x preind.inp >& loscar.log &

By default, the model output will be written to .dat files
in Loscar's root directory. See below for content of the files
and how to change the output directory.

(1b) Windoze

Freely available C Compiler (gcc) for Windows may be obtained 
here, for example:

http://gcc.gnu.org/install/binaries.html
http://www.bloodshed.net/devcpp.html

gzip 1.3.12, tar 1.13, make 3.81 can be obtained here:

http://gnuwin32.sourceforge.net/packages/gzip.htm
http://gnuwin32.sourceforge.net/packages/gtar.htm
http://gnuwin32.sourceforge.net/packages/make.htm

Installing Loscar under Windoze may require more work than under 
Linux or Mac OS. Detailed instruction for Win XP are given at 
the end of this file.

Once these tools are installed under Windows, installing and 
running Loscar is similar to the steps described under Linux 
(replace x.x.x by actual version number):

***for the Modern version:

gzip -d Loscar-x.x.x.tar.gz
tar -xvf Loscar-x.x.x.tar
cd Loscar-x.x.x
make loscar
loscar.exe preind.inp

***for the Paleo version:

gzip -d Loscar-x.x.x.tar.gz
tar -xvf Loscar-x.x.x.tar
cd Loscar-x.x.x
make loscar PALEO=1
loscar.exe prepetm.inp


(1c) Mac OS

Installing Xcode, for instance, should provide the necessary 
tools, including gcc. The steps to install and run Loscar are 
then analogous to the Linux steps described above.

Xcode may be found on the system install disk. Otherwise it 
can be downloaded. However, note that the file size for latest 
Xcode version may be very large (more info is given at the 
end of this file).

Note also that due to long directory names, it could be cumbersome 
to cd to the directory in which the Loscar package was saved.
On some systems the following may work. Open a terminal,
type 'cd ' (including the space) and drag the folder containing
the Loscar package into the terminal. The folder path name 
should appear in the terminal. Proceed as described under Linux.

(2) Model input

(2a) Input of parameter values (control file)

Loscar accepts a file name as command line argument to read
parameter values from that file which will control the run. 
The Loscar package includes the input files 'preind.inp'
and 'prepetm.inp', which contain parameter values for a 
preindustrial and pre-PETM model setup. The file can be used as 
a starting point to create your own control file. If no file 
name is given as command line argument, then Loscar runs using 
internal default standard parameter values.

(2b) Input of restart values (ystart)

The subdirectory /dat contains restart files y0preind-x.x.x.dat
and y0prepetm-x.x.x.dat. These files were created at the end of a 
Loscar run and can be used to provide initial (restart) values for 
all variables to be integrated (# of lines must equal the # of 
equations to be solved). Both the name of the restart file to be 
read at the beginning of the run and the name of the restart file 
to be written at the end of the run can be set in the control file 
(see above). If no file names are specified in the control file 
(or such lines are commented out by '#'), then Loscar runs using 
internal default start values.

*NOTE* A steady-state configuration/run (all variables remain 
constant during a run), requires both initial conditions AND 
parameter values to match. In addition to saving the restart file,
it is hence useful to also save the corresponding parameters
(see e.g. 'check/parms.out-x.x.x').

(2c) Input of fossil fuel emissions

The subdirectory /dat/Emss contains fossil fuel emission input 
files such as 1000_0500.dat (see Zeebe et al., Science, 2008). The 
name of the emission file to be read by Loscar can be set in the 
control file (see above). If no emission file name is specified in 
the control file (or such line is commented out by '#'), then 
Loscar runs normal.

(3) Model output

By default, the model output is written to .dat files in Loscar's
root directory as follows (OCN = tracer properties/concentrations 
in ocean boxes, ATM = atmosphere variable, SED = sediment variables):
      
FILE.dat  UNIT     VARIABLE
-----------------------------------------------------------------
tmv       (y)           time
tcb       (deg C)   OCN temperature
dic       (mmol/kg) OCN total dissolved inorganic carbon
alk       (mmol/kg) OCN total alkalinity
po4       (umol/kg) OCN phosphate
dox       (mol/m3)  OCN dissolved oxygen
dicc      (mmol/kg) OCN DIC-13
d13c      (per mil) OCN delta13C(DIC)
d13ca     (per mil) ATM delta13C(atmosphere)
pco2a     (ppmv)    ATM atmospheric pCO2
co3       (umol/kg) OCN carbonate ion concentration
ph        (-)       OCN pH (total scale)
pco2ocn   (uatm)    OCN ocean pCO2
omegaclc  (-)       OCN calcite saturation state
omegaarg  (-)       OCN aragonite saturation state
fca       (-)       SED calcite content Atlantic
fci       (-)       SED calcite content Indian
fcp       (-)       SED calcite content Pacific
fct       (-)       SED calcite content Tethys (PALEO only)
ccda      (m)       SED calcite compens. depth Atlantic
ccdi      (m)       SED calcite compens. depth Indian
ccdp      (m)       SED calcite compens. depth Pacific
ccdt      (m)       SED calcite compens. depth Tethys (PALEO only)
------------------------------------------------------------------

In ocean tracer files, rows represent time steps and columns 
represent ocean boxes in order LA,LI,LP,IA,II,IP,DA,DI,DP,H,
(LT,IT,DT); where L=Low, I=Intermediate, D=Deep, A=Atlantic, 
I=Indian, P=Pacific, H=High latitude box, (T=Tethys).

In sediment files (calcite content), rows represent time steps 
and columns represent sediment boxes at different water depths, 
starting with the shallowest and ending with the deepest sediment 
box.

In addition, several .out files are generated:

FILE     CONTENT
--------------------------------------------------------
parms    various parameters  
ystart   start values. #lines = #equations to solve
dsv      (m) depth of sediment boxes  
zv       (m) depth variable (continuous)
--------------------------------------------------------

To plot the results of a Loscar run, the MATLAB script 
'PlotLoscar.m' included in the package may be used. You may 
of course use a plot program of your choice and supply your 
own plot routine.

The shell script 'cleanup' can be used to clean the Loscar root 
directory. *NOTE* Under Windoze, rename the file 'cleanup.dos'
to 'cleanup.bat' before using.

CAUTION: all executables, data, object files etc. will be deleted!!!
Only source (.c,.h), input, and some other text files remain. 

The directory to which Loscar writes output may be changed as 
follows (Linux example). Simply run Loscar from a different directory 
by copying the required input files to that directory (e.g. control, 
restart, emission) and specify the path to the executable (in the 
Loscar root directory). Example:

cd ~/myloscar
/usr/Loscar-x.x.x/loscar.x preind.inp

Output files will be written to ~/myloscar
Don't forget to change the input file paths in the control file!

(4) Tests

The package includes a shell script 'runtest', which executes 
several test runs. Invoke the script by running e.g. 

***for the Modern version:

./runtest >& runtest.log &


***for the Paleo version:

./runtest PALEO >& runtest.log &

and examine runtest.log. Except for one run, no warnings or errors 
should appear (in the log file or terminal) and all runs should be 
completed indicated by the line: "LOSCAR V x.x.x  Done."
Check results for "Final Atm CO2" in check/test/TestResults-x.x.x.txt 

*NOTE* Rename 'runtest.dos' to 'runtest.bat' under Windows 
before use: 

runtest.bat > runtest.log 2>&1
or
runtest.bat PALEO > runtest.log 2>&1

(5) Subdirectories

'/check' contains check files for some parameters that can be 
compared to the outcome generated on different platforms. The 
check files were generated using default parameters and a 
fossil fuel emission scenario in which 1000 Pg C are emitted 
over 500 years (see Zeebe et al., Science, 2008). Input: 
dat/Emss/1000_0500.dat. Note that by default, ocean temperature 
sensitivity to pCO2 is OFF ('TSNS 0'). For exact results shown 
in Zeebe et al., set 'TSNS 1'.

*NOTE* Different OS/machines/optimizations produce slightly 
different outcome. In other words, the numbers in the check 
files will typically only agree with outcome generated on 
different platforms/optimizations to 8 significant digits
or so (for ocean tracers; for pCO2 occasionally less).

'/check/test' contains input files and results for the tests 
described under item 'Tests' above.

'/dat' contains preindustrial and pre-PETM restart files 
'y0preind-x.x.x.dat' 'y0prepetm-x.x.x.dat' (for corresponding 
parameter file, see 'check/parms.out-x.x.x').

'/dat/Emss' contains input files of fossil fuel emission 
scenarios (see Zeebe et al., Science, 2008).

'/docs' contains pdf copies of several papers relevant to 
model development and applications.


(6) Systems tested so far

Linux,  i486,    gcc 4.4.3
Linux,   x86_64, gcc 4.4.5 (compiler option: -m32)
Mac,    i686_64, gcc 4.0.1
Win XP,  x86,    gcc 3.4.2 (+GNU make for Win)


(42a) Mac Xcode


Issue (Nov 2011): latest Xcode 4 installer download is > 4 GB.
Note: Xcode 4.0.2 requires Mac OS X 10.6.6 or later.

One alternative option is to install an older, smaller and free 
version of Xcode. Xcode may be found on the system install disk.
Otherwise login to your apple developer account. If you don't 
have one, create it here:

http://developer.apple.com/programs/start/register/create.php

Search for Xcode (Apple seems to be constantly changing the
site). For the latest version they may require "Mac Developer 
Program member" credentials with annual fee! So look for an 
older (free) version of Xcode and download.

For instance, Xcode 3.1.x should install on Mac OS X 10.5.0 and 
higher. File size is typically less than 1 GB. If the installer 
allows customizing during install, you may be able to drop some 
packages but you will need Unix command line support.

At this point it seems that a functional gcc on Mac OS requires 
Xcode. As a result, Mac users have to download and install a
large (unnecessary) Apple package. Ironically, the original
(small) gcc package itself is GNU software!


(42b) Windoze XP Installation (example)

get C compiler from:
http://www.bloodshed.net/devcpp.html
go to download page 
download:
Dev-C++ 5.0 beta 9.2 (4.9.9.2) (9.0 MB) with Mingw/GCC 3.4.2
install to C:\Dev-Cpp\bin

Win XP - add 'C:\Dev-Cpp\bin' to path: 
go to 
> System Properties 
> Advanced
> Environment Variables
> System variables 
> Path (edit)
add: ;C:\Dev-Cpp\bin 
don't forget ';' to separate path entries!

get gzip from:
http://gnuwin32.sourceforge.net/packages/gzip.htm
download:
Complete package, except sources Setup
install to C:\Program Files\GnuWin32
Win XP - add to path (see above):
C:\Program Files\GnuWin32\bin

get gtar from:
http://gnuwin32.sourceforge.net/packages/gtar.htm
download:
Binaries Setup
install to C:\Program Files\GnuWin32
(already in path) 

get make from:
http://gnuwin32.sourceforge.net/packages/make.htm
download:
Complete package, except sources Setup
install to C:\Program Files\GnuWin32
(already in path) 
