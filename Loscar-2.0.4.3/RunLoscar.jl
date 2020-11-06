## alex cox 08/20/2020

# make and clean up the loscar files

# write the emission file format: start year, end year, GtC total
rm("dat/Emss/deccan_emss.dat")

# user input goes here [start, end, amount]
emis_co2 = [[390 390 115],[32 172 8350], [390 745 1241]];
# here ends user input
# change this to one column for emissions 08/28
emissions = zeros(1000);
times = (1:1000).*1000;

# find the years of emissions, and add it to an array
for i in 1:length(emis_co2[1,1])
    temp_length = emis_co2[i][2].+1-emis_co2[i][1]
    temp_emss = emis_co2[i][3] ./ (temp_length*1000);
    emissions[Int(emis_co2[i][1]):Int(emis_co2[i][2])] .+= temp_emss;
end

# for each new line of a dat file, put the year and sum of emissions
open("dat/Emss/deccan_emss.dat", "a+") do io
    for i in 1:1000
        write(io, string(times[i])*" "*string(emissions[i])*"\n") 
    end
end

# write the input file and define the values 

#=   ---------------------------------------------------------
#   RESTART   dat/y0prepetm-2.x.dat load restart file
#   SVSTART   dat/y0prepetm-2.x.dat save file for restart
#   EMSFILE  dat/Emss/1000_0500.dat load fossil fuel emission file
#   TSTART     0.        years     >= 0.0
#   TFINAL     200.e3    years     >  TSTART
#   CINP       1000.     Pg C      Carbon input (into atmosphere)
#   D13CIN     -55.      per mil   d13C Carbon input
#   TCIN0      0.        years     time start C input
#   TCSPAN     6.e3      years     time span C input (duration)
#   EPSLV      1.e-4     -         solver accuracy (1.e-4 ~sloppy)
#   FCSML      0.05      -         numerics: linear f_calcite drop
#                                  during dissolution if fc < FCSML.
#   KMAX       1000      -         approximate # output values.  
#                                   (range 2 to MAXSTP+1, see below)
#   TSNS       0         -         switch Temp sensitivity to pCO2 ON/OFF (1/0)
#   SCLIM      3.0       deg C     T-sensitivity (2.0-? C), requires TSNS 1
#   TEMP       25. 25. 25. 16. 16. 16. 12. 12. 12. 12. 18. 14. 12.           
#                        deg C     temp. boxes (3Low,3Int,3Deep,1High,3Tethys)
#   SAL        34.72 ... -         salinity boxes
#   PCO2SI     1000.0    ppmv      Final steady-state pCO2, balancing
#                                  long-term Si weath fluxes (> 200 ppmv)
#   FDICI      00.00     %         change init. DIC all boxes (range +- x%)
#   FALKI      00.00     %         change init. ALK all boxes (range +- x%)
#   FPO4I      00.00     %         change init. PO4 all boxes (range +- x%)
#   THC        25.       Sv        conveyor transport
#   FBIOL      0.80      -         Efficiency Low-Lat Corg Biopump 
#                                  (useful range 0.70 to 0.95)
#   CBIOH      1.8      mol C/m2/y Export High-Lat Corg Biopump
#   RRAIN      6.7       -         rain ratio (Corg/CaCO3)
#   FSHLF      4.5       -         raise shelf relative deep rain, see below 
#                                  (arbitrary factor, useful range: 1.0-6.0)
#   FINC       15.8e12  mol C/y    Initial CaCO3 riverine flux
#   CALC       20.0e-3  mol/kg     Seawater [Ca2+]
#   MAGN       30.0e-3  mol/kg     Seawater [Mg2+]
#   SULF       14.0e-3  mol/kg     Seawater [SO42-]
#   --------------------------------------------------------- =#

rm("deccan.inp")

RESTART = "dat/deccanrestart.dat";
SVSTART = "dat/deccanrestart.dat";

EMSFILE = "dat/Emss/deccan_emss.dat";

TSTART  = 0;
TFINAL  = 1000.e3;
#CINP    = 1000;
D13CIN  = -55;
#TCIN0   = 0;
#TCSPAN  = 6.e3;
EPSLV   = 1.e-4;
FCSML   = 0.05;
KMAX    = 5000;
TSNS    = 1;
SCLIM   = 3.0;
# these two need to be strings
TEMP    = "25, 25, 25, 16, 16, 16, 12, 12, 12, 12, 18, 14, 12";
SAL     = "34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72, 34.72";
PCO2SI  = 600.0;
FDICI   = 00.00;
FALKI   = 00.00;
FPO4I   = 00.00;
THC     = 25.;
FBIOL   = 0.80;
CBIOH   = 1.8;
RRAIN   = 6.7;
FSHLF   = 4.5;
FINC    = 15.83409493e12;
CALC    = 21.0e-3;
MAGN    = 42.0e-3;
SULF    = 14.0e-3;

# convert data into a .inp file -------------

# all the currently declared variables (including some Julia defaults which LOSCAR ignores)
inputvars = names(Main)[1:end-2];

# on each new line, add the variable name and its value.

open("deccan.inp", "a+") do io
    for i in 1:length(inputvars)
    write(io, string(inputvars[i])*" "*string(eval(inputvars[i]))*'\n')
    end
end;

# do the loscar run -------------------------

loscarrun = `./loscar.x deccan.inp`;

run(loscarrun);

# extract the required data files -----------

# FILE.dat  UNIT     VARIABLE
#= -----------------------------------------------------------------
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
------------------------------------------------------------------ =#

