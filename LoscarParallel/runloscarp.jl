## run loscar ac.gr@ 11/12
using DelimitedFiles
function runloscarp(timevals,CO2vals,Svals,co2doubling)
    # takes an array of timevals with associated CO2 and SO2 emissions,
    # runs loscar, and return the output times, pC02, and temperature 
    # given a user-specified CO2 doubling rate.
    #cd("Loscar-2.0.4.3")
    # remove the previous emissions file
    
    if isfile("LoscarParallel/deccan_CO2emss.dat")
        rm("LoscarParallel/deccan_CO2emss.dat");
    end
    if isfile("LoscarParallel/deccan_Semss.dat")
        rm("LoscarParallel/deccan_Semss.dat");
    end
    if isfile("deccan.inp")
        rm("deccan.inp")
    end
    # the values for the C and S array
    CO2vals = CO2vals;
    Svals = Svals;
    timevals = timevals;
    # the values for the input file ("deccan.inp")

    RESTART = "LoscarParallel/dat/deccanrestart.dat";
    SVSTART = "LoscarParallel/dat/deccanrestart.dat";

    EMSFILE = "LoscarParallel/deccan_CO2emss.dat";
    SEMSFILE = "LoscarParallel/deccan_Semss.dat";
    TSTART  = 0;
    TFINAL  = last(timevals);
    CINP    = 0;
    D13CIN  = -55;
    TCIN0   = 0;
    TCSPAN  = 0;
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


    inputstring = "RESTART  $RESTART
    EMSFILE $EMSFILE
    SEMSFILE $SEMSFILE
    TSTART  $TSTART
    TFINAL  $TFINAL
    CINP    $CINP
    D13CIN  $D13CIN
    TCIN0   $TCIN0
    TCSPAN  $TCSPAN
    EPSLV   $EPSLV
    FCSML   $FCSML
    KMAX    $KMAX
    TSNS    $TSNS
    SCLIM   $SCLIM
    TEMP    $TEMP
    SAL     $SAL
    PCO2SI  $PCO2SI
    FDICI   $FDICI
    FALKI   $FALKI
    FPO4I   $FPO4I
    THC     $THC
    FBIOL   $FBIOL
    CBIOH   $CBIOH
    RRAIN   $RRAIN
    FSHLF   $FSHLF
    FINC    $FINC
    CALC    $CALC
    MAGN    $MAGN
    SULF    $SULF";

    # write the three emissions files for inputs, CO2 and S emissions
    open("deccan.inp", "a+") do io
        write(io, inputstring)
    end;

    open("LoscarParallel/deccan_CO2emss.dat", "a+") do io
        for i in 1:length(CO2vals)
        write(io, string(timevals[i])*" "*string(CO2vals[i])*'\n')
        end
    end;

    open("LoscarParallel/deccan_Semss.dat", "a+") do io
        for i in 1:length(Svals)
        write(io, string(timevals[i])*" "*string(Svals[i])*'\n')
        end
    end;
    
    #do the make and run ***CHANGE THIS
    #loscarmake = `make loscar PALEO=1`

    #run(loscarmake)

    loscarrun = `./LoscarParallel/loscar.x deccan.inp`;

    run(loscarrun);

    # use the co2 doubling to turn pco2 into temperatures.
    co2doubling = co2doubling;

    time_vals = readdlm("LoscarParallel/tmv.dat", '\t', Float64, '\n')  
    pco2 = readdlm("LoscarParallel/pco2a.dat", '\t', Float64, '\n')
    temp = (pco2./600) .- 1;
    temp = co2doubling .*temp;

    #loscarcleanup = `./cleanup`;
    #run(loscarcleanup);

    return time_vals, pco2, temp

end