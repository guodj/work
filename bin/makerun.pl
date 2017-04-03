#!/usr/bin/perl -s

$maindir = "/home/guod/data/";
$datadir = "/raid3/Amie/Data/srcData/" if (!$datadir);
$bindir = "/home/guod/work/bin";

$help = $h;

$cputime = 23.5 if (!$cputime);
$nlats = 2 if (!$nlats);
$nlons = 2 if (!$nlons);
$dtout = "1800" if (!$dtout);
$amien = "none" if (!$amien);
$amies = "none" if (!$amies);
$startday   = $day   if (!$startday);
$startmonth = $month if (!$startmonth);
$startyear  = $year if (!$startyear);

$endday   = $day   if (!$endday);
$endmonth = $month if (!$endmonth);
$endyear  = $year  if (!$endyear);

$help = 1 if (!$startyear);
$help = 1 if (!$startmonth);
$help = 1 if (!$startday);

$f107set = 1 if ($f107);
$BuildImf = 1 if (length($imf) == 1);
$BuildPower = 1 if (length($power) == 1);
$BuildChamp = 1 if (length($champ) == 1);
$BuildGrace = 1 if (length($grace) == 1);
$BuildEUV   = 1 if (length($euv) == 1);

$l = &leapDay($startyear,3);
@theJulianDate = ( 0, 31, 59+$l, 90+$l, 120+$l, 151+$l, 
		   181+$l, 212+$l, 243+$l, 273+$l, 304+$l, 334+$l );

$startjday = $theJulianDate[$startmonth-1] + $startday;

$endjday   = $theJulianDate[$endmonth-1] + $endday;

$noymd = 1 if ($endjday != $startjday);

if ($endjday < $startjday) {
    if ($endyear > $startyear) {
	$endjday = $endjday + 365 + &leapDay($startyear,12);
    } else {
	print "Start Day after End Day!\n";
	$help = 1;
    }
}

$MaxjDay = 365 + &leapDay($startyear,12);

$dtsat = "60.0" if (!$dtsat);

if ($help) {
    print "\n";
    print "Usage:\n";
    print " makerun.pl -year=yyyy -month=mm -day=dd [-gitm or -swmf]\n";
    print "            -startyear=yyyy -endyear=yyyy\n";
    print "            -startmonth=mm -endmonth=mm\n";
    print "            -startday=dd -endday=dd\n";
    print "            -nlats=N -nlons=M [-restart]\n";
    print "            [-imf] [-power] \n";
    print "            e.g. \n";
    print "            makerun.pl -year=2003 -month=01 -startday=4 -endday=6 -gitm -power -imf -rm \n";
    print "                This will make UAM.in.Start, UAM.in.Restart, power.dat and imf.dat files\n";
    print "\n";
    print "            e.g. \n";
    print "            makerun.pl -year=2005 -startmonth=08 -endmonth=09 -startday=30 -endday=1 -nlats=12 -nlons=8 -gitm -power -imf -rm -dynamo -champ -grace -euv\n";
    print "\n";
    print "            -imf[=filename]  [-wind] (wind can be used for IMF)\n";
    print "            -champ[=filename]  (can build a champ satellite file for you\n";
    print "            -grace[=filename]  (can build a grace satellite file for you\n";
    print "            [-guvi]\n";
    print "            [-f107=value]\n";
    print "            -sats=file1,file2,file3,etc  [-buildsat] (see below)\n";
    print "            -mhdsats (searches for the standard MHD satellites)\n";
    print "            -dontmove (doesn't move the satellite files)\n";
    print "            -noymd (don't put year-month-day info into the inputfile\n";
    print "            -job [-cont] [-ncpus] \n";
    print "              (can build a job file \n";
    print "              cont adds a qsub to the end of the job file,\n";
    print "               ncpus is nlatsxnlons by default)\n";
    print "            [-power] (use NOAA HPI file)\n";
    print "            [-tar] (create a tarball of the files)\n";
    print "            [-cputime] (max cputime for the run, in hours)\n";
    print "            [-euv] (use fism euv data)\n\n";
    print "            -buildsat the code will create satellite files\n";
    print "               from the file1,file2,etc, which have to be 3-line files\n";
    print "               with longitude, latitude, and altitude as the 3-lines\n";
    print "               if -build is not specified, it assumes that the files are\n";
    print "               the actual satellite files\n\n";
    print "            -champ[=file] if the file is included, that is added as a satellite\n";
    print "               file.  If not, it builds a satellite file from the CHAMP database\n";
	
    exit(1);
}

# Determine the next day (this is the #ENDTIME in the UAM.in file)
$nextjday = $endjday + 1;
$nextyear = $endyear;
if ($endjday == $MaxjDay) {
    $nextyear++;
    $nextjday = $nextjday - $MaxjDay;
}
$nextjday = $nextjday - $MaxjDay if ($nextjday > $MaxjDay); 
($nextmonth,$nextday) = jday_to_md($nextyear,$nextjday);

if ($amie) {
    open(RUNAMIE,">runamie.sh");
    print RUNAMIE "#!/bin/sh\n";
    print RUNAMIE "\n";
}

for ($jDayCounter = $startjday; $jDayCounter <= $endjday; $jDayCounter++) {

    if ($restart || $jDayCounter != $startjday) {
	$restart = "T";
    } else {
	$restart = "F";
    }

    $year = $startyear;
    $jday = $jDayCounter;
    if ($jday > $MaxjDay) {
	$year++;
	$jday = $jday - $MaxjDay;
    }

    ($month,$day) = jday_to_md($year,$jday);

    if ($month < 10) {
	$month = $month + 0;
	$month = "0$month";
    }

    if ($day < 10) {
	$day = $day + 0;
	$day = "0$day";
    }

    print "Processing: $year $month $day\n";

    $f107 = getf107($year,$month,$day) if (!$f107set);
    print "Got f107 : $f107\n";

    $ymd = $year.$month.$day;
    $ym = $year.$month;

    $hpi = 20.0 if (!$power);

    if ($gitm) {

	$command = "pwd > remoteloc";
	system $command;

	if ($job) {
	
	    $ncpus = $nlons*$nlats if (!$ncpus);
	    $file = "job.$year$month$day";
	    print "Creating file $file\n";
	
	    open(JOB,">$file");

#	    if ($nasa) {
		$walltime = 24 if (!$walltime);
		$nodes = $nlats*$nlons/20;
		print JOB "#PBS -S /bin/csh\n";
		print JOB "#PBS -N gitm.$ym\n";
		print JOB "\n";
		print JOB "# set the number of CPU-s by changing select: nProc = select*mpiprocs\n";
		print JOB "# ncpus should be fixed to 12 (or 20) because there are 12 (or 20) cores on \n";
		print JOB "# each node. The faster nodes mean faster execution, but also more charged\n";
		print JOB "# to the account. \n";
		print JOB "\n";
		print JOB "# To run on the 20-core ivy bridge nodes\n";
		print JOB "#PBS -l select=".$nodes.":ncpus=20:model=ivy\n";
		print JOB "#PBS -l walltime=".$walltime.":00:00\n";
		print JOB "#PBS -j oe\n";
		print JOB "#PBS -m e\n";
		print JOB "\n";
		print JOB "cd \$PBS_O_WORKDIR\n";
		print JOB "\n";
		print JOB "module load comp-intel/2012.0.032 mpi-sgi/mpt.2.11r8\n";
		print JOB "\n";
#	    }
	    print JOB "\n";
	    print JOB "rm -f UAM.in ; ln -s UAM.in.Start UAM.in\n";
	    print JOB "#rm -f UAM.in ; ln -s UAM.in.Restart UAM.in\n";
#	    if ($nasa) {
		print JOB "mpiexec GITM.exe > log.$year$month$day\n";
#	    }
	    print JOB "\n";
#	    print JOB "if (-e GITM.DONE) then\n";
#	    print JOB "   cd UA\n";
#	    print JOB "   mv restartOUT/ restartOUT.$year$month$day.$job\n";
#	    print JOB "   mkdir restartOUT\n";
#	    print JOB "   ./pGITM\n" if (!$ejet);
#	    print JOB "   mv data data.$year$month$day.$job\n";
#	    print JOB "   mkdir data\n";
#	    print JOB "   rm -f restartIN; ln -s restartOUT.$year$month$day.$job restartIN\n";
#	    print JOB "   \n";
#	    print JOB "   cd ..\n";
#	    if(length($cont) == 1){
#		print JOB "   qsub job.$year$month$nextday\n";
#	    }
#	    print JOB "endif";
#	    print JOB "\n";
	
	    close(JOB);
	}
	$job = 0;
    }

    if ($BuildPower) {
	print "Grabbing hemispheric power file... ";
	$powerfile = "power_".$year.".txt";
	$command = "cp ".$maindir."HemisphericPower/".$powerfile." .";
	system $command;
	open(IN,"<$powerfile");
	open(OUT,">power$ymd.dat");
	$sjday = $jday;
	if ($year == 2006) {
	    $sjday = "0$jday" if ($jday<100);
	    $sjday = "0$jday" if ($jday<10);
	} else {
	    $sjday = " $jday" if ($jday<100);
	    $sjday = " $jday" if ($jday<10);
	}
	while(<IN>) {
	    if (/#/ || 
                /^\n/ || 
                /^\r\n/ || 
                /^$year$/ || 
                /^$year\r\n/ || 
                /^$year-$month-$day/ || 
                /$sjday\d\d\d\d\b/) {
		print OUT $_;
	    }
	}
	close(IN);
	close(OUT);
	$command = "/bin/rm -f $powerfile";
	system $command;

	$nLines = 0;
	open(IN,"power$ymd.dat");
	while (<IN>) {
	    $nLines++;
	}
	close(IN);

	print "lines read ... $nLines\n";

    }

    if ($BuildImf) {

	print "Creating IMF file... ";

	open(IDLFILE,">.idl_imf_input");

	print IDLFILE ".r cdf_to_mhd\n";
	if (!$wind) {
	    print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/data/ace/mag/level_2_cdaweb/mfi_h0/$year/ac_h0_mfi_".$year.$month.$day."_v*.cdf\n";
	} else {
	    print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/data/wind/mfi/mfi_h0/$year/wi_h0_mfi_".$year.$month.$day."_v*.cdf\n";
	}
	print IDLFILE "i\n";
	print IDLFILE "imf$ymd.dat\n";
	if (!$wind) {
	    print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/data/ace/swepam/level_2_cdaweb/swe_h0/$year/ac_h0_swe_".$year.$month.$day."_v*.cdf\n";
	} else {
	    print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/data/wind/swe/swe_h0/$year/wi_h0_swe_".$year.$month.$day."_v*.cdf\n";
	}
	print IDLFILE "d\n";

	close(IDLFILE);

	$command = "idl < .idl_imf_input >& .log_imf";
	system $command;

	if ($noymd) {
	    if ($yearmonth) {
		$imf = "imf$ym.dat";
	    } else {
		$imf = "imf.dat";
	    }
	} else {
	    $imf = "imf$ymd.dat";
	}

	$nLines = 0;
	open(IN,"imf$ymd.dat");
	while (<IN>) {
	    $nLines++;
	}
	close(IN);

	print "lines read ... $nLines\n";

    }

    if ($BuildChamp) {

	print "Creating CHAMP file... ";

	open(IDLFILE,">.idl_champ_input");

	print IDLFILE ".r thermo_convert_champfiles\n";
	print IDLFILE "\n";
	print IDLFILE "$year\n";
	print IDLFILE "$month\n";
	print IDLFILE "$day\n";

	close(IDLFILE);

	$command = "idl < .idl_champ_input >& .log_champ";
	system $command;

	if ($noymd) {
	    if ($yearmonth) {
		$champFile = "champ$ym.dat";
	    } else {
		$champFile = "champ.dat";
	    }
	} else {
	    $champFile = "champ$year$month$day.dat";
	}

	$nLines = 0;
	open(IN,"champ$ymd.dat");
	while (<IN>) {
	    $nLines++;
	}
	close(IN);

	print "lines read ... $nLines";

	$noChamp = 0;
	$noChamp = 1 if ($nLines < 1);

	print " yikes! no Champ file found!" if ($noChamp);
	print "\n";

    }

    if ($BuildGrace) {

	print "Creating Grace file... ";

	open(IDLFILE,">.idl_grace_input");

	print IDLFILE ".r thermo_convert_gracefiles\n";
	print IDLFILE "\n";
	print IDLFILE "$year\n";
	print IDLFILE "$month\n";
	print IDLFILE "$day\n";

	close(IDLFILE);

	$command = "idl < .idl_grace_input >& .log_grace";
	system $command;

	if ($noymd) {
	    if ($yearmonth) {
		$graceFile = "grace$ym.dat";
	    } else {
		$graceFile = "grace.dat";
	    }
	} else {
	    $graceFile = "grace$year$month$day.dat";
	}

	$nLines = 0;
	open(IN,"grace$ymd.dat");
	while (<IN>) {
	    $nLines++;
	}
	close(IN);

	print "lines read ... $nLines";

	$noGrace = 0;
	$noGrace = 1 if ($nLines < 1);

	print " yikes! no Grace file found!" if ($noGrace);
	print "\n";

    }

    if ($magnetometer || $magnetometers) {

	print "Downloading magnetometer data from SuperMag Site\n";

	if ($year >= 2004) {
	    $command = "$bindir/getSuperMAG-day.sh $ymd";
	}
	if ($year == 2003) {
	    $command = "$bindir/getSuperMAG-day-in-2003.sh $ymd";
	}
	if ($year < 2003) {
	    print "Can't make magnetometer file before 2003!\n";
	    exit(1);
	}

	system $command;

	print "Reformating the magnetometer file ... ";

	open(IDLFILE,">.idl_magnetometer_input");

	print IDLFILE ".r supermag_AMIE_main\n";
	print IDLFILE "$ymd.0000.supermag.txt\n";

	close(IDLFILE);

	$command = "idl < .idl_magnetometer_input >& .log_magnetometer";
	system $command;

    }

    if ($BuildEUV) {
	
	print "Grabbing EUV file\n";

	$command = "cp ".$maindir."FISM/BinnedFiles/".$year."/fismflux".$year.$month.$day.".dat .";
	system $command;
    }

    if ($pkr) {
	print "Creating PKR satellite files\n";
	makeFPIfiles("pkr",212.5,65.1,240.0);
    }

    if ($fyu) {
	print "Creating FYU satellite files\n";
	makeFPIfiles("fyu",214.73,66.56,240.0);
    }

    if ($pk1) {
	print "Creating PKR satellite file\n";
	makeSingleFile("pkr",212.5,65.1,240.0);
    }

    if ($msh) {
	print "Creating MSH satellite file\n";
	makeSingleFile("msh",288.51,42.62,300.0);
    }

    if ($stf) {
	print "Creating STF satellite file\n";
	makeSingleFile("stf",282.623,67.08,300.0);
    }

    if ($jic) {
	print "Creating Jicamarca satellite file\n";
	makeSingleFile("jic",283.13,-11.95,300.0);
    }

    if ($eis) {
	print "Creating Eiscat satellite file\n";
	makeSingleFile("eis",16.02,78.09,300.0);
    }

    # -------------------------------------------------------------
    # Build UAM.in file
    # -------------------------------------------------------------

    if ($gitm) {

	if ($jDayCounter == $startjday ||
	    $jDayCounter == $endjday) {

	    if ($nextmonth < 10) {
		$nextmonth = $nextmonth + 0;
		$nextmonth = "0$nextmonth";
	    }

	    if ($nextday < 10) {
		$nextday = $nextday + 0;
		$nextday = "0$nextday";
	    }

	    if ($startmonth < 10) {
		$startmonth = $startmonth + 0;
		$startmonth = "0$startmonth";
	    }

	    if ($startday < 10) {
		$startday = $startday + 0;
		$startday = "0$startday";
	    }

	    if ($startjday == $endjday) {
		$file = "UAM.in.$year$month$day";
	    } else {
		$file = "UAM.in.Start" if ($jDayCounter == $startjday);
		if ($jDayCounter == $endjday) {
		    $file = "UAM.in.Restart"; 
		    $restart = "T";
		}
	    }

	    #$file = $file.".$job" if ($job);
	    print "Creating file $file\n";
	    open(UAM,">$file");

	    if ($tar) {
		$tarfile = "uam.".$year.$month.$day.".tgz";
		@tarlist = $file;
	    }

	    print UAM "\n";
	    print UAM "#DEBUG\n";
	    print UAM "0	debug level\n";
	    print UAM "0	cpu to watch\n";
	    print UAM "10.0	dt between normal code output to stdout\n";
	    print UAM "F	usebarriers - forces the code to stop and wait more often\n";
	    print UAM "\n";
	    print UAM "#RESTART\n";
	    print UAM "$restart\n";
	    print UAM "\n";
	    print UAM "#GRID\n";
	    print UAM "$nlons	    lons\n";
	    print UAM "$nlats	    lats\n";
	    print UAM "-90.0	    minimum latitude to model\n";
	    print UAM "90.0	    maximum latitude to model\n";
	    print UAM "0.0	    start longitude to model\n";
	    print UAM "0.0	    end longitude to model\n";
	    print UAM "\n";
	    print UAM "#DIFFUSION\n";
	    print UAM "T\n";
	    print UAM "500.0	   Eddy Diffusion Coefficient (Should be about 37.5 for 1-D runs)\n";
	    print UAM "0.00100	   Total Eddy Diffusion applied at alts below this pressures level\n";
	    print UAM "0.00010	   No Eddy Diffusion at altitudes above this pressure level\n";
	    print UAM "\n";
	    print UAM "#PHOTOELECTRON\n";
	    print UAM "0.02	   Efficiency of photoelectron heating\n";
	    print UAM "\n";
	    print UAM "#NEUTRALHEATING\n";
	    print UAM "0.05	   Efficiency of photoelectron heating\n";
	    print UAM "\n";
	    print UAM "#THERMALCONDUCTION\n";
	    print UAM "3.6e-4	   Thermal conductivity (o2)\n";
	    print UAM "5.6e-4	   Thermal conductivity (o)\n";
	    print UAM "0.69	   Thermal conductivity (^s)\n";
	    print UAM "\n";
	    print UAM "#CPUTIMEMAX\n";
	    $c = $cputime*3600.0;
	    print UAM "$c	   Maximum amount of cputime to use before stopping the code\n";
	    print UAM "\n";
	    print UAM "#TIMESTART\n";
	    print UAM "$year        year\n";
	    print UAM "$startmonth        month\n";
	    print UAM "$startday        day\n";
	    print UAM "00        hour\n";
	    print UAM "00        minute\n";
	    print UAM "00        second\n";
	    print UAM "\n";
	    print UAM "#TIMEEND\n";
	    print UAM "$nextyear      year\n";
	    print UAM "$nextmonth        month\n";
	    print UAM "$nextday        day\n";
	    print UAM "00        hour\n";
	    print UAM "00        minute\n";
	    print UAM "00        second\n";
	    print UAM "\n";
	    print UAM "#CFL\n";
	    print UAM "0.70	 percentage of maximum allowable time-step to take\n";
	    print UAM "\n";
	    print UAM "#LIMITER\n";
	    print UAM "mc	only limiter available\n";
	    print UAM "1.75     Beta=1.6 seems to be more stable than 2.0\n";
	    print UAM "\n";
	    print UAM "#STATISTICALMODELSONLY\n";
	    print UAM "F	if you want to run with msis and iri only (i.e. not GITM)\n";
	    print UAM "1800.0	  time step to take if you run with msis and iri\n";
	    print UAM "\n";
	    print UAM "#LOGFILE\n";
	    print UAM "1.0		dt for output to a log file\n";
	    print UAM "\n";
	    print UAM "#SAVEPLOTS\n";
	    print UAM "7200.0		dt for writing restart files\n";
	    print UAM "1		how many output files do you want\n";
	    print UAM "3DALL		output style\n";
	    print UAM "$dtout		dt for output\n";
	    print UAM "\n";
	    print UAM "#ELECTRODYNAMICS\n";
	    print UAM "60.0		how often to update potential\n";
	    print UAM "60.0		how often to update aurora and euv\n";
	    print UAM "\n";
	    print UAM "#KP\n";
	    print UAM "1.0\n";
	    print UAM "\n";
	    print UAM "#ALTITUDE\n";
	    print UAM "100.0		minimum altitude to use\n";
	    print UAM "600.0		maximum altitude to use (ignored unless the following is F)\n";
	    print UAM "T		use stretched grid\n";
	    print UAM "\n";
	    print UAM "#INITIAL\n";
	    print UAM "T		initialize thermosphere using MSIS\n";
	    print UAM "T		initialize ionosphere using IRI\n";
	    print UAM "100.0		if msis is false, then this is the temperature at the base\n";
	    print UAM "1000.0		if msis is false, then this is the temperature at the top\n";
	    print UAM "5.0e17		if msis is false, then this is the N(species1) at the base\n";
	    print UAM "7.0e18		if msis is false, then this is the N(species2) at the base\n";
	    print UAM "3.0e19		if msis is false, then this is the N(species3) at the base\n";
	    print UAM "\n";
	    print UAM "#APEX\n";
	    print UAM "T		Use apex magnetic coordinate system\n";
	    print UAM "\n";
	    
	    print UAM "#AMIEFILES\n";
	    print UAM "$amien      northern hemisphere amie file\n";
	    print UAM "$amies      southern hemisphere amie file\n";
	    print UAM "\n";

	    if ($dynamo) {

		if ($gswm) {
		    print UAM "#TIDES\n";
		    print UAM "F               UseMSISFlat\n";
		    print UAM "F               UseMSISTides\n";
		    print UAM "T               UseGSWMTides\n";
		    print UAM "F               UseWACCMTides\n";
		    print UAM "\n";
		    print UAM "#GSWMCOMP\n";
		    print UAM "T               Diurnal Migrating\n";
		    print UAM "F               Diurnal NonMigrating\n";
		    print UAM "T               Semidiurnal Migrating\n";
		    print UAM "F               Semidiurnal NonMigrating\n";
		} else {
		    print UAM "#TIDES\n";
		    print UAM "F               UseMSISFlat\n";
		    print UAM "T               UseMSISTides\n";
		    print UAM "F               UseGSWMTides\n";
		    print UAM "F               UseWACCMTides\n";
		}

		print UAM "\n";
		print UAM "#DYNAMO\n";
		print UAM "T\n";
		print UAM "45.0            Latitude to start dynamo\n";
		print UAM "500             iterations to use for the solve\n";
		print UAM "1.0             minimum residual for the solver\n";
		print UAM "\n";

	    }

	    if ($f107) {
		print UAM "#F107\n";
		print UAM "$f107		f10.7\n";
		print UAM "$f107		f10.7 averaged over 81 days\n";
	    }
	    if (!$f107set) {
		print UAM "\n";
		print UAM "#NGDC_INDICES\n";
		print UAM "DataIn/f107.txt\n";
	    }
	    print UAM "\n";

	    if ($ovation) {
		print UAM "\n";
		print UAM "#OVATIONSME\n";
		print UAM "T         usediffuse\n";
		print UAM "T         usemono\n";
		print UAM "T         usewave\n";
		print UAM "F         useions (doesn't work yet)\n";
		print UAM "\n";

		$aefile = $maindir."/AE/sme_$year.dat";
		$onsetfile = $maindir."/AE/onsets_$year.dat";

		$command = "cp $aefile .";
		print $command."\n";
		system $command."\n";
		$command = "cp $onsetfile .";
		print $command;
		system $command;

		print UAM "#SME_INDICES\n";
		print UAM "sme_$year.dat\n";
		print UAM "onsets_$year.dat\n";
		print UAM "\n";
	
	    }

	    if ($hpi) {
		print UAM "#HPI\n";
		print UAM "$hpi		hemispheric power\n";
	    } else {
		if ($power) {
		    print UAM "#NOAAHPI_INDICES\n";
		    if ($noymd) {
			if ($yearmonth) {
			    print UAM "power$ym.dat\n";
			} else {
			    print UAM "power.dat\n";
			}
			push(@tarlist,"power.dat") if ($tar);
		    } else {
			print UAM "power$ymd.dat\n";
			push(@tarlist,"power$ymd.dat") if ($tar);
		    }
		}
	    }
	    print UAM "\n";

	    if ($by || $bz || $vx) { 
		
		print UAM "#SOLARWIND\n";
		print UAM "0.0		IMF Bx\n";
		print UAM "$by		IMF By\n";
		print UAM "$bz		IMF Bz\n";
		print UAM "$vx		Solar wind Vx\n";
		print UAM "\n";
	    }

	    if ($imf) {
		print UAM "#MHD_INDICES\n";
		print UAM "$imf\n";
		print UAM "\n";
		push(@tarlist,$imf) if ($tar);
	    }

	    if ($newell) {
		print UAM "#NEWELLAURORA\n";
		print UAM "T          UseNewellAurora\n";
		print UAM "T          UseNewellAveraged\n";
		print UAM "T          UseNewellMono\n";
		print UAM "T          UseNewellWave\n";
		print UAM "F          UseNewellRemoveSpikes\n";
		print UAM "T          DoNewellAverage\n";
		print UAM "\n";
	    }

	    if ($euv) {
		if (length($euv) == 1) {
		    if ($noymd) {
			if ($yearmonth) {
			    $fismfile = "fismflux$ym.dat";
			} else {
			    $fismfile = "fismflux.dat";
			}
		    } else {
			$fismfile = "fismflux$year$month$day.dat";
		    }
		} else {
		    $fismfile=$euv;
		}
		print UAM "#EUV_DATA\n";
		print UAM "T\n";
		print UAM $fismfile."\n";
		print UAM "\n";
		push(@tarlist,$fismfile) if ($tar);
	    }

#-------------------------------
# add some satellites
#-------------------------------

	    $nsats = 0;
	    if (($champ) && (!$noChamp)) {
		$nsats = 1;
	    }
	    if (($grace) && (!$noGrace)) {
		$nsats = $nsats + 1;
	    }
	    if ($guvi) {
		$nsats = $nsats + 1;
	    }
	    if ($fyu) {
		$nsats = $nsats + 5;
	    }
	    if ($pkr) {
		$nsats = $nsats + 5;
	    }
    
	    if ($pk1) { $nsats = $nsats + 1; }
	    if ($stf) { $nsats = $nsats + 1; }
	    if ($msh) { $nsats = $nsats + 1; }
	    if ($jic) { $nsats = $nsats + 1; }
	    if ($eis) { $nsats = $nsats + 1; }
    
	    if ($sats) {
		@sats = split(',',$sats);
		$nsats = $nsats + $#sats + 1;
	    }

	    if ($nsats) {
		print UAM "#SATELLITES\n";
		print UAM "$nsats\n";
		if (($champ) && (!$noChamp)) {
		    print UAM "$champFile\n";
		    print UAM "$dtsat\n";
		    push(@tarlist,$champFile) if ($tar);
		}
		if (($grace) && (!$noGrace)) {
		    print UAM "$graceFile\n";
		    print UAM "$dtsat\n";
		    push(@tarlist,$graceFile) if ($tar);
		}
		if ($guvi) {
		    print UAM "$guvi\n";
		    print UAM "$dtsat\n";
		    push(@tarlist,$guvi) if ($tar);
		}
		writeFPI(UAM,"pkr") if ($pkr);
		writeFPI(UAM,"fyu") if ($fyu);
		writeSingle(UAM,"pkr") if ($pk1);
		writeSingle(UAM,"stf") if ($stf);
		writeSingle(UAM,"msh") if ($msh);
		writeSingle(UAM,"jic") if ($jic);
		writeSingle(UAM,"eis") if ($eis);

		while ($file=pop(@sats)) {
		
		    if ($buildsat) {

			$outfile = $file.".$year$month$day";
		
			open(INFILE, "<$file");
			$lon = <INFILE>; $lon =~ s/\n//;
			$lat = <INFILE>; $lat =~ s/\n//;
			$alt = <INFILE>; $alt =~ s/\n//;
			close(INFILE);

			open(OUTFILE, ">$outfile");

			print OUTFILE "\n";
			print OUTFILE "#START\n";
			print OUTFILE "$year $month $day 0 0 0 0 $lon $lat $alt\n";
			print OUTFILE "$year $month $nextday 0 0 0 0 $lon $lat $alt\n";

			close(OUTFILE);
		
			$file = $outfile;

		    }

		    print UAM "$file\n";
		    print UAM "$dtsat\n";
		    push(@tarlist,$file) if ($tar);
	    
		}

		print UAM "\n";
	    }
    
	    print UAM "#THERMO\n";
	    print UAM "T		 UseSolarHeating\n";
	    print UAM "T		 UseJouleHeating\n";
	    print UAM "T		 UseAuroralHeating\n";
	    print UAM "T		 UseNOCooling\n";
	    print UAM "T		 UseOCooling\n";
	    print UAM "T		 UseConduction\n";
	    print UAM "T		 UseTurbulentConduction\n";
	    print UAM "F		 UseUpdatedTurbulentConduction\n";
	    print UAM "1.0		 EddyScalingFactor\n";
	    print UAM "\n";
	    print UAM "#FORCING\n";
	    print UAM "T		UsePressureGradient\n";
	    print UAM "T		UseIonDrag\n";
	    print UAM "T		UseNeutralDrag\n";
	    print UAM "T		UseViscosity\n";
	    print UAM "T		UseCoriolis\n";
	    print UAM "T		UseGravity\n";
	    print UAM "\n";
	    print UAM "#IONFORCING\n";
	    print UAM "T\n";
	    print UAM "T\n";
	    print UAM "T\n";
	    print UAM "T\n";
	    if ($dynamo) {
		print UAM "T        dynamo is on!\n";
		print UAM "\n";
		#print UAM "#STRETCH\n";
		#print UAM "0.0	! location of minimum grid spacing\n";
		#print UAM "0.7	! Amount of stretch 0 (none) to 1 (lots)\n";
		#print UAM "0.8 	! More control of stretch ( > 1 stretch less < 1 stretch more)\n";
		print UAM "\n";
	    } else {
		print UAM "F        dynamo is off!\n";
		print UAM "\n";
		print UAM "#STRETCH\n";
		print UAM "65.0 ! location of minimum grid spacing\n";
		print UAM "0.0  ! Amount of stretch 0 (none) to 1 (lots)\n";
		print UAM "1.0	! More control of stretch ( > 1 stretch less < 1 stretch more)\n";
		print UAM "\n";
	    }
	    print UAM "\n";
	    print UAM "#CHEMISTRY\n";
	    print UAM "T		UseIonChemistry\n";
	    print UAM "T		UseIonAdvection\n";
	    print UAM "T		UseNeutralChemistry\n";
	    print UAM "\n";
	    print UAM "#END\n";
	    print UAM "\n";
	    
	    close(UAM);

	}

    }

    if ($amie) {

	$ntimes = 86400.0/$dtout;

	$file = "AMIE.$ymd.in";
	print "Creating file $file\n";
	open(AMIE,">$file");

	print AMIE "\n";
	print AMIE "#DEBUG\n";
	print AMIE "0\n";
	print AMIE "\n";
	print AMIE "#STARTTIME\n";
	print AMIE "$year\n";
	print AMIE "$month\n";
	print AMIE "$day\n";
	print AMIE "00\n";
	print AMIE "00\n";
	print AMIE "00\n";
	print AMIE "\n";
	print AMIE "#DT\n";
	print AMIE "$dtout\n";
	print AMIE "0         - NEED TO CHANGE THIS when we get e-field data!\n";
	print AMIE "\n";
	print AMIE "#NTIMES\n";
	print AMIE "$ntimes\n";
	print AMIE "\n";
	print AMIE "#NORTH\n";
	if (!$south) {
	    print AMIE "T\n";
	} else {
	    print AMIE "F\n";
	}
	print AMIE "\n";
	print AMIE "#F107\n";
	print AMIE "$f107\n";
	print AMIE "\n";
	print AMIE "#BACKGROUND\n";
	print AMIE "$datadir\n";
	print AMIE "weimer96\n";
	print AMIE "ihp\n";
	print AMIE "kroehl\n";
	print AMIE "\n";
#	print AMIE "#NGDC_INDICES\n";
#	print AMIE "AMIE-$ymd.dat\n";
#	print AMIE "0      - this is zero because Eric put a delay in the data time stamp!!\n";
#	print AMIE "\n";
	print AMIE "#NOAA_INDICES\n";
	print AMIE "power.dat\n";
	print AMIE "\n";
	print AMIE "#MHD_INDICES\n";
	print AMIE "1\n";
	print AMIE "imf.dat\n";
	print AMIE "\n";
	print AMIE "#MAGNETOMETER\n";
	print AMIE "./mag$ymd.final\n";
	print AMIE "./master$ymd.txt\n";
	print AMIE "T\n";
	print AMIE "T\n";
	if (!$south) {
	    print AMIE "1.0\n";
	    print AMIE "1.0\n";
	    print AMIE "1.0\n";
	} else {
	    print AMIE "0.25\n";
	    print AMIE "0.25\n";
	    print AMIE "0.25\n";
	}
	print AMIE "\n";
	print AMIE "#OUTPUT\n";
	print AMIE "b$ymd.wm\n";
	print AMIE ".true.\n";
	print AMIE ".true.\n";
	print AMIE ".false.\n";
	print AMIE "\n";
	print AMIE "#END\n";
	print AMIE "\n";

	close(AMIE);

	print RUNAMIE "rm -f AMIE.in ; ln -s AMIE.$ymd.in AMIE.in\n";
	print RUNAMIE "./amie\n";
	print RUNAMIE "\n";

    }

}

if ($amie) {
    close(RUNAMIE);
    $command = "chmod a+x runamie.sh\n";
    system $command;
}


if ($noymd) {
    # There is a very good chance that we are doing multiple days and 
    # would like to concatinate the files....
    if ($rm) {
	$command = "$bindir/cat_sat_input.pl -rm";
    } else {
	$command = "$bindir/cat_sat_input.pl";
    }
    print "Executing : $command\n";
    system $command;

    if (($BuildPower) && ($yearmonth)) {
	$command = "mv power.dat power$ym.dat";
	print "Executing : $command\n";
	system $command;
    }

    if (($BuildImf) && ($yearmonth)) {
	$command = "mv imf.dat imf$ym.dat";
	print "Executing : $command\n";
	system $command;
    }

    if (($BuildEUV) && ($yearmonth)) {
	$command = "mv fismflux.dat fismflux$ym.dat";
	print "Executing : $command\n";
	system $command;
    }

    if (($BuildChamp) && ($yearmonth)) {
	$command = "mv champ.dat champ$ym.dat";
	print "Executing : $command\n";
	system $command;
    }

    if (($BuildGrace) && ($yearmonth)) {
	$command = "mv grace.dat grace$ym.dat";
	print "Executing : $command\n";
	system $command;
    }

}

if ($imf) {
    open(IDLFILE,">.idl_input");
    print IDLFILE ".r imf_plot\n";
    print IDLFILE "imf.dat\n";
    print IDLFILE "n\n";
    close(IDLFILE);
    $command = "idl < .idl_input";
    printf "Executing: $command\n";
    system $command;
}

if ($power) {
    open(IDLFILE,">.idl_input");
    print IDLFILE ".r power_plot\n";
    print IDLFILE "power.dat\n";
    print IDLFILE "power.ps\n";
    close(IDLFILE);
    $command = "idl < .idl_input";
    printf "Executing: $command\n";
    system $command;
}

exit(1);


if ($mhdsats) {

    print "Creating SAT files\n";

    open(IDLFILE,">.idl_sat_input");

    # cluster
    print IDLFILE ".r cdf_to_mhd\n";
    print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/istp/cluster/cl/sp/aux/$year/cl_sp_aux_".$ymd."_v*.cdf\n";
    print IDLFILE "s\n";
    print IDLFILE "cluster_$ymd.dat\n\n";
    
    # polar
    print IDLFILE ".r cdf_to_mhd\n";
    print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/istp/polar/mfe/$year/po_k0_mfe_".$ymd."_v*.cdf\n";
    print IDLFILE "s\n";
    print IDLFILE "polar_$ymd.dat\n\n";

    # goes 8
    print IDLFILE ".r cdf_to_mhd\n";
    print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/istp/goes/8_mag/$year/g8_k0_mag_".$ymd."_v*.cdf\n";
    print IDLFILE "s\n";
    print IDLFILE "goes08_$ymd.dat\n\n";

    # goes 9
    print IDLFILE ".r cdf_to_mhd\n";
    print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/istp/goes/9_mag/$year/g9_k0_mag_".$ymd."_v*.cdf\n";
    print IDLFILE "s\n";
    print IDLFILE "goes09_$ymd.dat\n\n";

    # goes 10
    print IDLFILE ".r cdf_to_mhd\n";
    print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/istp/goes/0_mag/$year/g0_k0_mag_".$ymd."_v*.cdf\n";
    print IDLFILE "s\n";
    print IDLFILE "goes10_$ymd.dat\n\n";

    # goes 11
    print IDLFILE ".r cdf_to_mhd\n";
    print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/istp/goes/11_mag/$year/goes11_k0_mag_".$ymd."_v*.cdf\n";
    print IDLFILE "s\n";
    print IDLFILE "goes11_$ymd.dat\n\n";

    # goes 12
    print IDLFILE ".r cdf_to_mhd\n";
    print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/istp/goes/12_mag/$year/goes12_k0_mag_".$ymd."_v*.cdf\n";
    print IDLFILE "s\n";
    print IDLFILE "goes12_$ymd.dat\n\n";

    if (!$wind) {
	# wind
	print IDLFILE ".r cdf_to_mhd\n";
	print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/istp/wind/mfi_h0/$year/wi_h0_mfi_".$ymd."_v*.cdf\n";
	print IDLFILE "s\n";
	print IDLFILE "wind_$ymd.dat\n\n";
    }

    # geotail
    print IDLFILE ".r cdf_to_mhd\n";
    print IDLFILE $maindir."cdaweb.gsfc.nasa.gov/pub/istp/geotail/mgf/$year/ge_k0_mgf_".$ymd."_v*.cdf\n";
    print IDLFILE "s\n";
    print IDLFILE "geotail_$ymd.dat\n\n";

    # LANL01
    print IDLFILE "cdf_file  = '".$maindir."cdaweb.gsfc.nasa.gov/pub/istp/lanl/01a_mpa/$year/a1_k0_mpa_".$ymd."*.cdf'\n";
    print IDLFILE "sat_name  = 'LANL-01'\n";
    print IDLFILE "file_name = 'lanl01_$ymd.dat'\n";
    print IDLFILE "orbitmaker_lanl_cdf, cdf_file, sat_name=sat_name, file_name=file_name\n";
    print IDLFILE "retall\n\n";

    # LANL02
    print IDLFILE "cdf_file  = '".$maindir."cdaweb.gsfc.nasa.gov/pub/istp/lanl/02a_mpa/$year/a2_k0_mpa_".$ymd."*.cdf'\n";
    print IDLFILE "sat_name  = 'LANL-02'\n";
    print IDLFILE "file_name = 'lanl02_$ymd.dat'\n";
    print IDLFILE "orbitmaker_lanl_cdf, cdf_file, sat_name=sat_name, file_name=file_name\n";
    print IDLFILE "retall\n\n";

    # LANL89
    print IDLFILE "cdf_file  = '".$maindir."cdaweb.gsfc.nasa.gov/pub/istp/lanl/89_mpa/$year/l9_k0_mpa_".$ymd."*.cdf'\n";
    print IDLFILE "sat_name  = 'LANL-89'\n";
    print IDLFILE "file_name = 'lanl89_$ymd.dat'\n";
    print IDLFILE "orbitmaker_lanl_cdf, cdf_file, sat_name=sat_name, file_name=file_name\n";
    print IDLFILE "retall\n\n";

    # LANL90
    print IDLFILE "cdf_file  = '".$maindir."cdaweb.gsfc.nasa.gov/pub/istp/lanl/90_mpa/$year/l0_k0_mpa_".$ymd."*.cdf'\n";
    print IDLFILE "sat_name  = 'LANL-90'\n";
    print IDLFILE "file_name = 'lanl90_$ymd.dat'\n";
    print IDLFILE "orbitmaker_lanl_cdf, cdf_file, sat_name=sat_name, file_name=file_name\n";
    print IDLFILE "retall\n\n";

    # LANL91
    print IDLFILE "cdf_file  = '".$maindir."cdaweb.gsfc.nasa.gov/pub/istp/lanl/91_mpa/$year/l1_k0_mpa_".$ymd."*.cdf'\n";
    print IDLFILE "sat_name  = 'LANL-91'\n";
    print IDLFILE "file_name = 'lanl91_$ymd.dat'\n";
    print IDLFILE "orbitmaker_lanl_cdf, cdf_file, sat_name=sat_name, file_name=file_name\n";
    print IDLFILE "retall\n\n";

    # LANL94
    print IDLFILE "cdf_file  = '".$maindir."cdaweb.gsfc.nasa.gov/pub/istp/lanl/94_mpa/$year/l4_k0_mpa_".$ymd."*.cdf'\n";
    print IDLFILE "sat_name  = 'LANL-94'\n";
    print IDLFILE "file_name = 'lanl94_$ymd.dat'\n";
    print IDLFILE "orbitmaker_lanl_cdf, cdf_file, sat_name=sat_name, file_name=file_name\n";
    print IDLFILE "retall\n\n";

    # LANL97
    print IDLFILE "cdf_file  = '".$maindir."cdaweb.gsfc.nasa.gov/pub/istp/lanl/97_mpa/$year/l7_k0_mpa_".$ymd."*.cdf'\n";
    print IDLFILE "sat_name  = 'LANL-97'\n";
    print IDLFILE "file_name = 'lanl97_$ymd.dat'\n";
    print IDLFILE "orbitmaker_lanl_cdf, cdf_file, sat_name=sat_name, file_name=file_name\n";
    print IDLFILE "retall\n\n";

    close(IDLFILE);

    $command = "idl < .idl_sat_input";
    system $command;

    @satfiles = glob("*".$ymd."*.dat");

    if ($#satfiles > 0) {
	$isat = 1;
	$isattotal = 0;
	while ($isattotal < $#satfiles) {
	    if (!($satfiles[$isattotal] =~ /imf/)) {
		$isat++;
	    }
	    $isattotal++;
	}

	if (!(-d "sat_files_$ymd") && !$dontmove) {
	    print "sat_files directory doesn't exist - creating\n";
	    $command = "mkdir sat_files_$ymd";
	    system $command;
	}

	print "A total of $isat satellite files were found\n";
	if ($noymd) {
	    open(SATFILE,">PARAM.sats");
	} else {
	    open(SATFILE,">PARAM.sats.$ymd");
	}
	print SATFILE "\n";
	print SATFILE "#SATELLITE\n";
	print SATFILE "$isat\n";
	while (@satfiles) {
	    $satfile = shift @satfiles;
	    if (!($satfile =~ /imf/)) {
		print SATFILE "MHD file RAY\n";
		print SATFILE "-1\n";
		print SATFILE "5.0\n";
		if ($dontmove) {
		    print SATFILE "$satfile\n";
		} else {
		    print SATFILE "sat_files_$ymd/$satfile\n";
		    $command = "mv $satfile sat_files_$ymd";
		    system $command;
		}
	    }
	}
	print SATFILE "\n";
	close(SATFILE);

    }
}

if (length($guvi) == 1) {

    $guvi = "guvi_$year$month$day.dat";

}

if ($swmf) {

    if ($noymd) {
	if ($restart =~ /T/) {
	    open(PARAM,">PARAM.in.Restart");
	} else {
	    open(PARAM,">PARAM.in.Start");
	}
    } else {
	open(PARAM,">PARAM.in.$year$month$day");
    }

    if ($tar) {
	$tarfile = "param.".$year.$month.$day.".tgz";
	@tarlist = "PARAM.in.".$year.$month.$day;
    }

    print PARAM "#ECHO\n";
    print PARAM "T\n";
    print PARAM "\n";
    print PARAM "#DESCRIPTION\n";
    print PARAM "$ymd Run\n";
    print PARAM "\n";
    if ($restart=~/T/) {

	print PARAM "#INCLUDE\n";
	print PARAM "RESTART.in\n";
	print PARAM "\n";

    } else {

	print PARAM "#TIMEACCURATE\n";
	print PARAM "F			DoTimeAccurate\n";
	print PARAM "\n";
	print PARAM "#STARTTIME\n";
	print PARAM "$year		year\n";
	print PARAM "$month		month\n";
	print PARAM "$day		day\n";
	print PARAM "00		hour\n";
	print PARAM "00		minute\n";
	print PARAM "00		second\n";
	print PARAM "0.0 		FracSecond\n";
	print PARAM "\n";
	print PARAM "#COMPONENT\n";
	print PARAM "IM                      NameComp\n";
	print PARAM "F                       UseComp\n";
	print PARAM "\n";
	print PARAM "#COUPLE2\n";
	print PARAM "GM			NameComp1\n";
	print PARAM "IE			NameComp2\n";
	print PARAM "10                      DnCouple\n";
	print PARAM "-1.0                     DtCouple\n";
	print PARAM "\n";
	print PARAM "#BEGIN_COMP GM ---------------------------------------------------------------\n";
	print PARAM "\n";
	print PARAM "#INCLUDE\n";
	print PARAM "PARAM.gm.1st\n";
	print PARAM "\n";
	print PARAM "#INCLUDE\n";
	if ($noymd) {
	    print PARAM "PARAM.sats\n";
	} else {
	    print PARAM "PARAM.sats.$ymd\n";
	}
	print PARAM "\n";
	print PARAM "#UPSTREAM_INPUT_FILE\n";
	print PARAM "T                       UseUpstreamInputFile\n";
	print PARAM "$imf                    UpstreamFileName\n";
	print PARAM "0.0                     Satellite_Y_Pos\n";
	print PARAM "0.0                     Satellite_Z_Pos\n";
	print PARAM "\n";
	print PARAM "#END_COMP GM -----------------------------------------------------------------\n";
	print PARAM "\n";
	print PARAM "#BEGIN_COMP IE ---------------------------------------------------------------\n";
	print PARAM "\n";
	print PARAM "#IONOSPHERE\n";
	print PARAM "5                       TypeConductanceModel\n";
	print PARAM "F                       UseFullCurrent\n";
	print PARAM "F	 		UseFakeRegion2\n";
	print PARAM "$f107                   F107Flux\n";
	print PARAM "1.0                     StarLightPedConductance\n";
	print PARAM "0.25                    PolarCapPedConductance\n";
	print PARAM "\n";
	print PARAM "#SAVEPLOT\n";
	print PARAM "1			nPlotFile\n";
	print PARAM "min idl                 StringPlot\n";
	print PARAM "100                     DnSavePlot\n";
	print PARAM "-1.0			DtSavePlot\n";
	print PARAM "\n";
	print PARAM "#SPS\n";
	print PARAM "T\n";
	print PARAM "\n";
	print PARAM "#DEBUG\n";
	print PARAM "0\n";
	print PARAM "0\n";
	print PARAM "\n";
	print PARAM "#END_COMP IE -----------------------------------------------------------------\n";
	print PARAM "\n";
	print PARAM "#STOP\n";
	print PARAM "700                     MaxIter\n";
	print PARAM "-1.			TimeMax\n";
	print PARAM "\n";
	print PARAM "#RUN	!!!!!!!!!!!!!!!!!!!!!!\n";
	print PARAM "\n";
	print PARAM "#BEGIN_COMP GM ---------------------------------------------------------------\n";
	print PARAM "\n";
	print PARAM "#INCLUDE\n";
	print PARAM "PARAM.gm.2nd\n";
	print PARAM "\n";
	print PARAM "#END_COMP GM -----------------------------------------------------------------\n";
	print PARAM "\n";
	print PARAM "#STOP\n";
	print PARAM "1500                    MaxIter\n";
	print PARAM "-1.			TimeMax\n";
	print PARAM "\n";
	print PARAM "#RUN	!!!!!!!!!!!!!!!!!!!!!\n";
	print PARAM "\n";

    }

    print PARAM "#INCLUDE\n";
    print PARAM "PARAM.sw.ta\n";
    print PARAM "\n";
    print PARAM "#BEGIN_COMP GM ---------------------------------------------------------------\n";
    print PARAM "\n";

    if ($restart=~/T/) {

	print PARAM "#INCLUDE\n";
	print PARAM "GM/restartIN/restart.H\n";
	print PARAM "\n";
    }

    print PARAM "#INCLUDE\n";
    print PARAM "PARAM.gm.ta\n";
    print PARAM "\n";

    if ($restart=~/T/) {

	print PARAM "#OUTERBOUNDARY\n";
	print PARAM "outflow			TypeBC_I(east_)\n";
	print PARAM "vary			TypeBC_I(west_)\n";
	print PARAM "float			TypeBC_I(south_)\n";
	print PARAM "float			TypeBC_I(north_)\n";
	print PARAM "float			TypeBC_I(bot_)\n";
	print PARAM "float			TypeBC_I(top_)\n";
	print PARAM "\n";
	print PARAM "#INNERBOUNDARY\n";
	print PARAM "ionosphere		\n";
	print PARAM "\n";
	print PARAM "#INCLUDE\n";
	if ($noymd) {
	    print PARAM "PARAM.sats\n";
	} else {
	    print PARAM "PARAM.sats.$ymd\n";
	}
	print PARAM "\n";
	print PARAM "#UPSTREAM_INPUT_FILE\n";
	print PARAM "T                       UseUpstreamInputFile\n";
	print PARAM "$imf                    UpstreamFileName\n";
	print PARAM "0.0                     Satellite_Y_Pos\n";
	print PARAM "0.0                     Satellite_Z_Pos\n";
	print PARAM "\n";

    }

    print PARAM "#END_COMP GM -----------------------------------------------------------------\n";
    print PARAM "\n";
    print PARAM "#BEGIN_COMP IM ---------------------------------------------------------------\n";
    print PARAM "\n";
    print PARAM "#RESTART\n";
    print PARAM "$restart\n";
    print PARAM "\n";
    print PARAM "#SAVEPLOT\n";
    print PARAM "1\n";
    print PARAM "3d max tec\n";
    print PARAM "-1\n";
    print PARAM "1800.0\n";
    print PARAM "\n";
    print PARAM "#END_COMP IM -----------------------------------------------------------------\n";
    print PARAM "\n";
    print PARAM "#BEGIN_COMP IE ---------------------------------------------------------------\n";
    print PARAM "\n";
    print PARAM "#IONOSPHERE\n";
    print PARAM "5                       TypeConductanceModel\n";
    print PARAM "F                       UseFullCurrent\n";
    print PARAM "F			UseFakeRegion2\n";
    print PARAM "$f107                   F107Flux\n";
    print PARAM "1.0                     StarLightPedConductance\n";
    print PARAM "0.25                    PolarCapPedConductance\n";
    print PARAM "\n";
    print PARAM "#SAVEPLOT\n";
    print PARAM "1			nPlotFile\n";
    print PARAM "aur idl                 StringPlot\n";
    print PARAM "-1                      DnSavePlot\n";
    print PARAM "60.0			DtSavePlot\n";
    print PARAM "\n";
    print PARAM "#SPS\n";
    print PARAM "T\n";
    print PARAM "\n";
    print PARAM "#DEBUG\n";
    print PARAM "0\n";
    print PARAM "0\n";
    print PARAM "\n";
    print PARAM "#END_COMP IE -----------------------------------------------------------------\n";
    print PARAM "\n";
    print PARAM "#STOP\n";
    print PARAM "-1			MaxIter\n";
    print PARAM "86400.0			TimeMax\n";
    print PARAM "\n";
    print PARAM "#END\n";

    close(PARAM);

}

if ($tar) {

    $command = "tar -cvzf $tarfile ";
    
    while ($file=pop(@tarlist)) {
	$command = $command.$file." ";
    }

    print $command,"\n";
    system $command;

}


exit(1);

#===============================================================

sub getf107{

    my $file = "DataIn/f107.txt";

    open(F107, $file);

    my $yyyy = shift(@_);
    my $mm   = shift(@_)+0;
    my $dd   = shift(@_)+0;

    if ($mm < 10) {
	$mm = "0".$mm;
    }

    if ($dd < 10) {
	$dd = "0".$dd;
    }

    while(<F107>) {
	if (/$yyyy.$mm.$dd.\d\d.\d\d.(\d\d...)/ ||
	    /$yyyy.$mm.$dd.\d\d.\d\d..(\d\d...)/ ||
	    /$yyyy.$mm.$dd.\d\d.\d\d...(\d\d...)/) {
	    $f107 = $1;
	}
    }

    close(F107);
    return int($f107);

}


sub makeFPIfiles{

    my $name = shift(@_);
    my $lon  = shift(@_)+0.0;
    my $lat  = shift(@_)+0.0;
    my $alt  = shift(@_)+0.0;

    my $latn = $lat + $alt/111.0;
    my $lats = $lat - $alt/111.0;

    my $lone = $lon + $alt/(111.0*cos($lat*3.1415/180.0));
    my $lonw = $lon - $alt/(111.0*cos($lat*3.1415/180.0));

    print "$name$ymd.dat\n";

    # Base station
    open(FPI, ">$name$ymd.dat");
    print FPI "\n";
    print FPI "#START\n";
    print FPI "$year $month $day 0 0 0 0 $lon $lat $alt\n";
    print FPI "$year $month $nextday 0 0 0 0 $lon $lat $alt\n";
    close(FPI);

    # North station
    open(FPI, ">$name"."n$ymd.dat");
    print FPI "\n";
    print FPI "#START\n";
    print FPI "$year $month $day 0 0 0 0 $lon $latn $alt\n";
    print FPI "$year $month $nextday 0 0 0 0 $lon $latn $alt\n";
    close(FPI);

    # South station
    open(FPI, ">$name"."s$ymd.dat");
    print FPI "\n";
    print FPI "#START\n";
    print FPI "$year $month $day 0 0 0 0 $lon $lats $alt\n";
    print FPI "$year $month $nextday 0 0 0 0 $lon $lats $alt\n";
    close(FPI);

    # East station
    open(FPI, ">$name"."e$ymd.dat");
    print FPI "\n";
    print FPI "#START\n";
    print FPI "$year $month $day 0 0 0 0 $lone $lat $alt\n";
    print FPI "$year $month $nextday 0 0 0 0 $lone $lat $alt\n";
    close(FPI);

    # West station
    open(FPI, ">$name"."w$ymd.dat");
    print FPI "\n";
    print FPI "#START\n";
    print FPI "$year $month $day 0 0 0 0 $lonw $lat $alt\n";
    print FPI "$year $month $nextday 0 0 0 0 $lonw $lat $alt\n";
    close(FPI);

}

sub makeSingleFile{

    my $name = shift(@_);
    my $lon  = shift(@_)+0.0;
    my $lat  = shift(@_)+0.0;
    my $alt  = shift(@_)+0.0;

    print "$name.dat\n";

    # Base station
    open(FPI, ">$name.dat");
    print FPI "\n";
    print FPI "#START\n";
    print FPI "$startyear $startmonth $startday 0 0 0 0 $lon $lat $alt\n";
    print FPI "$endyear $endmonth $nextday 0 0 0 0 $lon $lat $alt\n";
    close(FPI);

}

sub writeFPI{

    my $file = shift(@_);
    my $base = shift(@_);

    $f = $base;
    $f = $f.$ymd if (!$noymd);
    $f = $f.".dat";
    print $file $f."\n";
    print $file "$dtsat\n";
    push(@tarlist,$f) if ($tar);

    $f = $base."n";
    $f = $f.$ymd if (!$noymd);
    $f = $f.".dat";
    print $file $f."\n";
    print $file "$dtsat\n";
    push(@tarlist,$f) if ($tar);

    $f = $base."s";
    $f = $f.$ymd if (!$noymd);
    $f = $f.".dat";
    print $file $f."\n";
    print $file "$dtsat\n";
    push(@tarlist,$f) if ($tar);

    $f = $base."e";
    $f = $f.$ymd if (!$noymd);
    $f = $f.".dat";
    print $file $f."\n";
    print $file "$dtsat\n";
    push(@tarlist,$f) if ($tar);

    $f = $base."w";
    $f = $f.$ymd if (!$noymd);
    $f = $f.".dat";
    print $file $f."\n";
    print $file "$dtsat\n";
    push(@tarlist,$f) if ($tar);

}

sub writeSingle{

    my $file = shift(@_);
    my $base = shift(@_);

    $f = $base;
    $f = $f.$ymd if (!$noymd);
    $f = $f.".dat";
    print $file $f."\n";
    print $file "$dtsat\n";
    push(@tarlist,$f) if ($tar);

}

#************************************************************************
#****   Return 1 if we are after the leap day in a leap year.       *****
#************************************************************************

sub leapDay                  
{                            

    if ($_[0] % 4 != 0) {
	return(0);
    }

    if (!($_[0] % 100)) {             # years that are multiples of 100
				      # are not leap years
	if ($_[0] % 400) {            # unless they are multiples of 400
	    return(0);
	}
    }

    if ($_[1] < 2) {
	return(0);
    } elsif (($_[1] == 2) && ($_[1] < 29)) {
	return(0);
    } else {
	return(1);
    }
}

sub jday_to_md {

    $y = $_[0];
    $j = $_[1];

    $m = 0;
    for ($mc=0;$mc<12;$mc++) {
	$m = $mc if ($j > $theJulianDate[$mc]);
    }
    $d = $j - $theJulianDate[$m];
    $m++;

    return ($m,$d);

}
jday_to_md($year,$jday);
