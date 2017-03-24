#!/usr/bin/perl -s

sub catfiles {

    my $OutFile = shift(@_);
    print "OutFile : $OutFile\n";

    open(OUTPUT, ">$OutFile");

    $first = 1;
    $last = "crap\n";
    while (@_) {
	my $file = shift @_;
	print "processing : $file\n";
	open(INPUT, "<$file");
	$start = 0;
	$prestart = 0;
	while (<INPUT>) {
	    if (/START/) { 
		$start = 1;
		print OUTPUT $_ if ($first);
	    } else {
		if ($first || $start) {
		    print OUTPUT $_ if (!/$last/);
#		    print "-->".$_,"   ".$last;
		    $last = $_;
#		    print length($last),"\n";
		    if (length($last) < 5) {$last = "craa\n";}
		}
		if ($prestart == 2) { $start = 1 };
		if ($prestart == 1) { $prestart = 2 };
		if (/Normalizing/) {$prestart = 1};
	    }
	}
	$first = 0;
	close(INPUT);
	if ($rm) {
	    $command = "/bin/rm -f $file";
	    print $command."\n";
	    system $command;
	}
    }

    close(OUTPUT);

}

if ($#ARGV < 0) {

    # We have to figure out all of the files to process!

    # Let's make some assumptions:
    # 1. All files end with YYYYMMDD.dat
    # 2. We should drop the _ in the sat_yyyymmdd.dat and output sat.dat

    @AllFiles = <*.dat>;

    $IsBad = 0;
    $IsGood = 0;

    while (!$IsGood) {
	$OldFile = shift @AllFiles;
	print "there ->$OldFile<-\n";
	if ($OldFile =~ /(.*)(\d\d\d\d\d\d\d\d).dat/) {
	    $_ = $1;
	    s/_\z//;
	    $Header = $_;
	    push @FileList, $_.".dat";
	    push @FileList, $OldFile;
	    $IsGood = 1;
	}
    }

    print $Header,"\n";

    while (@AllFiles) {

	$IsGood = 0;
	$IsBad = 0;
	while (!$IsGood && !$IsBad) {
	    $File = shift @AllFiles;
	    if ($File =~ /(.*)(\d\d\d\d\d\d\d\d).dat/) {
		$IsGood = 1;
	    }
	    $IsBad = 1 if (length($File) == 0);
	}

	if ($File =~ /$Header(\d\d\d\d\d\d\d\d)/) {
	    print "Match! $File\n";
	    push @FileList, $File;
	} else {
	    catfiles(@FileList);
	    print "New Header $File!\n";
	    if ($File =~ /(.*)(\d\d\d\d\d\d\d\d).dat/) {
		$_ = $1;
		s/_\z//;
		$Header = $_;
		@FileList = $_.".dat";
		push @FileList, $File;
	    } else {
		print "File : $File doesn't fit regular expression!\n";
		exit(1) if (!IsBad);
	    }
	}

    }

    catfiles(@FileList) if (!$IsBad);

} else {

    while (@ARGV) {
	$file = shift @ARGV;
	if ($file =~/-out\w*=(\w*)/ || $file =~/-file\w*=(\w*)/) { 
	    $Output = $1; 
	}
	push @FileList, $file if (!($file =~/-/));
    }
    
    if (!$Output) {
	$Output = "out";
    }

    @FinalFileList = $Output;

    while (@FileList) {
	$file = shift @FileList;
	push @FinalFileList, $file;
    }

    catfiles(@FinalFileList);

}
