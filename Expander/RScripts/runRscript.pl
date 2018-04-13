
use strict;

my $rExeFile = $ARGV[0];    # R exe file location
my $scriptFile = $ARGV[1];  # R script to be run by this perl script
my $dummy = $ARGV[2];
my $envDir = $ARGV[3];
my $logFile = $ARGV[4];

my $debug = 1;

### Check usage:
$#ARGV == 4 or die "Usage: runRscript.pl  R_location  R_script_path DUMMY File_dir log_file_path\n";

my $ret = system("$rExeFile --no-save $dummy $envDir < $scriptFile > $logFile");

#my $ret = system("$rExeFile --no-save < $scriptFile");

#my $ret = system("$rExeFile CMD BATCH $scriptFile");

#my $ret = system("$rExeFile  --no-save  < $scriptFile > $logFile");

if ($debug){
	#print "return value: $ret\n";
}

#exit $ret;

###################################################################################
# An example call for this script:
# perl runRscript.pl  "E:/seagull/R-2.8.1/bin/R.exe" "tryArgs.R" "dummy" "rLog.txt"
###################################################################################
