#!/user/bin/perl
use strict;
use warnings;
use File::Path;
use Cwd;

#my $workdir = "/afs/cern.ch/work/k/kreis/likelihood/UserCode/SusyAnalysis/RA2b/Statistics/toys_8TeV/LikelihoodBuilder";
my $workdir = cwd();

my $dir = $ARGV[0];

open( my $flist, "<", "$dir/$dir.dat") or die "Can't open input file!";

mkpath($dir.'/submission_scripts');
mkpath($dir.'/log_files');
mkpath($dir.'/dat_files');
mkpath($dir.'/root_files');

my $subDirectory = $dir."/submission_scripts";
my $datDirectory = $dir."/dat_files";
my $rootDirectory = $dir."/root_files";
my $tmpDirectory = "/tmp/ra2b_kreis";

while(<$flist>) {
  
  if(/($dir)(\/)(\S+)/) {
    
    my $inputDirectory = $3;

    open(my $fsub, ">", "$subDirectory/submit_$inputDirectory.sh")  or die  "Can't open output file!";
    
    print $fsub "#! /bin/bash -f\n";
    print $fsub "export BATCH_DEBUG=true\n";
    print $fsub "if [ ! -d \"$tmpDirectory\" ]; then mkdir $tmpDirectory; fi;\n";
    print $fsub "cd $workdir\n";
    print $fsub "source /afs/cern.ch/sw/lcg/contrib/gcc/4.3/x86_64-slc5/setup.sh\n";
    print $fsub "source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.33.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh\n";
    print $fsub "root -q -l -b callRunFit.C'(\"$workdir/$dir/$inputDirectory/\",\"$tmpDirectory/\",\"".$dir."_$inputDirectory\")'\n";
    print $fsub "if [ -e \"$tmpDirectory/dat_".$dir."_$inputDirectory.dat\" ]; then mv $tmpDirectory/dat_".$dir."_$inputDirectory.dat $workdir/$datDirectory/dat_".$dir."_$inputDirectory.dat; fi\n";
    print $fsub "if [ -e \"$tmpDirectory/likelihood_".$dir."_$inputDirectory.root\" ]; then mv $tmpDirectory/likelihood_".$dir."_$inputDirectory.root $workdir/$rootDirectory/likelihood_".$dir."_$inputDirectory.root; fi\n";

    close $fsub;
    
    #system("bsub -q 1nh -N -oo $dir/log_files/log_$inputDirectory.txt < $subDirectory/submit_$inputDirectory.sh");
    system("bsub -q 8nh -R \"rusage[mem=2048]\" -M 2000000 -N -oo $dir/log_files/log_$inputDirectory.txt < $subDirectory/submit_$inputDirectory.sh");
        
  }
  
  
  open(my $fdel, ">", "$dir/clear.sh") or die "Can't open clear file";
  print $fdel "rm -r dat_files\n";  
  print $fdel "rm -r root_files\n";  
  print $fdel "rm -r log_files\n";  
  print $fdel "rm -r submission_scripts\n";  
  close $fdel;
  
}
