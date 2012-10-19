#!/user/bin/perl
use strict;
use warnings;

opendir my($dh), $ARGV[0] or die  "Couldn't open dir '$ARGV[0]': $!";
my @list = grep {!/^\./ }readdir $dh;
closedir $dh;

open(my $fout, ">", "$ARGV[0]/$ARGV[0].dat") or die "Couldn't open output : $!";

foreach(@list) {
  
  my $file = $ARGV[0]."/".$_;
  
  print "Translating $file...\n";
  system("perl translate.pl $file");
  
  
  if($file =~ /(\S*)(.dat)/) {
    #system("cp signalModelFile.dat $1/sig1/setupSignalNominal.dat");
    system("cp setupSignalNominal.dat $1/sig1/setupSignalNominal.dat");
    print $fout "$1\n";
  }
}

close $fout;
