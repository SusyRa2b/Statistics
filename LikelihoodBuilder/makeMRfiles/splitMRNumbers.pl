#!/user/bin/perl
use strict;
use warnings;

#ugly code -- assumes structure is 
#nameA
#number_1
#...
#number_nPerFile
#nameB
#...

open(my $fin,  "<", "$ARGV[0]") or die "Can't open input file.";

my $nPerFile = 13;

my $btag = 0;
my $ht = 0;
my $met = 0;

my $line = 0;
my $fout;

my $firstFile = 1;

while(<$fin>) {
  
  if(($line%($nPerFile+1)) == 0) {
    if(/(\d)(B)(\d)(HT)(\d)(MET)/) {
      $btag = $1;
      $ht = $3;
      $met = $5;
      print "$btag $ht $met\n";
      if($firstFile != 1) {
	close $fout;
	$firstFile = 0;
      }
      open($fout, ">", "smmctest$ht$met$btag.txt");
    }
    else {die "expected new name.\n"; }
  }
  elsif(/(\S+)/) {
    print $fout "$1\n";
  }
  
  $line += 1;
  
}
