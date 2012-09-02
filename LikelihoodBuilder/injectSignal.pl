#!/user/bin/perl
use strict;
use warnings;

my $inputCrossSection = $ARGV[0];
my $m0 = 850;
my $m12 = 300;


my $dir="";
$dir = "injected\_$m0\_$m12\_$inputCrossSection";
mkdir $dir;

my $lumi = 0;
#get the luminosity
open(my $setupfin,  "<", "setupFile.dat") or die "Can't open signalModelFile.";
while(<$setupfin>) {
  if(/(Luminosity)(\s*)(\S*)/) {
    $lumi = $3;
  }
}
print "Luminosity = $lumi\n";

my $signalLine = "";
open(my $signalfin,  "<", "signalModelFile.dat") or die "Can't open signalModelFile.";

while(<$signalfin>) {
  $signalLine = $_;
  if($signalLine =~ /(\d*)(\s*)(\d*)/) {
    if( $1 eq $m0 && $3 eq $m12 ) { last; }
  }
}

#print "$signalLine\n";

my @signalLineValues = split(' ', $signalLine);

my $numberOfEvents = $signalLineValues[2];
my $nBins = scalar(@signalLineValues);
$nBins -= 3;
$nBins = $nBins/6;

print "m0 = $m0, m12 = $m12, Ntot = $numberOfEvents, Nbins = $nBins \n";

open(my $binsfin,  "<", "binFilesFile.dat") or die "Can't open binFilesFile.";

my $bins = "";
while(<$binsfin>) {
  if(/(\S*)(\s*)(\S*)/) {
    $bins = $bins . " $3";
  }
}

my @nBinFilesValues = split(' ', $bins);
my $nBinFilesValuesSize = scalar(@nBinFilesValues);

if($nBins ne $nBinFilesValuesSize) { die "Inconsistency in number of bins!"; }

my $binNumber = -1;

my $weight = $inputCrossSection*$lumi/$numberOfEvents;
print "weight = $weight\n";

foreach(@nBinFilesValues) {
  
  $binNumber += 1;
  
  my $thisFileName = $_;
  
  print "bin $binNumber, $thisFileName\n";
  
  open(my $thisFile, "<", "originalInputs/".$thisFileName) or die "Can't open $thisFileName.";
  open(my $thisFileInjected, ">", "$dir/$thisFileName") or die "Can't open $thisFileName.";

  while(<$thisFile>) {

    my $thisLine = $_;
    my $num = 0;

    if(/(zeroLeptonCount)(\s*)(\S*)/) {
      $num = $3;
      $num += @signalLineValues[3+$binNumber*6]*$weight;
      print $thisFileInjected "$1$2$num\n";
    } elsif (/(oneMuonCount)(\s*)(\S*)/) {
      $num = $3;
      $num += @signalLineValues[3+$binNumber*6+4]*$weight;
      print $thisFileInjected "$1$2$num\n";
    } elsif (/(oneElectronCount)(\s*)(\S*)/) {
      $num = $3;
      $num += @signalLineValues[3+$binNumber*6+4]*$weight;
      print $thisFileInjected "$1$2$num\n";
    } elsif (/(zeroLeptonLowDeltaPhiNCount)(\s*)(\S*)/) {
      $num = $3;
      $num += @signalLineValues[3+$binNumber*6+2]*$weight;
      print $thisFileInjected "$1$2$num\n";
    } else {
      print $thisFileInjected "$thisLine";
    }

  }


  close $thisFile;
  close $thisFileInjected;

}


#remember to close files

