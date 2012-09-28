#!/user/bin/perl
use strict;
use warnings;

print "WARNING: xsec is hardcoded!\n";

open(my $fin,  "<",  "$ARGV[0]") or  die "Can't open input file.";
open(my $fout,  ">", "signalModelFile.dat") or die "Can't open output file.";

my $nHTbins;
my $nMETbins;
my $nBTAGbins = 3;

my $sigtot = 0;

while(<$fin>) {
  if(/(Using HT bins:)(\s*)(\S*)(\s*)(Using MET bins:)(\s*)(\S*)/) {
      my $HTbit = $3;
      my $METbit = $7;
      
      $nHTbins = $HTbit =~ tr/-//;
      $nMETbins = $METbit =~ tr/-//;
      
      print "Number of HT bins: $nHTbins\n"; 
      print "Number of MET bins: $nMETbins\n";
      print "Number of BTAG bins: $nBTAGbins\n";
    }

}
close $fin;



open($fin,  "<",  "$ARGV[0]") or  die "Can't open  input file.";
while(<$fin>) {
  
  if(/(Using HT bins:)(\s*)(\S*)(\s*)(Using MET bins:)(\s*)(\S*)/) {
    #print $fout "$_";
  }
  else{
    
    my $line = $_;
    my @linea = split(" ",$line);
    my $lineasize = scalar @linea;
    
    $lineasize = $lineasize - 3;
    
    my $numNumbers = 9*$nHTbins*$nMETbins*$nBTAGbins;
    if($lineasize != $numNumbers) { die "Inconsistency!"; }
    
    print $fout "$linea[0] $linea[1] $linea[2] ";
    
    for(my $i=0;  $i<= $nHTbins*$nMETbins*$nBTAGbins; $i++) {
      #IN:  sig sl  ldp
      #OUT: sig ldp sl
      
      #my $xsec = 4.41988000000000034e-02;#850 600 7TeV
      my $xsec = 9.66802999999999968e-02; #850 600 8TeV
      my $weight = 1.0/0.5/$xsec;

      my $sig = int($linea[3*$i+3]*$weight+0.5);
      $sigtot += $sig; 
      #my $sigunrounded = $linea[3*$i+3]*$weight;
      #print "$sigunrounded -> $sig\n";
      
      my $sigr = sqrt($sig);
      my $sl = int($linea[3*$i+4]*$weight+0.5);
      my $slr = sqrt($sl);
      my $ldp = int($linea[3*$i+5]*$weight+0.5);
      my $ldpr = sqrt($ldp); 

      print $fout "$sig $sigr $ldp $ldpr $sl $slr ";
      
    }
    print $fout "\n";
    
    print "sigtot = $sigtot\n";

  }
}

close $fin;
