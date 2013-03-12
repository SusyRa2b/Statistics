#!/user/bin/perl
use strict;
use warnings;

open(my $fin,  "<", "$ARGV[0]") or die "Can't open input file.";
open(my $fout, ">", "saveGenerate.dat") or die "Can't open output file.";

while(<$fin>)
  {

    if(/(\S+)(,)(\s)(\S+)(\s+)(met,ht,nbjet bin \()(\d)(,)(\d)(,)(\d)(\))(.+)(npass=)(\s*)(\S+)(\s)(\+\/-)(\s*)(\S+)/){

      my $region = $1;
      my $sample = $4;
      my $met = $7;
      my $ht = $9;
      my $btag = $11;
      my $val = $16;
      my $err = $20;
      #print $fout "$region $sample $met $ht $btag $val $err\n";
      $met += 1;
      $ht += 1;
      $btag += 1;
      print $fout "$region"."_$sample"."_M$met"."_H$ht"."_$btag"."b $val\n";

    }
    elsif(/(trigeff)(\S+)(\s+)(:)(\s+)(\S+)/)
      {
	print $fout "$1$2 $6\n";
      }


  }
