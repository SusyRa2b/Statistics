#/user/bin/perl
use strict;
use warnings;

open (my $fin, "<", "listMR.txt") or die "Can't open input file.";
open (my $fout, ">", "saveInputMR.dat") or die "Can't open input file.";


while(<$fin>)
  {

    if(/(\S+)/)
      {
	my $fileName = $1;

	if(/(\d)(\d)(\d)/)
	  {

	    my $ht = $1;
	    my $met = $2;
	    my $btag = $3;
	    
	    $ht+=1;
	    $met+=1;


	    #get trigger efficiency
	    my $trigeff = 0;
	    open (my $ft, "<", "$ARGV[0]/bin_M$met"."_H$ht"."_$btag"."b.dat") or die "can't open bin file.";
	    while(<$ft>)
	      {
		if(/(oneLeptonTriggerEfficiency)(\s+)(\S+)/)
		  {
		    $trigeff = $3;
		  }
	      }
	    close $ft;
	    print "$trigeff\n";

	    open (my $fc, "<", "$fileName");
	    
	    my $line = 0;
	    while(<$fc>)
	      {
		if(/(\S+)/)
		  {
		    my $value = $1/$trigeff;

		    my $binName = "M$met"."_H$ht"."_$btag"."b";

		    $line++;
		    if($line==1){ print $fout "N_MRzeroLepton_TopWJets_$binName $value\n"; }
		    if($line==2){ print $fout "N_oneTightMu_Theta1_TopWJets_$binName $value\n"; }
		    if($line==3){ print $fout "N_oneLooseLep_Theta1_TopWJets_$binName $value\n"; }
		    if($line==4){ print $fout "N_oneTightMu_Theta2_TopWJets_$binName $value\n"; }
		    if($line==5){ print $fout "N_oneLooseLep_Theta2_TopWJets_$binName $value\n"; }
		    if($line==6){ print $fout "N_oneTightMu_Theta3_TopWJets_$binName $value\n"; }
		    if($line==7){ print $fout "N_oneLooseLep_Theta3_TopWJets_$binName $value\n"; }
		    if($line==8){ print $fout "N_oneTightMu_Theta4_TopWJets_$binName $value\n"; }
		    if($line==9){ print $fout "N_oneLooseLep_Theta4_TopWJets_$binName $value\n"; }
		    if($line==10){ print $fout "N_oneTightMu_Theta5_TopWJets_$binName $value\n"; }
		    if($line==11){ print $fout "N_oneLooseLep_Theta5_TopWJets_$binName $value\n"; }
		    if($line==12){ print $fout "N_twoTightMu_TopWJets_$binName $value\n"; }
		    if($line==13){ print $fout "N_twoLooseLep_TopWJets_$binName $value\n"; }
		  }
	      }	
	  }
      }
  }





