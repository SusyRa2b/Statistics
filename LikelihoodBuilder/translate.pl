#!/user/bin/perl
use strict;
use warnings;

my $dir="";
if($ARGV[0] =~ /(\S*)(.dat)/) {
  $dir = $1;
  mkdir $dir;
  mkdir "$dir/sig1";
}
else{
  die "Error making directory.";
}

open(my $fin,  "<", "$ARGV[0]") or die "Can't open input file.";
open(my $globalout,  ">", "$dir/setupFile.dat") or die "Can't open setupFile file.";
open(my $binlistout,  ">", "$dir/binFilesFile.dat") or die "Can't open binFilesFile file.";

my $nHTbins;
my $nMETbins;
my $nBTAGbins = 3;

my @binFileNames;
my @binFileHandles;
my %binFileHash;

my @ldpMCSum;

my $linecount = 0;
while(<$fin>) {
  my $line = $_;
  $linecount = $linecount+1;
  
  #with first line, get binning and make output files
  if($linecount eq 1) {
    if($line =~ /(Using HT bins:)(\s*)(\S*)(\s*)(Using MET bins:)(\s*)(\S*)/) {
      my $HTbit = $3;
      my $METbit = $7;
      
      $nHTbins = $HTbit =~ tr/-//;
      $nMETbins = $METbit =~ tr/-//;
      
      #print "Number of HT bins: $nHTbins\n"; 
      #print "Number of MET bins: $nMETbins\n";
      #print "Number of BTAG bins: $nBTAGbins\n";
    }
    else {die "Can't understand first line.";}
    
    #make array of bin names, and more
    for(my $i=1; $i<=$nMETbins; $i++) {
      for(my $j=1; $j<=$nHTbins; $j++) {
	for(my $k=1; $k<=$nBTAGbins; $k++) {
	  my $tempBinFile = "M$i\_H$j\_$k";
	  $tempBinFile = $tempBinFile."b";
	  push(@binFileNames, $tempBinFile);
	  push(@ldpMCSum,0);
	}
      }
    }
    
    #make output files
    my $tempCount=0;
    foreach(@binFileNames) {
      #print "Opening file for $_\n";
      my $bfn = $_;
      my $handle;
      open($handle, ">", "$dir/bin_$bfn.dat");
      push @binFileHandles, $handle;
      $binFileHash{$bfn} = $tempCount;
      $tempCount = $tempCount+1;
      print $binlistout "bin$tempCount bin_$bfn.dat\n";
      #print $binlistout "bin_$bfn bin_$bfn.dat\n";
    }
    
    #FOLLOWING LINES ARE FOR TESTING
    #print hash -- e.g. M4_H1_1b key gives 36 value
    #while ((my $key, my $value) = each %binFileHash)
    #  {
    #	print "$key key gives $binFileHash{$key} value\n";
    #  }
    #print {$binFileHandles[$binFileHash{"M4_H1_1b"}]} "testing";
    
  }#end of first line check
  else{
    
    if($line =~ /(\S*)(\s*)(\S*)/) {
      my $fullOKAname = $1;
      my $value = $3;
      
      if($fullOKAname =~ /(N_0lep)(_M)(\d)(_H)(\d)(_)(\d)(b)/) {
      	my $bin = "M$3_H$5_$7"."b";
	print {$binFileHandles[$binFileHash{$bin}]} "zeroLeptonCount $value\n";
      }
      
      elsif($fullOKAname =~ /(N_1lep)(_M)(\d)(_H)(\d)(_)(\d)(b)/) {
      	my $bin = "M$3_H$5_$7"."b";
	my $half = $value/2.0;
	print {$binFileHandles[$binFileHash{$bin}]} "oneLeptonCount $value\n";
	print {$binFileHandles[$binFileHash{$bin}]} "oneMuonCount $half\n";
	print {$binFileHandles[$binFileHash{$bin}]} "oneElectronCount $half\n";
      }
      
      elsif($fullOKAname =~ /(N_ldp)(_M)(\d)(_H)(\d)(_)(\d)(b)/) {
      	my $bin = "M$3_H$5_$7"."b";
	print {$binFileHandles[$binFileHash{$bin}]} "zeroLeptonLowDeltaPhiNCount $value\n";
      }
      
      elsif($fullOKAname =~ /(R_lsb)(_H)(\d)(_)(1b)(_*)(\S*)/) {#Only take 1b since it's independent of b-tagging!
	my $dim = "H$3";
	my $err = "$7";
	foreach(@binFileNames) {
	  my $binFileName = $_;
	  if($binFileName =~ /$dim/) {
	    if($err ne "err") {
	      print {$binFileHandles[$binFileHash{$binFileName}]} "deltaPhiNRatioName $dim\n";
	      print {$binFileHandles[$binFileHash{$binFileName}]} "deltaPhiNRatio $value\n";
	    }
	    else{
	      print {$binFileHandles[$binFileHash{$binFileName}]} "deltaPhiNRatioError $value\n";
	      print {$binFileHandles[$binFileHash{$binFileName}]} "lowDeltaPhiNScalingName $dim\n";
	    }
	  }
	}
      }

      elsif($fullOKAname =~ /(ttwj_mc_ldpover0lep_ratio)(_M)(\d)(_H)(\d)(_)(\d)(b)(_*)(\S*)/) {
	my $binFileName = "M$3_H$5_$7b";
	my $err = $10;
	if($err ne "err") {
	  print {$binFileHandles[$binFileHash{$binFileName}]} "topWJetsLowDeltaPhiNOverZeroLeptonRatioMC $value\n";
	}
	else {
	  #currently not used
	  #print {$binFileHandles[$binFileHash{$binFileName}]} "topWJetsLowDeltaPhiNOverZeroLeptonRatioMCError $value\n";
	}
      }

      elsif($fullOKAname =~ /(znn_mc_ldpover0lep_ratio)(_M)(\d)(_H)(\d)(_)(\d)(b)(_*)(\S*)/) {
	my $binFileName = "M$3_H$5_$7b";
	my $err = $10;
	if($err ne "err") {
	  print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMC $value\n";
	}
	else {
	  #currently not used
	  #print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoNuNuLowDeltaPhiNOverZeroLeptonRatioMCError $value\n";
	}
      }
      
      elsif($fullOKAname =~ /(N_Zee)(_M)(\d)(_H)(\d)/) {
	my $dim = "M$3_H$5";
	foreach(@binFileNames) {
	  my $binFileName = $_;
	  if($binFileName =~ /$dim/) {
	    print {$binFileHandles[$binFileHash{$binFileName}]} "diElectronCountName $dim\n";
	    print {$binFileHandles[$binFileHash{$binFileName}]} "diElectronCount $value\n";
	  }
	}
      }
      
      elsif($fullOKAname =~ /(N_Zmm)(_M)(\d)(_H)(\d)/) {
	my $dim = "M$3_H$5";
	foreach(@binFileNames) {
	  my $binFileName = $_;
	  if($binFileName =~ /$dim/) {
	    print {$binFileHandles[$binFileHash{$binFileName}]} "diMuonCountName $dim\n";
	    print {$binFileHandles[$binFileHash{$binFileName}]} "diMuonCount $value\n";
	  }
	}
      }
      
      elsif($fullOKAname =~ /(N_ttbarsingletopzjetsmc_ldp_)(\S*)/) {
	my $binFileName = $2;
	$ldpMCSum[$binFileHash{$binFileName}] += $value;
      }
      elsif($fullOKAname =~ /(N_WJmc_ldp_)(\S*)/) {
        my $binFileName = $2;
        $ldpMCSum[$binFileHash{$binFileName}] += $value;
      }
      elsif($fullOKAname =~ /(N_Znnmc_ldp_)(\S*)/) {
	my $binFileName = $2;
	$ldpMCSum[$binFileHash{$binFileName}] += $value;
      }
      
      elsif($fullOKAname =~ /(acc_Zee)(_M)(\d)(_*)(\S*)/) {
	my $dim = "M$3";
	my $err = "$5";
	foreach(@binFileNames) {
	  my $binFileName = $_;
	  if($binFileName =~ /$dim/) {
	    if($err ne "err") {
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoeeAcceptanceName $dim\n";
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoeeAcceptance $value\n";
	    }
	    else{
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoeeAcceptanceError $value\n";
	    }
	  }
	}
      }
      
      elsif($fullOKAname =~ /(acc_Zmm)(_M)(\d)(_*)(\S*)/) {
	my $dim = "M$3";
	my $err = "$5";
	foreach(@binFileNames) {
	  my $binFileName = $_;
	  if($binFileName =~ /$dim/) {
	    if($err ne "err") {
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtomumuAcceptanceName $dim\n";
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtomumuAcceptance $value\n";
	    }
	    else{
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtomumuAcceptanceError $value\n";
	    }
	  }
	}
      }
      
      elsif($fullOKAname =~ /(Z_ee_eff_err)/) {
	print $globalout "ZtoeeEfficiencyError $value\n";
      }
      elsif($fullOKAname =~ /(Z_ee_eff)/) {#order important
	print $globalout "ZtoeeEfficiency $value\n";
      }
      
      elsif($fullOKAname =~ /(Z_mm_eff_err)/) {
	print $globalout "ZtomumuEfficiencyError $value\n";
      }
      elsif($fullOKAname =~ /(Z_mm_eff)/) {#order important
	print $globalout "ZtomumuEfficiency $value\n";
      }
      
      elsif($fullOKAname =~ /(knn_)(\d)(b)(_*)(\S*)/) {
	my $dim = "$2b";
	my $err = "$5";
	foreach(@binFileNames) {
	  my $binFileName = $_;
	  if($binFileName =~ /$dim/) {
	    if($err ne "err") {
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoNuNubTagScalingName $dim\n";
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoNuNubTagScaling $value\n";
	    }
	    else{
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoNuNubTagScalingError $value\n";
	    }
	  }
	}
      }
      
      elsif($fullOKAname =~ /(Z_ee_pur_err)/) {
	print $globalout "ZtoeePurityError $value\n";
      }
      elsif($fullOKAname =~ /(Z_ee_pur)/) {#order important
        print $globalout "ZtoeePurity $value\n";
      }
      
      elsif($fullOKAname =~ /(Z_mm_pur_err)/) {
        print $globalout "ZtomumuPurityError $value\n";
      }
      elsif($fullOKAname =~ /(Z_mm_pur)/) {#order important
        print $globalout "ZtomumuPurity $value\n";
      }
      
      elsif($fullOKAname =~ /(sf_mc_err)/) {
	print $globalout "MCUncertainty $value\n";
      }
      
      elsif($fullOKAname =~ /(sf_qcd)(_M)(\d)(_H)(\d)(_)(\d)(b)(_*)(\S*)/) {
	my $binFileName = "M$3_H$5_$7b";
	my $err = $10;
	if($err ne "err") {
	  print {$binFileHandles[$binFileHash{$binFileName}]} "qcdClosure $value\n";
	}
	else {
	  print {$binFileHandles[$binFileHash{$binFileName}]} "qcdClosureError $value\n";
	}
      }
      
      elsif($fullOKAname =~ /(sf_ttwj)(_M)(\d)(_H)(\d)(_)(\d)(b)(_*)(\S*)/) {
	my $binFileName = "M$3_H$5_$7b";
	my $err = $10;
	if($err ne "err") {
	  print {$binFileHandles[$binFileHash{$binFileName}]} "topWJetsClosure $value\n";
	}
	else {
	  print {$binFileHandles[$binFileHash{$binFileName}]} "topWJetsClosureError $value\n";
	}
      }
      
      # elsif($fullOKAname =~ /(sf_ee_)(\d)(b)(_*)(\S*)/) {
      # 	my $dim = "$2b";
      # 	my $err = "$5";
      # 	foreach(@binFileNames) {
      # 	  my $binFileName = $_;
      # 	  if($binFileName =~ /$dim/) {
      # 	    if($err ne "err") {
      # 	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoeeSystematicName $dim\n";
      # 	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoeeSystematic $value\n";
      # 	    }
      # 	    else{
      # 	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoeeSystematicError $value\n";
      # 	    }
      # 	  }
      # 	}
      # }
      
      # elsif($fullOKAname =~ /(sf_mm_)(\d)(b)(_*)(\S*)/) {
      # 	my $dim = "$2b";
      # 	my $err = "$5";
      # 	foreach(@binFileNames) {
      # 	  my $binFileName = $_;
      # 	  if($binFileName =~ /$dim/) {
      # 	    if($err ne "err") {
      # 	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtomumuSystematicName $dim\n";
      # 	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtomumuSystematic $value\n";
      # 	    }
      # 	    else{
      # 	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtomumuSystematicError $value\n";
      # 	    }
      # 	  }
      # 	}
      # }

      #BEN FIXME -- temporary until we get the right systematics from OAK
      elsif($fullOKAname =~ /(sf_ee_)(1)(b)(_*)(\S*)/) {
	my $dim = "$2b";
	my $err = "$5";
	foreach(@binFileNames) {
	  my $binFileName = $_;
	  #if($binFileName =~ /$dim/) {
	    if($err ne "err") {
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoeeSystematicName $dim\n";
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoeeSystematic $value\n";
	    }
	    else{
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtoeeSystematicError $value\n";
	    }
	  #}
	}
      }
      
      elsif($fullOKAname =~ /(sf_mm_)(1)(b)(_*)(\S*)/) {
	my $dim = "$2b";
	my $err = "$5";
	foreach(@binFileNames) {
	  my $binFileName = $_;
	  #if($binFileName =~ /$dim/) {
	    if($err ne "err") {
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtomumuSystematicName $dim\n";
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtomumuSystematic $value\n";
	    }
	    else{
	      print {$binFileHandles[$binFileHash{$binFileName}]} "ZtomumuSystematicError $value\n";
	    }
	  #}
	}
      }
      
      else {
	print "Doing nothing with line: $line";
      }
      
    }#line to translate (Name Value)
   
    else {
      print "Doing nothing with line: $line";
    }
    
  }#not the first line
}

#print the values that are sums
foreach(@binFileNames) {
  my $binFileName = $_;
  my $value = $ldpMCSum[$binFileHash{$binFileName}];
  print {$binFileHandles[$binFileHash{$binFileName}]} "zeroLeptonLowDeltaPhiNMC $value\n";
}


#Efficiencies, not in OKA's code yet
foreach(@binFileNames) {
  my $binFileName = $_;
  print {$binFileHandles[$binFileHash{$binFileName}]} "zeroLeptonTriggerEfficiency 1.0\n";
  print {$binFileHandles[$binFileHash{$binFileName}]} "zeroLeptonTriggerEfficiencyError 0.01\n";
  print {$binFileHandles[$binFileHash{$binFileName}]} "zeroLeptonLowDeltaPhiNTriggerEfficiency 1.0\n";
  print {$binFileHandles[$binFileHash{$binFileName}]} "zeroLeptonLowDeltaPhiNTriggerEfficiencyError 0.01\n";
  print {$binFileHandles[$binFileHash{$binFileName}]} "oneLeptonTriggerEfficiency 1.0\n";
  print {$binFileHandles[$binFileHash{$binFileName}]} "oneLeptonTriggerEfficiencyError 0.01\n";
  print {$binFileHandles[$binFileHash{$binFileName}]} "oneElectronTriggerEfficiency 1.0\n";
  print {$binFileHandles[$binFileHash{$binFileName}]} "oneElectronTriggerEfficiencyError 0.01\n";
  print {$binFileHandles[$binFileHash{$binFileName}]} "oneMuonTriggerEfficiency 1.0\n";
  print {$binFileHandles[$binFileHash{$binFileName}]} "oneMuonTriggerEfficiencyError 0.01\n";
}

print $globalout "ZtollOverZtoNuNuRatio 0.168067227\n";
print $globalout "Luminosity 15.0\n";


close $globalout;
close $binlistout;
foreach(@binFileHandles){
  close $_;
}

print "Translated files are in directory $dir.\n";
