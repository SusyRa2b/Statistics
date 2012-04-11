
void runLimitWideSB_withKristen(char* cut , char* model , int m0, int m12, bool isMeasured, bool isHybrid, bool doToys = true)
{
  //gSystem->CompileMacro("RooBetaPdf.cxx","kO") ;
  //gSystem->CompileMacro("RooGammaPdf.cxx","kO") ;
  //gSystem->CompileMacro("RooBetaPdfWithPoissonGenerator.cxx","kO") ;
  gSystem->CompileMacro("DebuggingToyMCSampler.cxx","kO") ;
  gSystem->CompileMacro("HybridToyMCSampler.cxx","kO") ;
  gSystem->CompileMacro("ProfileLikelihoodTestStat_New.cxx","kO") ;
  gSystem->CompileMacro("profileLikelihoodLimit.C","kO") ;
  gSystem->CompileMacro("ra2bRoostatsClass7_noNSig_wideSideband_withKristen_moreObservables.c","kO") ;
  TString outputfile("/tmp/ra2b/ws_");
  if(isMeasured) outputfile+="measured_";
  else outputfile+="expected_";
  if(isHybrid) outputfile+="hybrid_";
  else outputfile+="frequentist_";
  outputfile+=cut;
  outputfile+="_";
  outputfile+=model;
  outputfile+="_";
  outputfile+=m0;
  outputfile+="_";
  outputfile+=m12;
  outputfile+=".root";
  ra2bRoostatsClass7 ra ;
  TString inputfile("input-files/byhand-data-");
  if(TString(cut).Contains("ge1bLoose")) inputfile+="1BL";
  if(TString(cut).Contains("ge1bTight")) inputfile+="1BT";
  if(TString(cut).Contains("ge2bLoose")) inputfile+="2BL";
  if(TString(cut).Contains("ge2bTight")) inputfile+="2BT";
  if(TString(cut).Contains("ge3bLoose")) inputfile+="3B";
  if(!isMeasured) inputfile+="-expected";
  inputfile+="-tc.txt";
  TString systematicsfile("/afs/cern.ch/user/j/joshmt/public/RA2bFall2011/signalSyst.");
  systematicsfile+=model;
  systematicsfile+=".";
  systematicsfile+=cut;
  systematicsfile+="WideSB.dat";
  ra.initialize(inputfile,systematicsfile,m0,m12,false,1.0,outputfile,false,false,false);
  TString name(model);
  name+="_";
  name+=cut;
  name+="_kristen";
  if(doToys){
    if(isMeasured) {
      if(isHybrid) profileLikelihoodLimit(outputfile.Data(),"ws","SbModel", "BModel","ra2b_observed_rds",name.Data(),true,
					  m0,m12,1000,10000,10,false,true);
      else profileLikelihoodLimit(outputfile.Data(),"ws","SbModel", "BModel","ra2b_observed_rds",name.Data(),true,
				  m0,m12,1000,10000,10,true,false);
    }
    else {
      if(isHybrid) profileLikelihoodLimit(outputfile.Data(),"ws","SbModel", "BModel","ra2b_observed_rds",name.Data(),false,
					  m0,m12,1000,10000,10,false,true);
      else profileLikelihoodLimit(outputfile.Data(),"ws","SbModel", "BModel","ra2b_observed_rds",name.Data(),false,
				  m0,m12,1000,10000,10,true,false);
    }
  }
  else {
    if(isMeasured) profileLikelihoodLimit(outputfile.Data(),"ws","SbModel", "BModel","ra2b_observed_rds",name.Data(),true,
					  m0,m12,1000,10000,10,false,false);
    else profileLikelihoodLimit(outputfile.Data(),"ws","SbModel", "BModel","ra2b_observed_rds",name.Data(),false,
				m0,m12,1000,10000,10,false,false);
  }
}
