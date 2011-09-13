#include "Riostream.h"

void runcls_ln_ge1btight_mSugra_7_ra2b(int m0, int m12, double nevt) {

  if (nevt==0) return ; 

  gROOT->LoadMacro("ra2bRoostatsClass7.c+") ;

  ra2bRoostatsClass7 ra2b(1,0,2) ; 
  
  ra2b.initialize("an-11-257-v3-files/input-files/byhand-data-ge1b-tight.txt",
                  "an-11-257-v3-files/input-files/signalSyst.mSUGRAtanb40.ge1bTight.1143invpb.dat", m0, m12, false, 1.0 ) ;
  
  gROOT->LoadMacro("RA2bHypoTestInvDemo.c+") ;

  stringstream M0; M0 << m0 ;
  stringstream M12; M12 << m12 ;

  TString fname = "partial_results/my_output_file_";
  fname += M0.str();
  fname += "_";
  fname += M12.str();
  fname += ".root";

  RA2bHypoTestInvDemo ("ws.root", "ws", "SbModel", "BModel", 
		       "ra2b_observed_rds", 0, 3, true, 1, nevt, nevt, 1000,
		       m0, m12, fname);

}



