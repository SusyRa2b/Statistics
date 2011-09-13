#include "Riostream.h"

void runcls_ln_ge1btight_t1bbbb_7_ra2b(int mgl, int mlsp) {

  gROOT->LoadMacro("ra2bRoostatsClass7.c+") ;

  ra2bRoostatsClass7 ra2b(1,0,2) ; 
  
  ra2b.initialize("an-11-257-v3-files/input-files/byhand-data-ge1b-tight.txt",
                  "an-11-257-v3-files/input-files/signalSyst.T1bbbb.ge1bTight.dat", mgl, mlsp, true, 1.0 ) ;
  
  gROOT->LoadMacro("RA2bHypoTestInvDemo.c+") ;

  stringstream Mgl; Mgl << mgl ;
  stringstream Mlsp; Mlsp << mlsp ;

  TString fname = "my_output_file_";
  fname += Mgl.str();
  fname += "_";
  fname += Mlsp.str();
  fname += ".root";

  RA2bHypoTestInvDemo ("ws.root", "ws", "SbModel", "BModel", 
		       "ra2b_observed_rds", 0, 3, true, 9, 5, 37, 1000,
		       mgl, mlsp, fname);

}



