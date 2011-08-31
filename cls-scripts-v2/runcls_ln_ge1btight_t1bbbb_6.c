

void runcls_ln_ge1btight_t1bbbb_6() {

  gROOT->LoadMacro("ra2bRoostatsClass6.c+") ;

  ra2bRoostatsClass6 ra2b(1,0,2) ; 
  
  ra2b.initialize("an-11-257-v3-files/input-files/byhand-data-ge1b-tight.txt") ;
  ra2b.setSusyScanPoint("an-11-257-v3-files/input-files/signalSyst.T1bbbb.ge1bTight.dat", 1000, 800, true, 1.0 ) ;
  
  gROOT->LoadMacro("StandardHypoTestInvDemo.C+") ;
  
  StandardHypoTestInvDemo ("ws.root", "ws", "SbModel", "BModel", 
			   "ra2b_observed_rds", 0, 3, true, 4, 20, 50, 1000);

}



