
void make_newmcinp_t1bbbb_ws_3B( double mgl=600., double mlsp=400. ) {

  gROOT->LoadMacro("ra2bRoostatsClass7.c+") ;

  ra2bRoostatsClass7 ra2b(1,0) ;

  ra2b.initialize("paper-2011-files/input-files/likelihood-input-from-mc-HT400-SIGMET250-nb3.txt",
                  "paper-2011-files/input-files/signalSyst.T1bbbb.ge3bLoose.dat", mgl, mlsp, true, 1.0 ) ;

  gSystem->Exec("mkdir -p output-files") ;



  char wsfilename[10000] ;
  sprintf( wsfilename, "output-files/ws-newmcinp-t1bbbb-3B-mgl%.0f-mlsp%.0f.root", mgl, mlsp ) ;

  char command[10000] ;
  sprintf( command, "mv %s %s-old", wsfilename, wsfilename ) ;
  gSystem->Exec( command ) ;
  sprintf( command, "mv ws.root %s", wsfilename ) ;
  gSystem->Exec( command ) ;


}



