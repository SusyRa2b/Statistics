
void make_newfit_t1bbbb_ws_1BL( double mgl=600., double mlsp=400. ) {

  gROOT->LoadMacro("ra2bRoostatsClass9.c+") ;

  ra2bRoostatsClass9 ra2b(0,1) ;

  ra2b.initialize("paper-2011-files/input-files/likelihood-newfit-input-from-mc-HT400-SIGMET250.txt",
                  "paper-2011-files/input-files/likelihood-newfit-t1bbbb-systfile-HT400-SIGMET250.txt", mgl, mlsp, true, 1.0 ) ;

  char wsfilename[10000] ;
  sprintf( wsfilename, "output-files/ws-newfit-t1bbbb-1BL-mgl%.0f-mlsp%.0f.root", mgl, mlsp ) ;

  gSystem->Exec("mkdir -p output-files") ;

  char command[10000] ;
  sprintf( command, "mv %s %s-old", wsfilename, wsfilename ) ;
  gSystem->Exec( command ) ;
  sprintf( command, "mv ws.root %s", wsfilename ) ;
  gSystem->Exec( command ) ;

}



