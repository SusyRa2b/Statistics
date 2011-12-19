
void make_lm9_ws_3B_sbvars() {

  gROOT->LoadMacro("ra2bRoostatsClass7.c+") ;

  ra2bRoostatsClass7 ra2b(0,0) ;

  ra2b.initialize("paper-2011-files/input-files/byhand-data-3B.txt",
                  "paper-2011-files/input-files/signalSyst.LM9.ge3bLoose.dat", 0., 0., false, 1.0 ) ;

  gSystem->Exec("mkdir -p output-files") ;
  gSystem->Exec("mv output-files/ws-lm9-3B-sbvars.root output-files/ws-lm9-3B-sbvars.root-old") ;
  gSystem->Exec("mv ws.root output-files/ws-lm9-3B-sbvars.root") ;

}



