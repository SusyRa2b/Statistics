
void make_lm9_ws_1BT_sbvars() {

  gROOT->LoadMacro("ra2bRoostatsClass7.c+") ;

  ra2bRoostatsClass7 ra2b(0,0) ;

  ra2b.initialize("paper-2011-files/input-files/byhand-data-1BT.txt",
                  "paper-2011-files/input-files/signalSyst.LM9.ge1bTight.dat", 0., 0., false, 1.0 ) ;

  gSystem->Exec("mkdir -p output-files") ;
  gSystem->Exec("mv output-files/ws-lm9-1BT-sbvars.root output-files/ws-lm9-1BT-sbvars.root-old") ;
  gSystem->Exec("mv ws.root output-files/ws-lm9-1BT-sbvars.root") ;

}



