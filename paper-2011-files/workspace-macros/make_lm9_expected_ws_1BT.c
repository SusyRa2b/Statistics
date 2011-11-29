
void make_lm9_expected_ws_1BT() {

  gROOT->LoadMacro("ra2bRoostatsClass7.c+") ;

  ra2bRoostatsClass7 ra2b(1,0) ;

  ra2b.initialize("paper-2011-files/input-files/byhand-data-1BT-expected.txt",
                  "paper-2011-files/input-files/signalSyst.LM9.ge1bTight.dat", 0., 0., false, 1.0 ) ;

  gSystem->Exec("mkdir -p output-files") ;
  gSystem->Exec("mv output-files/expected-ws-lm9-1BT.root output-files/expected-ws-lm9-1BT.root-old") ;
  gSystem->Exec("mv ws.root output-files/expected-ws-lm9-1BT.root") ;

}



