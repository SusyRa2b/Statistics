
void make_lm9_ws_2BT() {

  gROOT->LoadMacro("ra2bRoostatsClass7.c+") ;

  ra2bRoostatsClass7 ra2b(1,0) ;

  ra2b.initialize("paper-2011-files/input-files/byhand-data-2BT.txt",
                  "paper-2011-files/input-files/signalSyst.LM9.ge2bTight.dat", 0., 0., false, 1.0 ) ;

  gSystem->Exec("mkdir -p output-files") ;
  gSystem->Exec("mv output-files/ws-lm9-2BT.root output-files/ws-lm9-2BT.root-old") ;
  gSystem->Exec("mv ws.root output-files/ws-lm9-2BT.root") ;

}



