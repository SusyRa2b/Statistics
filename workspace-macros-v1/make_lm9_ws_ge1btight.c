
void make_lm9_ws_ge1btight() {

  gROOT->LoadMacro("ra2bRoostatsClass7.c+") ;

  ra2bRoostatsClass7 ra2b(1,0,2) ; 

  ra2b.initialize("an-11-257-v3-files/input-files/byhand-data-ge1b-tight.txt",
                  "an-11-257-v3-files/input-files/signalSyst.LM9.ge1bTight.1143invpb.dat", 0., 0., true, 1.0 ) ;

  gSystem->Exec("mv ws-lm9-ge1btight.root ws-lm9-ge1btight.root-old") ;
  gSystem->Exec("mv ws.root ws-lm9-ge1btight.root") ;

}



