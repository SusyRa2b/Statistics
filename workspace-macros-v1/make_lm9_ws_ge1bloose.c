
void make_lm9_ws_ge1bloose() {

  gROOT->LoadMacro("ra2bRoostatsClass7.c+") ;

  ra2bRoostatsClass7 ra2b(1,0) ;

  ra2b.initialize("an-11-257-v3-files/input-files/byhand-data-ge1b-loose-newformat.txt",
                  "an-11-257-v3-files/input-files/signalSyst.LM9.ge1bLoose.1143invpb.dat", 0., 0., true, 1.0 ) ;

  gSystem->Exec("mkdir -p output-files") ;
  gSystem->Exec("mv output-files/ws-lm9-ge1bloose.root output-files/ws-lm9-ge1bloose.root-old") ;
  gSystem->Exec("mv ws.root output-files/ws-lm9-ge1bloose.root") ;

}



