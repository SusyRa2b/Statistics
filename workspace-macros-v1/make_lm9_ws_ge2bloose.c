
void make_lm9_ws_ge2bloose() {

  gROOT->LoadMacro("ra2bRoostatsClass7.c+") ;

  ra2bRoostatsClass7 ra2b(1,0) ;

  ra2b.initialize("an-11-257-v3-files/input-files/byhand-data-ge2b-loose-newformat.txt",
                  "an-11-257-v3-files/input-files/signalSyst.LM9.ge2bLoose.1143invpb.dat", 0., 0., true, 1.0 ) ;

  gSystem->Exec("mkdir -p output-files") ;
  gSystem->Exec("mv output-files/ws-lm9-ge2bloose.root output-files/ws-lm9-ge2bloose.root-old") ;
  gSystem->Exec("mv ws.root output-files/ws-lm9-ge2bloose.root") ;

}



