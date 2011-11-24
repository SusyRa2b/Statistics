
void make_mc_smsusy_lm9_ws_ge1bloose() {

  gROOT->LoadMacro("ra2bRoostatsClass7.c+") ;

  ra2bRoostatsClass7 ra2b(1,0,1) ;

  ra2b.initialize("an-11-257-v3-files/input-files/mc-inputs-SM+lm9-loose-sel-newformat.txt",
                  "an-11-257-v3-files/input-files/signalSyst.LM9.ge1bLoose.1143invpb.dat", 0., 0., true, 1.0 ) ;

  gSystem->Exec("mkdir -p output-files") ;
  gSystem->Exec("mv output-files/ws-mc-smsusy-lm9-ge1bloose.root output-files/ws-mc-smsusy-lm9-ge1bloose.root-old") ;
  gSystem->Exec("mv ws.root output-files/ws-mc-smsusy-lm9-ge1bloose.root") ;

}



