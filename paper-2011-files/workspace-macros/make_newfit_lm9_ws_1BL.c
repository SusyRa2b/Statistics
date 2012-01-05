
void make_newfit_lm9_ws_1BL() {

  gROOT->LoadMacro("ra2bRoostatsClass9.c+") ;

  ra2bRoostatsClass9 ra2b(0,1) ;

  ra2b.initialize("paper-2011-files/input-files/likelihood-newfit-input-from-mc-HT400-SIGMET250.txt",
                  "paper-2011-files/input-files/likelihood-newfit-susyinput-LM9-HT400-SIGMET250.txt", 0., 0., false, 1.0 ) ;

  gSystem->Exec("mkdir -p output-files") ;
  gSystem->Exec("mv output-files/ws-newfit-lm9-1BL.root output-files/ws-newfit-lm9-1BL.root-old") ;
  gSystem->Exec("mv ws.root output-files/ws-newfit-lm9-1BL.root") ;

}



