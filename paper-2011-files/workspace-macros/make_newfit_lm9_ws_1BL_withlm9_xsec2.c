
void make_newfit_lm9_ws_1BL_withlm9_xsec2() {

  gROOT->LoadMacro("ra2bRoostatsClass9.c+") ;

  ra2bRoostatsClass9 ra2b(0,1) ;

  ra2b.initialize("paper-2011-files/input-files/likelihood-newfit-input-from-mc-HT400-SIGMET250-withLM9-xSec2.0.txt",
                  "paper-2011-files/input-files/likelihood-newfit-susyinput-LM9-HT400-SIGMET250.txt",
                  0., 0., false, 1.0,
                  "paper-2011-files/input-files/likelihood-newfit-syst-deff_dbtageff-LM9-HT400-SIGMET250.txt"
                 ) ;

  gSystem->Exec("mkdir -p output-files") ;
  gSystem->Exec("mv output-files/ws-newfit-lm9-1BL-withLM9-xSec2.0.root output-files/ws-newfit-lm9-1BL-withLM9-xSec2.0.root-old") ;
  gSystem->Exec("mv ws.root output-files/ws-newfit-lm9-1BL-withLM9-xSec2.0.root") ;

}



