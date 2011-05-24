

   void dofit_sigvars_ge1() {

      gROOT->LoadMacro("ra2bRoostatsClass.c+") ;

      ra2bRoostatsClass rfit(1) ;
      rfit.initialize("ra2b-roostats-v4-input-2011b-ge1.txt") ;
      rfit.doFit() ;


   }
