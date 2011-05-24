

   void dotoy_sigvars_mc() {

      gROOT->LoadMacro("ra2bRoostatsClass.c+") ;

      ra2bRoostatsClass rgen(1) ;
      rgen.initialize("ra2b-roostats-v4-mcinput-2011b-ge1.txt") ;
      rgen.generateToyDatasetsFromLikelihood("toy-datasets-2011b-ge1-sigvars-mc.root",1000) ;


      ra2bRoostatsClass rfit(1) ;
      rfit.initialize("ra2b-roostats-v4-mcinput-2011b-ge1.txt") ;
      rfit.doToyStudy("toy-datasets-2011b-ge1-sigvars-mc.root","toy-fits-2011b-ge1-sigvars-mc.root",0,1000) ;


   }
