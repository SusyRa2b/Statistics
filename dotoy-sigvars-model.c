

   void dotoy_sigvars_model() {

      gROOT->LoadMacro("ra2bRoostatsClass.c+") ;

      ra2bRoostatsClass rgen(1) ;
      rgen.initialize("ra2b-roostats-v4-mcinput-2011-ge2.txt") ;
      rgen.generateToyDatasetsFromNVals("toy-datasets-2011-ge2-sigvars-model.root",1000) ;


      ra2bRoostatsClass rfit(1) ;
      rfit.initialize("ra2b-roostats-v4-mcinput-2011-ge2.txt") ;
      rfit.doToyStudy("toy-datasets-2011-ge2-sigvars-model.root","toy-fits-2011-ge2-sigvars-model.root",0,1000) ;


   }
