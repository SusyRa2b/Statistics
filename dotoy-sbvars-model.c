

   void dotoy_sbvars_model() {

      gROOT->LoadMacro("ra2bRoostatsClass.c+") ;

      ra2bRoostatsClass rgen(0) ;
      rgen.initialize("ra2b-roostats-v4-mcinput-2011-ge1.txt") ;
      rgen.generateToyDatasetsFromNVals("toy-datasets-2011-ge1-sbvars-model.root",1000) ;


      ra2bRoostatsClass rfit(0) ;
      rfit.initialize("ra2b-roostats-v4-mcinput-2011-ge1.txt") ;
      rfit.doToyStudy("toy-datasets-2011-ge1-sbvars-model.root","toy-fits-2011-ge1-sbvars-model.root",0,1000) ;


   }
