

   void dotoy_sigvars_pure() {

      gROOT->LoadMacro("ra2bRoostatsClass.c+") ;

      ra2bRoostatsClass rgen(1) ;
      rgen.initialize("ra2b-roostats-v4-pureinput-2011-ge1.txt") ;
      rgen.generateToyDatasetsFromLikelihood("toy-datasets-2011-ge1-sigvars-pure.root",1000) ;


      ra2bRoostatsClass rfit(1) ;
      rfit.initialize("ra2b-roostats-v4-pureinput-2011-ge1.txt") ;
      rfit.doToyStudy("toy-datasets-2011-ge1-sigvars-pure.root","toy-fits-2011-ge1-sigvars-pure.root",0,1000) ;


   }
