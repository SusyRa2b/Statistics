

   void dotoy_sbvars() {

      gROOT->LoadMacro("ra2bRoostatsClass.c+") ;

      ra2bRoostatsClass rgen(0) ;
      rgen.initialize("ra2b-roostats-v4-input-2011-ge1.txt") ;
      rgen.generateToyDatasets("toy-datasets-2011-ge1-sbvars.root",1000) ;


      ra2bRoostatsClass rfit(0) ;
      rfit.initialize("ra2b-roostats-v4-input-2011-ge1.txt") ;
      rfit.doToyStudy("toy-datasets-2011-ge1-sbvars.root","toy-fits-2011-ge1-sbvars.root",0,1000) ;


   }
