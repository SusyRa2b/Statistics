

#include "TTree.h"
#include "TFile.h"
#include <iostream>

  using namespace std ;




    void toyreadtest( const char* fname="ra2b-roostats-v4-input-2011-ge1.txt-toy.root" ) {

       TFile intoyfile( fname, "READ" ) ;
       gDirectory->ls() ;
       TTree* toytree = (TTree*) gDirectory->FindObjectAny("likelihoodData") ;

       if ( toytree == 0 ) {
          printf(" \n\n *** Can't find TTree likelihoodData in file %s\n", fname ) ;
       } else {
          toytree->Print("toponly") ;
       }

       double toyNsig ;
       TBranch* b_toyNsig ;

       toytree->SetBranchAddress("Nsig", &toyNsig, &b_toyNsig ) ;

       for ( int dsi=0; dsi<10; dsi++ ) {
          toytree->GetEntry( dsi ) ;
          printf("  toy dataset %3d : Nsig = %5.0f\n", dsi, toyNsig ) ;
       }


    }




