

   void run_likelihood_closure( const char* ws_filebase = "ws-met3-ht3-v2",
                                int nMetBins = 3,
                                int nHtBins = 3 ) {

       gROOT -> LoadMacro("ws_constrained_profile3D.c+") ;

       char ws_filename[1000] ;
       sprintf( ws_filename, "rootfiles/%s.root", ws_filebase ) ;

       char command[1000] ;
       sprintf( command, "mkdir -p rootfiles/likelihood-closure-%s", ws_filebase ) ;
       gSystem -> Exec( command ) ;

       for ( int mbi=1; mbi<=nMetBins; mbi++ ) {
          for ( int hbi=1; hbi<=nMetBins; hbi++ ) {
             for ( int bbi=1; bbi<=nMetBins; bbi++ ) {
                char new_poi_name[1000] ;
                sprintf( new_poi_name, "n_M%d_H%d_%db", mbi, hbi, bbi ) ;
                char ignore_observable_name[1000] ;
                sprintf( ignore_observable_name, "N_0lep_M%d_H%d_%db", mbi, hbi, bbi ) ;
                ws_constrained_profile3D( ws_filename, new_poi_name, ignore_observable_name, 16 ) ;
                char new_rootfile_name[1000] ;
                sprintf( new_rootfile_name, "rootfiles/likelihood-closure-%s/likelihood-closure-%s-%s.root", ws_filebase, ws_filebase, new_poi_name ) ;
                sprintf( command, "mv cscan.root %s", new_rootfile_name ) ;
                gSystem -> Exec( command ) ;
                gROOT->Reset() ;
             } // bbi.
          } // hbi.
       } // mbi.

   }
