

   void draw_likelihood_closure( const char* ws_filebase = "ws-met3-ht3-v2",
                                int nMetBins = 3,
                                int nHtBins = 3 ) {


       char command[1000] ;
       sprintf( command, "mkdir -p rootfiles/likelihood-closure-%s", ws_filebase ) ;
       gSystem -> Exec( command ) ;

       for ( int bbi=1; bbi<=nMetBins; bbi++ ) {
          for ( int mbi=1; mbi<=nMetBins; mbi++ ) {
             for ( int hbi=1; hbi<=nMetBins; hbi++ ) {

                char new_poi_name[1000] ;
                sprintf( new_poi_name, "n_M%d_H%d_%db", mbi, hbi, bbi ) ;
                char rootfile_name[1000] ;
                sprintf( rootfile_name, "rootfiles/likelihood-closure-%s/likelihood-closure-%s-%s.root", ws_filebase, ws_filebase, new_poi_name ) ;

                TFile f( rootfile_name ) ;
                char hname[1000] ;
                sprintf( hname, "hscanout_%s", new_poi_name ) ;
                TH1F* hp = (TH1F*) f.Get( hname ) ;
                if ( hp == 0x0 ) {
                   printf("\n\n *** histogram %s missing from file %s\n\n\n", hname, rootfile_name ) ;
                   return ;
                }

                double obsVal = hp->GetBinContent(1) ;
                double modelPlus1sd  = hp->GetBinContent(2) ;
                double modelVal      = hp->GetBinContent(3) ;
                double modelMinus1sd = hp->GetBinContent(4) ;

                printf("   met=%d, ht=%d, nb=%d :  obs=%6.0f,    model=%6.1f    [%.1f, %.1f]\n", mbi, hbi, bbi, obsVal, modelVal,
                    modelMinus1sd, modelPlus1sd ) ;


             } // hbi.
          } // mbi.
       } // bbi.

   }
