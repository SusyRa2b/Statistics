
#include "TChain.h"
#include "TH2F.h"
#include <fstream>



    void MakeDataFile( const char* outputfilename = "datfiles/data-vals-unblind.dat" ) {

  //-- met4-ht4-v15
       const int nBinsBjets = 3 ;
       const int nBinsMET   = 4 ;
       const int nBinsHT    = 4 ;
       float Mbins[nBinsMET+1] = {125.,150.,250.,350.,99999.};
       float Hbins[nBinsHT+1] = {400.,500.,800.,1000.,99999.};

       bool binIsBlind[nBinsMET][nBinsHT][nBinsBjets] ;
       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
             for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
                binIsBlind[mbi][hbi][bbi] = false ;
             } // bbi.
          } // hbi.
       } // mbi.

     //--- nB=2, highest met bin
       binIsBlind[3][0][1] = false ;
       binIsBlind[3][1][1] = false ;
       binIsBlind[3][2][1] = false ;
       binIsBlind[3][3][1] = false ;

     //--- nB=3, top 3 met bins
       binIsBlind[1][0][2] = false ;
       binIsBlind[1][1][2] = false ;
       binIsBlind[1][2][2] = false ;
       binIsBlind[1][3][2] = false ;

       binIsBlind[2][0][2] = false ;
       binIsBlind[2][1][2] = false ;
       binIsBlind[2][2][2] = false ;
       binIsBlind[2][3][2] = false ;

       binIsBlind[3][0][2] = false ;
       binIsBlind[3][1][2] = false ;
       binIsBlind[3][2][2] = false ;
       binIsBlind[3][3][2] = false ;

       TChain datach("tree") ;

       datach.Add( "filesHCP_53_v6/MET_2012A.root" ) ;
       datach.Add( "filesHCP_53_v6/MET_2012B.root" ) ;
       datach.Add( "filesHCP_53_v6/MET_2012C_pr.root" ) ;
       datach.Add( "filesHCP_53_v6/MET_2012C_rr.root" ) ;
       datach.Add( "filesHCP_53_v6/HT_2012A.root" ) ;
       datach.Add( "filesHCP_53_v6/HTMHT_2012B.root" ) ;
       datach.Add( "filesHCP_53_v6/HTMHT_2012C_pr.root" ) ;
       datach.Add( "filesHCP_53_v6/HTMHT_2012C_rr.root" ) ;
       datach.Add( "filesHCP_53_v6/JetHT_2012B.root" ) ;
       datach.Add( "filesHCP_53_v6/JetHT_2012C_pr.root" ) ;
       datach.Add( "filesHCP_53_v6/JetHT_2012C_rr.root" ) ;

       TH2F* h_ldp[nBinsBjets] ;
       TH2F* h_zl [nBinsBjets] ;
       TH2F* h_sl [nBinsBjets] ; 
  

       for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

          char hname[100] ;
          char htitle[100] ;

          sprintf( hname, "h_zl_%db", bbi+1 ) ;
          sprintf( htitle, "Data, ZL, nB=%d", bbi+1 ) ;
          h_zl [bbi] = new TH2F( hname, htitle, nBinsMET, Mbins, nBinsHT, Hbins ) ;

          sprintf( hname, "h_ldp_%db", bbi+1 ) ;
          sprintf( htitle, "Data, LDP, nB=%d", bbi+1 ) ;
          h_ldp[bbi] = new TH2F( hname, htitle, nBinsMET, Mbins, nBinsHT, Hbins ) ;

          sprintf( hname, "h_sl_%db", bbi+1 ) ;
          sprintf( htitle, "Data, SL, nB=%d", bbi+1 ) ;
          h_sl [bbi] = new TH2F( hname, htitle, nBinsMET, Mbins, nBinsHT, Hbins ) ;

       } // bbi.

       const int nJetsCut = 3 ;     // #jets >= nJetsCut
       double minLeadJetPt = 70. ;
       double min3rdJetPt = 50. ;

       char commoncuts[10000] ;
       sprintf( commoncuts, "maxChNMultDiff<40&&pfOcaloMET<2.0&&nJets>=%d&&(pt_1st_leadJet>%.0f&&pt_2nd_leadJet>%.0f&&pt_3rd_leadJet>%.0f)&&nB>0&&MET>125&&HT>400&&passedTrigger==1",
             nJetsCut, minLeadJetPt, minLeadJetPt, min3rdJetPt ) ;

       char basecuts_ldp[10000] ;
       sprintf( basecuts_ldp, "%s&&minDelPhiN<=4&&(nMu==0&&nEl==0)&&nIsoTrk==0", commoncuts ) ;

       char basecuts_sl[10000] ;
       sprintf( basecuts_sl, "%s&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&MT<100", commoncuts ) ;

       char basecuts_zl[10000] ;
       sprintf( basecuts_zl, "%s&&minDelPhiN>4&&(nMu==0&&nEl==0)&&nIsoTrk==0", commoncuts ) ;

       char bcutstring[3][100] = { "nB==1", "nB==2", "nB>=3" } ;

       for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

          char arg1[1000] ;
          char cuts[10000] ;

          sprintf( arg1, "HT:MET>>h_ldp_%db", bbi+1 ) ;
          sprintf( cuts, "(%s)&&(%s)", basecuts_ldp, bcutstring[bbi] ) ;
          printf( "%s : %s\n", arg1, cuts ) ;
          datach.Draw( arg1, cuts ) ;

          sprintf( arg1, "HT:MET>>h_sl_%db", bbi+1 ) ;
          sprintf( cuts, "(%s)&&(%s)", basecuts_sl, bcutstring[bbi] ) ;
          printf( "%s : %s\n", arg1, cuts ) ;
          datach.Draw( arg1, cuts ) ;

          sprintf( arg1, "HT:MET>>h_zl_%db", bbi+1 ) ;
          sprintf( cuts, "(%s)&&(%s)", basecuts_zl, bcutstring[bbi] ) ;
          printf( "%s : %s\n", arg1, cuts ) ;
          datach.Draw( arg1, cuts ) ;

       } // bbi.

     //--- make sure we are blind.
       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
             for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
                if ( binIsBlind[mbi][hbi][bbi] ) {
                   h_zl[bbi] -> SetBinContent( mbi+1, hbi+1, 0. ) ;
                }
             } // bbi.
          } // hbi.
       } // mbi.






     //--- print tables.

       printf("\n                     LDP, nB=1                                        LDP, nB=2                                         LDP, nB>=3\n" ) ;
       printf("            H1      H2      H3      H4      |                H1      H2      H3      H4      |                 H1      H2      H3      H4     |\n") ;
       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
             printf("   M%d  :", mbi+1 ) ;
             for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
                printf(" %6.0f ", h_ldp[bbi]->GetBinContent( mbi+1, hbi+1 ) ) ;
             } // hbi
             printf("    |    ") ;
          } // bbi
          printf("\n") ;
       } // mbi

       printf("\n\n") ;
       printf("\n                     SL , nB=1                                        SL , nB=2                                         SL , nB>=3\n" ) ;
       printf("            H1      H2      H3      H4      |                H1      H2      H3      H4      |                 H1      H2      H3      H4     |\n") ;
       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
             printf("   M%d  :", mbi+1 ) ;
             for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
                printf(" %6.0f ", h_sl[bbi]->GetBinContent( mbi+1, hbi+1 ) ) ;
             } // hbi
             printf("    |    ") ;
          } // bbi
          printf("\n") ;
       } // mbi

       printf("\n\n") ;
       printf("\n                     ZL , nB=1                                        ZL , nB=2                                         ZL , nB>=3\n" ) ;
       printf("            H1      H2      H3      H4      |                H1      H2      H3      H4      |                 H1      H2      H3      H4     |\n") ;
       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
             printf("   M%d  :", mbi+1 ) ;
             for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
                if ( binIsBlind[mbi][hbi][bbi] ) {
                   printf("   blind" ) ;
                } else {
                   printf(" %6.0f ", h_zl[bbi]->GetBinContent( mbi+1, hbi+1 ) ) ;
                }
             } // hbi
             printf("    |    ") ;
          } // bbi
          printf("\n") ;
       } // mbi
       printf("\n\n\n") ;




       //--- write results to file.

       FILE* fp = fopen( outputfilename, "w" ) ;

       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
             for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
                fprintf( fp, "N_0lep_M%d_H%d_%db    %g\n", mbi+1, hbi+1, bbi+1, h_zl[bbi]->GetBinContent( mbi+1, hbi+1 ) ) ;
             } // bbi
          } // hbi
       } // bbi

       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
             for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
                fprintf( fp, "N_1lep_M%d_H%d_%db    %g\n", mbi+1, hbi+1, bbi+1, h_sl[bbi]->GetBinContent( mbi+1, hbi+1 ) ) ;
             } // bbi
          } // hbi
       } // bbi

       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
             for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
                fprintf( fp, "N_ldp_M%d_H%d_%db    %g\n", mbi+1, hbi+1, bbi+1, h_ldp[bbi]->GetBinContent( mbi+1, hbi+1 ) ) ;
             } // bbi
          } // hbi
       } // bbi

       fclose( fp ) ;
       printf("\n\n Wrote results to %s\n\n", outputfilename ) ;

    } // MakeDataFile





