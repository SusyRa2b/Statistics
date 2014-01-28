// create the "SUSY cocktail" MC from the following samples:
//
// 0: HH -> bbbb
// 1: HH -> WWbb
// 2: HH -> tautaubb
//
// weighting them by reconstruction efficiency and product BF's
//


void GenerateSusyCocktail(int version) {

  gROOT->Reset() ;
  gStyle->SetOptStat(0) ;

  const int nBinsMET  = 3 ;
  const int nBinsHT   = 4 ;
  const int nBinsBtag = 3 ;

  TString sMbins[nBinsMET];
  TString sHbins[nBinsHT];
  TString sBbins[3] = {"_1b","_2b","_3b"};
     
  for (int i = 0 ; i < nBinsMET ; i++) {
    TString base = "_M";
    stringstream sbin;
    sbin << i+1;
    base += sbin.str();
    sMbins[i] = base;
  }
      
  for (int j = 0 ; j < nBinsHT ; j++) {
    TString base = "_H";
    stringstream sbin;
    sbin << j+1;
    base += sbin.str();
    sHbins[j] = base;
  }

  const int nPoints(15) ;
  const int nSamp(3) ;
  const double lumi = 19399. ;

  // samples:
  TString Samples[nSamp] = {"bbbb","WWbb","2taubb"} ; 

  // set here the (hardcoded) "product BF's"
  const double prodBF[nSamp] = { 0.3147,
				 0.2592,
				 0.0690 } ;

  double ZlepEff[nSamp] ;
  double ZlepTot[nSamp] ;
  double ExpZlepTot[nSamp] ;
  double weightFactor[nSamp] ;

  const int ArraySize = 291 ;
  double ArrayContent[nSamp][ArraySize] ;
  double ArrayCocktail[ArraySize] ;

  // as usual, hardcode the predicted Xsec:
  TH1F *Xsec_p = new TH1F("Xsec_p","",15,137.5,512.5);

  Xsec_p->SetBinContent(1,1.876);
  Xsec_p->SetBinContent(2,1.027);
  Xsec_p->SetBinContent(3,0.608);
  Xsec_p->SetBinContent(4,0.377);
  Xsec_p->SetBinContent(5,0.244);
  Xsec_p->SetBinContent(6,0.162);
  Xsec_p->SetBinContent(7,0.111);
  Xsec_p->SetBinContent(8,0.0779);
  Xsec_p->SetBinContent(9,0.0553);
  Xsec_p->SetBinContent(10,0.0401);
  Xsec_p->SetBinContent(11,0.0294);
  Xsec_p->SetBinContent(12,0.0218);
  Xsec_p->SetBinContent(13,0.0163);
  Xsec_p->SetBinContent(14,0.0123);
  Xsec_p->SetBinContent(15,0.0094);

  // Owen style histograms
  const int nbins = nBinsMET*(nBinsHT+1) + 1 ;

  TH1F *hOS_0lep_1b[nSamp][nPoints] ;
  TH1F *hOS_0lep_2b[nSamp][nPoints] ;
  TH1F *hOS_0lep_3b[nSamp][nPoints] ;

  TH1F *hOS_1lepSig_1b[nSamp][nPoints] ;
  TH1F *hOS_1lepSig_2b[nSamp][nPoints] ;
  TH1F *hOS_1lepSig_3b[nSamp][nPoints] ;

  TH1F *hOS_1lep_1b[nSamp][nPoints] ;
  TH1F *hOS_1lep_2b[nSamp][nPoints] ;
  TH1F *hOS_1lep_3b[nSamp][nPoints] ;

  for ( int i = 0 ; i < nPoints ; i++ ) {
    for ( int j = 0 ; j < nSamp ; j++ ) {
      
      TString Mass = "_" ;
      Mass += 125+25*i ;
      
      hOS_0lep_1b[j][i] = new TH1F("hOS_0lep_1b_"+Samples[j]+Mass,"hOS_0lep_1b",nbins,0.5,nbins+0.5);
      hOS_0lep_2b[j][i] = new TH1F("hOS_0lep_2b_"+Samples[j]+Mass,"hOS_0lep_2b",nbins,0.5,nbins+0.5);
      hOS_0lep_3b[j][i] = new TH1F("hOS_0lep_3b_"+Samples[j]+Mass,"hOS_0lep_3b",nbins,0.5,nbins+0.5);
      
      hOS_1lepSig_1b[j][i] = new TH1F("hOS_1lepSig_1b_"+Samples[j]+Mass,"hOS_1lepSig_1b",nbins,0.5,nbins+0.5);
      hOS_1lepSig_2b[j][i] = new TH1F("hOS_1lepSig_2b_"+Samples[j]+Mass,"hOS_1lepSig_2b",nbins,0.5,nbins+0.5);
      hOS_1lepSig_3b[j][i] = new TH1F("hOS_1lepSig_3b_"+Samples[j]+Mass,"hOS_1lepSig_3b",nbins,0.5,nbins+0.5);
      
      hOS_1lep_1b[j][i] = new TH1F("hOS_1lep_1b_"+Samples[j]+Mass,"hOS_1lep_1b",nbins,0.5,nbins+0.5);
      hOS_1lep_2b[j][i] = new TH1F("hOS_1lep_2b_"+Samples[j]+Mass,"hOS_1lep_2b",nbins,0.5,nbins+0.5);
      hOS_1lep_3b[j][i] = new TH1F("hOS_1lep_3b_"+Samples[j]+Mass,"hOS_1lep_3b",nbins,0.5,nbins+0.5);
      
    }
  }
  
  
  // get efficiencies from file
  TH1F *heff[nSamp] ;
  TFile *effFile[nSamp] ;

  for ( int samp = 0 ; samp < nSamp ; samp++ ) {
    
    TString inFile = "cocktail_dir/efficiency_TChiHH_" ;
    inFile += Samples[samp] ;
    inFile += "_v" ;
    inFile += version ;
    inFile += ".root" ;

    effFile[samp] = new TFile(inFile);
    heff[samp] = (TH1F*)effFile[samp]->Get("heff");

  }
    
  TString outFile = "cocktail_dir/TChiHH_cocktail_v" ;
  outFile += version ;
  outFile += ".dat" ;

  ofstream outStream ;
  outStream.open(outFile) ;

  bool filledSamp[3] = {false, false, false} ;
  
  // loop on Mass points
  for ( int ChiMass = 150 ; ChiMass <= 500 ; ChiMass += 25 ) {
      
    int bin = ChiMass/25 - 5 ;
    
    // initialize the "cocktail" array
    ArrayCocktail[0] = ChiMass ;
    ArrayCocktail[1] = 0 ;

    // the cocktail is set up to contain the expected number of events
    // for the given lumi and efficiency. To get the efficiency right, the
    // number of generated events must correspond to the total number of 
    // TChiHH events (before any Higgs decay and any efficiency):

    ArrayCocktail[2] = lumi * Xsec_p->GetBinContent(bin) ; 
    
    for (int i = 3; i < 147; ++ i) {
      ArrayCocktail[i] = 0. ;
    }
    for (int i = 147; i < ArraySize; ++ i) {
      ArrayCocktail[i] = 0.00001 ;
    }
    
    for ( int samp = 0 ; samp < nSamp ; samp++ ) {
      
      // first, compute the expected 0-lep yield
      ExpZlepTot[samp] = lumi * prodBF[samp] * Xsec_p->GetBinContent(bin) * heff[samp]->GetBinContent(bin) ;
      
      // then read in the input file
      TString inFile = "cocktail_dir/TChiHH_" ;
      inFile += Samples[samp] ;
      inFile += "-v" ;
      inFile += version ;
      inFile += ".dat" ;
      
      ifstream infp ;
      infp.open(inFile) ;
      
      while ( infp.good() ) {
	
	for (int i = 0; infp && i < ArraySize; ++ i) {
	  infp >> ArrayContent[samp][i];
	}
	
	if ( ArrayContent[samp][0] != ChiMass || filledSamp[samp] ) continue ;
	
	ZlepTot[samp] = 0. ;
	  
	for (int i = 3; i < 39 ; i++ ) {
	  ZlepTot[samp] += ArrayContent[samp][i] ;
	}
	
	// weight factor
	weightFactor[samp] = ExpZlepTot[samp] / ZlepTot[samp] ;
	cout << "\n\ndebugging: point (" << ChiMass << ",0) ; sample = " << Samples[samp] << " ; weightFactor = " << weightFactor[samp] << endl ;
	
	// rescale yields by the weight factor (relative errors stay the same) 
	for (int i = 3; i < 147; ++ i) {
	  ArrayContent[samp][i] = ArrayContent[samp][i] * weightFactor[samp] ;
	  ArrayCocktail[i] += ArrayContent[samp][i] ;
	}
	
	// fill the array of errors with the absolute errors
	for (int i = 147; i < ArraySize; ++ i) {
	  ArrayContent[samp][i] = ArrayContent[samp][i] * weightFactor[samp] ;
	  ArrayCocktail[i] = sqrt( pow(ArrayCocktail[i],2) + pow(ArrayContent[samp][i],2) ) ;
	}
	
	if ( ChiMass == 500 ) { 
	  filledSamp[samp] = true ;
	  cout << "Setting filledSamp to true!" << endl ;
	}
	
	// fill Owen style histograms
	for ( int i = 0 ; i < 3 ; i++ ) {
	  for ( int j = 0 ; j < 4 ; j++ ) {
	    
	    int binIndex = 1 + 5*i + j + 1 ;
	    
	    TString binLabel_0lep    = "0 lep";
	    TString binLabel_1lepSig = "1 lepSig";
	    TString binLabel_1lep    = "1 lep";
	    
	    binLabel_0lep    += sMbins[i]+sHbins[j] ;
	    binLabel_1lepSig += sMbins[i]+sHbins[j] ;
	    binLabel_1lep    += sMbins[i]+sHbins[j] ;
	    
	    
	    hOS_0lep_1b[samp][bin-1]->SetBinContent(binIndex,ArrayContent[samp][3+3*j+12*i]) ;
	    hOS_0lep_2b[samp][bin-1]->SetBinContent(binIndex,ArrayContent[samp][4+3*j+12*i]) ;
	    hOS_0lep_3b[samp][bin-1]->SetBinContent(binIndex,ArrayContent[samp][5+3*j+12*i]) ;
	    
	    hOS_1lepSig_1b[samp][bin-1]->SetBinContent(binIndex,ArrayContent[samp][39+3*j+12*i]) ;
	    hOS_1lepSig_2b[samp][bin-1]->SetBinContent(binIndex,ArrayContent[samp][40+3*j+12*i]) ;
	    hOS_1lepSig_3b[samp][bin-1]->SetBinContent(binIndex,ArrayContent[samp][41+3*j+12*i]) ;
	      
	    hOS_1lep_1b[samp][bin-1]->SetBinContent(binIndex,ArrayContent[samp][75+3*j+12*i]) ;
	    hOS_1lep_2b[samp][bin-1]->SetBinContent(binIndex,ArrayContent[samp][76+3*j+12*i]) ;
	    hOS_1lep_3b[samp][bin-1]->SetBinContent(binIndex,ArrayContent[samp][77+3*j+12*i]) ;
	    
	  }
	}
	
      }

    }  // loop on nSamp


    // dump to file new cocktail input file

    for (int i = 0; i < ArraySize; ++ i) {
      outStream << ArrayCocktail[i] << " " ;  
    }
    outStream << endl ;

  } // loop on ChiMass

  outStream.close() ;


  // draw "Owen style" histograms of the SUSY cocktail

  for ( int ChiMass = 150 ; ChiMass <= 500 ; ChiMass += 25 ) {

    int bin = ChiMass/25 - 5 ;

    THStack *hOSS_0lep_1b, *hOSS_0lep_2b, *hOSS_0lep_3b ;
    THStack *hOSS_1lepSig_1b, *hOSS_1lepSig_2b, *hOSS_1lepSig_3b ;
    THStack *hOSS_1lep_1b, *hOSS_1lep_2b, *hOSS_1lep_3b ;

    hOSS_0lep_1b = new THStack("hOSS_0lep_1b","" );
    hOSS_0lep_2b = new THStack("hOSS_0lep_2b","" );
    hOSS_0lep_3b = new THStack("hOSS_0lep_3b","" );

    hOSS_1lepSig_1b = new THStack("hOSS_1lepSig_1b","" );
    hOSS_1lepSig_2b = new THStack("hOSS_1lepSig_2b","" );
    hOSS_1lepSig_3b = new THStack("hOSS_1lepSig_3b","" );

    hOSS_1lep_1b = new THStack("hOSS_1lep_1b","" );
    hOSS_1lep_2b = new THStack("hOSS_1lep_2b","" );
    hOSS_1lep_3b = new THStack("hOSS_1lep_3b","" );


    // style

    hOS_0lep_1b[0][bin-1]->SetFillColor(kMagenta);
    hOS_0lep_2b[0][bin-1]->SetFillColor(kMagenta);
    hOS_0lep_3b[0][bin-1]->SetFillColor(kMagenta);
    
    hOS_1lepSig_1b[0][bin-1]->SetFillColor(kMagenta);
    hOS_1lepSig_2b[0][bin-1]->SetFillColor(kMagenta);
    hOS_1lepSig_3b[0][bin-1]->SetFillColor(kMagenta);
    
    hOS_1lep_1b[0][bin-1]->SetFillColor(kMagenta);
    hOS_1lep_2b[0][bin-1]->SetFillColor(kMagenta);
    hOS_1lep_3b[0][bin-1]->SetFillColor(kMagenta);
    

    hOS_0lep_1b[1][bin-1]->SetFillColor(kBlue+2);
    hOS_0lep_2b[1][bin-1]->SetFillColor(kBlue+2);
    hOS_0lep_3b[1][bin-1]->SetFillColor(kBlue+2);
    
    hOS_1lepSig_1b[1][bin-1]->SetFillColor(kBlue+2);
    hOS_1lepSig_2b[1][bin-1]->SetFillColor(kBlue+2);
    hOS_1lepSig_3b[1][bin-1]->SetFillColor(kBlue+2);
    
    hOS_1lep_1b[1][bin-1]->SetFillColor(kBlue+2);
    hOS_1lep_2b[1][bin-1]->SetFillColor(kBlue+2);
    hOS_1lep_3b[1][bin-1]->SetFillColor(kBlue+2);
    

    hOS_0lep_1b[2][bin-1]->SetFillColor(kGreen+2);
    hOS_0lep_2b[2][bin-1]->SetFillColor(kGreen+2);
    hOS_0lep_3b[2][bin-1]->SetFillColor(kGreen+2);
    
    hOS_1lepSig_1b[2][bin-1]->SetFillColor(kGreen+2);
    hOS_1lepSig_2b[2][bin-1]->SetFillColor(kGreen+2);
    hOS_1lepSig_3b[2][bin-1]->SetFillColor(kGreen+2);
    
    hOS_1lep_1b[2][bin-1]->SetFillColor(kGreen+2);
    hOS_1lep_2b[2][bin-1]->SetFillColor(kGreen+2);
    hOS_1lep_3b[2][bin-1]->SetFillColor(kGreen+2);
    
    for ( int samp = 0 ; samp < nSamp ; samp++ ) {

      hOSS_0lep_1b->Add(hOS_0lep_1b[samp][bin-1]) ;
      hOSS_0lep_2b->Add(hOS_0lep_2b[samp][bin-1]) ;
      hOSS_0lep_3b->Add(hOS_0lep_3b[samp][bin-1]) ;

      hOSS_1lepSig_1b->Add(hOS_1lepSig_1b[samp][bin-1]) ;
      hOSS_1lepSig_2b->Add(hOS_1lepSig_2b[samp][bin-1]) ;
      hOSS_1lepSig_3b->Add(hOS_1lepSig_3b[samp][bin-1]) ;

      hOSS_1lep_1b->Add(hOS_1lep_1b[samp][bin-1]) ;
      hOSS_1lep_2b->Add(hOS_1lep_2b[samp][bin-1]) ;
      hOSS_1lep_3b->Add(hOS_1lep_3b[samp][bin-1]) ;

    }


    // cout the amount of signal for each dataset

    double Tot0lep(0), Tot1lepSig(0), Tot1lep(0) ;
    
    for ( int i = 0 ; i < nBinsMET ; i++ ) {
      for ( int j = 0 ; j < nBinsHT ; j++ ) {

	int binIndex = 1 + 5*i + j + 1 ;

	Tot0lep += hOS_0lep_1b[0][bin-1]->GetBinContent(binIndex) ;
	Tot0lep += hOS_0lep_1b[1][bin-1]->GetBinContent(binIndex) ;
	Tot0lep += hOS_0lep_1b[2][bin-1]->GetBinContent(binIndex) ;

	Tot0lep += hOS_0lep_2b[0][bin-1]->GetBinContent(binIndex) ;
	Tot0lep += hOS_0lep_2b[1][bin-1]->GetBinContent(binIndex) ;
	Tot0lep += hOS_0lep_2b[2][bin-1]->GetBinContent(binIndex) ;

	Tot0lep += hOS_0lep_3b[0][bin-1]->GetBinContent(binIndex) ;
	Tot0lep += hOS_0lep_3b[1][bin-1]->GetBinContent(binIndex) ;
	Tot0lep += hOS_0lep_3b[2][bin-1]->GetBinContent(binIndex) ;


	Tot1lepSig += hOS_1lepSig_1b[0][bin-1]->GetBinContent(binIndex) ;
	Tot1lepSig += hOS_1lepSig_1b[1][bin-1]->GetBinContent(binIndex) ;
	Tot1lepSig += hOS_1lepSig_1b[2][bin-1]->GetBinContent(binIndex) ;

	Tot1lepSig += hOS_1lepSig_2b[0][bin-1]->GetBinContent(binIndex) ;
	Tot1lepSig += hOS_1lepSig_2b[1][bin-1]->GetBinContent(binIndex) ;
	Tot1lepSig += hOS_1lepSig_2b[2][bin-1]->GetBinContent(binIndex) ;

	Tot1lepSig += hOS_1lepSig_3b[0][bin-1]->GetBinContent(binIndex) ;
	Tot1lepSig += hOS_1lepSig_3b[1][bin-1]->GetBinContent(binIndex) ;
	Tot1lepSig += hOS_1lepSig_3b[2][bin-1]->GetBinContent(binIndex) ;


	Tot1lep += hOS_1lep_1b[0][bin-1]->GetBinContent(binIndex) ;
	Tot1lep += hOS_1lep_1b[1][bin-1]->GetBinContent(binIndex) ;
	Tot1lep += hOS_1lep_1b[2][bin-1]->GetBinContent(binIndex) ;

	Tot1lep += hOS_1lep_2b[0][bin-1]->GetBinContent(binIndex) ;
	Tot1lep += hOS_1lep_2b[1][bin-1]->GetBinContent(binIndex) ;
	Tot1lep += hOS_1lep_2b[2][bin-1]->GetBinContent(binIndex) ;

	Tot1lep += hOS_1lep_3b[0][bin-1]->GetBinContent(binIndex) ;
	Tot1lep += hOS_1lep_3b[1][bin-1]->GetBinContent(binIndex) ;
	Tot1lep += hOS_1lep_3b[2][bin-1]->GetBinContent(binIndex) ;

      }
    }

    cout << "\n\n\n Point (" << ChiMass << ",0):\n" << endl ; 
    cout << " Total 0-lep yield    = " << Tot0lep << " events," << endl ;
    cout << " Total 1-lepSig yield = " << Tot1lepSig << " events," << endl ;
    cout << " Total 1-lep yield    = " << Tot1lep << " events," << endl ;
    

    TCanvas *c0 = new TCanvas("c0","c0",1000,700);
    c0->Divide(3,3);

    c0->cd(1);
    hOSS_0lep_1b->Draw();
    c0->cd(2);
    hOSS_0lep_2b->Draw();
    c0->cd(3);
    hOSS_0lep_3b->Draw();

    c0->cd(4);
    hOSS_1lepSig_1b->Draw();
    c0->cd(5);
    hOSS_1lepSig_2b->Draw();
    c0->cd(6);
    hOSS_1lepSig_3b->Draw();

    c0->cd(7);
    hOSS_1lep_1b->Draw();
    c0->cd(8);
    hOSS_1lep_2b->Draw();
    c0->cd(9);
    hOSS_1lep_3b->Draw();
    c0->cd(0);

    TString outPlot = "cocktail_plots/OSplot_" ;
    outPlot += ChiMass ;
    outPlot += "_v" ;
    outPlot += version ;
    outPlot += ".gif" ;

    c0->SaveAs(outPlot) ;

  }

  return ;

}
