TGraph* getReferenceLimit(double scale, TH2D* theLimit)
{
  TFile *f_xsec = new TFile("reference_xSec.root");
  TH1D* referenceXSec = (TH1D*)f_xsec->Get("gluino");
  TAxis* refAxis = referenceXSec->GetXaxis();
  TString name("limitCurve");
  TGraph* theCurve = new TGraph();
  int nPoints = 0;
  TAxis* xAxis = theLimit->GetXaxis();
  TAxis* yAxis = theLimit->GetYaxis();
  int refBin = 1;
  for ( int iX = 1 ; iX <= theLimit->GetNbinsX() ; iX++ ) {
    double gluinoMass =xAxis->GetBinLowEdge(iX);
    for(;refAxis->GetBinLowEdge(refBin)<gluinoMass; refBin++);
    double refXSec = scale * referenceXSec->GetBinContent(refBin);
    for ( int iY = 1 ; iY <= theLimit->GetNbinsY() ; iY++ ) {
      double lspMass = yAxis->GetBinLowEdge(iY);
      if(lspMass+400 >= gluinoMass){
	theCurve->SetPoint(nPoints,gluinoMass,yAxis->GetBinLowEdge(iY-1));
	nPoints++;
	break;
      }
      double testXSec = theLimit->GetBinContent(iX,iY);
      if(testXSec < refXSec ) continue;
      if(iY == 1) {
	theCurve->SetPoint(nPoints,gluinoMass,lspMass);
	nPoints++;
	break;
      }
      else{
	//linear Approximation
	double lastXSec = theLimit->GetBinContent(iX,iY-1);
	double lastLspMass = yAxis->GetBinLowEdge(iY-1);
	double linLspMass = (lspMass - lastLspMass)/(testXSec-lastXSec) * (refXSec-lastXSec) + lastLspMass;
	theCurve->SetPoint(nPoints,gluinoMass,linLspMass);
	nPoints++;
	break;
      }
    }
  }
  return theCurve;
}

CombineExpectedT1tttt(char* likelihood,char* limit){

  char* thisTitle;

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".4g");
  gStyle->SetPadGridX(1) ;
  gStyle->SetPadGridY(1) ;

  const double Lumi(4674.);

  TTree tree1BL_expected("tree1BL_expected","tree1BL_expected");
  TTree tree1BT_expected("tree1BT_expected","tree1BT_expected");
  TTree tree2BL_expected("tree2BL_expected","tree2BL_expected");
  TTree tree2BT_expected("tree2BT_expected","tree2BT_expected");
  TTree  tree3B_expected( "tree3B_expected", "tree3B_expected");

  TTree tree1BL_observed("tree1BL_observed","tree1BL_observed");
  TTree tree1BT_observed("tree1BT_observed","tree1BT_observed");
  TTree tree2BL_observed("tree2BL_observed","tree2BL_observed");
  TTree tree2BT_observed("tree2BT_observed","tree2BT_observed");
  TTree  tree3B_observed( "tree3B_observed", "tree3B_observed");

  TString branchDescriptor("mGluino:mLSP:asymptoticLimit:toyLimit:toyLimit_err:expectedToyLimit:expectedToyLimit_errHigh:expectedToyLimit_errLow:CLsplusb_toyLimit:CLsplusb_toyLimit_err:CLsplusb_expectedToyLimit:CLsplusb_expectedToyLimit_errHigh:CLsplusb_expectedToyLimit_errLow:hybrid_toyLimit:hybrid_toyLimit_err:hybrid_expectedToyLimit:hybrid_expectedToyLimit_errHigh:hybrid_expectedToyLimit_errLow:hybrid_CLsplusb_toyLimit:hybrid_CLsplusb_toyLimit_err:hybrid_CLsplusb_expectedToyLimit:hybrid_CLsplusb_expectedToyLimit_errHigh:hybrid_CLsplusb_expectedToyLimit_errLow");

  tree1BL_expected.ReadFile(TString(likelihood)+TString(".expected.1BL.dat"),branchDescriptor);
  tree1BT_expected.ReadFile(TString(likelihood)+TString(".expected.1BT.dat"),branchDescriptor);
  tree2BL_expected.ReadFile(TString(likelihood)+TString(".expected.2BL.dat"),branchDescriptor);
  tree2BT_expected.ReadFile(TString(likelihood)+TString(".expected.2BT.dat"),branchDescriptor);
   tree3B_expected.ReadFile(TString(likelihood)+TString(".expected.3B.dat" ),branchDescriptor);

  tree1BL_observed.ReadFile(TString(likelihood)+TString(".measured.1BL.dat"),branchDescriptor);
  tree1BT_observed.ReadFile(TString(likelihood)+TString(".measured.1BT.dat"),branchDescriptor);
  tree2BL_observed.ReadFile(TString(likelihood)+TString(".measured.2BL.dat"),branchDescriptor);
  tree2BT_observed.ReadFile(TString(likelihood)+TString(".measured.2BT.dat"),branchDescriptor);
   tree3B_observed.ReadFile(TString(likelihood)+TString(".measured.3B.dat" ),branchDescriptor);

  TH2D *h2d_1BL_expected = new TH2D("h2d_1BL_expected","1BL, expected UL",31,450,1225,31,50,825);
  TH2D *h2d_1BT_expected = new TH2D("h2d_1BT_expected","1BT, expected UL",31,450,1225,31,50,825);
  TH2D *h2d_2BL_expected = new TH2D("h2d_2BL_expected","2BL, expected UL",31,450,1225,31,50,825);
  TH2D *h2d_2BT_expected = new TH2D("h2d_2BT_expected","2BT, expected UL",31,450,1225,31,50,825);
  TH2D  *h2d_3B_expected = new TH2D( "h2d_3B_expected", "3B, expected UL",31,450,1225,31,50,825);

  tree1BL_expected.Draw("mLSP:mGluino>>+h2d_1BL_expected",limit,"goff");
  tree1BT_expected.Draw("mLSP:mGluino>>+h2d_1BT_expected",limit,"goff");
  tree2BL_expected.Draw("mLSP:mGluino>>+h2d_2BL_expected",limit,"goff");
  tree2BT_expected.Draw("mLSP:mGluino>>+h2d_2BT_expected",limit,"goff");
   tree3B_expected.Draw( "mLSP:mGluino>>+h2d_3B_expected",limit,"goff");

  TH2D *h2d_1BL_observed = new TH2D("h2d_1BL_observed","1BL, observed UL",31,450,1225,31,50,825);
  TH2D *h2d_1BT_observed = new TH2D("h2d_1BT_observed","1BT, observed UL",31,450,1225,31,50,825);
  TH2D *h2d_2BL_observed = new TH2D("h2d_2BL_observed","2BL, observed UL",31,450,1225,31,50,825);
  TH2D *h2d_2BT_observed = new TH2D("h2d_2BT_observed","2BT, observed UL",31,450,1225,31,50,825);
  TH2D  *h2d_3B_observed = new TH2D( "h2d_3B_observed", "3B, observed UL",31,450,1225,31,50,825);

  tree1BL_observed.Draw("mLSP:mGluino>>+h2d_1BL_observed",limit,"goff");
  tree1BT_observed.Draw("mLSP:mGluino>>+h2d_1BT_observed",limit,"goff");
  tree2BL_observed.Draw("mLSP:mGluino>>+h2d_2BL_observed",limit,"goff");
  tree2BT_observed.Draw("mLSP:mGluino>>+h2d_2BT_observed",limit,"goff");
   tree3B_observed.Draw( "mLSP:mGluino>>+h2d_3B_observed",limit,"goff");

  TCanvas *c0 = new TCanvas("c0", "",30,98,1100,500);

  c0->SetBottomMargin(0.1);
  c0->SetBorderSize(2);
  c0->SetFrameFillColor(0);
  c0->SetFillColor(10) ;

  // draw histograms

  h2d_1BL_expected->SetMaximum(600);
  h2d_1BT_expected->SetMaximum(20);
  h2d_2BL_expected->SetMaximum(240);
  h2d_2BT_expected->SetMaximum(100);  
  h2d_3B_expected->SetMaximum(100);  

  h2d_1BL_expected->SetMinimum(0);
  h2d_1BT_expected->SetMinimum(0);
  h2d_2BL_expected->SetMinimum(0);
  h2d_2BT_expected->SetMinimum(0);  
  h2d_3B_expected->SetMinimum(0);  

  h2d_1BL_expected->SetTitle("1BL - expected signal yield UL");
  h2d_1BL_expected->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_1BL_expected->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_1BL_expected->Draw("textcolz");
  c0->cd(0);
  TString plotTitle("plots/1BL_nevtUL_expected.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 
  
  h2d_1BT_expected->SetTitle("1BT - expected signal yield UL");
  h2d_1BT_expected->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_1BT_expected->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_1BT_expected->Draw("textcolz");
  c0->cd(0);
  TString plotTitle("plots/1BT_nevtUL_expected.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 


  h2d_2BL_expected->SetTitle("2BL - expected signal yield UL");
  h2d_2BL_expected->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_2BL_expected->GetYaxis()->SetTitle("LSP mass (GeV)");  
  h2d_2BL_expected->Draw("textcolz");
  c0->cd(0);
  TString plotTitle("plots/2BL_nevtUL_expected.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 


  h2d_2BT_expected->SetTitle("2BT - expected signal yield UL");
  h2d_2BT_expected->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_2BT_expected->GetYaxis()->SetTitle("LSP mass (GeV)");  
  h2d_2BT_expected->Draw("textcolz");
  c0->cd(0);
  TString plotTitle("plots/2BT_nevtUL_expected.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 


  h2d_3B_expected->SetTitle("3B - expected signal yield UL");
  h2d_3B_expected->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_3B_expected->GetYaxis()->SetTitle("LSP mass (GeV)");  
  h2d_3B_expected->Draw("textcolz");
  c0->cd(0);
  TString plotTitle("plots/3B_nevtUL_expected.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  h2d_1BL_observed->SetMaximum(600);
  h2d_1BT_observed->SetMaximum(20);
  h2d_2BL_observed->SetMaximum(240);
  h2d_2BT_observed->SetMaximum(100);  
  h2d_3B_observed->SetMaximum(100);  

  h2d_1BL_observed->SetMinimum(0);
  h2d_1BT_observed->SetMinimum(0);
  h2d_2BL_observed->SetMinimum(0);
  h2d_2BT_observed->SetMinimum(0);  
  h2d_3B_observed->SetMinimum(0);  

  h2d_1BL_observed->SetTitle("1BL - observed signal yield UL");
  h2d_1BL_observed->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_1BL_observed->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_1BL_observed->Draw("textcolz");
  c0->cd(0);
  TString plotTitle("plots/1BL_nevtUL_observed.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 
  
  h2d_1BT_observed->SetTitle("1BT - observed signal yield UL");
  h2d_1BT_observed->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_1BT_observed->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_1BT_observed->Draw("textcolz");
  c0->cd(0);
  TString plotTitle("plots/1BT_nevtUL_observed.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 


  h2d_2BL_observed->SetTitle("2BL - observed signal yield UL");
  h2d_2BL_observed->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_2BL_observed->GetYaxis()->SetTitle("LSP mass (GeV)");  
  h2d_2BL_observed->Draw("textcolz");
  c0->cd(0);
  TString plotTitle("plots/2BL_nevtUL_observed.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 


  h2d_2BT_observed->SetTitle("2BT - observed signal yield UL");
  h2d_2BT_observed->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_2BT_observed->GetYaxis()->SetTitle("LSP mass (GeV)");  
  h2d_2BT_observed->Draw("textcolz");
  c0->cd(0);
  TString plotTitle("plots/2BT_nevtUL_observed.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 


  h2d_3B_observed->SetTitle("3B - observed signal yield UL");
  h2d_3B_observed->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_3B_observed->GetYaxis()->SetTitle("LSP mass (GeV)");  
  h2d_3B_observed->Draw("textcolz");
  c0->cd(0);
  TString plotTitle("plots/3B_nevtUL_observed.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  // now compute the xs UL's

  TTree tree1BL_efficiency("tree1BL_efficiency","tree1BL_efficiency");
  TTree tree1BT_efficiency("tree1BT_efficiency","tree1BT_efficiency");
  TTree tree2BL_efficiency("tree2BL_efficiency","tree2BL_efficiency");
  TTree tree2BT_efficiency("tree2BT_efficiency","tree2BT_efficiency");
  TTree  tree3B_efficiency( "tree3B_efficiency", "tree3B_efficiency");
  
  TString efficiencyBranchDescriptor("mGluino:mLSP:nEntries:sig:sb:sigsl:sbsl:sigldp:sbldp:eps_sig:eps_sb:eps_sigsl:eps_sbsl:eps_sigldp:eps_sbldp:syst_sig:syst_sb:syst_sigsl:syst_sbsl:syst_sigldp:syst_sbldp");
  
  tree1BL_efficiency.ReadFile("/afs/cern.ch/user/j/joshmt/public/RA2bFall2011/signalSyst.T1tttt.ge1bLoose.dat",efficiencyBranchDescriptor);
  tree1BT_efficiency.ReadFile("/afs/cern.ch/user/j/joshmt/public/RA2bFall2011/signalSyst.T1tttt.ge1bTight.dat",efficiencyBranchDescriptor);
  tree2BL_efficiency.ReadFile("/afs/cern.ch/user/j/joshmt/public/RA2bFall2011/signalSyst.T1tttt.ge2bLoose.dat",efficiencyBranchDescriptor);
  tree2BT_efficiency.ReadFile("/afs/cern.ch/user/j/joshmt/public/RA2bFall2011/signalSyst.T1tttt.ge2bTight.dat",efficiencyBranchDescriptor);
   tree3B_efficiency.ReadFile("/afs/cern.ch/user/j/joshmt/public/RA2bFall2011/signalSyst.T1tttt.ge3bLoose.dat",efficiencyBranchDescriptor);
  
  TH2D *eps_1BL = new TH2D("eps_1BL","efficiency - 1BL",31,450,1225,31,50,825);
  TH2D *eps_1BT = new TH2D("eps_1BT","efficiency - 1BT",31,450,1225,31,50,825);
  TH2D *eps_2BL = new TH2D("eps_2BL","efficiency - 2BL",31,450,1225,31,50,825);
  TH2D *eps_2BT = new TH2D("eps_2BT","efficiency - 2BT",31,450,1225,31,50,825);
  TH2D *eps_3B  = new TH2D("eps_3B" ,"efficiency - 3B", 31,450,1225,31,50,825);
  
  tree1BL_efficiency.Draw("mLSP:mGluino>>+eps_1BL","sig*eps_sig/nEntries","goff");
  tree1BT_efficiency.Draw("mLSP:mGluino>>+eps_1BT","sig*eps_sig/nEntries","goff");
  tree2BL_efficiency.Draw("mLSP:mGluino>>+eps_2BL","sig*eps_sig/nEntries","goff");
  tree2BT_efficiency.Draw("mLSP:mGluino>>+eps_2BT","sig*eps_sig/nEntries","goff");
   tree3B_efficiency.Draw("mLSP:mGluino>>+eps_3B" ,"sig*eps_sig/nEntries","goff");

  TH2D *h2d_1BL_xsUL_expected = new TH2D("h2d_1BL_xsUL_expected","1BL, expected xs UL (pb)",31,450,1225,31,50,825);
  TH2D *h2d_1BT_xsUL_expected = new TH2D("h2d_1BT_xsUL_expected","1BT, expected xs UL (pb)",31,450,1225,31,50,825);
  TH2D *h2d_2BL_xsUL_expected = new TH2D("h2d_2BL_xsUL_expected","2BL, expected xs UL (pb)",31,450,1225,31,50,825);
  TH2D *h2d_2BT_xsUL_expected = new TH2D("h2d_2BT_xsUL_expected","2BT, expected xs UL (pb)",31,450,1225,31,50,825);
  TH2D  *h2d_3B_xsUL_expected = new TH2D( "h2d_3B_xsUL_expected", "3B, expected xs UL (pb)",31,450,1225,31,50,825);

  TH2D *h2d_1BL_xsUL_observed = new TH2D("h2d_1BL_xsUL_observed","1BL, observed xs UL (pb)",31,450,1225,31,50,825);
  TH2D *h2d_1BT_xsUL_observed = new TH2D("h2d_1BT_xsUL_observed","1BT, observed xs UL (pb)",31,450,1225,31,50,825);
  TH2D *h2d_2BL_xsUL_observed = new TH2D("h2d_2BL_xsUL_observed","2BL, observed xs UL (pb)",31,450,1225,31,50,825);
  TH2D *h2d_2BT_xsUL_observed = new TH2D("h2d_2BT_xsUL_observed","2BT, observed xs UL (pb)",31,450,1225,31,50,825);
  TH2D  *h2d_3B_xsUL_observed = new TH2D( "h2d_3B_xsUL_observed", "3B, observed xs UL (pb)",31,450,1225,31,50,825);

  TH2D *h2d_bestSel  = new TH2D("h2d_bestSel" ,"Best selection - expected from DD bkg",31,450,1225,31,50,825);
  TH2D *h2d_bestUL   = new TH2D("h2d_bestUL"  ,"Best expected xs UL (pb)",31,450,1225,31,50,825);
  TH2D *h2d_bestUL_observed   = new TH2D("h2d_bestUL_observed" ,"observed xs UL (pb) using best expected selection",31,450,1225,31,50,825);
  TH2D *h2d_bestSel_observed  = new TH2D("h2d_bestSel_observed" ,"Best selection - observed",31,450,1225,31,50,825);
  TH2D *h2d_bestUL_observed_noexpected   = new TH2D("h2d_bestUL_observed_noexpected" ,"observed xs UL (pb) using best observed selection",31,450,1225,31,50,825);
  
  //TFile *f_eps = new TFile("SavedEfficiencies.root");
  //eps_1BL = (TH2D*)f_eps->Get("eps_1BL");
  //eps_1BT = (TH2D*)f_eps->Get("eps_1BT");
  //eps_2BL = (TH2D*)f_eps->Get("eps_2BL");
  //eps_2BT = (TH2D*)f_eps->Get("eps_2BT");
  //eps_3B  = (TH2D*)f_eps->Get("eps_3B");

  double NevtUL(0.), eps(1.), xsUL(0.);

  for ( int i = 1 ; i <= 31 ; i++ ) {
    for ( int j = 1 ; j <= 31 ; j++ ) {

      if ( j > i ) continue ;
      
      // 1b-loose
      NevtUL = h2d_1BL_expected->GetBinContent(i,j);
      eps = eps_1BL->GetBinContent(i,j);
      if (eps > 0 ) {
	xsUL = NevtUL/(eps*Lumi);
	h2d_1BL_xsUL_expected->SetBinContent(i,j,xsUL);
      }

      // 1b-tight
      NevtUL = h2d_1BT_expected->GetBinContent(i,j);
      eps = eps_1BT->GetBinContent(i,j);
      if (eps > 0 ) {
	xsUL = NevtUL/(eps*Lumi);
	h2d_1BT_xsUL_expected->SetBinContent(i,j,xsUL);
      }      

      // 2b-loose
      NevtUL = h2d_2BL_expected->GetBinContent(i,j);
      eps = eps_2BL->GetBinContent(i,j);
      if (eps > 0 ) {
	xsUL = NevtUL/(eps*Lumi);
	h2d_2BL_xsUL_expected->SetBinContent(i,j,xsUL);
      }

      // 2b-tight
      NevtUL = h2d_2BT_expected->GetBinContent(i,j);
      eps = eps_2BT->GetBinContent(i,j);
      if (eps > 0 ) {
	xsUL = NevtUL/(eps*Lumi); 
	h2d_2BT_xsUL_expected->SetBinContent(i,j,xsUL);
      }

      // 3b
      NevtUL = h2d_3B_expected->GetBinContent(i,j);
      eps = eps_3B->GetBinContent(i,j);
      if (eps > 0 /*&& NevtUL > 5*/) {       // remove points with unrealistic Nevt UL
	xsUL = NevtUL/(eps*Lumi); 
	h2d_3B_xsUL_expected->SetBinContent(i,j,xsUL);
      }
    }
  }

  for ( int i = 1 ; i <= 31 ; i++ ) {
    for ( int j = 1 ; j <= 31 ; j++ ) {

      if ( j > i ) continue ;
      
      // 1b-loose
      NevtUL = h2d_1BL_observed->GetBinContent(i,j);
      eps = eps_1BL->GetBinContent(i,j);
      if (eps > 0 ) {
	xsUL = NevtUL/(eps*Lumi);
	h2d_1BL_xsUL_observed->SetBinContent(i,j,xsUL);
      }

      // 1b-tight
      NevtUL = h2d_1BT_observed->GetBinContent(i,j);
      eps = eps_1BT->GetBinContent(i,j);
      if (eps > 0 ) {
	xsUL = NevtUL/(eps*Lumi);
	h2d_1BT_xsUL_observed->SetBinContent(i,j,xsUL);
      }      

      // 2b-loose
      NevtUL = h2d_2BL_observed->GetBinContent(i,j);
      eps = eps_2BL->GetBinContent(i,j);
      if (eps > 0 ) {
	xsUL = NevtUL/(eps*Lumi);
	h2d_2BL_xsUL_observed->SetBinContent(i,j,xsUL);
      }

      // 2b-tight
      NevtUL = h2d_2BT_observed->GetBinContent(i,j);
      eps = eps_2BT->GetBinContent(i,j);
      if (eps > 0 ) {
	xsUL = NevtUL/(eps*Lumi); 
	h2d_2BT_xsUL_observed->SetBinContent(i,j,xsUL);
      }

      // 3b
      NevtUL = h2d_3B_observed->GetBinContent(i,j);
      eps = eps_3B->GetBinContent(i,j);
      if (eps > 0 /*&& NevtUL > 5*/) {       // remove points with unrealistic Nevt UL
	xsUL = NevtUL/(eps*Lumi); 
	h2d_3B_xsUL_observed->SetBinContent(i,j,xsUL);
      }
    }
  }


  // one more (separate) loop to find the best selection

  Int_t whichSel(-9);
  Double_t bestUL(-9.);
  Double_t bestUL_observed(0.);

  for ( int i = 1 ; i <= 31 ; i++ ) {
    for ( int j = 1 ; j <= 31 ; j++ ) {

      if ( j > i ) continue ;

      bestUL = 9999;

      whichSel = 1;
      if ( h2d_1BL_xsUL_expected->GetBinContent(i,j) > 0 ) bestUL = h2d_1BL_xsUL_expected->GetBinContent(i,j);
      
      if ( h2d_1BT_xsUL_expected->GetBinContent(i,j) > 0 && h2d_1BT_xsUL_expected->GetBinContent(i,j) < bestUL ) { 
	whichSel = 2;
	bestUL = h2d_1BT_xsUL_expected->GetBinContent(i,j);
	bestUL_observed = h2d_1BT_xsUL_observed->GetBinContent(i,j);
      }
      
      if ( h2d_2BL_xsUL_expected->GetBinContent(i,j) > 0 && h2d_2BL_xsUL_expected->GetBinContent(i,j) < bestUL ) { 
	whichSel = 3;
	bestUL = h2d_2BL_xsUL_expected->GetBinContent(i,j);
	bestUL_observed = h2d_2BL_xsUL_observed->GetBinContent(i,j);
      }
      
      if ( h2d_2BT_xsUL_expected->GetBinContent(i,j) > 0 && h2d_2BT_xsUL_expected->GetBinContent(i,j) < bestUL ) { 
	whichSel = 4;
	bestUL = h2d_2BT_xsUL_expected->GetBinContent(i,j);
	bestUL_observed = h2d_2BT_xsUL_observed->GetBinContent(i,j);
      }
      
      if ( h2d_3B_xsUL_expected->GetBinContent(i,j) > 0 && h2d_3B_xsUL_expected->GetBinContent(i,j) < bestUL ) { 
	whichSel = 5;
	bestUL = h2d_3B_xsUL_expected->GetBinContent(i,j);
	bestUL_observed = h2d_3B_xsUL_observed->GetBinContent(i,j);
      }

      // fill "best" histograms

      h2d_bestSel->SetBinContent(i,j,whichSel);
      h2d_bestUL->SetBinContent(i,j,bestUL);
      h2d_bestUL_observed->SetBinContent(i,j,bestUL_observed);

    }
  }

  for ( int i = 1 ; i <= 31 ; i++ ) {
    for ( int j = 1 ; j <= 31 ; j++ ) {

      if ( j > i ) continue ;

      bestUL = 9999;

      whichSel = 1;
      if ( h2d_1BL_xsUL_observed->GetBinContent(i,j) > 0 ) bestUL = h2d_1BL_xsUL_observed->GetBinContent(i,j);
      
      if ( h2d_1BT_xsUL_observed->GetBinContent(i,j) > 0 && h2d_1BT_xsUL_observed->GetBinContent(i,j) < bestUL ) { 
	whichSel = 2;
	bestUL = h2d_1BT_xsUL_observed->GetBinContent(i,j);
      }
      
      if ( h2d_2BL_xsUL_observed->GetBinContent(i,j) > 0 && h2d_2BL_xsUL_observed->GetBinContent(i,j) < bestUL ) { 
	whichSel = 3;
	bestUL = h2d_2BL_xsUL_observed->GetBinContent(i,j);
      }
      
      if ( h2d_2BT_xsUL_observed->GetBinContent(i,j) > 0 && h2d_2BT_xsUL_observed->GetBinContent(i,j) < bestUL ) { 
	whichSel = 4;
	bestUL = h2d_2BT_xsUL_observed->GetBinContent(i,j);
      }
      
      if ( h2d_3B_xsUL_observed->GetBinContent(i,j) > 0 && h2d_3B_xsUL_observed->GetBinContent(i,j) < bestUL ) { 
	whichSel = 5;
	bestUL = h2d_3B_xsUL_observed->GetBinContent(i,j);
      }

      // fill "best" histograms

      h2d_bestSel_observed->SetBinContent(i,j,whichSel);
      h2d_bestUL_observed_noexpected->SetBinContent(i,j,bestUL);

    }
  }
 
  TGraph* refCurve;
  TGraph* thirdRefCurve;
  TGraph* threeRefCurve;

  h2d_1BL_xsUL_expected->SetMaximum(150.);
  h2d_1BT_xsUL_expected->SetMaximum(150.);
  h2d_2BL_xsUL_expected->SetMaximum(150.);
  h2d_2BT_xsUL_expected->SetMaximum(150.);
  h2d_3B_xsUL_expected->SetMaximum(150.);
  h2d_bestUL->SetMaximum(60);

  h2d_1BL_xsUL_expected->SetMinimum(0.005);
  h2d_1BT_xsUL_expected->SetMinimum(0.005);
  h2d_2BL_xsUL_expected->SetMinimum(0.005);
  h2d_2BT_xsUL_expected->SetMinimum(0.005);
  h2d_3B_xsUL_expected->SetMinimum(0.005);
  h2d_bestUL->SetMinimum(0.005);

  c0->SetLogz(1);

  h2d_1BL_xsUL_expected->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_1BL_xsUL_expected->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_1BL_xsUL_expected->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_1BL_xsUL_expected );
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_1BL_xsUL_expected );
  threeRefCurve = getReferenceLimit(3.0 , h2d_1BL_xsUL_expected );
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  TString plotTitle("plots/1BL_xsUL_expected.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve;

  h2d_1BT_xsUL_expected->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_1BT_xsUL_expected->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_1BT_xsUL_expected->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_1BT_xsUL_expected );
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_1BT_xsUL_expected );
  threeRefCurve = getReferenceLimit(3.0 , h2d_1BT_xsUL_expected );
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  plotTitle = TString("plots/1BT_xsUL_expected.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle);

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve; 

  h2d_2BL_xsUL_expected->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_2BL_xsUL_expected->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_2BL_xsUL_expected->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_2BL_xsUL_expected );
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_2BL_xsUL_expected );
  threeRefCurve = getReferenceLimit(3.0 , h2d_2BL_xsUL_expected );
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  plotTitle = TString("plots/2BL_xsUL_expected.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve;

  h2d_2BT_xsUL_expected->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_2BT_xsUL_expected->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_2BT_xsUL_expected->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_2BT_xsUL_expected );
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_2BT_xsUL_expected );
  threeRefCurve = getReferenceLimit(3.0 , h2d_2BT_xsUL_expected );
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  plotTitle = TString("plots/2BT_xsUL_expected.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve;

  h2d_3B_xsUL_expected->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_3B_xsUL_expected->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_3B_xsUL_expected->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_3B_xsUL_expected );
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_3B_xsUL_expected );
  threeRefCurve = getReferenceLimit(3.0 , h2d_3B_xsUL_expected );
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  plotTitle = TString("plots/3B_xsUL_expected.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve;

  h2d_1BL_xsUL_observed->SetMaximum(150.);
  h2d_1BT_xsUL_observed->SetMaximum(150.);
  h2d_2BL_xsUL_observed->SetMaximum(150.);
  h2d_2BT_xsUL_observed->SetMaximum(150.);
  h2d_3B_xsUL_observed->SetMaximum(150.);
  h2d_bestUL_observed->SetMaximum(60);
  h2d_bestUL_observed_noexpected->SetMaximum(60);

  h2d_1BL_xsUL_observed->SetMinimum(0.005);
  h2d_1BT_xsUL_observed->SetMinimum(0.005);
  h2d_2BL_xsUL_observed->SetMinimum(0.005);
  h2d_2BT_xsUL_observed->SetMinimum(0.005);
  h2d_3B_xsUL_observed->SetMinimum(0.005);
  h2d_bestUL_observed->SetMinimum(0.005);
  h2d_bestUL_observed_noexpected->SetMinimum(0.005);

  h2d_1BL_xsUL_observed->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_1BL_xsUL_observed->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_1BL_xsUL_observed->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_1BL_xsUL_observed );
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_1BL_xsUL_observed );
  threeRefCurve = getReferenceLimit(3.0 , h2d_1BL_xsUL_observed );
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  TString plotTitle("plots/1BL_xsUL_observed.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve;

  h2d_1BT_xsUL_observed->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_1BT_xsUL_observed->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_1BT_xsUL_observed->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_1BT_xsUL_observed );
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_1BT_xsUL_observed );
  threeRefCurve = getReferenceLimit(3.0 , h2d_1BT_xsUL_observed );
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  plotTitle = TString("plots/1BT_xsUL_observed.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve;

  h2d_2BL_xsUL_observed->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_2BL_xsUL_observed->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_2BL_xsUL_observed->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_2BL_xsUL_observed );
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_2BL_xsUL_observed );
  threeRefCurve = getReferenceLimit(3.0 , h2d_2BL_xsUL_observed );
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  plotTitle = TString("plots/2BL_xsUL_observed.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve;

  h2d_2BT_xsUL_observed->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_2BT_xsUL_observed->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_2BT_xsUL_observed->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_2BT_xsUL_observed );
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_2BT_xsUL_observed );
  threeRefCurve = getReferenceLimit(3.0 , h2d_2BT_xsUL_observed );
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  plotTitle = TString("plots/2BT_xsUL_observed.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve;

  h2d_3B_xsUL_observed->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_3B_xsUL_observed->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_3B_xsUL_observed->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_3B_xsUL_observed );
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_3B_xsUL_observed );
  threeRefCurve = getReferenceLimit(3.0 , h2d_3B_xsUL_observed );
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  plotTitle = TString("plots/3B_xsUL_observed.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve;

  // "best" histograms

  h2d_bestUL->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_bestUL->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_bestUL->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_bestUL);
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_bestUL);
  threeRefCurve = getReferenceLimit(3.0 , h2d_bestUL);
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  plotTitle = TString("plots/bestUL_expected.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve;

  h2d_bestUL_observed->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_bestUL_observed->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_bestUL_observed->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_bestUL_observed);
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_bestUL_observed);
  threeRefCurve = getReferenceLimit(3.0 , h2d_bestUL_observed);
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  plotTitle = TString("plots/bestUL_observed.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve;

  h2d_bestUL_observed_noexpected->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_bestUL_observed_noexpected->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_bestUL_observed_noexpected->Draw("colztext");

  refCurve = getReferenceLimit(1.0 , h2d_bestUL_observed_noexpected);
  thirdRefCurve = getReferenceLimit(1.0/3.0 , h2d_bestUL_observed_noexpected);
  threeRefCurve = getReferenceLimit(3.0 , h2d_bestUL_observed_noexpected);
  refCurve->SetLineWidth(2);
  thirdRefCurve->SetLineWidth(2);
  threeRefCurve->SetLineWidth(2);
  thirdRefCurve->SetLineStyle(4);
  threeRefCurve->SetLineStyle(5);
  refCurve->Draw("C");
  thirdRefCurve->Draw("C");
  threeRefCurve->Draw("C");

  c0->cd(0);
  plotTitle = TString("plots/bestUL_observed_noexpected.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  delete refCurve;
  delete thirdRefCurve;
  delete threeRefCurve;

  c0->SetLogz(0);

  h2d_bestSel->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_bestSel->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_bestSel->Draw("colztext");
  c0->cd(0);
  plotTitle = TString("plots/bestSel.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 

  h2d_bestSel_observed->GetXaxis()->SetTitle("gluino mass (GeV)");
  h2d_bestSel_observed->GetYaxis()->SetTitle("LSP mass (GeV)");
  h2d_bestSel_observed->Draw("colztext");
  c0->cd(0);
  plotTitle = TString("plots/bestSel_observed.");
  plotTitle+=likelihood;   plotTitle+=".";   plotTitle+=limit;
  plotTitle+=".pdf";
  thisTitle = plotTitle.Data() ; c0->SaveAs(thisTitle); 


  // save all relevant histograms to file

  plotTitle = TString(likelihood); plotTitle+=".";   plotTitle+=limit; plotTitle+="."; 
  plotTitle+="xsULs.root";
  thisTitle = plotTitle.Data();

  TFile *f_xsUL = new TFile(thisTitle,"RECREATE");
  f_xsUL->cd();

  h2d_1BL_xsUL_expected->Write();
  h2d_1BT_xsUL_expected->Write();
  h2d_2BL_xsUL_expected->Write();
  h2d_2BT_xsUL_expected->Write();
  h2d_3B_xsUL_expected->Write();

  h2d_1BL_xsUL_observed->Write();
  h2d_1BT_xsUL_observed->Write();
  h2d_2BL_xsUL_observed->Write();
  h2d_2BT_xsUL_observed->Write();
  h2d_3B_xsUL_observed->Write();

  h2d_bestSel->Write();
  h2d_bestUL->Write();
  h2d_bestUL_observed->Write();

  h2d_bestSel_observed->Write();
  h2d_bestUL_observed_noexpected->Write();


  f_xsUL->Close();
  
}
