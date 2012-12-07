#define reducedTree_cxx
#include "reducedTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

double reducedTree::getWeight3()
{
  fChain->GetEntry(0);
  return weight3;
}

void reducedTree::Loop(char* outfile, Bool_t doQCD, Bool_t doBlind)
{
cout << "inside loop with doQCD = " << doQCD << endl;
//   In a ROOT session, you can do:
//      Root > .L reducedTree.C
//      Root > reducedTree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;



   //First get trigger efficiency histograms from Harold
   // File is from /afs/cern.ch/user/n/nguyenh/public/RA2bTrigger2012/trigEff_coarse_bin_newTriggerScheme_ecalLaserFilter_v3.root
   TFile trigeff = TFile("trigEff_coarse_bin_newTriggerScheme_ecalLaserFilter_v3.root");
   trigeff.GetObject("hnum_ht400to500_ele_abc", hnum_ht400to500_ele);
   trigeff.GetObject("hnum_ht500to800_ele_abc", hnum_ht500to800_ele);
   trigeff.GetObject("hnum_ht800to1000_ele_abc", hnum_ht800to1000_ele);
   trigeff.GetObject("hnum_ht1000toInf_ele_abc", hnum_ht1000toInf_ele);
   trigeff.GetObject("hnum_ht400to500_mu_abc", hnum_ht400to500_mu);
   trigeff.GetObject("hnum_ht500to800_mu_abc", hnum_ht500to800_mu);
   trigeff.GetObject("hnum_ht800to1000_mu_abc", hnum_ht800to1000_mu);
   trigeff.GetObject("hnum_ht1000toInf_mu_abc", hnum_ht1000toInf_mu);
   trigeff.GetObject("hden_ht400to500_ele_abc", hden_ht400to500_ele);
   trigeff.GetObject("hden_ht500to800_ele_abc", hden_ht500to800_ele);
   trigeff.GetObject("hden_ht800to1000_ele_abc", hden_ht800to1000_ele);
   trigeff.GetObject("hden_ht1000toInf_ele_abc", hden_ht1000toInf_ele);
   trigeff.GetObject("hden_ht400to500_mu_abc", hden_ht400to500_mu);
   trigeff.GetObject("hden_ht500to800_mu_abc", hden_ht500to800_mu);
   trigeff.GetObject("hden_ht800to1000_mu_abc", hden_ht800to1000_mu);
   trigeff.GetObject("hden_ht1000toInf_mu_abc", hden_ht1000toInf_mu);
   // now 0L efficiency is fake met efficiency for QCD
   trigeff.GetObject("hnum_ht400to500_0L_abc_qcd", hnum_ht400to500_0L);
   trigeff.GetObject("hnum_ht500to800_0L_abc_qcd", hnum_ht500to800_0L);
   trigeff.GetObject("hnum_ht800to1000_0L_abc_qcd", hnum_ht800to1000_0L);
   trigeff.GetObject("hnum_ht1000toInf_0L_abc_qcd", hnum_ht1000toInf_0L);
   trigeff.GetObject("hden_ht400to500_0L_abc_qcd", hden_ht400to500_0L);
   trigeff.GetObject("hden_ht500to800_0L_abc_qcd", hden_ht500to800_0L);
   trigeff.GetObject("hden_ht800to1000_0L_abc_qcd", hden_ht800to1000_0L);
   trigeff.GetObject("hden_ht1000toInf_0L_abc_qcd", hden_ht1000toInf_0L);
   hnum_ht400to500_mu->Divide(hden_ht400to500_mu);
   hnum_ht500to800_mu->Divide(hden_ht500to800_mu);
   hnum_ht800to1000_mu->Divide(hden_ht800to1000_mu);
   hnum_ht1000toInf_mu->Divide(hden_ht1000toInf_mu);
   hnum_ht400to500_ele->Divide(hden_ht400to500_ele);
   hnum_ht500to800_ele->Divide(hden_ht500to800_ele);
   hnum_ht800to1000_ele->Divide(hden_ht800to1000_ele);
   hnum_ht1000toInf_ele->Divide(hden_ht1000toInf_ele);
   hnum_ht400to500_0L->Divide(hden_ht400to500_0L);
   hnum_ht500to800_0L->Divide(hden_ht500to800_0L);
   hnum_ht800to1000_0L->Divide(hden_ht800to1000_0L);
   hnum_ht1000toInf_0L->Divide(hden_ht1000toInf_0L);

   TFile* outFile = new TFile(outfile,"RECREATE");  
   outFile->cd();					    
   

   TTree *newtree = new TTree("tree","RA2b selected events");
   int nBCSVM50;
   float mgluino = -1;
   float mlsp = -1;
   float trigWeight = -1;
   float pfOcaloMET = -999;
   newtree->Branch("nJets",&njets,"nJets/i");
   newtree->Branch("HT",&HT,"HT/f");
   newtree->Branch("MET",&MET,"MET/f");
   newtree->Branch("nB",&nBCSVM50,"nB/i");
   newtree->Branch("run",&runNumber,"run/i");
   newtree->Branch("event",&eventNumber,"event/i");
   newtree->Branch("lumi",&lumiSection,"lumi/i");
   newtree->Branch("nPV",&nGoodPV,"nPV/i");
   newtree->Branch("nMu",&nMuons,"nMu/i");
   newtree->Branch("nEl",&nElectrons,"nEl/i");
   newtree->Branch("minDelPhiN",&minDeltaPhiN_asin,"minDelPhiN/f");
   newtree->Branch("minDelPhi",&minDeltaPhi,"minDelPhi/f");
   newtree->Branch("pt_1st_leadJet",&jetpt1,"pt_1st_leadJet/f");
   newtree->Branch("pt_2nd_leadJet",&jetpt2,"pt_2nd_leadJet/f");
   newtree->Branch("pt_3rd_leadJet",&jetpt3,"pt_3rd_leadJet/f");
   newtree->Branch("mgluino",&mgluino,"mgluino/f");
   newtree->Branch("mlsp",&mlsp,"mlsp/f");
   newtree->Branch("weightPU",&PUweight,"weightPU/f");
   newtree->Branch("weightPUSystVar",&PUweightSystVar,"weightPUSystVar/f");
   newtree->Branch("prob1", &prob1, "prob1/f"); 
   newtree->Branch("prob2", &prob2, "prob2/f"); 
   newtree->Branch("probge3", &probge3, "probge3/f"); 
   newtree->Branch("passedTrigger",&cutTrigger2,"passedTrigger/O");
   newtree->Branch("trigWeight", &trigWeight, "trigWeight/f");
   newtree->Branch("qEl", &elecharge1, "qEl/I");
   newtree->Branch("qMu", &muoncharge1, "qMu/I");
   newtree->Branch("MT", &MT_Wlep, "MT/f");
   newtree->Branch("pfOcaloMET", &pfOcaloMET, "pfOcaloMET/f");
   newtree->Branch("spher", &transverseSphericity_jets, "spher/f");
   newtree->Branch("spherMET", &transverseSphericity_jetsMet, "spherMET/f");
   newtree->Branch("spherMETLep", &transverseSphericity_jetsMetLeptons, "spherMETLep/f");
   newtree->Branch("spherLep", &transverseSphericity_jetsLeptons, "spherLep/f");
   newtree->Branch("nIsoTrk",        &nIsoTracks15_005_03,   "nIsoTrk/I");
   newtree->Branch("maxChNMultDiff", &maxTOBTECjetDeltaMult, "maxChNMultDiff/I");
   newtree->Branch("chMultTOBTECjet",&TOBTECjetChMult,       "chMultTOBTECjet/I");

   Long64_t nentries = fChain->GetEntriesFast();
   bool filled = false;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<1000;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (cutPV!=1) continue; 
      if (HT<150) continue;
      if (MET<50) continue;
      if (njets<3) continue;
      if (csctighthaloFilter!=1) continue;
      //if (eenoiseFilter!=1) continue;	    
      if (greedymuonFilter!=1) continue;       
      if (hbhenoiseFilter!=1) continue;      
      if (inconsistentmuonFilter!=1) continue;      
      if (ra2ecaltpFilter!=1) continue;      
      if (scrapingvetoFilter!=1) continue;	
      if (trackingfailureFilter!=1) continue;	  
//also filtering on hcallaser, eebadsc and PBNRcode
      if (badjetFilter!=1) continue;
      if (passCleaning!=1) continue;
      if (buggyEvent) continue; 			     
      //if (nbjets<1) continue;
      //Int_t           nElectrons;
      //Int_t	      nMuons;
      //Float_t	      minDeltaPhiN_asin;
      //Float_t         jetpt1;
      //Float_t         jetpt2;
      //Float_t         jetpt3;
      if (nbjets>=3 && bjetpt3>50) nBCSVM50 = nbjets;
      else if (nbjets>=2 && bjetpt2>50) nBCSVM50 = 2;
      else if (nbjets>=1 && bjetpt1>50) nBCSVM50 = 1;
      else nBCSVM50 = 0;
      
      if (nMuons==1&&nElectrons==0) trigWeight = getTrigWeightMu();
      if (nMuons==0&&nElectrons==1) trigWeight = getTrigWeightEl();
      if (nMuons==0&&nElectrons==0 &&!doQCD) trigWeight = getTrigWeightMu();
      if (nMuons==0&&nElectrons==0 && doQCD) trigWeight = getTrigWeight0L();
      
      mgluino = m0;
      mlsp = m12;

      pfOcaloMET = MET/caloMET;


      //lastly check for blinding
      if(doBlind) {
        // blind only high signal bins
//	if(HT>400&&MET>350&&njets>=3&&nBCSVM50>=2&&nElectrons==0&&nMuons==0&&minDeltaPhiN_asin>4) continue;   
//	if(HT>400&&MET>150&&njets>=3&&nBCSVM50>=3&&nElectrons==0&&nMuons==0&&minDeltaPhiN_asin>4) continue;
	trigWeight = 1.0;
      }
    if(!filled) {cout << "filling tree" << endl;  filled=true;}
      newtree->Fill();
   } 
  cout << "writing and closing output file" << endl;
  outFile->Write();
  outFile->Close();
}


float reducedTree::getTrigWeightMu() {

  float Mbins[5] = {125.,150.,250.,350.,99999.};    
  float Hbins[5] = {400.,500.,800.,1000.,99999.};    

  float toreturnold = -1;
  float toreturnnew = -1;

  if (MET>125 && HT>400) {
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnnew = hnum_ht400to500_mu->GetBinContent( hnum_ht400to500_mu->GetXaxis()->FindBin(MET) );    
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnnew =  hnum_ht500to800_mu->GetBinContent( hnum_ht500to800_mu->GetXaxis()->FindBin(MET) );    
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnnew =  hnum_ht800to1000_mu->GetBinContent( hnum_ht800to1000_mu->GetXaxis()->FindBin(MET) );  
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnnew =  hnum_ht1000toInf_mu->GetBinContent( hnum_ht1000toInf_mu->GetXaxis()->FindBin(MET) );	 
  }


/* Old trigger efficiencies without HTMHT650
  if ( (MET>Mbins[0]) && (MET<Mbins[1]) ){
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) return 0.962;
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) return 0.950;
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) return 0.886;
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) return 0.842;
  }
  if ( (MET>Mbins[1]) && (MET<Mbins[2]) ){
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) return 0.995;
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) return 0.991;
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) return 0.991;
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) return 0.963;
  }
  if ( (MET>Mbins[2]) && (MET<Mbins[3]) ){
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) return 1.000;
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) return 0.997;
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) return 1.000;
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) return 1.000;
  }
  if ( (MET>Mbins[3]) && (MET<Mbins[4]) ){
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) return 1.000;
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) return 1.000;
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) return 1.000;
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) return 1.000;
  }
*/
  if ( (MET>Mbins[0]) && (MET<=Mbins[1]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnold =  0.962;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnold =  0.964;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnold =  0.943;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnold =  0.947;
  }
  if ( (MET>Mbins[1]) && (MET<=Mbins[2]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnold =  0.995;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnold =  0.994;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnold =  1.000;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnold =  0.929;
  }
  if ( (MET>Mbins[2]) && (MET<=Mbins[3]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnold =  1.000;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnold =  1.000;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnold =  1.000;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnold =  1.000;
  }
  if ( (MET>Mbins[3]) && (MET<=Mbins[4]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnold =  1.000;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnold =  1.000;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnold =  1.000;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnold =  1.000;
  }

//cout << "1mu: HT = " << HT << " MET " << MET << "   old value: " << toreturnold << "  new value: " << toreturnnew << endl;
//cout << "1mu eff to return = " << toreturnnew << " from MET bin " <<  hnum_ht400to500_mu->GetXaxis()->FindBin(MET) << " with MET = " << MET << " and HT = " << HT << endl;

return toreturnnew;

  return -1.0;
}

float reducedTree::getTrigWeightEl() {
  float Mbins[5] = {125.,150.,250.,350.,99999.};    
  float Hbins[5] = {400.,500.,800.,1000.,99999.};    
  float toreturnold = -1;
  float toreturnnew = -1;

  if (MET>125 && HT>400) {
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) toreturnnew =  hnum_ht400to500_ele->GetBinContent( hnum_ht400to500_ele->GetXaxis()->FindBin(MET) );    
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) toreturnnew =  hnum_ht500to800_ele->GetBinContent( hnum_ht500to800_ele->GetXaxis()->FindBin(MET) );    
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) toreturnnew =  hnum_ht800to1000_ele->GetBinContent( hnum_ht800to1000_ele->GetXaxis()->FindBin(MET) );  
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) toreturnnew =  hnum_ht1000toInf_ele->GetBinContent( hnum_ht1000toInf_ele->GetXaxis()->FindBin(MET) );	 
  }


/* eff without PFHT650
  if ( (MET>Mbins[0]) && (MET<Mbins[1]) ){
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) return 0.792;
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) return 0.748;
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) return 0.794;
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) return 0.733;
  }
  if ( (MET>Mbins[1]) && (MET<Mbins[2]) ){
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) return 0.964;
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) return 0.953;
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) return 0.964;
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) return 0.907;
  }
  if ( (MET>Mbins[2]) && (MET<Mbins[3]) ){
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) return 1.000;
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) return 0.996;
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) return 0.976;
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) return 1.000;
  }
  if ( (MET>Mbins[3]) && (MET<Mbins[4]) ){
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) return 1.000;
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) return 1.000;
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) return 1.000;
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) return 1.000;
  }
  */
    if ( (MET>Mbins[0]) && (MET<=Mbins[1]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnold =  0.833;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnold =  0.918;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnold =  0.941;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnold =  0.933;
  }
  if ( (MET>Mbins[1]) && (MET<=Mbins[2]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnold =  0.975;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnold =  0.986;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnold =  1.000;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnold =  0.929;
  }
  if ( (MET>Mbins[2]) && (MET<=Mbins[3]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnold =  1.000;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnold =  0.996;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnold =  0.976;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnold =  1.000;
  }
  if ( (MET>Mbins[3]) && (MET<=Mbins[4]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnold =  1.000;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnold =  1.000;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnold =  1.000;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnold =  1.000;
  }

//  cout << "1el: HT = " << HT << " MET " << MET << "   old value: " << toreturnold << "  new value: " << toreturnnew << endl;
//cout << "1el eff to return = " << toreturnnew << " from MET bin " <<  hnum_ht400to500_ele->GetXaxis()->FindBin(MET) << " with MET = " << MET << " and HT = " << HT << endl;

  return toreturnnew;
  
  return -1.0;
}

float reducedTree::getTrigWeight0L() {
  float Mbins[5] = {125.,150.,250.,350.,99999.};    
  float Hbins[5] = {400.,500.,800.,1000.,99999.};    
  float toreturnold = -1;
  float toreturnnew = -1;
  float toreturnfakemet = -1;

  if (MET>125 && HT>400) {
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) toreturnnew =  hnum_ht400to500_0L->GetBinContent( hnum_ht400to500_0L->GetXaxis()->FindBin(MET) );    
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) toreturnnew =  hnum_ht500to800_0L->GetBinContent( hnum_ht500to800_0L->GetXaxis()->FindBin(MET) );    
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) {
      toreturnnew =  hnum_ht800to1000_0L->GetBinContent( hnum_ht800to1000_0L->GetXaxis()->FindBin(MET) );  
      // for coarse binning, set HT3 MET4 to 100% by hand because no events
      //if (MET>350) toreturnnew = 1.0;
      // for variable binning set two highest HT bins to 100% (file binning not good)
      if (MET>125) toreturnnew = 1.0;
    }
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) {
      toreturnnew =  hnum_ht1000toInf_0L->GetBinContent( hnum_ht1000toInf_0L->GetXaxis()->FindBin(MET) );    
      // for variable binning set two highest HT bins to 100% (file binning not good)
      if (MET>125) toreturnnew = 1.0;
    }
  }

/* eff without PFHT650
  if ( (MET>Mbins[0]) && (MET<Mbins[1]) ){
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) return 0.696;
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) return 0.750;
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) return 0.492;
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) return 0.483;
  }
  if ( (MET>Mbins[1]) && (MET<Mbins[2]) ){
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) return 0.895;
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) return 0.864;
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) return 0.751;
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) return 0.667;
  }
  if ( (MET>Mbins[2]) && (MET<Mbins[3]) ){
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) return 1.000;
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) return 1.000;
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) return 0.936;
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) return 0.867;
  }
  if ( (MET>Mbins[3]) && (MET<Mbins[4]) ){
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) return 1.000;
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) return 1.000;
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) return 0.992;
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) return 0.980;
  }
  */
  
   if ( (MET>Mbins[0]) && (MET<=Mbins[1]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnold =  0.696;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnold =  0.800;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnold =  1.000;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnold =  1.000;
  }
  if ( (MET>Mbins[1]) && (MET<=Mbins[2]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnold =  0.895;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnold =  0.864;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnold =  1.000;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnold =  1.000;
  }
  if ( (MET>Mbins[2]) && (MET<=Mbins[3]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnold =  1.000;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnold =  1.000;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnold =  1.000;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnold =  1.000;
  }
  if ( (MET>Mbins[3]) && (MET<=Mbins[4]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnold =  1.000;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnold =  1.000;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnold =  1.000;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnold =  1.000;
  }
  
//   cout << "0L: HT = " << HT << " MET " << MET << "   old value: " << toreturnold << "  new value: " << toreturnnew << endl;


  if (MET>Mbins[0] && MET<=Mbins[2] && HT>400) {
    if ( (HT>Hbins[0]) && (HT<Hbins[1]) ) toreturnfakemet =  hnum_ht400to500_0L->GetBinContent( hnum_ht400to500_0L->GetXaxis()->FindBin(MET) );    
    if ( (HT>Hbins[1]) && (HT<Hbins[2]) ) toreturnfakemet =  hnum_ht500to800_0L->GetBinContent( hnum_ht500to800_0L->GetXaxis()->FindBin(MET) );    
    if ( (HT>Hbins[2]) && (HT<Hbins[3]) ) toreturnfakemet =  hnum_ht800to1000_0L->GetBinContent( hnum_ht800to1000_0L->GetXaxis()->FindBin(MET) );  
    if ( (HT>Hbins[3]) && (HT<Hbins[4]) ) toreturnfakemet =  hnum_ht1000toInf_0L->GetBinContent( hnum_ht1000toInf_0L->GetXaxis()->FindBin(MET) );	 
  }
  if ( (MET>Mbins[2]) && (MET<=Mbins[3]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnfakemet =  1.0;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnfakemet =  1.0;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnfakemet =  1.0;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnfakemet =  1.0;
  }
  if ( (MET>Mbins[3]) && (MET<=Mbins[4]) ){
    if ( (HT>Hbins[0]) && (HT<=Hbins[1]) ) toreturnfakemet =  1.0;
    if ( (HT>Hbins[1]) && (HT<=Hbins[2]) ) toreturnfakemet =  1.0;
    if ( (HT>Hbins[2]) && (HT<=Hbins[3]) ) toreturnfakemet =  1.0;
    if ( (HT>Hbins[3]) && (HT<=Hbins[4]) ) toreturnfakemet =  1.0;
  }

//cout << "0L eff to return = " << toreturnfakemet << " from MET bin " <<  hnum_ht400to500_0L->GetXaxis()->FindBin(MET) << " with MET = " << MET << " and HT = " << HT << endl;

  return toreturnfakemet;
  
  return -1.0;
}



