#define SmallTree_cxx
#include "SmallTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void SmallTree::Loop(TH2F *histo, int si, int k)
{
using namespace std;
//   In a ROOT session, you can do:
//      Root > .L SmallTree.C
//      Root > SmallTree t
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

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

float weightTree = fChain->GetWeight();
//cout << "tree weight = " << weightTree << endl;
//cout << "nentries = " << nentries << endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<50;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) { cout << "breaking from no tree" << endl; break;}
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      if ( !(nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)) ) continue;
      if(si==0 && !(minDelPhiN>4&&nMu==0&&nEl==0)) continue;
      if(si==1 && !(minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )) ) continue;
      if(si==2 && !(minDelPhiN<4&&nMu==0&&nEl==0 ) ) continue;
      if(k==0 && !(nB==1)) continue;
      if(k==1 && !(nB==2)) continue;
      if(k==2 && !(nB>=3)) continue;
      
      histo->Fill(MET,HT,weightTree*weightPU);
      //histo->Fill(MET,HT,weightTree);

   }
}
