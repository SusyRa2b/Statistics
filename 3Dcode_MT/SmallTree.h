//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 13 15:46:57 2012 by ROOT version 5.30/02
// from TTree tree/RA2b selected events
// found on file: files15fb_8TeV/TT.root
//////////////////////////////////////////////////////////

#ifndef SmallTree_h
#define SmallTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class SmallTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          nJets;
   Float_t         HT;
   Float_t         MET;
   UInt_t          nB;
   UInt_t          run;
   UInt_t          event;
   UInt_t          lumi;
   UInt_t          nMu;
   UInt_t          nEl;
   Float_t         minDelPhiN;
   Float_t         pt_1st_leadJet;
   Float_t         pt_2nd_leadJet;
   Float_t         pt_3rd_leadJet;
   Float_t         mgluino;
   Float_t         mlsp;
   Float_t         weightPU;
   Float_t         MT;

   // List of branches
   TBranch        *b_nJets;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_nB;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_nEl;   //!
   TBranch        *b_minDelPhiN;   //!
   TBranch        *b_pt_1st_leadJet;   //!
   TBranch        *b_pt_2nd_leadJet;   //!
   TBranch        *b_pt_3rd_leadJet;   //!
   TBranch        *b_mgluino;   //!
   TBranch        *b_mlsp;   //!
   TBranch        *b_weightPU;   //!
   TBranch        *b_MT;   //!


   SmallTree(TTree *tree=0);
   virtual ~SmallTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TH2F *histo, int si, int k);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SmallTree_cxx
SmallTree::SmallTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
   cout << "have tree = 0" << endl;
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("files15fb_8TeV/TT.root");
      if (!f || !f->IsOpen()) {
      cout << "getting tree from SmallTree.h" << endl;
         f = new TFile("files15fb_8TeV/TT.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

SmallTree::~SmallTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SmallTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SmallTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void SmallTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("nB", &nB, &b_nB);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("nEl", &nEl, &b_nEl);
   fChain->SetBranchAddress("minDelPhiN", &minDelPhiN, &b_minDelPhiN);
   fChain->SetBranchAddress("pt_1st_leadJet", &pt_1st_leadJet, &b_pt_1st_leadJet);
   fChain->SetBranchAddress("pt_2nd_leadJet", &pt_2nd_leadJet, &b_pt_2nd_leadJet);
   fChain->SetBranchAddress("pt_3rd_leadJet", &pt_3rd_leadJet, &b_pt_3rd_leadJet);
   fChain->SetBranchAddress("mgluino", &mgluino, &b_mgluino);
   fChain->SetBranchAddress("mlsp", &mlsp, &b_mlsp);
   fChain->SetBranchAddress("weightPU", &weightPU, &b_weightPU);
   fChain->SetBranchAddress("MT", &MT, &b_MT);
   Notify();
}

Bool_t SmallTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SmallTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SmallTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   if ( entry > 0 ) return 1 ;
   return 1;
}
#endif // #ifdef SmallTree_cxx
