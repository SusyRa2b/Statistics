#include <iostream>

#include "TTree.h"
#include "TH1D.h"

void drawCompare(TString fileName)
{

  TString branchDescriptor = "nameOAK/C:nameLB/C:valLB/D:valOAK:diff:percentDiff";
  
  TTree* fitTree = new TTree("fitTree", "fitTree");

  fitTree->ReadFile(fileName, branchDescriptor);

  TH1D* hPercentDiff = new TH1D("hPercentDiff", "Percent difference", 100, -1, 1);

    
  fitTree->Project("hPercentDiff", "percentDiff");


  hPercentDiff->Draw();

  double pdiff;
  fitTree->SetBranchAddress("percentDiff", &pdiff);

  double maxpdiff = 0;

  unsigned int ntot = fitTree->GetEntries();
  for(unsigned int i =0; i<ntot; i++)
    {
      fitTree->GetEvent(i);

      if( fabs(pdiff)>maxpdiff ) maxpdiff = fabs(pdiff);

    }
  cout << "Max diff = " << maxpdiff << endl;
  

}
