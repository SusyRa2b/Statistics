// @(#)root/roostats:$Id: DebuggingToyMCSampler.h,v 1.2 2011/11/30 13:55:58 owen Exp $
// Author: Sven Kreiss and Kyle Cranmer    June 2010
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
// Additions and modifications by Mario Pelliccioni
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_DebuggingToyMCSampler
#define ROOSTATS_DebuggingToyMCSampler

//_________________________________________________
/*
BEGIN_HTML
<p>
DebuggingToyMCSampler is an implementation of the TestStatSampler interface.
It generates Toy Monte Carlo for a given parameter point and evaluates a
TestStatistic.
</p>

<p>
For parallel runs, DebuggingToyMCSampler can be given an instance of ProofConfig
and then run in parallel using proof or proof-lite. Internally, it uses
ToyMCStudy with the RooStudyManager.
</p>
END_HTML
*/
//

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

#include <vector>
#include <sstream>

#include "RooStats/TestStatSampler.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/TestStatistic.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProofConfig.h"
#include "RooStats/ToyMCSampler.h"

#include "RooWorkspace.h"
#include "RooMsgService.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"

#include "RooDataSet.h"

namespace RooStats {

  // only used inside DebuggingToyMCSampler, ie "private" in the cxx file
  //class NuisanceParametersSampler;
  
  class DebuggingToyMCSampler: public ToyMCSampler {
    
  public:
  DebuggingToyMCSampler() :
    ToyMCSampler(),fFitToData(NULL),fDebuggingData(NULL)
      {
      }
  DebuggingToyMCSampler(TestStatistic &ts, Int_t ntoys) :
    ToyMCSampler(ts, ntoys),fFitToData(NULL),fDebuggingData(NULL)
      {
      }
    
    virtual ~DebuggingToyMCSampler();
    
    //virtual RooAbsData* GenerateToyData(RooArgSet& paramPoint) const {
    //  if(fExpectedNuisancePar) oocoutE((TObject*)NULL,InputArguments) << "ToyMCSampler: using expected nuisance parameters but ignoring weight. Use GetSamplingDistribution(paramPoint, weight) instead." << endl;
    //  double weight;
    //  return GenerateToyData(paramPoint, weight);
    //}
    using ToyMCSampler::GenerateToyData;
    RooAbsData* GenerateToyData(RooArgSet& paramPoint, double& weight, RooAbsPdf& pdf) const;
    
    void SetFitToData(RooAbsData* d) { fFitToData = d; }
    void SetDebuggingData(RooAbsData* d) { fDebuggingData = d; }
    
    RooAbsData* fFitToData ;
    RooAbsData* fDebuggingData ;
    
  protected:
    ClassDef(DebuggingToyMCSampler,3) // A simple implementation of the ToyMCSampler interface
      };
}


#endif
