// @(#)root/roostats:$Id: HybridToyMCSampler.h,v 1.1 2011/12/16 11:08:25 winstrom Exp $
// Author: Sven Kreiss and Kyle Cranmer    June 2010
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
// Additions and modifications by Lucas Winstrom and Mario Pelliccioni
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_HybridToyMCSampler
#define ROOSTATS_HybridToyMCSampler

//_________________________________________________
/*
BEGIN_HTML
<p>
HybridToyMCSampler is an implementation of the TestStatSampler interface.
It generates Toy Monte Carlo for a given parameter point and evaluates a
TestStatistic.
</p>

<p>
For parallel runs, HybridToyMCSampler can be given an instance of ProofConfig
and then run in parallel using proof or proof-lite. Internally, it uses
ToyMCStudy with the RooStudyManager.
</p>
END_HTML
*/
//

//#ifndef ROOT_Rtypes
//#include "Rtypes.h"
//#endif
//
//#include <vector>
//#include <sstream>

//#include "RooStats/TestStatSampler.h"
//#include "RooStats/SamplingDistribution.h"
//#include "RooStats/TestStatistic.h"
//#include "RooStats/ModelConfig.h"
//#include "RooStats/ProofConfig.h"
#include "RooStats/ToyMCSampler.h"

//#include "RooWorkspace.h"
//#include "RooMsgService.h"
//#include "RooAbsPdf.h"
//#include "RooRealVar.h"
//
//#include "RooDataSet.h"

namespace RooStats {

// only used inside HybridToyMCSampler, ie "private" in the cxx file
class NuisanceParametersSamplerAndFitter;

class HybridToyMCSampler: public ToyMCSampler {

   public:
      HybridToyMCSampler() :
      ToyMCSampler(),fFitToData(NULL),fDebuggingData(NULL),fSampledNuisances(NULL)
      {
	fNuisanceParametersSamplerAndFitter = NULL;
      }
      HybridToyMCSampler(TestStatistic &ts, Int_t ntoys) :
      ToyMCSampler(ts, ntoys),fFitToData(NULL),fDebuggingData(NULL),fSampledNuisances(NULL)
      {
	fNuisanceParametersSamplerAndFitter = NULL;
      }

      virtual ~HybridToyMCSampler();

      //virtual SamplingDistribution* GetSamplingDistribution(RooArgSet& paramPoint);
      //virtual SamplingDistribution* GetSamplingDistributionSingleWorker(RooArgSet& paramPoint);
      using ToyMCSampler::GenerateToyData;
      //virtual RooAbsData* GenerateToyData(RooArgSet& paramPoint) const {
      //   if(fExpectedNuisancePar) oocoutE((TObject*)NULL,InputArguments) << "HybridToyMCSampler: using expected nuisance parameters but ignoring weight. Use GetSamplingDistribution(paramPoint, weight) instead." << endl;
      //   double weight;
      //   return GenerateToyData(paramPoint, weight);
      //}
      virtual RooAbsData* GenerateToyData(RooArgSet& paramPoint, double& weight, RooAbsPdf& pdf ) const;
      void SetFitToData(RooAbsData* d) { fFitToData = d; }
      void SetDebuggingData(RooAbsData* d) { fDebuggingData = d; }
      void SetSampledNuisances(RooArgSet* s) { fSampledNuisances = s; }

   protected:

      RooAbsData* fFitToData ;
      RooAbsData* fDebuggingData ;
      RooArgSet*  fSampledNuisances ;
      mutable NuisanceParametersSamplerAndFitter *fNuisanceParametersSamplerAndFitter; //!

   protected:
   ClassDef(HybridToyMCSampler,3) // A simple implementation of the TestStatSampler interface
};
}


#endif
