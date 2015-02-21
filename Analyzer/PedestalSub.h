#ifndef PedestalSub_h
#define PedestalSub_h 1

#include <typeinfo>

#include <TH1.h>
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"

#include "inverseGaussCDF.hh"
#include "sampleQuantile.hh"

class PedestalSub
{
 public:
  enum Method { DoNothing, AvgWithThresh, AvgWithoutThresh, AvgWithThreshNoPedSub, Percentile };

  PedestalSub();
  ~PedestalSub();
  
  void Init(Method method, float threshold, float quantile);
  
  // This is the CMSSW Implementation of the apply function
  //void apply(const CaloSamples & cs, const std::vector<int> & capidvec, const HcalCalibrations & calibs, std::vector<double> & correctedOutput) const;
  // This is the edited implementation for our standalone test code
  
  void Calculate(const std::vector<double> & inputCharge, const std::vector<double> & inputPedestal, std::vector<double> & corrCharge) const;

  Method fMethod;
  float fThreshold;
  float fQuantile;
  float fNoiseCorr;
  float fNoisePara; // need to figure out how to get this value for reak

 private:
  
};

#endif 
