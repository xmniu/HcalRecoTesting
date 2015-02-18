#ifndef PedestalSub_h
#define PedestalSub_h 1

#include <typeinfo>

#include <TH1.h>
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"

class PedestalSub
{
 public:
  PedestalSub();
  ~PedestalSub();
  // This is the CMSSW Implementation of the apply function
  //void apply(const CaloSamples & cs, const std::vector<int> & capidvec, const HcalCalibrations & calibs, std::vector<double> & correctedOutput) const;
  // This is the edited implementation for our standalone test code

  void Calculate(const std::vector<double> & inputCharge, const std::vector<double> & inputPedestal, const std::vector<double> & inputGain) const;  
  
 private:
  
};

#endif 
