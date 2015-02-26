#ifndef TimeSlewPar_h
#define TimeSlewPar_h 1

#include <typeinfo>
#include <algorithm>
#include <vector>

// #include "DataFormats/HcalDetId/interface/HcalDetId.h"
//#include "HcalPulseShapes.h"
#include "PedestalSub.h"
// #include "CalibFormats/HcalObjects/interface/HcalCoder.h"
// #include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
//#include "HybridMinimizer.h"

#include <TMinuit.h>

//#include <TH1.h>
#include "Minuit2/FCNBase.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"



class TimeSlewPar 
{
public:
    TimeSlewPar();
    ~TimeSlewPar();
    TF1 *slewFit;
    // This is the CMSSW Implementation of the apply function
    //void apply(const CaloSamples & cs, const std::vector<int> & capidvec, const HcalCalibrations & calibs, std::vector<double> & correctedOutput) const;
    // This is the edited implementation for our standalone test code
    void getParameters(const std::vector<double> & inputCharge, const std::vector<double> & inputPedestal, double & RatioTS54, double & TimeSlew, double & Pulse, TF1 *slewFit, PedestalSub pedSubFxn_) const;

private:

};

#endif // HLTAnalyzer_h
