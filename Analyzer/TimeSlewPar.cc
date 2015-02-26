#include <iostream>
#include <cmath>
#include <climits>
#include "TimeSlewPar.h"
//#include "isFinite.h"


TimeSlewPar::TimeSlewPar() {
}

TimeSlewPar::~TimeSlewPar() { 
}

void TimeSlewPar::getParameters(const std::vector<double> & inputCharge, const std::vector<double> & inputPedestal, double & RatioTS54, double & TimeSlew, double & Pulse, TF1 *slewFit, PedestalSub pedSubFxn_) const
{

// Pedestal/Baseline(2.7fC) Subtraction Method
///*
  int Counter = 0;
  double TotalPedestal = 0.;
  double AveragePedestal = 0.;
  double TS[10] = {0.};

  for(int i = 0; i < 10; i++){
//  No Pedestal Subtraction
//  TS[i] = inputCharge[i];
//  Pedestal Subtraction
    TS[i] = inputCharge[i] - inputPedestal[i];
    if(i != 4 && i != 5){
      TotalPedestal += TMath::Min(TS[i],2.7);
      Counter++;
    }
  }

  AveragePedestal = TotalPedestal / Counter;

//Parametrize before Baselien Subtraction
//  RatioTS54 = TS[5] / TS[4];
//  Pulse = TS[4];
//Parametrize after Baseline Subtraction
  RatioTS54 = (TS[5] - AveragePedestal) / (TS[4] - AveragePedestal);
  Pulse = TS[4] - AveragePedestal;
  TimeSlew = -slewFit->GetX(RatioTS54);
//*/

// Quantile Subtraction Method
/*
  std::vector<double> corrCharge = {};
  pedSubFxn_.Calculate(inputCharge, inputPedestal, corrCharge);
  RatioTS54 = corrCharge[5] / corrCharge[4];
  Pulse = corrCharge[4];
   TimeSlew = -slewFit->GetX(RatioTS54);
*/
}
