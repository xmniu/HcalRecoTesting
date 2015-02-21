#include <iostream>
#include <cmath>
#include <climits>
#include "PedestalSub.h"

using namespace std;

PedestalSub::PedestalSub() {
}

PedestalSub::~PedestalSub() { 
}

void PedestalSub::Init(Method method=AvgWithThresh, float threshold=0.0, float quantile=0.0) {
  fMethod=method;
  fThreshold=threshold;
  fQuantile=quantile;
  
  if ( (fMethod==DoNothing||fMethod==AvgWithoutThresh||fMethod==Percentile)&&(fThreshold!=0.0) ) {
    cout << "You are almost certainly doing something wrong. Check your PedestalSub::Calculate parameters, you're passing a threshold for a fMethod that doesn't use them!" << endl;
  }
  else if ( (fMethod!=Percentile)&&(fQuantile!=0.0) ) {
    cout << "You are almost certainly doing something wrong. Check your PedestalSub::Calculate parameters, you're passing a quantile for a fMethod that doesn't use them!" << endl;
  }
  else if ( (fMethod==AvgWithThresh||fMethod==AvgWithThreshNoPedSub)&&(fThreshold==0) ) {
    cout << "You are almost certainly doing something wrong. Check your PedestalSub::Calculate parameters, you're using 0.0 as your threshold!" << endl;
  }

  if (fMethod==Percentile) {
    //need to fix these parameters;

    fNoisePara=0.6; //fC, but only from pedestal
    fNoiseCorr=-inverseGaussCDF(fQuantile);

  }
  
}


void PedestalSub::Calculate(const std::vector<double> & inputCharge, const std::vector<double> & inputPedestal, std::vector<double> & corrCharge) const
{

  Float_t baseline=0;

  if (fMethod==DoNothing) { 
    baseline=0;
    for (Int_t i=0; i<10; i++) {
      corrCharge.push_back(inputCharge[i]-inputPedestal[i]);
    }  
  }
  else if (fMethod==AvgWithThresh) {
    for (Int_t i=0; i<10; i++) {
      if (i==4||i==5) continue;
      if ( (inputCharge[i]-inputPedestal[i])<fThreshold) {
	baseline+=(inputCharge[i]-inputPedestal[i]);
      }
      else baseline+=fThreshold;
    }
    baseline/=8;
    for (Int_t i=0; i<10; i++) {
      corrCharge.push_back(inputCharge[i]-inputPedestal[i]-baseline);
    }
  }
  else if (fMethod==AvgWithoutThresh) {
    for (Int_t i=0; i<10; i++) {
      if (i==4||i==5) continue;
      baseline+=(inputCharge[i]-inputPedestal[i]);
    }
    baseline/=8;
    for (Int_t i=0; i<10; i++) {
      corrCharge.push_back(inputCharge[i]-inputPedestal[i]-baseline);
    }
  }
  if (fMethod==AvgWithThreshNoPedSub) {
    for (Int_t i=0; i<10; i++) {
      if (i==4||i==5) continue;
      if ( (inputCharge[i])<fThreshold) {
	baseline+=(inputCharge[i]);
      }
      else baseline+=fThreshold;
    }
    baseline/=8;
    for (Int_t i=0; i<10; i++) {
      corrCharge.push_back(inputCharge[i]-baseline);
    }
  }
  else if (fMethod==Percentile) {

    baseline=sampleQuantile<10>(&inputCharge[0],fQuantile);
    baseline+=fNoisePara*fNoiseCorr;

    for (Int_t i=0; i<10; i++) {
      corrCharge.push_back(inputCharge[i]-baseline);
    }
  }
  
}  
