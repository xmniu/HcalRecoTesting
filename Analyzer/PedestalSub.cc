#include <iostream>
#include <cmath>
#include <climits>
#include "PedestalSub.h"

using namespace std;

PedestalSub::PedestalSub() {
}

PedestalSub::~PedestalSub() { 
}

void PedestalSub::Calculate(const std::vector<double> & inputCharge, const std::vector<double> & inputPedestal, const std::vector<double> & inputGain) const
{

  Float_t corrCharge[10];

  Int_t iDump=0;

  for (Int_t i=0; i<10; i++) {
    corrCharge[i]=inputCharge[i]-inputPedestal[i]; //is this the "database pedestal" i need to subtract?
    if (inputCharge[i]<0) {
      cout << "wtfInputCharge! " << inputCharge[i] << endl;
      iDump=1;
    }
    //if (corrCharge[i]<0) cout << "wtfCorrCharge! " << corrCharge[i] << endl;
  }

  if (iDump==1) {

    for (Int_t i=0; i<10; i++) {
      cout << inputCharge[i] << ", ";
    }
    cout << endl;

  }

}  
