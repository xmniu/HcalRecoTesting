#include <iostream>
#include <cmath>
#include <climits>
#include "HLTv2.h"

using namespace std;

HLTv2::HLTv2() {
}

HLTv2::~HLTv2() { 
}

void HLTv2::Init(HcalTimeSlew::ParaSource tsParam, HcalTimeSlew::BiasSetting bias, NegStrategy nStrat, PedestalSub pedSubFxn_) {
  fTimeSlew=tsParam;
  fTimeSlewBias=bias;
  fNegStrat=nStrat;
  fPedestalSubFxn_=pedSubFxn_;
}

void HLTv2::apply(const std::vector<double> & inputCharge, const std::vector<double> & inputPedestal, std::vector<double> & HLTOutput) const {

  std::vector<double> corrCharge;

  fPedestalSubFxn_.Calculate(inputCharge, inputPedestal, corrCharge);
  
  // Iteration one assuming time slew                                                                                                                                                                       
  Float_t tsShift=HcalTimeSlew::delay(corrCharge[4],fTimeSlew,fTimeSlewBias); 
  
  Float_t fracL_prev=0;
  getLandauFrac(tsShift-25, tsShift, fracL_prev);
  Float_t fracL_intime=0;
  getLandauFrac(tsShift,tsShift+25,fracL_intime);
  Float_t fracL_next=0;
  getLandauFrac(tsShift+25,tsShift+50,fracL_next);
  
  Float_t ch3 = ( (fracL_intime*fracL_intime-fracL_prev*fracL_next)*corrCharge[3] - fracL_intime*fracL_next*corrCharge[4] + fracL_next*fracL_next*corrCharge[5] ) /( fracL_intime*(fracL_intime*fracL_intime-2*fracL_prev*fracL_next ) );
  Float_t ch4 =  (fracL_prev*corrCharge[3] - fracL_intime*corrCharge[4] + fracL_next*corrCharge[5])/( 2*fracL_prev*fracL_next-fracL_intime*fracL_intime);
  Float_t ch5 = (fracL_prev*fracL_prev*corrCharge[3]+(fracL_intime*fracL_intime-fracL_prev*fracL_next)*corrCharge[5]-fracL_next*fracL_intime*corrCharge[4]) / (fracL_intime*(fracL_intime*fracL_intime-2*fracL_prev*fracL_next ) );

  HLTOutput.clear();
  HLTOutput.push_back(ch3);
  HLTOutput.push_back(ch4);
  HLTOutput.push_back(ch5);

}

// Landau function integrated in 1 ns intervals
//Landau pulse shape from https://indico.cern.ch/event/345283/contribution/3/material/slides/0.pdf
//Landau turn on by default at left edge of time slice 
// normalized to 1 on [0,10000]
void HLTv2::getLandauFrac(Float_t tStart, Float_t tEnd, Float_t &sum) const{
  // can be further optimized to reduce computational time
  if (tStart==0 && tEnd==25) {
    sum=0.613668;
    return;
  }
  else if (tStart==0 && tEnd==50) {
    sum=0.853668;
    return;
  }

  Float_t landauFrac[100] = { 7.6377e-05, 0.000342278, 0.00111827, 0.00283151, 0.00583919, 0.0102101, 0.0156382, 0.0215401, 0.0272533, 0.0322193, 0.036083, 0.0387059, 0.0401229, 0.0404802, 0.0399773, 0.0388241, 0.0372142, 0.0353118, 0.0332481, 0.0311228, 0.029008, 0.0269535, 0.0249918, 0.023142, 0.0214139, 0.0198106, 0.0183306, 0.0169693, 0.0157203, 0.0145765, 0.0135299, 0.0125729, 0.0116979, 0.0108976, 0.0101654, 0.00949491, 0.00888051, 0.00831695, 0.00779945, 0.00732371, 0.00688584, 0.00648234, 0.00611004, 0.0057661, 0.00544794, 0.00515327, 0.00488, 0.00462627, 0.00439037, 0.0041708, 0.00396618, 0.00377526, 0.00359691, 0.00343013, 0.00327399, 0.00312765, 0.00299034, 0.00286138, 0.00274014, 0.00262604, 0.00251856, 0.00241722, 0.00232157, 0.00223122, 0.0021458, 0.00206497, 0.00198842, 0.00191587, 0.00184705, 0.00178172, 0.00171965, 0.00166064, 0.0016045, 0.00155105, 0.00150013, 0.00145159, 0.00140527, 0.00136107, 0.00131884, 0.00127848, 0.00123989, 0.00120296, 0.0011676, 0.00113373, 0.00110126, 0.00107013, 0.00104026, 0.00101159, 0.000984047, 0.000957586, 0.000932148, 0.000907681, 0.00088414, 0.000861478, 0.000839653, 0.000818625, 0.000798357, 0.000778814, 0.000759962 };

  sum=0;
  for (Int_t i=int(ceil(tStart)); i<int(ceil(tEnd)); i++) {
    if (i<0) sum+=0;
    else sum+=landauFrac[i];
  }
  return;

}
