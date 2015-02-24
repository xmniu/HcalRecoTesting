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
  
  Float_t tsShift3=HcalTimeSlew::delay(corrCharge[3],HcalTimeSlew::MCShift,fTimeSlewBias); 
  Float_t tsShift4=HcalTimeSlew::delay(corrCharge[4],HcalTimeSlew::MCShift,fTimeSlewBias); 
  Float_t tsShift5=HcalTimeSlew::delay(corrCharge[5],HcalTimeSlew::MCShift,fTimeSlewBias); 
  
  Float_t i3=0;
  getLandauFrac(tsShift3,tsShift3+25,i3);
  Float_t n3=0;
  getLandauFrac(tsShift3+25,tsShift3+50,n3);

  Float_t i4=0;
  getLandauFrac(tsShift4,tsShift4+25,i4);
  Float_t n4=0;
  getLandauFrac(tsShift4+25,tsShift4+50,n4);

  Float_t i5=0;
  getLandauFrac(tsShift5,tsShift5+25,i5);
  Float_t n5=0;
  getLandauFrac(tsShift5+25,tsShift5+50,n5);

  Float_t ch3=corrCharge[3]/i3;
  Float_t ch4=(i3*corrCharge[4]-n3*corrCharge[3])/(i3*i4);
  Float_t ch5=(n3*n4*corrCharge[3]-i3*n4*corrCharge[4]+i3*i4*corrCharge[5])/(i3*i4*i5);

  HLTOutput.clear();
  HLTOutput.push_back(ch3);
  HLTOutput.push_back(ch4);
  HLTOutput.push_back(ch5);

}

void HLTv2::applyXM(const std::vector<double> & inputCharge, const std::vector<double> & inputPedestal, std::vector<double> & HLTOutput) const {

  std::vector<double> corrCharge;
  fPedestalSubFxn_.Calculate(inputCharge, inputPedestal, corrCharge);

  double TS35[3];
  double TS46[3];
  double TS57[3];
  PulseFraction(corrCharge[3], TS35);
  PulseFraction(corrCharge[4], TS46);
  PulseFraction(corrCharge[5], TS57);

  double a3[3] = {TS35[0], TS35[1], TS35[2]};
  double b3[3] = {0., TS46[0], TS46[1]};
  double c3[3] = {0., 0., TS57[0]};
  double d3[3] = {corrCharge[3], corrCharge[4], corrCharge[5]};

  double deno3 = Det3(a3, b3, c3);

  double A3 = Det3(d3, b3, c3) / deno3;
  double A4 = Det3(a3, d3, c3) / deno3;
  double A5 = Det3(a3, b3, d3) / deno3;

  HLTOutput.clear();
  HLTOutput.push_back(A3);
  HLTOutput.push_back(A4);
  HLTOutput.push_back(A5);

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

void HLTv2::PulseFraction(Double_t fC, Double_t *TS46) const{

  //static Double_t TS3par[3] = {0.44, -18.6, 5.136}; //Gaussian parameters: norm, mean, sigma for the TS3 fraction                    
  static Double_t TS4par[3] = {0.71, -5.17, 12.23}; //Gaussian parameters: norm, mean, sigma for the TS4 fraction                      
  static Double_t TS5par[3] = {0.258, 0.0178, 4.786e-4}; // pol2 parameters for the TS5 fraction                                       
  static Double_t TS6par[4] = {0.06391, 0.002737, 8.396e-05, 1.475e-06};// pol3 parameters for the TS6 fraction                        

  Double_t tslew = HcalTimeSlew::delay(fC,HcalTimeSlew::MC,fTimeSlewBias);

  TS46[0] = TS4par[0] * TMath::Gaus(tslew,TS4par[1],TS4par[2]); // fraction of pulse in the TS4                                        
  TS46[1] = TS5par[0] + TS5par[1]*tslew + TS5par[2]*tslew*tslew; // fraction of pulse in the T5S                                       
  TS46[2] = TS6par[0] + TS6par[1]*tslew + TS6par[2]*tslew*tslew + TS6par[3]*tslew*tslew*tslew; //fraction of pulse in the TS6          

  return;
}

double HLTv2::Det2(double *b, double *c) const{
  return b[1]*c[2]-b[2]*c[1];
}

double HLTv2::Det3(double *a, double *b, double *c) const{
  return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]);
}
