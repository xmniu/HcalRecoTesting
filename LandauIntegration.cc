#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>
#include <TH1F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <Rtypes.h>

#include <TMath.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAttLine.h>
#include <TPaveText.h>
#include <TColor.h>

#include "TTree.h"
#include "TCanvas.h"

#include "HttStyles.h"

#endif

void LandauIntegration() {

  TCanvas *c1 = MakeCanvas("c1", "c1", 800, 600);

  //Float_t intCharge[10] = {4.84858, 2.05987, 4.97552, 2.38107, 17.8341, 13.292, 7.22145, 4.62927, 7.10693, 3.18308};

  //Landau pulse shape from https://indico.cern.ch/event/345283/contribution/3/material/slides/0.pdf
  //Landau turn on by default at left edge of time slice 

  Float_t timeSlewShift=10;

  // normalized to: cout << pulse_shape->Integral(-25,10000) << endl;
  TF1* pulse_shape = new TF1("pulse_shape","[0]*TMath::Landau((x+[1]),[2],[3])",-50,100);
  pulse_shape->SetParameter(0,1/4.458);
  pulse_shape->SetParameter(1,timeSlewShift);
  pulse_shape->SetParameter(2,14.36);
  pulse_shape->SetParameter(3,4.46);

  pulse_shape->SetTitle("");
  pulse_shape->GetXaxis()->SetTitle("Time");
  pulse_shape->GetYaxis()->SetTitle("Charge");

  pulse_shape->Draw("");

  TF1* ts0 = new TF1("pulse_shape","[0]*TMath::Landau((x+[1]),[2],[3])",-25,0);
  TF1* ts1 = new TF1("pulse_shape","[0]*TMath::Landau((x+[1]),[2],[3])",0,25);
  TF1* ts2 = new TF1("pulse_shape","[0]*TMath::Landau((x+[1]),[2],[3])",25,50);
  TF1* ts3 = new TF1("pulse_shape","[0]*TMath::Landau((x+[1]),[2],[3])",50,75);

  ts0->SetParameter(0,1/4.458);
  ts0->SetParameter(1,timeSlewShift);
  ts0->SetParameter(2,14.36);
  ts0->SetParameter(3,4.46);

  ts1->SetParameter(0,1/4.458);
  ts1->SetParameter(1,timeSlewShift);
  ts1->SetParameter(2,14.36);
  ts1->SetParameter(3,4.46);

  ts2->SetParameter(0,1/4.458);
  ts2->SetParameter(1,timeSlewShift);
  ts2->SetParameter(2,14.36);
  ts2->SetParameter(3,4.46);

  ts3->SetParameter(0,1/4.458);
  ts3->SetParameter(1,timeSlewShift);
  ts3->SetParameter(2,14.36);
  ts3->SetParameter(3,4.46);

  pulse_shape->SetLineColor(kBlack);
  pulse_shape->SetLineWidth(3);
  ts0->SetLineColor(kBlack);
  ts1->SetLineColor(kBlack);
  ts2->SetLineColor(kBlack);
  ts3->SetLineColor(kBlack);
  ts0->SetFillColor(kBlack);

  ts0->SetFillColor(kYellow-7);
  //ts0->SetLineColor(kYellow-7);
  ts0->SetFillStyle(1001);
  //ts1->SetLineColor(kGreen-7);
  ts1->SetFillColor(kGreen-7);
  ts1->SetFillStyle(1001);
  //ts2->SetLineColor(kCyan-9);
  ts2->SetFillColor(kCyan-9);
  ts2->SetFillStyle(1001);
  //ts3->SetLineColor(kViolet-9);
  ts3->SetFillColor(kViolet-9);
  ts3->SetFillStyle(1001);

  ts0->Draw("same");
  ts1->Draw("same");
  ts2->Draw("same");
  ts3->Draw("same");

  pulse_shape->Draw("same l");

  c1->SaveAs("landau_with_slew_shift.png");

  /*  cout << "{ ";
  for (Int_t i=0; i<99; i++) {
    cout << pulse_shape->Integral(i,i+1) << ", ";

  }
  cout << "}";*/
  /*
  Float_t landauFrac[100] = { 7.6377e-05, 0.000342278, 0.00111827, 0.00283151, 0.00583919, 0.0102101, 0.0156382, 0.0215401, 0.0272533, 0.0322193, 0.036083, 0.0387059, 0.0401229, 0.0404802, 0.0399773, 0.0388241, 0.0372142, 0.0353118, 0.0332481, 0.0311228, 0.029008, 0.0269535, 0.0249918, 0.023142, 0.0214139, 0.0198106, 0.0183306, 0.0169693, 0.0157203, 0.0145765, 0.0135299, 0.0125729, 0.0116979, 0.0108976, 0.0101654, 0.00949491, 0.00888051, 0.00831695, 0.00779945, 0.00732371, 0.00688584, 0.00648234, 0.00611004, 0.0057661, 0.00544794, 0.00515327, 0.00488, 0.00462627, 0.00439037, 0.0041708, 0.00396618, 0.00377526, 0.00359691, 0.00343013, 0.00327399, 0.00312765, 0.00299034, 0.00286138, 0.00274014, 0.00262604, 0.00251856, 0.00241722, 0.00232157, 0.00223122, 0.0021458, 0.00206497, 0.00198842, 0.00191587, 0.00184705, 0.00178172, 0.00171965, 0.00166064, 0.0016045, 0.00155105, 0.00150013, 0.00145159, 0.00140527, 0.00136107, 0.00131884, 0.00127848, 0.00123989, 0.00120296, 0.0011676, 0.00113373, 0.00110126, 0.00107013, 0.00104026, 0.00101159, 0.000984047, 0.000957586, 0.000932148, 0.000907681, 0.00088414, 0.000861478, 0.000839653, 0.000818625, 0.000798357, 0.000778814,0.000759962 };

  Float_t tStart=0.5, tEnd=25.5;
  
  cout << std::int(std::ceil(tStart)) << ", " << std::int(std::ceil(tEnd)) << endl;
  Float_t sum=0;
  for (Int_t i=std::int(std::ceil(tStart)); i<std::int(std::ceil(tEnd)); i++) {
    sum+=landauFrac[i];
  }
  cout <<  sum << endl;

  cout << pulse_shape->Integral(0,25) << endl;
  cout << pulse_shape->Integral(25,50) << endl;*/
/*
  //Float_t totalCharge=pulse_shape->Integral(-25,10000);
  //cout << totalCharge << endl;
  // fraction of Landau pulse falling in "in time" time slice
  Float_t fracLandauPREV = pulse_shape->Integral(-25,0)/totalCharge;
  Float_t fracLandauIT = pulse_shape->Integral(0,25)/totalCharge;
  // fraction of Landau pulse falling in "next" time slice
  Float_t fracLandauNEXT = pulse_shape->Integral(25,50)/totalCharge;
  // fraction of Landau pulse falling in "in time"+"next" time slice
  Float_t fracLandauITNEXT = pulse_shape->Integral(0,50)/totalCharge;
  
  cout << "Integrated charge assuming one pulse at T3: " << endl;
  cout << " No correction: " << intCharge[3] << endl;
  cout << " Contribution in TS4: " << intCharge[3]*fracLandauNEXT/fracLandauIT << endl;
  Float_t contribPREV=intCharge[3]*fracLandauNEXT/fracLandauIT;

  cout << "Integrated charge assuming one pulse at T6: " << endl;
  cout << " No correction: " << intCharge[6] << endl;
  cout << " Expected contribution in T5: " << intCharge[6]*fracLandauPREV/fracLandauIT << endl;
  Float_t contribNN=intCharge[6]*fracLandauPREV/fracLandauIT;

  cout << "Integrated charge in T5: " << endl;
  cout << " No correction: " << intCharge[5] << endl;
  cout << " With Subtraction from T6: " <<  intCharge[5]-contribNN << endl;

  cout << "Integrated charge in T4: " << endl;
  cout << " No correction: " << intCharge[4] << endl;
  cout << " With subtraction from T3: " << intCharge[4]-contribPREV << endl;

  cout << "Total charge from T4: " << intCharge[4]+intCharge[5]-contribNN-contribPREV << endl;

  cout << "Frac IT/NEXT: " << fracLandauIT/fracLandauNEXT << endl;
  cout << "Calculated: " << (intCharge[4]-contribPREV)/(intCharge[5]-contribNN) << endl;


}

double getLandauFrac(Float_t tStart, Float_t tEnd) const{
*/
}
