#ifndef Analysis_H
#define Analysis_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMinuit.h"

// Include all of the other classes
#include "HcalPulseShapes.h"
#include "HcalPulseShape.h"
#include "HybridMinimizer.h"
#include "PulseShapeFitOOTPileupCorrection.h"
#include "HLTv2.h"
#include "TimeSlewPar.h"

#include "inverseGaussCDF.hh"
#include "sampleQuantile.hh"
#include "PedestalSub.h"

#include "analysistree.h"

using namespace std;
const double EL_CHARGE=1.60217657*pow(10,-19);
const double FEMTO=pow(10,-15);

double x[10],y[10],errory[10];

class Analysis : public analysistree
{
 public:

  enum HcalRegion {All=0, Barrel=1, Endcap=2};

  string Input_File;
  string Output_File;
  string Plot_Dir;
  int Entries;
  int Region;
  int Condition;
  int Baseline;
  int Time_Slew;
  int Neg_Charges;
  float Threshold;
  float Quantile;

//For time slew parametrization only
  TH2F *TimeSlewPulse_All;
  TH2F *TimeSlewPulse_HB;
  TH2F *TimeSlewPulse_HE;

  TF1 *slewFit;
  TF1 *timeslewFit;
  TF1 *chtimeslewFit;

  TH2F *TimeSlewPulse[59][72][3];
  TH2F *Par0[3];
  TH2F *Par1[3];
  TH2F *Par2[3];

  Analysis(TTree *tree);
  ~Analysis();

  void Init(char* paramfile);
  void DefineHistograms();
  void TSP();
  void Process();
  void MakeCutflow();
  void FillHistograms();
  void Finish();

  void MakePedestalPlots();
  void DoHlt();
  void MakeTimeSlewPlots();

  void useMethod2(){psFitOOTpuCorr_ = std::auto_ptr<PulseShapeFitOOTPileupCorrection>(new PulseShapeFitOOTPileupCorrection()); }
  std::auto_ptr<PedestalSub> pedSubFxn_= std::auto_ptr<PedestalSub>(new PedestalSub());
         
 private:
  TFile *fout;
  TFile *fTSP;
  std::auto_ptr<PulseShapeFitOOTPileupCorrection> psFitOOTpuCorr_= std::auto_ptr<PulseShapeFitOOTPileupCorrection>(new PulseShapeFitOOTPileupCorrection());
  std::auto_ptr<HLTv2> hltv2_= std::auto_ptr<HLTv2>(new HLTv2());
  std::auto_ptr<TimeSlewPar> TimeSlewParameters = std::auto_ptr<TimeSlewPar>(new TimeSlewPar());
  HcalPulseShapes theHcalPulseShapes_;

};
#endif // Analysis_H 
