#ifndef Analysis_H
#define Analysis_H

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
#include "HLTAnalyzer.h"
#include "HLTv2.h"
#include "PedestalSub.h"

#include <string>
#include <vector>
#include <math.h>

#include "analysistree.h"

using namespace std;
const double EL_CHARGE=1.60217657*pow(10,-19);
const double FEMTO=pow(10,-15);

double x[10],y[10],errory[10];

class Analysis : public analysistree
{
 public:

  string Input_File;
  string Output_File;

  int nevents =0;
  
  TH2F *RatioPulse;
  TH2F *TimeSlewPulse;
  
  TH1F *Norm0;
  TH1F *Norm1;
  TH1F *Norm2;
  
  //    TH1F *Ped;
  //    TH1F *Time;
  //    TH1F *Chi2;
  TF1 *slewFit;
  TF1 *logtimeslewFit;
  TF1 *exptimeslewFit;
  //Variables

  //Histograms
  TH1D *NUMBER_TS_ABOVE_THR_HB;
  TH1D *NUMBER_TS_ABOVE_THR_HE;
  
  TH1D *CHARGE_TSTOT_HB_FIT;
  TH1D *CHARGE_TSTOT_HE_FIT;
  TH1D *PULSE_ARRIVAL_HB_FIT;
  TH1D *PULSE_ARRIVAL_HE_FIT;

  TH1D *hEdepDist_all;
  TH1D *hEdepDist_not3456;
  TH1D *hEdepDist_not345;
  TH1D *hEdepDist_not45;
  TH1D *hEdepDist_least;
  TH1D *hEdepDist_least4;
  TH1D *hEdepDist_least_not3456;
  TH1D *hEdepDist_least_not345;
  TH1D *hEdepDist_least_not45;
  TH1D *hEdepDist_least4_not3456;
  TH1D *hEdepDist_least4_not345;
  TH1D *hEdepDist_least4_not45;
  
  //TH2D *hCharge_Method2_v_HLT;
  //TH2D *hCharge_Method2_v_JAY;
  //TH2D *hCharge_HLT_v_JAY;

  std::vector<TH1D*> vHistPed;
  std::vector<TH1D*> vHistVal;

  TH1D* hPedSub1;
  TH1D* hPedSub2;
  TH1D* hPedDiff;

  TH2D *h45vHLT0;
  TH2D *h45vHLT1;
  TH2D *h45vHLT2;
  TProfile *p45vHLT0;
  TProfile *p45vHLT1;
  TProfile *p45vHLT2;

  TH2D *hM2vHLT0;
  TH2D *hM2vHLT1;
  TH2D *hM2vHLT2;
  TProfile *pM2vHLT0;
  TProfile *pM2vHLT1;
  TProfile *pM2vHLT2;

  //TProfile *hHLTResolution;
  //TProfile *hJayResolution;
 
  Analysis(TTree *tree);
  ~Analysis();

  void Init(char* paramfile);
  void DefineHistograms();
  void Process();
  void MakeCutflow();
  void FillHistograms();
  void Finish();

  void MakePedestalPlots(int* n);
  void useMethod2(){psFitOOTpuCorr_ = std::auto_ptr<PulseShapeFitOOTPileupCorrection>(new PulseShapeFitOOTPileupCorrection()); }
         
 private:
  TFile *fout;
  std::auto_ptr<PulseShapeFitOOTPileupCorrection> psFitOOTpuCorr_= std::auto_ptr<PulseShapeFitOOTPileupCorrection>(new PulseShapeFitOOTPileupCorrection());
  std::auto_ptr<HLTAnalyzer> hltThing_= std::auto_ptr<HLTAnalyzer>(new HLTAnalyzer());
  std::auto_ptr<HLTv2> hltv2_= std::auto_ptr<HLTv2>(new HLTv2());
  std::auto_ptr<PedestalSub> pedSubFxn_= std::auto_ptr<PedestalSub>(new PedestalSub());
  HcalPulseShapes theHcalPulseShapes_;
  
   
};
#endif // Analysis_H 
