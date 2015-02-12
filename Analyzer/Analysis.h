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
//#include "TMinuit.h"

// Include all of the other classes
#include "HcalPulseShapes.h"
#include "HcalPulseShape.h"
#include "HybridMinimizer.h"
#include "PulseShapeFitOOTPileupCorrection.h"
#include "HLTAnalyzer.h"

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
  
  TH2D *hCharge_Method2_v_HLT;
  
  TProfile *hHLTResolution;
 
  Analysis(TTree *tree);
  ~Analysis();

  void Init(char* paramfile);
  void DefineHistograms();
  void Process();
  void MakeCutflow();
  void FillHistograms();
  void Finish();
 // void useMethod2(){psFitOOTpuCorr_ = std::auto_ptr<PulseShapeFitOOTPileupCorrection>(new PulseShapeFitOOTPileupCorrection()); }
         
 private:
  TFile *fout;
  std::auto_ptr<PulseShapeFitOOTPileupCorrection> psFitOOTpuCorr_= std::auto_ptr<PulseShapeFitOOTPileupCorrection>(new PulseShapeFitOOTPileupCorrection());
  std::auto_ptr<HLTAnalyzer> hltThing_= std::auto_ptr<HLTAnalyzer>(new HLTAnalyzer());
  HcalPulseShapes theHcalPulseShapes_;
  
   
};
#endif // Analysis_H 
