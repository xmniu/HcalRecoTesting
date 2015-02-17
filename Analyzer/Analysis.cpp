#define Analysis_cxx

#include "Analysis.h"
#include "readparameters/readparameters.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>

using namespace std;

double round_nplaces(double value, int to){
  int places = 1, whole = value;
  for(int i = 0; i < to; i++) places *= 10;
  value -= whole; //leave decimals
  value *= places; //0.1234 -> 123.4
  value = round(value);//123.4 -> 123
  value /= places; //123 -> .123
  value += whole; //bring the whole value back
  return value;
}

string int2string(int i){
  stringstream ss;
  string ret;
  ss<<i;
  ss>>ret;
  return ret;
}

// Ugly hack to apply energy corrections to some HB- cells
double eCorr(int ieta, int iphi, double energy) {
// return energy correction factor for HBM channels 
// iphi=6 ieta=(-1,-15) and iphi=32 ieta=(-1,-7)
// I.Vodopianov 28 Feb. 2011
  static const float low32[7]  = {0.741,0.721,0.730,0.698,0.708,0.751,0.861};
  static const float high32[7] = {0.973,0.925,0.900,0.897,0.950,0.935,1};
  static const float low6[15]  = {0.635,0.623,0.670,0.633,0.644,0.648,0.600,
				  0.570,0.595,0.554,0.505,0.513,0.515,0.561,0.579};
  static const float high6[15] = {0.875,0.937,0.942,0.900,0.922,0.925,0.901,
				  0.850,0.852,0.818,0.731,0.717,0.782,0.853,0.778};
  
  double slope, mid, en;
  double corr = 1.0;

  if (!(iphi==6 && ieta<0 && ieta>-16) && !(iphi==32 && ieta<0 && ieta>-8)) 
    return corr;

  int jeta = -ieta-1;
  double xeta = (double) ieta;
  if (energy > 0.) en=energy;
  else en = 0.;

  if (iphi == 32) {
    slope = 0.2272;
    mid = 17.14 + 0.7147*xeta;
    if (en > 100.) corr = high32[jeta];
    else corr = low32[jeta]+(high32[jeta]-low32[jeta])/(1.0+exp(-(en-mid)*slope));
  }
  else if (iphi == 6) {
    slope = 0.1956;
    mid = 15.96 + 0.3075*xeta;
    if (en > 100.0) corr = high6[jeta];
    else corr = low6[jeta]+(high6[jeta]-low6[jeta])/(1.0+exp(-(en-mid)*slope));
  }

  return corr;
}

int main(int argc, char **argv)
{
  int ret=0;
  if (argc!=2) {
    cerr<<"Usage: ./Analysis <paramfile>"<<endl;
    ret=1;
  } else {

    readparameters rp(argv[1]);
    //TChain* ch = new TChain("HcalNoiseTree");
    TChain* ch = new TChain("ExportTree/HcalTree");

    string filelistname;

    filelistname=rp.get<string>((string("in_filelist")).c_str());
    string line;
    ifstream filelist(filelistname.c_str());
    if (filelist.fail()) { //catch
      cerr << "\nERROR: Could not open " << filelistname << endl;
      exit(1);
    }
    while (getline(filelist,line)) {
      ch->Add(line.c_str());
    }

  Analysis Ana25ns((TTree*) ch);

  Ana25ns.Init(argv[1]);
  Ana25ns.DefineHistograms();
  Ana25ns.Process();
  Ana25ns.Finish();
  }
 return ret;
}

Analysis::~Analysis() {
}
Analysis::Analysis(TTree *tree):analysistree(tree){};

void Analysis::Init(char* paramfile)
{
  try {
    readparameters rp(paramfile);
    try {Output_File=rp.get<string>("Output_File");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
  } 
  catch (exception& e) {cerr<<e.what()<<endl;} 
  return; 
}

void Analysis::Process() {
  if (fChain == 0) return;

  int nentries = fChain->GetEntries();

  cout<<"Number of Entries to Process: "<<nentries<<endl<<endl;

  int nbytes = 0, nb = 0;

  for (int jentry=0; jentry<nentries;jentry++) {
    if(nevents%10==0) cout<<" "<<nevents<<"\t out of  "<<nentries<<"\t have already been processed ("<<round_nplaces((double)nevents/nentries*100,1)<<"%/100%)"<<endl;

    nb = fChain->GetEntry(jentry);   nbytes += nb;

    nevents++;
    MakeCutflow();
    FillHistograms();
  }

  TCanvas *c_tstot_hb_fit = new TCanvas("c_tstot_hb_fit");
  CHARGE_TSTOT_HB_FIT->GetXaxis()->SetTitle("Total Charge from Fit [fC]");
  CHARGE_TSTOT_HB_FIT->GetXaxis()->SetTitleSize(0.05);
  CHARGE_TSTOT_HB_FIT->GetYaxis()->SetTitle("Hits");
  CHARGE_TSTOT_HB_FIT->GetYaxis()->SetTitleSize(0.05);
  CHARGE_TSTOT_HB_FIT->Draw();
  c_tstot_hb_fit->SetLogy();
  c_tstot_hb_fit->SaveAs("CHARGE_TSTOT_HB_FIT.png");

  TCanvas *c_tstot_he_fit = new TCanvas("c_tstot_he_fit");
  CHARGE_TSTOT_HE_FIT->GetXaxis()->SetTitle("Total Charge from Fit [fC]");
  CHARGE_TSTOT_HE_FIT->GetXaxis()->SetTitleSize(0.05);
  CHARGE_TSTOT_HE_FIT->GetYaxis()->SetTitle("Hits");
  CHARGE_TSTOT_HE_FIT->GetYaxis()->SetTitleSize(0.05);
  CHARGE_TSTOT_HE_FIT->Draw();
  c_tstot_he_fit->SetLogy();
  c_tstot_he_fit->SaveAs("CHARGE_TSTOT_HE_FIT.png");

  TCanvas *c_t_hb_fit = new TCanvas("c_t_hb_fit");
  PULSE_ARRIVAL_HB_FIT->GetXaxis()->SetTitle("Arrival Time from Fit [ns]");
  PULSE_ARRIVAL_HB_FIT->GetXaxis()->SetTitleSize(0.05);
  PULSE_ARRIVAL_HB_FIT->GetYaxis()->SetTitle("Hits");
  PULSE_ARRIVAL_HB_FIT->GetYaxis()->SetTitleSize(0.05);
  PULSE_ARRIVAL_HB_FIT->Draw();
  //c_t_hb_fit->SetLogy();
  c_t_hb_fit->SaveAs("PULSE_ARRIVAL_HB_FIT.png");

  TCanvas *c_t_he_fit = new TCanvas("c_t_he_fit");
  PULSE_ARRIVAL_HE_FIT->GetXaxis()->SetTitle("Arrival Time from Fit [ns]");
  PULSE_ARRIVAL_HE_FIT->GetXaxis()->SetTitleSize(0.05);
  PULSE_ARRIVAL_HE_FIT->GetYaxis()->SetTitle("Hits");
  PULSE_ARRIVAL_HE_FIT->GetYaxis()->SetTitleSize(0.05);
  PULSE_ARRIVAL_HE_FIT->Draw();
 // c_t_he_fit->SetLogy();
  c_t_he_fit->SaveAs("PULSE_ARRIVAL_HE_FIT.png");
  
  TCanvas *c_Jay_profile = new TCanvas("c_jay_profile");
  hJayResolution->GetXaxis()->SetTitle("M.2 Energy [GeV]");
  hJayResolution->GetXaxis()->SetTitleSize(0.05);
  hJayResolution->GetYaxis()->SetTitle("(M.2 - HLT v2)/M2 Energy");
  hJayResolution->GetYaxis()->SetTitleSize(0.05);
  hJayResolution->Draw();
 // c_t_he_fit->SetLogy();
  c_Jay_profile->SaveAs("hJayResolution.png");

  TCanvas *c_hlt_profile = new TCanvas("c_hlt_profile");
  hHLTResolution->GetXaxis()->SetTitle("M.2 Energy [GeV]");
  hHLTResolution->GetXaxis()->SetTitleSize(0.05);
  hHLTResolution->GetYaxis()->SetTitle("(M.2 - HLT)/M2 Energy");
  hHLTResolution->GetYaxis()->SetTitleSize(0.05);
  hHLTResolution->Draw();
 // c_t_he_fit->SetLogy();
  c_hlt_profile->SaveAs("hHLTResolution.png");
  
  TCanvas *c_hlt_v_m2 = new TCanvas("c_hlt_v_m2");
  hCharge_Method2_v_HLT->GetXaxis()->SetTitle("Energy of M.2 Rechit [GeV]");
  hCharge_Method2_v_HLT->GetXaxis()->SetTitleSize(0.05);
  hCharge_Method2_v_HLT->GetYaxis()->SetTitle("Energy of HLT Rechit [GeV]");
  hCharge_Method2_v_HLT->GetYaxis()->SetTitleSize(0.05);
  hCharge_Method2_v_HLT->Draw("colz");
 // c_t_he_fit->SetLogy();
  c_hlt_v_m2->SaveAs("hCharge_Method2_v_HLT.png");

  TCanvas *c_jay_v_m2 = new TCanvas("c_jay_v_m2");
  hCharge_Method2_v_JAY->GetXaxis()->SetTitle("Energy of M.2 Rechit [GeV]");
  hCharge_Method2_v_JAY->GetXaxis()->SetTitleSize(0.05);
  hCharge_Method2_v_JAY->GetYaxis()->SetTitle("Energy of HLT v2 Rechit [GeV]");
  hCharge_Method2_v_JAY->GetYaxis()->SetTitleSize(0.05);
  hCharge_Method2_v_JAY->Draw("colz");
  // c_t_he_fit->SetLogy();
  c_jay_v_m2->SaveAs("hCharge_Method2_v_jay.png");

  TCanvas *c_jay_v_hlt = new TCanvas("c_jay_v_hlt");
  hCharge_HLT_v_JAY->GetXaxis()->SetTitle("Energy of HLT Rechit [GeV]");
  hCharge_HLT_v_JAY->GetXaxis()->SetTitleSize(0.05);
  hCharge_HLT_v_JAY->GetYaxis()->SetTitle("Energy of HLT v2 Rechit [GeV]");
  hCharge_HLT_v_JAY->GetYaxis()->SetTitleSize(0.05);
  hCharge_HLT_v_JAY->Draw("colz");
  // c_t_he_fit->SetLogy();
  c_jay_v_hlt->SaveAs("hCharge_hlt_v_jay.png");

 }
 
void Analysis::DefineHistograms()
{
  fout = new TFile(Output_File.c_str(), "RECREATE");
  
  // ---------------Plots for HLT ------------
  RatioPulse = new TH2F("RatioPulse","TS5/TS4 vs TS45",20,0.,2.,100,0.,500.);
  TimeSlewPulse = new TH2F("TimeSlewPulse","Time Slew vs TS45",25,-14.5,10.5,100,0.,500.);
  
  Norm0 = new TH1F("fC0","Amplitude in Pulse ealier [fC]",100,0.,500.);
  Norm1 = new TH1F("fC1","Amplitude in in-time Pulse [fC]",100,0.,500.);
  Norm2 = new TH1F("fC2","Amplitude in Pulse later [fC]",100,0.,500.);
  
  slewFit = new TF1("slewFit","pol4*expo(5)",-10.,14.);
  slewFit->SetParameters(1.07618e-02,-4.19145e-06,2.70310e-05,-8.71584e-08,1.86597e-07,3.59216e+00,-1.02057e-01);
    
  logtimeslewFit = new TF1("logtimeslewFit", "[0]+[1]*TMath::Log(x+[2])",0.,500.);
  logtimeslewFit->SetParameters(3.89838e+01, -6.93560e+00, 8.52052e+01);
  
  exptimeslewFit = new TF1("exptimeslewFit", "[0]+[1]*TMath::Exp([2]*x)",0.,500.);
  exptimeslewFit->SetParameters(-2.69330e+00, 1.09162e+01, -7.60722e-03);
  
  // Output Plots
  hHLTResolution=new TProfile("hHLTResolution","",20,0,100,-1.0,1.0);
  hJayResolution=new TProfile("hJayResolution","",20,0,100,-1.0,1.0);

  hCharge_Method2_v_HLT=new TH2D("hCharge_Method2_v_HLT","",50,0,250,50,0,250);

  hCharge_Method2_v_JAY=new TH2D("hCharge_Method2_v_jay","",50,0,250,50,0,250);
  hCharge_HLT_v_JAY=new TH2D("hCharge_Method2_v_jay","",50,0,250,50,0,250);
  
  PULSE_ARRIVAL_HB_FIT=new TH1D("PULSE_ARRIVAL_HB_FIT","",50,-30.0,30.0);
  PULSE_ARRIVAL_HE_FIT=new TH1D("PULSE_ARRIVAL_HE_FIT","",50,-30.0,30.0);

  CHARGE_TSTOT_HB_FIT=new TH1D("CHARGE_TSTOT_HB_FIT","",150,0,1500);
  CHARGE_TSTOT_HE_FIT=new TH1D("CHARGE_TSTOT_HE_FIT","",150,0,1500);

}

void Analysis::MakeCutflow() 
{

  for (int j = 0; j < (int)PulseCount; j++) {
    
    //=========================================================================      
    // These are the values we should set for Method 2 config with the python files
    // Don't currently have a setup to input with python but we shouldn't have to
    // change these values for our tests
    
    // --------------------------------------------------------------------
    bool iPedestalConstraint = true;
    bool iTimeConstraint = true;
    bool iAddPulseJitter = false;
    bool iUnConstrainedFit = false;
    bool iApplyTimeSlew = true;
    double iTS4Min = 5.;
    double iTS4Max = 500.;
    double iPulseJitter = 1.;
    double iTimeMean = -5.5;
    double iTimeSig = 5.;
    double iPedMean = 0.;
    double iPedSig = 0.5;
    double iNoise = 1.;
    double iTMin = -18.;
    double iTMax = 7.;
    double its3Chi2 = 5.;
    double its4Chi2 = 15.;
    double its345Chi2 = 100.;
    double iChargeThreshold = 6.;
    int iFitTimes = 1;
    
    //========================================================================

    // Set the Method 2 Parameters here
    psFitOOTpuCorr_->setPUParams(iPedestalConstraint,iTimeConstraint,iAddPulseJitter,iUnConstrainedFit,iApplyTimeSlew,
                  iTS4Min, iTS4Max, iPulseJitter,iTimeMean,iTimeSig,iPedMean,iPedSig,iNoise,iTMin,iTMax,its3Chi2,its4Chi2,its345Chi2,
                  iChargeThreshold,HcalTimeSlew::Medium, iFitTimes);

                  
    // Now set the Pulse shape type
    psFitOOTpuCorr_->setPulseShapeTemplate(theHcalPulseShapes_.getShape(105));

    // void PulseShapeFitOOTPileupCorrection::apply(const CaloSamples & cs, const std::vector<int> & capidvec, const HcalCalibrations & calibs, std::vector<double> & correctedOutput)
    // Changing to take the inputs (vectors) Charge and Pedestal and correctedOutput vector for the moment
    // and will remain this way unless we change the ntuple

    std::vector<double> correctedOutput, HLTOutput, Jay1Output, Jay2Output;
    std::vector<double> inputCaloSample, inputPedestal;
    std::vector<double> inputGain;

    for(int i = 0; i < 10; ++i) {

      // Note: In CMSSW the "Charge" vector is not already pedestal subtracted, unlike here
      // so I add the pedestal back to the charge so we can keep the same CMSSW implementation
      inputCaloSample.push_back(Charge[j][i]+Pedestal[j][i]);
      inputPedestal.push_back(Pedestal[j][i]);
      inputGain.push_back(Gain[j][i]);
    }
    
    // Begin Method 2
    psFitOOTpuCorr_->apply(inputCaloSample,inputPedestal,inputGain,correctedOutput);

    // Begin Xinmei's implementation of HLT
    double RatioTS54, TimeSlew, Pulse = 0.;
    hltThing_->apply(inputCaloSample,inputPedestal,inputGain,HLTOutput, RatioTS54, TimeSlew, Pulse, slewFit);

    // Begin Jay's implementation of HLT
    hltv2_->applyOnce(inputCaloSample,inputPedestal,inputGain,Jay1Output, RatioTS54, TimeSlew, Pulse, slewFit);

    hltv2_->applyOnceWithTS(inputCaloSample,inputPedestal,inputGain,Jay2Output, RatioTS54, TimeSlew, Pulse, slewFit);

    
    RatioPulse->Fill(RatioTS54, Pulse);
    TimeSlewPulse->Fill(TimeSlew, Pulse);
    Norm0->Fill(HLTOutput.at(0));
    Norm1->Fill(HLTOutput.at(1));
    Norm2->Fill(HLTOutput.at(2));

    // ------------ Fill our comparison histograms -------------------
    // Make sure that the output vectors have non-zero number of entries
    if(correctedOutput.size() > 1 && Jay2Output.size() > 1){
      
      hCharge_Method2_v_HLT->SetBinContent(correctedOutput.at(0)*0.2+1, HLTOutput.at(1)*Gain[j][0]*0.2+1, 1+hCharge_Method2_v_HLT->GetBinContent(correctedOutput.at(0)*0.2+1, HLTOutput.at(1)*Gain[j][0]*0.2+1));
      
      hCharge_Method2_v_JAY->SetBinContent(correctedOutput.at(0)*0.2+1, Jay2Output.at(1)*Gain[j][0]*0.2+1, 1+hCharge_Method2_v_JAY->GetBinContent(correctedOutput.at(0)*0.2+1, Jay2Output.at(1)*Gain[j][0]*0.2+1));
      
      hCharge_HLT_v_JAY->SetBinContent(HLTOutput.at(1)*Gain[j][0]*0.2+1, Jay2Output.at(1)*Gain[j][0]*0.2+1, 1+hCharge_HLT_v_JAY->GetBinContent(HLTOutput.at(1)*Gain[j][0]*0.2+1, Jay2Output.at(1)*Gain[j][0]*0.2+1));

      // Fill the TProfile with the % difference of energies
      double resolution = 0;
      correctedOutput.at(0) > 0 ? resolution = ( correctedOutput.at(0) - Jay2Output.at(1)*Gain[j][0] )/correctedOutput.at(0) : resolution = 0 ;
      hJayResolution->Fill(correctedOutput.at(0), resolution, 1);

      // Should do something with the correctedOutput vector here, such as fill a histogram...
      if(IEta[j] < 16 && correctedOutput.size() > 1) {
	if(correctedOutput.at(1) > -99.){
	  CHARGE_TSTOT_HB_FIT->Fill(correctedOutput.at(0));
	  PULSE_ARRIVAL_HB_FIT->Fill(correctedOutput.at(1));
	}
      } else if(IEta[j] >= 16 && correctedOutput.size() > 1){
	CHARGE_TSTOT_HE_FIT->Fill(correctedOutput.at(0));
	PULSE_ARRIVAL_HE_FIT->Fill(correctedOutput.at(1));
      }
    }
  }//PulseCount

}// End MakeCutflow Function

void Analysis::FillHistograms()
{
}

void Analysis::Finish()
{
  gStyle->SetOptFit(1);
  TimeSlewPulse->Draw("BOX");
  TimeSlewPulse->ProfileY("Y Profile",1,-1,"do")->Fit("pol1");
  TimeSlewPulse->ProfileY("X Profile",1,-1,"do")->Fit("pol1");

  fout->cd();
  fout->Write();
  fout->Close();
}
