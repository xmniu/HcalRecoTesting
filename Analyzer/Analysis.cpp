#define Analysis_cxx

#include "Analysis.h"
#include "readparameters/readparameters.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include <sys/stat.h>

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
    try {Entries=rp.get<int>("Entries");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Plot_Dir=rp.get<string>("Plot_Dir");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Region=rp.get<int>("Region");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Baseline=rp.get<int>("Baseline");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Time_Slew=rp.get<int>("Time_Slew");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Neg_Charges=rp.get<int>("Neg_Charges");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Threshold=rp.get<float>("Threshold");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Quantile=rp.get<float>("Quantile");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
  } 
  catch (exception& e) {cerr<<e.what()<<endl;} 

  cout << "Running on ";
  if (Entries==-1) cout << "all events." << endl;
  else cout << "first " << Entries << " events. " << endl;

  cout << "Output ROOT file: " << Output_File << endl;

  cout << "Using channels in ";
  if (Region==All) { cout << "all HCAL regions. " << endl; }
  else if (Region==Barrel) { cout << "HCAL barrel. " << endl; }
  else if (Region==Endcap) { cout << "HCAL endcap. " << endl; }

  int check=mkdir(Plot_Dir.c_str(),755);
  if (!check) {
    cout << "Saving files to: " << Plot_Dir << endl;
  }
  else {
    cout << "Couldn't make plot directory. Not running :(" << endl;
    exit(1);
  }

  if (Baseline==PedestalSub::DoNothing) {
    cout << "Pedestal subtraction only." << endl;
  }
  else if (Baseline==PedestalSub::AvgWithThresh) {
    cout << "Pedestal subtraction, average baseline subtraction with threshold: " << Threshold << endl;
  }
  else if (Baseline==PedestalSub::AvgWithoutThresh) {
    cout << "Pedestal subtraction, average baseline subtraction with no threshold. " << endl;
  }
  else if (Baseline==PedestalSub::AvgWithThreshNoPedSub) {
    cout << "Average baseline+pedestal subtraction with threshold: " << Threshold << endl;
  }
  else if (Baseline==PedestalSub::Percentile) {
    cout << "Percentile-based pedestal subtraction";
      if (Quantile<0 || Quantile>1) {
	cout << endl << "Quantile value out of range. Not running." << endl;
	exit(1);
      }
      else  {
	cout << "with quantile value: " << Quantile << endl;
      }
  }

  if (Time_Slew==HcalTimeSlew::TestStand) cout << "Using test stand medium WP time slew parameterization." << endl;
  else cout << "Sorry, you asked for a time slew parameterization I don't have implemented yet. Using test stand medium WP time slew parameterization." << endl;

  if (Neg_Charges==HLTv2::DoNothing) cout << "Not requiring positive charge outputs." << endl;
  else cout << "Sorry, you asked me to require positive charge outputs, but I don't have that implemented yet." << endl;

  pedSubFxn_->Init(((PedestalSub::Method)Baseline), Threshold, Quantile);

  return; 
}

void Analysis::Process() {
  if (fChain == 0) return;

  if (Entries==-1) Entries=fChain->GetEntries();

  int nbytes = 0, nb = 0;

  for (int jentry=0; jentry<Entries;jentry++) {
    //if(nevents%10==0) cout<<" "<<nevents<<"\t out of  "<<nentries<<"\t have already been processed ("<<round_nplaces((double)nevents/nentries*100,1)<<"%/100%)"<<endl;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    nevents++;
    //DoHltTests();
    MakeCutflow();
    //MakePedestalPlots();
    //FillHistograms();
  }
 
}
 
void Analysis::DefineHistograms()
{
  fout = new TFile(Output_File.c_str(), "RECREATE");
  
  Int_t xBins=100;
  Int_t xMin=0;
  Int_t xMax=100;

  Int_t xBins2=50;
  Int_t xMin2=0;
  Int_t xMax2=100;

  h45vHLT = new TH2D("h45vHLT0", "", xBins,xMin,xMax,xBins,xMin,xMax);
  p45vHLT = new TProfile("p45vHLT0", "", xBins2,xMin2,xMax2,-10,10);
  hM2vHLT = new TH2D("hM2vHLT0", "", xBins,xMin2,xMax2,xBins,xMin2,xMax2);
  pM2vHLT = new TProfile("pM2vHLT0", "", xBins2,xMin2,xMax2,-10,10);

  hPedSub = new TH1D("hPedSub", "", 80,-1,7);

  a3 = new TH1D("a3", "", 200, -10, 10);
  a4 = new TH1D("a4", "", 300, -25, 50);
  a5 = new TH1D("a5", "", 200, -5, 15);

  a4v3 = new TH2D("a4v3", "", 300, -25, 50, 300, -25, 50);
  a4v5 = new TH2D("a4v5", "", 300, -25, 50, 300, -25, 50);
  a5v3 = new TH2D("a5v3", "", 300, -25, 50, 300, -25, 50);

}

void Analysis::MakeCutflow() 
{

  for (int j = 0; j < (int)PulseCount; j++) {

    //if (IEta[j]>16) continue;
    //if (IEta[j]<17) continue;

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
    
    std::vector<double> correctedOutput, hltOutput;
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

    HcalTimeSlew::ParaSource useTS=(HcalTimeSlew::ParaSource)Time_Slew;
    HLTv2::NegStrategy useNS=(HLTv2::NegStrategy)Neg_Charges;

    // Begin Jay's implementation of HLT
    hltv2_->apply(inputCaloSample,inputPedestal,hltOutput,useTS,HcalTimeSlew::Medium,useNS,*pedSubFxn_);

    // ------------ Fill our comparison histograms -------------------
    // Make sure that the output vectors have non-zero number of entries

    if (correctedOutput.size() > 1) {
      if(hltOutput.size() > 1){
	h45vHLT->Fill( (Charge[j][4]), hltOutput.at(1),1);
	hM2vHLT->Fill( correctedOutput.at(0)/Gain[j][0], hltOutput.at(1),1);
	p45vHLT->Fill((Charge[j][4]), -(hltOutput.at(1)-(Charge[j][4]))/((Charge[j][4])),1);
	pM2vHLT->Fill( correctedOutput.at(0)/Gain[j][0], -(hltOutput.at(1)*Gain[j][0]-(correctedOutput.at(0)))/((correctedOutput.at(0))),1);

	a3->Fill(hltOutput.at(0));
	a4->Fill(hltOutput.at(1));
	a5->Fill(hltOutput.at(2));

	a4v3->Fill(hltOutput.at(1), hltOutput.at(0));
	a4v5->Fill(hltOutput.at(1), hltOutput.at(2));
	a5v3->Fill(hltOutput.at(2), hltOutput.at(0));

      }
    }
    
  }//PulseCount

}// End MakeCutflow Function

void Analysis::FillHistograms()
{
}

void Analysis::MakePedestalPlots() {

  /*  for (int j = 0; j < (int)PulseCount; j++) {
    std::vector<double> inputCaloSample, inputPedestal;
    std::vector<double> inputGain;
  }

  TCanvas *c1 = new TCanvas("c1");
  gStyle->SetOptStat    (0);

  hPedSub->GetXaxis()->SetTitle("Baseline [fC]");
  hPedSub->GetXaxis()->SetTitleSize(0.05);
  hPedSub->GetYaxis()->SetTitle("Counts");
  hPedSub->GetYaxis()->SetTitleSize(0.05);
  hPedSub->Draw("");
  c1->SaveAs(TString(Plot_Dir.c_str())+"/baseline.png");
  */

}

void Analysis::Finish()
{
  
  TCanvas *c1 = new TCanvas("c1");
  gStyle->SetOptStat    (0);
  
  h45vHLT->GetXaxis()->SetTitle("Charge TS4 [fC]");
  h45vHLT->GetXaxis()->SetTitleSize(0.05);
  h45vHLT->GetYaxis()->SetTitle("Charge from HLT [fC]");
  h45vHLT->GetYaxis()->SetTitleSize(0.05);
  h45vHLT->Draw("colz");
  c1->SaveAs(TString(Plot_Dir.c_str())+"/h4vHLT.png");
  
  p45vHLT->GetXaxis()->SetTitle("Charge in TS4 [fC]");
  p45vHLT->GetXaxis()->SetTitleSize(0.05);
  p45vHLT->GetYaxis()->SetTitle("(HLT - E4)/E4");
  p45vHLT->GetYaxis()->SetTitleSize(0.05);
  p45vHLT->GetYaxis()->SetRangeUser(-1,1);
  p45vHLT->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/p4vHLT.png");
  
  hM2vHLT->GetXaxis()->SetTitle("M2 Charge [fC]");
  hM2vHLT->GetXaxis()->SetTitleSize(0.05);
  hM2vHLT->GetYaxis()->SetTitle("Charge from HLT [fC]");
  hM2vHLT->GetYaxis()->SetTitleSize(0.05);
  hM2vHLT->Draw("colz");
  c1->SaveAs(TString(Plot_Dir.c_str())+"/hM2vHLT.png");
  
  pM2vHLT->GetXaxis()->SetTitle("M2 Charge [fC]");
  pM2vHLT->GetXaxis()->SetTitleSize(0.05);
  pM2vHLT->GetYaxis()->SetTitle("(HLT - M2)/(M2)");
  pM2vHLT->GetYaxis()->SetTitleSize(0.05);
  pM2vHLT->GetYaxis()->SetRangeUser(-1,1);
  pM2vHLT->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/pM2vHLT.png");

  a3->GetXaxis()->SetTitle("A3 [fC]");
  a3->GetXaxis()->SetTitleSize(0.05);
  a3->GetYaxis()->SetTitle("Counts");
  a3->GetYaxis()->SetTitleSize(0.05);
  a3->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a3.png");

  a4->GetXaxis()->SetTitle("A4 [fC]");
  a4->GetXaxis()->SetTitleSize(0.05);
  a4->GetYaxis()->SetTitle("Counts");
  a4->GetYaxis()->SetTitleSize(0.05);
  a4->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a4.png");

  a5->GetXaxis()->SetTitle("A5 [fC]");
  a5->GetXaxis()->SetTitleSize(0.05);
  a5->GetYaxis()->SetTitle("Counts");
  a5->GetYaxis()->SetTitleSize(0.05);
  a5->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a5.png");

  a4v3->GetXaxis()->SetTitle("A4 [fC]");
  a4v3->GetXaxis()->SetTitleSize(0.05);
  a4v3->GetYaxis()->SetTitle("A3 [fC]");
  a4v3->GetYaxis()->SetTitleSize(0.05);
  a4v3->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a4v3.png");

  a4v5->GetXaxis()->SetTitle("A4 [fC]");
  a4v5->GetXaxis()->SetTitleSize(0.05);
  a4v5->GetYaxis()->SetTitle("A5 [fC]");
  a4v5->GetYaxis()->SetTitleSize(0.05);
  a4v5->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a4v5.png");

  a5v3->GetXaxis()->SetTitle("A5 [fC]");
  a5v3->GetXaxis()->SetTitleSize(0.05);
  a5v3->GetYaxis()->SetTitle("A3 [fC]");
  a5v3->GetYaxis()->SetTitleSize(0.05);
  a5v3->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a5v3.png");
  
  fout->cd();
  fout->Write();
  fout->Close();
}
