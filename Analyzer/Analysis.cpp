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

  int nSeq[5]={0,0,0,0,0};

  for (int jentry=0; jentry<nentries;jentry++) {
    //if(nevents%10==0) cout<<" "<<nevents<<"\t out of  "<<nentries<<"\t have already been processed ("<<round_nplaces((double)nevents/nentries*100,1)<<"%/100%)"<<endl;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    nevents++;
    MakeCutflow();
    //MakePedestalPlots(nSeq);
    FillHistograms();
  }

  /*  TCanvas *c_tstot_hb_fit = new TCanvas("c_tstot_hb_fit");
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
  c_t_he_fit->SaveAs("PULSE_ARRIVAL_HE_FIT.png");*/
  
  TCanvas *c1 = new TCanvas("c1");
  gStyle->SetOptStat    (0);
  /*
  hPedSub1->GetXaxis()->SetTitle("Baseline Energy M45 [GeV]");
  hPedSub1->GetXaxis()->SetTitleSize(0.05);
  hPedSub1->GetYaxis()->SetTitle("Counts");
  hPedSub1->GetYaxis()->SetTitleSize(0.05);
  //c1->SetLogy();
  hPedSub1->Draw("");
  c1->SaveAs("pedSub1.png");

  hPedSub2->GetXaxis()->SetTitle("Baseline Energy M012 [GeV]");
  hPedSub2->GetXaxis()->SetTitleSize(0.05);
  hPedSub2->GetYaxis()->SetTitle("Counts");
  hPedSub2->GetYaxis()->SetTitleSize(0.05);
  hPedSub2->Draw("");
  c1->SaveAs("pedSub2.png");

  hPedDiff->GetXaxis()->SetTitle("Diff. BE (M45-M012) [GeV]");
  hPedDiff->GetXaxis()->SetTitleSize(0.05);
  hPedDiff->GetYaxis()->SetTitle("Counts");
  hPedDiff->GetYaxis()->SetTitleSize(0.05);
  hPedDiff->Draw("");
  c1->SaveAs("pedDiff.png");
  */
  h45vHLT0->GetXaxis()->SetTitle("Energy Deposition TS4+5  [GeV]");
  h45vHLT0->GetXaxis()->SetTitleSize(0.05);
  h45vHLT0->GetYaxis()->SetTitle("Energy Deposition from HLT m0 [GeV]");
  h45vHLT0->GetYaxis()->SetTitleSize(0.05);
  h45vHLT0->Draw("colz");
  c1->SaveAs("h45vHLT0.png");

  h45vHLT1->GetXaxis()->SetTitle("Energy Deposition TS4+5  [GeV]");
  h45vHLT1->GetXaxis()->SetTitleSize(0.05);
  h45vHLT1->GetYaxis()->SetTitle("Energy Deposition from HLT m1 [GeV]");
  h45vHLT1->GetYaxis()->SetTitleSize(0.05);
  h45vHLT1->Draw("colz");
  c1->SaveAs("h45vHLT1.png");

  h45vHLT2->GetXaxis()->SetTitle("Energy Deposition TS4+5  [GeV]");
  h45vHLT2->GetXaxis()->SetTitleSize(0.05);
  h45vHLT2->GetYaxis()->SetTitle("Energy Deposition from HLT m2 [GeV]");
  h45vHLT2->GetYaxis()->SetTitleSize(0.05);
  h45vHLT2->Draw("colz");
  c1->SaveAs("h45vHLT2.png");

  p45vHLT0->GetXaxis()->SetTitle("Energy Deposition TS4+5 [GeV]");
  p45vHLT0->GetXaxis()->SetTitleSize(0.05);
  p45vHLT0->GetYaxis()->SetTitle("(E4+5 - HLT)/E4+5");
  p45vHLT0->GetYaxis()->SetTitleSize(0.05);
  p45vHLT0->Draw();
  c1->SaveAs("p45vHLT0.png");

  p45vHLT1->GetXaxis()->SetTitle("Energy Deposition 4+5 [GeV]");
  p45vHLT1->GetXaxis()->SetTitleSize(0.05);
  p45vHLT1->GetYaxis()->SetTitle("(E4+5 - HLT)/E4+5");
  p45vHLT1->GetYaxis()->SetTitleSize(0.05);
  p45vHLT1->Draw();
  c1->SaveAs("p45vHLT1.png");

  p45vHLT2->GetXaxis()->SetTitle("Energy Deposition 4+5 [GeV]");
  p45vHLT2->GetXaxis()->SetTitleSize(0.05);
  p45vHLT2->GetYaxis()->SetTitle("(E4+5 - HLT)/E4+5");
  p45vHLT2->GetYaxis()->SetTitleSize(0.05);
  p45vHLT2->Draw();
  c1->SaveAs("p45vHLT2.png");

  hM2vHLT0->GetXaxis()->SetTitle("M2 Energy [GeV]");
  hM2vHLT0->GetXaxis()->SetTitleSize(0.05);
  hM2vHLT0->GetYaxis()->SetTitle("Energy Deposition from HLT m0 [GeV]");
  hM2vHLT0->GetYaxis()->SetTitleSize(0.05);
  hM2vHLT0->Draw("colz");
  c1->SaveAs("hM2vHLT0.png");

  hM2vHLT1->GetXaxis()->SetTitle("M2 Energy [GeV]");
  hM2vHLT1->GetXaxis()->SetTitleSize(0.05);
  hM2vHLT1->GetYaxis()->SetTitle("Energy Deposition from HLT m1 [GeV]");
  hM2vHLT1->GetYaxis()->SetTitleSize(0.05);
  hM2vHLT1->Draw("colz");
  c1->SaveAs("hM2vHLT1.png");

  hM2vHLT2->GetXaxis()->SetTitle("M2 Energy [GeV]");
  hM2vHLT2->GetXaxis()->SetTitleSize(0.05);
  hM2vHLT2->GetYaxis()->SetTitle("Energy Deposition from HLT m2 [GeV]");
  hM2vHLT2->GetYaxis()->SetTitleSize(0.05);
  hM2vHLT2->Draw("colz");
  c1->SaveAs("hhM2vHLT2.png");

  pM2vHLT0->GetXaxis()->SetTitle("M2 Energy [GeV]");
  pM2vHLT0->GetXaxis()->SetTitleSize(0.05);
  pM2vHLT0->GetYaxis()->SetTitle("(M2 - HLT)/(M2)");
  pM2vHLT0->GetYaxis()->SetTitleSize(0.05);
  pM2vHLT0->Draw();
  c1->SaveAs("pM2vHLT0.png");

  pM2vHLT1->GetXaxis()->SetTitle("M2 Energy [GeV]");
  pM2vHLT1->GetXaxis()->SetTitleSize(0.05);
  pM2vHLT1->GetYaxis()->SetTitle("(M2 - HLT)/(M2)");
  pM2vHLT1->GetYaxis()->SetTitleSize(0.05);
  pM2vHLT1->Draw();
  c1->SaveAs("pM2vHLT1.png");

  pM2vHLT2->GetXaxis()->SetTitle("M2 Energy [GeV]");
  pM2vHLT2->GetXaxis()->SetTitleSize(0.05);
  pM2vHLT2->GetYaxis()->SetTitle("(M2 - HLT)/(M2)");
  pM2vHLT2->GetYaxis()->SetTitleSize(0.05);
  pM2vHLT2->Draw();
  c1->SaveAs("pM2vHLT2.png");

  /*
  TH1D *hEdepDist_all;
  TH1D *hEdepDist_not3456;
  TH1D *hEdepDist_not345;
  TH1D *hEdepDist_least;
  TH1D *hEdepDist_least4;
  TH1D *hEdepDist_least_not3456;
  TH1D *hEdepDist_least_not345;
  TH1D *hEdepDist_least4_not3456;
  TH1D *hEdepDist_least4_not345;*/
    
  /*  TCanvas *cEdep = new TCanvas("cEdep");
  cEdep->SetLogy();

  //nothing excluded
  hEdepDist_all->GetXaxis()->SetTitle("Energy deposition per TS [GeV]");
  hEdepDist_all->GetXaxis()->SetTitleSize(0.05);
  hEdepDist_all->GetYaxis()->SetTitle("Counts");
  hEdepDist_all->GetYaxis()->SetTitleSize(0.05);
  //hEdepDist->GetYaxis()->SetRangeUser(0,100);
  hEdepDist_all->Draw();
  cEdep->SaveAs("edep_all.png");

  //excluding TS3,4,5,6
  hEdepDist_not3456->GetXaxis()->SetTitle("Energy deposition per TS, not TS3/4/5/6 [GeV]");
  hEdepDist_not3456->GetXaxis()->SetTitleSize(0.05);
  hEdepDist_not3456->GetYaxis()->SetTitle("Counts");
  hEdepDist_not3456->GetYaxis()->SetTitleSize(0.05);
  //hEdepDist->GetYaxis()->SetRangeUser(0,100);
  hEdepDist_not3456->Draw();
  cEdep->SaveAs("edep_not3456.png");

  //excluding TS3,4,5
  hEdepDist_not345->GetXaxis()->SetTitle("Energy deposition per TS, not TS3/4/5 [GeV]");
  hEdepDist_not345->GetXaxis()->SetTitleSize(0.05);
  hEdepDist_not345->GetYaxis()->SetTitle("Counts");
  hEdepDist_not345->GetYaxis()->SetTitleSize(0.05);
  //hEdepDist->GetYaxis()->SetRangeUser(0,100);
  hEdepDist_not345->Draw();
  cEdep->SaveAs("edep_not345.png");

  //excluding TS4,5
  hEdepDist_not45->GetXaxis()->SetTitle("Energy deposition per TS, not TS4/5 [GeV]");
  hEdepDist_not45->GetXaxis()->SetTitleSize(0.05);
  hEdepDist_not45->GetYaxis()->SetTitle("Counts");
  hEdepDist_not45->GetYaxis()->SetTitleSize(0.05);
  //hEdepDist->GetYaxis()->SetRangeUser(0,100);
  hEdepDist_not45->Draw();
  cEdep->SaveAs("edep_not45.png");

  //only least TS per 10 TS's
  hEdepDist_least->GetXaxis()->SetTitle("Energy deposition for least TS of 10 [GeV]");
  hEdepDist_least->GetXaxis()->SetTitleSize(0.05);
  hEdepDist_least->GetYaxis()->SetTitle("Counts");
  hEdepDist_least->GetYaxis()->SetTitleSize(0.05);
  //hEdepDist->GetYaxis()->SetRangeUser(0,100);
  hEdepDist_least->Draw();
  cEdep->SaveAs("edep_least.png");

  //only least 4 TS per 10 TS's
  hEdepDist_least4->GetXaxis()->SetTitle("Energy deposition for least 4 TS of 10 [GeV]");
  hEdepDist_least4->GetXaxis()->SetTitleSize(0.05);
  hEdepDist_least4->GetYaxis()->SetTitle("Counts");
  hEdepDist_least4->GetYaxis()->SetTitleSize(0.05);
  //hEdepDist->GetYaxis()->SetRangeUser(0,100);
  hEdepDist_least4->Draw();
  cEdep->SaveAs("edep_least4.png");

  //only least TS per 10 TS's excluding TS3,4,5,6
  hEdepDist_least_not3456->GetXaxis()->SetTitle("Energy deposition for least TS of 10, not TS3/4/5/6 [GeV]");
  hEdepDist_least_not3456->GetXaxis()->SetTitleSize(0.05);
  hEdepDist_least_not3456->GetYaxis()->SetTitle("Counts");
  hEdepDist_least_not3456->GetYaxis()->SetTitleSize(0.05);
  //hEdepDist->GetYaxis()->SetRangeUser(0,100);
  hEdepDist_least_not3456->Draw();
  cEdep->SaveAs("edep_least_not3456.png");

  //only least TS per 10 TS's excluding TS3,4,5
  hEdepDist_least_not345->GetXaxis()->SetTitle("Energy deposition for least TS of 10, not TS3/4/5 [GeV]");
  hEdepDist_least_not345->GetXaxis()->SetTitleSize(0.05);
  hEdepDist_least_not345->GetYaxis()->SetTitle("Counts");
  hEdepDist_least_not345->GetYaxis()->SetTitleSize(0.05);
  //hEdepDist->GetYaxis()->SetRangeUser(0,100);
  hEdepDist_least_not345->Draw();
  cEdep->SaveAs("edep_least_not345.png");

  //only least TS per 10 TS's excluding TS4,5
  hEdepDist_least_not45->GetXaxis()->SetTitle("Energy deposition for least TS of 10, not TS4/5 [GeV]");
  hEdepDist_least_not45->GetXaxis()->SetTitleSize(0.05);
  hEdepDist_least_not45->GetYaxis()->SetTitle("Counts");
  hEdepDist_least_not45->GetYaxis()->SetTitleSize(0.05);
  //hEdepDist->GetYaxis()->SetRangeUser(0,100);
  hEdepDist_least_not45->Draw();
  cEdep->SaveAs("edep_least_not45.png");

  //only least 4 TS per 10 TS's excluding TS3,4,5,6
  hEdepDist_least4_not3456->GetXaxis()->SetTitle("Energy deposition for least 4 TS of 10, not TS3/4/5/6 [GeV]");
  hEdepDist_least4_not3456->GetXaxis()->SetTitleSize(0.05);
  hEdepDist_least4_not3456->GetYaxis()->SetTitle("Counts");
  hEdepDist_least4_not3456->GetYaxis()->SetTitleSize(0.05);
  //hEdepDist->GetYaxis()->SetRangeUser(0,100);
  hEdepDist_least4_not3456->Draw();
  cEdep->SaveAs("edep_least4_not3456.png");

  //only least 4 TS per 10 TS's excluding TS3,4,5
  hEdepDist_least4_not345->GetXaxis()->SetTitle("Energy deposition for least 4 TS of 10, not TS3/4/5 [GeV]");
  hEdepDist_least4_not345->GetXaxis()->SetTitleSize(0.05);
  hEdepDist_least4_not345->GetYaxis()->SetTitle("Counts");
  hEdepDist_least4_not345->GetYaxis()->SetTitleSize(0.05);
  //hEdepDist->GetYaxis()->SetRangeUser(0,100);
  hEdepDist_least4_not345->Draw();
  cEdep->SaveAs("edep_least4_not345.png");

  //only least 4 TS per 10 TS's excluding TS4,5
  hEdepDist_least4_not45->GetXaxis()->SetTitle("Energy deposition for least 4 TS of 10, not TS4/5 [GeV]");
  hEdepDist_least4_not45->GetXaxis()->SetTitleSize(0.05);
  hEdepDist_least4_not45->GetYaxis()->SetTitle("Counts");
  hEdepDist_least4_not45->GetYaxis()->SetTitleSize(0.05);
  //hEdepDist->GetYaxis()->SetRangeUser(0,100);
  hEdepDist_least4_not45->Draw();
  cEdep->SaveAs("edep_least4_not45.png");
  */
 }
 
void Analysis::DefineHistograms()
{
  fout = new TFile(Output_File.c_str(), "RECREATE");
  
  // ---------------Plots for HLT ------------
  //RatioPulse = new TH2F("RatioPulse","TS5/TS4 vs TS45",20,0.,2.,100,0.,500.);
  //TimeSlewPulse = new TH2F("TimeSlewPulse","Time Slew vs TS45",25,-14.5,10.5,100,0.,500.);
  
  //Norm0 = new TH1F("fC0","Amplitude in Pulse ealier [fC]",100,0.,500.);
  //Norm1 = new TH1F("fC1","Amplitude in in-time Pulse [fC]",100,0.,500.);
  //Norm2 = new TH1F("fC2","Amplitude in Pulse later [fC]",100,0.,500.);
  
  //slewFit = new TF1("slewFit","pol4*expo(5)",-10.,14.);
  //slewFit->SetParameters(1.07618e-02,-4.19145e-06,2.70310e-05,-8.71584e-08,1.86597e-07,3.59216e+00,-1.02057e-01);
    
  //logtimeslewFit = new TF1("logtimeslewFit", "[0]+[1]*TMath::Log(x+[2])",0.,500.);
  //logtimeslewFit->SetParameters(3.89838e+01, -6.93560e+00, 8.52052e+01);
  
  //exptimeslewFit = new TF1("exptimeslewFit", "[0]+[1]*TMath::Exp([2]*x)",0.,500.);
  //exptimeslewFit->SetParameters(-2.69330e+00, 1.09162e+01, -7.60722e-03);
  
  // Output Plots
  //hHLTResolution=new TProfile("hHLTResolution","",20,0,100,-1.0,1.0);
  //hJayResolution=new TProfile("hJayResolution","",20,0,100,-1.0,1.0);

  Int_t xBins=100;
  Int_t xMin=0;
  Int_t xMax=10;

  //TH2D *h45vHLT0;
  h45vHLT0 = new TH2D("h45vHLT0", "", xBins,xMin,xMax,xBins,xMin,xMax);
  //TH2D *h45vHLT1;
  h45vHLT1 = new TH2D("h45vHLT1", "", xBins,xMin,xMax,xBins,xMin,xMax);
  //TH2D *h45vHLT2;
  h45vHLT2 = new TH2D("h45vHLT2", "", xBins,xMin,xMax,xBins,xMin,xMax);

  //TProfile *p45vHLT0;
  p45vHLT0 = new TProfile("p45vHLT0", "", xBins/10,xMin,xMax,-10,10);
  //TProfile *p45vHLT1;
  p45vHLT1 = new TProfile("p45vHLT1", "", xBins/10,xMin,xMax,-10,10);
  //TProfile *p45vHLT2;
  p45vHLT2 = new TProfile("p45vHLT2", "", xBins/10,xMin,xMax,-10,10);

  //TH2D *hM2vHLT0;
  hM2vHLT0 = new TH2D("hM2vHLT0", "", xBins,xMin,xMax,xBins,xMin,xMax);
  //TH2D *hM2vHLT1;
  hM2vHLT1 = new TH2D("hM2vHLT1", "", xBins,xMin,xMax,xBins,xMin,xMax);
  //TH2D *hM2vHLT2;
  hM2vHLT2 = new TH2D("hM2vHLT2", "", xBins,xMin,xMax,xBins,xMin,xMax);

  //TProfile *pM2vHLT0;
  pM2vHLT0 = new TProfile("pM2vHLT0", "", xBins/10,xMin,xMax,-10,10);
  //TProfile *pM2vHLT1;
  pM2vHLT1 = new TProfile("pM2vHLT1", "", xBins/10,xMin,xMax,-10,10);
  //TProfile *pM2vHLT2;
  pM2vHLT2 = new TProfile("pM2vHLT2", "", xBins/10,xMin,xMax,-10,10);

  //hCharge_Method2_v_HLT=new TH2D("hCharge_Method2_v_HLT","",50,0,250,50,0,250);

  //hCharge_Method2_v_JAY=new TH2D("hCharge_Method2_v_jay","",50,0,250,50,0,250);
  //hCharge_HLT_v_JAY=new TH2D("hCharge_Method2_v_jay","",50,0,250,50,0,250);
  
  PULSE_ARRIVAL_HB_FIT=new TH1D("PULSE_ARRIVAL_HB_FIT","",50,-30.0,30.0);
  PULSE_ARRIVAL_HE_FIT=new TH1D("PULSE_ARRIVAL_HE_FIT","",50,-30.0,30.0);

  CHARGE_TSTOT_HB_FIT=new TH1D("CHARGE_TSTOT_HB_FIT","",150,0,1500);
  CHARGE_TSTOT_HE_FIT=new TH1D("CHARGE_TSTOT_HE_FIT","",150,0,1500);

  hPedSub1 = new TH1D("hPedSub1", "", 100,-1,1);
  hPedSub2 = new TH1D("hPedSub2", "", 100,-1,1);
  hPedDiff = new TH1D("hPedDiff", "", 100,-1,1);

  /*  TH1D *hEdepDist_all;
  TH1D *hEdepDist_not3456;
  TH1D *hEdepDist_not345;
  TH1D *hEdepDist_least;
  TH1D *hEdepDist_least4;
  TH1D *hEdepDist_least_not3456;
  TH1D *hEdepDist_least_not345;
  TH1D *hEdepDist_least4_not3456;
  TH1D *hEdepDist_least4_not345;*/

  hEdepDist_all            =new TH1D("hEdepDist_all", "",50,0,200);
  hEdepDist_not3456        =new TH1D("hEdepDist_not3456", "",50,0,100);
  hEdepDist_not345         =new TH1D("hEdepDist_not345", "",50,0,100);
  hEdepDist_not45          =new TH1D("hEdepDist_not45", "",50,0,100);
  hEdepDist_least          =new TH1D("hEdepDist_least", "",100,-5,5);
  hEdepDist_least4         =new TH1D("hEdepDist_least4", "",100,-5,15);
  hEdepDist_least_not3456  =new TH1D("hEdepDist_least_not3456", "",100,-5,5);
  hEdepDist_least_not345   =new TH1D("hEdepDist_least_not345", "",100,-5,5);
  hEdepDist_least_not45    =new TH1D("hEdepDist_least_not45", "",100,-5,5);
  hEdepDist_least4_not3456 =new TH1D("hEdepDist_least4_not3456", "",100,-5,15);
  hEdepDist_least4_not345  =new TH1D("hEdepDist_least4_not345", "",100,-5,15);
  hEdepDist_least4_not45  =new TH1D("hEdepDist_least4_not45", "",100,-5,15);

  char hname[50];
  for (Int_t i=0; i<10; i++) {
    sprintf(hname,"hPed_%i",i);
    vHistPed.push_back(new TH1D(hname,"",100,-10,10));
    sprintf(hname,"hVal_%i",i);
    vHistVal.push_back(new TH1D(hname,"",100,-10,10));
  }
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

    std::vector<double> correctedOutput, HLTOutput, Jay0Output, Jay1Output, Jay2Output;
    std::vector<double> inputCaloSample, inputPedestal;
    std::vector<double> inputGain;

    // ignore pulses of less than 1 GeV in TS 4+5
    if ( (Charge[j][4]+Charge[j][5])<5 ) continue;
    
    for(int i = 0; i < 10; ++i) {
      
      // Note: In CMSSW the "Charge" vector is not already pedestal subtracted, unlike here
      // so I add the pedestal back to the charge so we can keep the same CMSSW implementation
      inputCaloSample.push_back(Charge[j][i]+Pedestal[j][i]);
      inputPedestal.push_back(Pedestal[j][i]);
      inputGain.push_back(Gain[j][i]);

      //hEdepDist->Fill( (Pedestal[j][i]) );
    }

    // Begin Method 2
    psFitOOTpuCorr_->apply(inputCaloSample,inputPedestal,inputGain,correctedOutput);

    // Begin Xinmei's implementation of HLT
    //double RatioTS54, TimeSlew, Pulse = 0.;
    //hltThing_->apply(inputCaloSample,inputPedestal,inputGain,HLTOutput, RatioTS54, TimeSlew, Pulse, slewFit);

    // Begin Jay's implementation of HLT
    hltv2_->applyOnce(inputCaloSample,inputPedestal,inputGain,Jay0Output);
    hltv2_->applyOnceL4_45(inputCaloSample,inputPedestal,inputGain,Jay1Output);
    hltv2_->applyOnce012(inputCaloSample,inputPedestal,inputGain,Jay2Output);
    
    //RatioPulse->Fill(RatioTS54, Pulse);
    //TimeSlewPulse->Fill(TimeSlew, Pulse);
    //Norm0->Fill(Jay1Output.at(0));
    //Norm1->Fill(Jay1Output.at(1));
    //Norm2->Fill(Jay1Output.at(2));

    // ------------ Fill our comparison histograms -------------------
    // Make sure that the output vectors have non-zero number of entries

    //TH2D *h45vHLT0;
    //TH2D *h45vHLT1;
    //TH2D *h45vHLT2;
    if (correctedOutput.size() > 1) {
      if(Jay0Output.size() > 1){
	h45vHLT0->Fill( (Charge[j][4]+Charge[j][5])*Gain[j][4], Jay0Output.at(1)*Gain[j][4],1);
	hM2vHLT0->Fill( correctedOutput.at(0), Jay0Output.at(1)*Gain[j][4],1);
	p45vHLT0->Fill((Charge[j][4]+Charge[j][5])*Gain[j][4], (Jay0Output.at(1)-(Charge[j][4]+Charge[j][5]))/((Charge[j][4]+Charge[j][5])),1);
	pM2vHLT0->Fill( correctedOutput.at(0), (Jay0Output.at(1)*Gain[j][4]-(correctedOutput.at(0)))/((correctedOutput.at(0))),1);
      }      
      if(Jay1Output.size() > 1){
	h45vHLT1->Fill( (Charge[j][4]+Charge[j][5])*Gain[j][4], Jay1Output.at(1)*Gain[j][4],1);
	hM2vHLT1->Fill( correctedOutput.at(0), Jay1Output.at(1)*Gain[j][4],1);
	p45vHLT1->Fill((Charge[j][4]+Charge[j][5])*Gain[j][4], (Jay1Output.at(1)-(Charge[j][4]+Charge[j][5]))/((Charge[j][4]+Charge[j][5])),1);
	pM2vHLT1->Fill( correctedOutput.at(0), (Jay1Output.at(1)*Gain[j][4]-(correctedOutput.at(0)))/((correctedOutput.at(0))),1);
      }      

      if(Jay2Output.size() > 1){
	h45vHLT2->Fill( (Charge[j][4]+Charge[j][5])*Gain[j][4], Jay2Output.at(1)*Gain[j][4],1);
	hM2vHLT2->Fill( correctedOutput.at(0), Jay2Output.at(1)*Gain[j][4],1);
	p45vHLT2->Fill((Charge[j][4]+Charge[j][5])*Gain[j][4], (Jay2Output.at(1)-(Charge[j][4]+Charge[j][5]))/((Charge[j][4]+Charge[j][5])),1);
	pM2vHLT2->Fill( correctedOutput.at(0), (Jay2Output.at(1)*Gain[j][4]-(correctedOutput.at(0)))/((correctedOutput.at(0))),1);
      }      
      
    }

      /*      // Fill the TProfile with the % difference of energies
      double resolution = 0;
      correctedOutput.at(0) > 0 ? resolution = ( correctedOutput.at(0) - Jay2Output.at(1)*Gain[j][0] )/correctedOutput.at(0) : resolution = 0 ;
      //hJayResolution->Fill(correctedOutput.at(0), resolution, 1);

      // Should do something with the correctedOutput vector here, such as fill a histogram...
      if(IEta[j] < 16 && correctedOutput.size() > 1) {
	if(correctedOutput.at(1) > -99.){
	  CHARGE_TSTOT_HB_FIT->Fill(correctedOutput.at(0));
	  PULSE_ARRIVAL_HB_FIT->Fill(correctedOutput.at(1));
	}
      } else if(IEta[j] >= 16 && correctedOutput.size() > 1){
	CHARGE_TSTOT_HE_FIT->Fill(correctedOutput.at(0));
	PULSE_ARRIVAL_HE_FIT->Fill(correctedOutput.at(1));
	}*/
      //}
  }//PulseCount

}// End MakeCutflow Function

void Analysis::FillHistograms()
{
}

void Analysis::MakePedestalPlots(int *n) {

  *(n+0)+=PulseCount*10;

  for (int j = 0; j < (int)PulseCount; j++) {
    
    std::vector<double> inputCaloSample, inputPedestal;
    std::vector<double> inputGain;

    /*    
    TH1D *hEdepDist_least;
    TH1D *hEdepDist_least4;
    TH1D *hEdepDist_least_not3456;
    TH1D *hEdepDist_least_not345;
    TH1D *hEdepDist_least4_not3456;
    TH1D *hEdepDist_least4_not345;*/

    double all[10];
    double not45[8];
    double not345[7];
    double not3456[6];

    int n=0, m=0, p=0;
    for(int i = 0; i < 10; ++i) {
      
      // Note: In CMSSW the "Charge" vector is not already pedestal subtracted, unlike here
      // so I add the pedestal back to the charge so we can keep the same CMSSW implementation
      inputCaloSample.push_back(Charge[j][i]+Pedestal[j][i]);
      inputPedestal.push_back(Pedestal[j][i]);
      inputGain.push_back(Gain[j][i]);

      vHistPed[i]->Fill(Pedestal[j][i]);
      vHistVal[i]->Fill(Pedestal[j][i]+Charge[j][i]);

      Float_t tempC=Charge[j][i];//+Pedestal[j][i];

      hEdepDist_all->Fill( tempC*Gain[j][i] );

      all[i]=tempC*Gain[j][i];
      
      if (i!=3&&i!=4&&i!=5&&i!=6) {
	hEdepDist_not3456->Fill( tempC*Gain[j][i] );
	not3456[n]=tempC*Gain[j][i];
	n++;
      }
      if (i!=4&&i!=5) {
	hEdepDist_not45->Fill( tempC*Gain[j][i] );
	not45[p]=tempC*Gain[j][i];
	p++;
      }
      if (i!=3&&i!=4&&i!=5) {
	hEdepDist_not345->Fill( tempC*Gain[j][i] );
	not345[m]=tempC*Gain[j][i];
	m++;
      }

      //if (Charge[j][i]+Pedestal[j][i]<0) *(n+1)+=1;
      //if (Charge[j][i]<0) *(n+2)+=1;
      
    }

    //std::sort(std::begin(all), std::end(all));
    std::sort(std::begin(not45), std::end(not45));
    std::sort(std::begin(not345), std::end(not345));
    std::sort(std::begin(not3456), std::end(not3456));
    
    //hEdepDist_least->Fill( all[0] );
    hEdepDist_least_not45->Fill( not45[0] );
    hEdepDist_least_not345->Fill( not345[0] );
    hEdepDist_least_not3456->Fill( not3456[0] );

    Float_t ped45=0, ped012=0;

    for (int i=0; i<4; i++) { ped45+=not45[i]; }
    for (int i=0; i<3; i++) { ped012+=all[i]; }

    ped45=ped45/4;
    ped012=ped012/3;
    
    /*    for (int i=0; i<4; i++) {
      hEdepDist_least4->Fill( all[i] );
      hEdepDist_least4_not45->Fill( not45[i] );
      hEdepDist_least4_not345->Fill( not345[i] );
      hEdepDist_least4_not3456->Fill( not3456[i] );
      }*/

    for (int i=0; i<10; i++) {

      hPedSub1->Fill( ped45 );
      hPedSub2->Fill( ped012 );
      hPedDiff->Fill( ped45 - ped012 );

    }


    //cout << Charge[j][0] << " " << Charge[j][1] << " " << Charge[j][2] << " " << Charge[j][3] << endl;
    //pedSubFxn_->Calculate(inputCaloSample,inputPedestal,inputGain);
    
  }


}

void Analysis::Finish()
{
  /*  gStyle->SetOptFit(1);
  TimeSlewPulse->Draw("BOX");
  TimeSlewPulse->ProfileY("Y Profile",1,-1,"do")->Fit("pol1");
  TimeSlewPulse->ProfileY("X Profile",1,-1,"do")->Fit("pol1");
  */
  fout->cd();
  fout->Write();
  fout->Close();
}
