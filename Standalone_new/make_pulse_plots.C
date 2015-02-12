#include <iostream>
#include <string>
#include <sstream>
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TH2F.h"

string int2string(int i) {
  stringstream ss;
  string ret;
  ss << i;
  ss >> ret;
  return ret;
}

// bool isSingleTS(int ts, float digi[10]){
//   
// }

void make_pulse_plots(){
//   example of macro to read data from an ascii file and
//   create a root file with an histogram and an ntuple.
 
   gROOT->Reset();
   
   gStyle->SetOptStat(1);

   // Open the text file
   FILE *fp = fopen("log.txt","r");
   FILE *fp2 = fopen("failures_HI_data.txt","r");

   Float_t chi2,time1,energy1, edm;
   Float_t time2, energy2, time3, energy3, pedestal;
   Int_t ncols, ncols2, ncols3, ncols4, ncols5, ncols6;
   
   Float_t chi22,time21,energy21, edm2;
   Float_t time22, energy22, time23, energy23, pedestal2;
   Int_t ncols7, ncols8, ncols9, ncols10, ncols11, ncols12;
   
   Int_t nlines = 0;
   Float_t digi[10], pulse_1[10], pulse_2[10], pulse_3[10];
   Float_t digi2[10], pulse2_1[10], pulse2_2[10], pulse2_3[10];
   
   // Start a new root file
   TFile *f = new TFile("basic_HI.root","RECREATE");
   TCanvas *lC0 = new TCanvas("canvas","canvas",800,600); 

int nHistos = 20;

  TH1F **hPulse1 = new TH1F*[nHistos];
  TH1F **hPulse2 = new TH1F*[nHistos];
  TH1F **hPulse3 = new TH1F*[nHistos];
  TH1F **hData = new TH1F*[nHistos];
  TH1F **hPedestals = new TH1F*[nHistos];
  TH1F **hSum = new TH1F*[nHistos];
  
  TH1F **hPulse2_1 = new TH1F*[nHistos];
  TH1F **hPulse2_2 = new TH1F*[nHistos];
  TH1F **hPulse2_3 = new TH1F*[nHistos];
  TH1F **hData2 = new TH1F*[nHistos];
  TH1F **hPedestals2 = new TH1F*[nHistos];
  TH1F **hSum2 = new TH1F*[nHistos];
  
  TH2F *hEvsTime = new TH2F("hEvsTime","hEvsTime", 100, 0, 100, 100, -50, 50);
  TH1F *hFailureStats = new TH1F("hFailureStats","hFailureStats",6,0,6);
  
  // Some plots of EDM and energy of failed rechits
  TH2F *hEdmVsEnergy = new TH2F("hEdmVsEnergy","hEdmVsEnergy",200, 0, 200, 50, 0, 50);
  TH1F *hEnergyFailures = new TH1F("hEnergyFailures","hEnergyFailures", 50, 0, 100);
  TH1F *hEdmFailures = new TH1F("hEdmFailures","hEdmFailures", 10000, 0, 0.01);
  TH1F *hTimeFailures = new TH1F("hTimeFailures","hTimeFailures",100, -50, 50);
  
  TH1F *hEnergyNoFail = new TH1F("hEnergyNoFail","hEnergyNoFail", 50, 0, 100);
  TH1F *hEdmNoFail = new TH1F("hEdmNoFail","hEdmNoFail", 10000, 0, 0.01);
  TH1F *hTimeNoFail = new TH1F("hTimeNoFail","hTimeNoFail",100, -50, 50);
  
  TH2F *hEnergyCorrelation = new TH2F("hEnergyCorrelation","hEnergyCorrelation",200, 0, 200, 200, 0, 200);
  
  TH1F *hTiming_1 = new TH1F("hTiming_1","hTiming_13.5mean", 50, -40, 10);
  TH1F *hTiming_2 = new TH1F("hTiming_2","hTiming_14mean", 50, -40, 10);
  TH1F *hEnergy = new TH1F("hEnergy","hEnergy", 1000, -0.5, 0.5);
  TH1F *hEnergyCheck = new TH1F("hEnergyCheck","hEnergyCheck", 100, -0.5, 20);
  TH1F *hTimingOfZeros = new TH1F("hTimingOfZeros","hTimingOfZeros", 100, -100, 75);
  TH1F *hEnergyOtherPulsesLate = new TH1F("hEnergyOtherPulsesLate","hEnergyOtherPulsesLate", 100, -0.5, 20);
  TH1F *hEnergyOtherPulsesEarly = new TH1F("hEnergyOtherPulsesEarly","hEnergyOtherPulsesEarly", 100, -0.5, 20);
  
  TH1F *hTimingOfEarlyPulse = new TH1F("hTimingOfEarlyPulse","hTimingOfEarlyPulse", 100, -100, 75);
  TH1F *hTimingOfLatePulse = new TH1F("hTimingOfLatePulse","hTimingOfLatePulse", 100, -100, 75);
  TH1F *hRatioTS5toTS4 = new TH1F("hRatioTS5toTS4","hRatioTS5toTS4",100, 0, 10);  
  TH2F *hRatioToFindSpike = new TH2F("hRatioToFindSpike","hRatioToFindSpike", 20, 0, 1, 20, 0, 1);
  
  
  // loop through the possible layers & initialize histograms for HE
  for(int i1 = 0; i1 < nHistos; i1++) {
    std::stringstream pST; pST << "Event_" << i1+1;
    hPulse1[i1] = new TH1F(("Pulse_1_" + pST.str()).c_str(),("Pulse_1_" + pST.str()).c_str(),10,-100,150);
    hPulse2[i1] = new TH1F(("Pulse_2_" + pST.str()).c_str(),("Pulse_2_" + pST.str()).c_str(),10,-100,150);
    hPulse3[i1] = new TH1F(("Pulse_3_" + pST.str()).c_str(),("Pulse_3_" + pST.str()).c_str(),10,-100,150);
    hData[i1] = new TH1F(("Data_" + pST.str()).c_str(),("Data_" + pST.str()).c_str(),10,-100,150);
    hPedestals[i1] = new TH1F(("Peds_" + pST.str()).c_str(),("Peds_" + pST.str()).c_str(),10,-100,150);
    hSum[i1] = new TH1F(("Sum_" + pST.str()).c_str(),("Sum_" + pST.str()).c_str(),10,-100,150);
    
    hPulse2_1[i1] = new TH1F(("Small_Pulse_1_" + pST.str()).c_str(),("Small_Pulse_1_" + pST.str()).c_str(),10,-100,150);
    hPulse2_2[i1] = new TH1F(("Small_Pulse_2_" + pST.str()).c_str(),("Small_Pulse_2_" + pST.str()).c_str(),10,-100,150);
    hPulse2_3[i1] = new TH1F(("Small_Pulse_3_" + pST.str()).c_str(),("Small_Pulse_3_" + pST.str()).c_str(),10,-100,150);
    hData2[i1] = new TH1F(("Small_Data_" + pST.str()).c_str(),("Small_Data_" + pST.str()).c_str(),10,-100,150);
    hPedestals2[i1] = new TH1F(("Small_Peds_" + pST.str()).c_str(),("Small_Peds_" + pST.str()).c_str(),10,-100,150);
    hSum2[i1] = new TH1F(("Sum_" + pST.str()).c_str(),("Sum_" + pST.str()).c_str(),10,-100,150);
  }

   Int_t fillArray = 0;
   
   TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","chi2:time:energy:time2:energy2:time3:energy3:pedestal");

   int    nFailed  = 0;
   double nFailedE = 0;
   
   float rechitTime;
   float rechitE;
   
   int code1=0;
   int code2=0;
   int code3=0;
   int code4=0;
   
   while (1) {
      ncols = fscanf(fp,"%f %f %f %f %f %f %f %f ",&chi2, &edm, &time1, &energy1, &time2, &energy2, &time3, &energy3);
      ncols2 = fscanf(fp,"%f %f %f %f %f %f %f %f %f %f",  &digi[0], &digi[1], &digi[2], &digi[3], &digi[4], &digi[5], &digi[6], &digi[7], &digi[8], &digi[9]);
      ncols3 = fscanf(fp,"%f %f %f %f %f %f %f %f %f %f", &pulse_1[0], &pulse_1[1], &pulse_1[2], &pulse_1[3], &pulse_1[4], &pulse_1[5], &pulse_1[6], &pulse_1[7], &pulse_1[8], &pulse_1[9]);
      ncols4 = fscanf(fp, "%f %f %f %f %f %f %f %f %f %f", &pulse_2[0], &pulse_2[1], &pulse_2[2], &pulse_2[3], &pulse_2[4], &pulse_2[5], &pulse_2[6], &pulse_2[7], &pulse_2[8], &pulse_2[9]);
      ncols5 = fscanf(fp, "%f %f %f %f %f %f %f %f %f %f", &pulse_3[0], &pulse_3[1], &pulse_3[2], &pulse_3[3], &pulse_3[4], &pulse_3[5], &pulse_3[6], &pulse_3[7], &pulse_3[8], &pulse_3[9]);
      ncols6 = fscanf(fp, "%f", &pedestal);
      
      ncols7 = fscanf(fp2,"%f %f %f %f %f %f %f %f",&chi22, &edm2, &time21, &energy21, &time22, &energy22, &time23, &energy23);
      ncols8 = fscanf(fp2,"%f %f %f %f %f %f %f %f %f %f",  &digi2[0], &digi2[1], &digi2[2], &digi2[3], &digi2[4], &digi2[5], &digi2[6], &digi2[7], &digi2[8], &digi2[9]);
      ncols9 = fscanf(fp2,"%f %f %f %f %f %f %f %f %f %f", &pulse2_1[0], &pulse2_1[1], &pulse2_1[2], &pulse2_1[3], &pulse2_1[4], &pulse2_1[5], &pulse2_1[6], &pulse2_1[7], &pulse2_1[8], &pulse2_1[9]);
      ncols10 = fscanf(fp2, "%f %f %f %f %f %f %f %f %f %f", &pulse2_2[0], &pulse2_2[1], &pulse2_2[2], &pulse2_2[3], &pulse2_2[4], &pulse2_2[5], &pulse2_2[6], &pulse2_2[7], &pulse2_2[8], &pulse2_2[9]);
      ncols11 = fscanf(fp2, "%f %f %f %f %f %f %f %f %f %f", &pulse2_3[0], &pulse2_3[1], &pulse2_3[2], &pulse2_3[3], &pulse2_3[4], &pulse2_3[5], &pulse2_3[6], &pulse2_3[7], &pulse2_3[8], &pulse2_3[9]);
      ncols12 = fscanf(fp2, "%f", &pedestal2);
      
      if (ncols < 0 || ncols2 < 0 || ncols3 < 0 || ncols4 < 0 || ncols5 < 0 || ncols6 < 0 || ncols7 < 0 || ncols8 < 0 || ncols9 < 0 || ncols10 < 0 || ncols11 < 0 || ncols12 < 0) break;    
      if (time == -999. || time == -9999.) break;
      if (nlines < 5)  printf("x=%i, y=%8f, z=%8f\n",nlines,time1,energy1);
        
    //  if(nlines > 80000) break;
      
      
      hTiming_1->Fill(time1);
      hTiming_2->Fill(time21);
      
      if(chi2==1) code1++;
      if(chi2==2) code2++;
      if(chi2==3) code3++;
      if(chi2==4) code4++;
      
      hEnergy->Fill(energy1); 
      if(energy1 < 0.0001) {
        hEnergyCheck->Fill(digi[4]);
        hTimingOfZeros->Fill(time1);
        hTimingOfLatePulse->Fill(time2);
        hTimingOfEarlyPulse->Fill(time3);
        hEnergyOtherPulsesLate->Fill(energy2);
        hEnergyOtherPulsesEarly->Fill(energy3);
      }
      
      hEdmVsEnergy->SetBinContent(energy1+1, edm+1, 1+hEdmVsEnergy->GetBinContent(energy1+1, edm+1));
      
      
      
      if(chi2 >0){
        double approxEnergy = digi[4] + digi[5];
        hEnergyFailures->Fill(energy1);
        hEdmFailures->Fill(edm);
        hEnergyCorrelation->SetBinContent(approxEnergy+1, energy1+1, 1+hEnergyCorrelation->GetBinContent(approxEnergy+1, energy1+1));
        hTimeFailures->Fill(time1);
      }
      if(chi2 == 0){
        hEnergyNoFail->Fill(energy1);
        hEdmNoFail->Fill(edm);
        hTimeNoFail->Fill(time1);
      }
      
      
      // The key to fit failures: 
      /*
      status = 1    : Covariance was made pos defined
      status = 2    : Hesse is invalid
      status = 3    : Edm is above max 
      status = 4    : Reached call limit
      status = 5    : Any other failure 
      */
      
      hFailureStats->Fill(chi2);
      rechitTime = time1;
      rechitE = energy1;
      
      
      // Search for the single-TS events
      
      
      if(chi2 > 0 && digi[4] > digi[5] && digi[4] > digi[3]) hRatioToFindSpike->Fill(fabs(digi[5]*20/digi[4]), fabs(digi[3]*20/digi[4]), 1+hRatioToFindSpike->GetBinContent(fabs(digi[5]*20/digi[4]), fabs(digi[3]*20/digi[4])));
      
      if(fillArray<nHistos) {
       // if(chi2 == 0 ){
          //printf("chi2 big window=%8f, chi2 small window=%8f\n",chi2,chi22);
          
          // find the average pulse energy
          ++nFailed;
          nFailedE+=rechitE;
          
         // if(digi[4] > digi[5] && digi[4] > digi[3]) hRatioToFindSpike->Fill(fabs(digi[5]*20/digi[4]), fabs(digi[3]*20/digi[4]), 1+hRatioToFindSpike->GetBinContent(fabs(digi[5]*20/digi[4]), fabs(digi[3]*20/digi[4])));
          
           if(digi[4] != 0) hRatioTS5toTS4->Fill(fabs(digi[5]/digi[4]));
          
          for(int j = 0; j < 10; j++){
            hPulse1[fillArray]->SetBinContent(j+1, pulse_1[j]);
            hPulse2[fillArray]->SetBinContent(j+1, pulse_2[j]);
            hPulse3[fillArray]->SetBinContent(j+1, pulse_3[j]);
            hData[fillArray]->SetBinContent(j+1, digi[j]);
            hPedestals[fillArray]->SetBinContent(j+1, pedestal);
            hSum[fillArray]->SetBinContent(j+1, pulse_1[j]+pulse_2[j]+pulse_3[j]+pedestal);
            
            
            hPulse2_1[fillArray]->SetBinContent(j+1, pulse2_1[j]);
            hPulse2_2[fillArray]->SetBinContent(j+1, pulse2_2[j]);
            hPulse2_3[fillArray]->SetBinContent(j+1, pulse2_3[j]);
            hData2[fillArray]->SetBinContent(j+1, digi2[j]);
            hPedestals2[fillArray]->SetBinContent(j+1, pedestal2);
            hSum2[fillArray]->SetBinContent(j+1, pulse2_1[j]+pulse2_2[j]+pulse2_3[j]+pedestal2);
          }
          fillArray++;
         // printf("Even # %i time1 = %8f energy1 = %8f time2 = %8f energy2 = %8f time3 = %8f energy3 = %8f\n", fillArray, time1, energy1, time2, energy2, time3, energy3);
          //printf("Even # %i time1 = %8f energy1 = %8f time2 = %8f energy2 = %8f time3 = %8f energy3 = %8f\n", fillArray, time21, energy21, time22, energy22, time23, energy23);
        //}
      }
      hEvsTime->SetBinContent(energy1, time1+50, 1+hEvsTime->GetBinContent(energy1, time1+50));
      ntuple->Fill(chi2,edm,time1,energy1,time2, energy2, time3, energy3, pedestal);
     // printf("on line %i\n", nlines);
      nlines++;
   }

   printf(" found %d points\n",nlines);
   printf("Average Energy of Failed Fits = %4f\n", nFailedE/nFailed);
   printf("Number with Error codes: 1 has %i, 2 has %i, 3 has %i, 4 has %i\n", code1, code2, code3, code4);
   printf("Percentage of Fits with Failure status: %4f\n", (code1+code2+code3+code4)*100/nlines);
   
   // now we want to overlay the histograms & save as
  for(int i = 0; i < nHistos; i++){
    
hData[i]->SetFillColor(15);
    hData[i]->SetFillStyle(3002);
    hData[i]->SetLineColor(1);
    hData[i]->SetLineWidth(2);
    hPulse1[i]->SetLineColor(2);
    hPulse1[i]->SetLineWidth(2);
    hPulse2[i]->SetLineColor(3);
    hPulse2[i]->SetLineWidth(2);
    hPulse3[i]->SetLineColor(4);
    hPulse3[i]->SetLineWidth(2);
    hPedestals[i]->SetLineColor(5);
    hPedestals[i]->SetLineWidth(2);
    hSum[i]->SetLineColor(9);
    hSum[i]->SetLineWidth(2);
    
    hData[i]->SetTitle(("Event " + int2string(i+1)).c_str());
    hData[i]->GetXaxis()->SetTitle("Time [ns]");
    hData[i]->GetYaxis()->SetTitle("Amplitude [GeV]");

    hData[i]->Draw();
    hPedestals[i]->Draw("same");
    hSum[i]->Draw("same");
    hPulse1[i]->Draw("same");
    hPulse2[i]->Draw("same");
    hPulse3[i]->Draw("same");
 
    
    TLegend* leg = new TLegend(0.6, 0.65, 0.9, 0.9);
   // SetLegendStyle(leg);
    leg->AddEntry(hPulse1[i], "In-Time Pulse" , "f" );
    leg->AddEntry(hPulse2[i], "Late OOT Pulse" , "f" );
    leg->AddEntry(hPulse3[i], "Early OOT Pulse" , "f" );
    leg->AddEntry(hPedestals[i], "Pedestal" , "f" );
    leg->AddEntry(hSum[i], "Sum of Pulses + Ped" , "f" );
    leg->AddEntry(hData[i], "Data" , "f" );
    leg->Draw();
    
    //leg->Draw();
    
    lC0->SaveAs(("Plots/Failed_EventNo_"+int2string(i+1)+".png").c_str());
   // lC0->Clear();
  }
  
  lC0->Clear();
  hEnergyCorrelation->SetTitle("Energy Comparison, Failed Fit from Minuit");
  hEnergyCorrelation->GetXaxis()->SetTitle("Energy TS4 + TS5 [GeV]");
  hEnergyCorrelation->GetYaxis()->SetTitle("Energy from Fit [GeV]");
  hEnergyCorrelation->Draw("colz");
  lC0->SaveAs("Plots/Failed_Energy_Correlation.png");
  
  
//   for(int i = 0; i < nHistos; i++){
//     
//     hData2[i]->SetLineColor(1);
//     hData2[i]->SetLineWidth(2);
//     hPulse2_1[i]->SetLineColor(2);
//     hPulse2_1[i]->SetLineWidth(2);
//     hPulse2_2[i]->SetLineColor(3);
//     hPulse2_2[i]->SetLineWidth(2);
//     hPulse2_3[i]->SetLineColor(4);
//     hPulse2_3[i]->SetLineWidth(2);
//     hPedestals2[i]->SetLineColor(5);
//     hPedestals2[i]->SetLineWidth(2);
//     hSum2[i]->SetLineColor(9);
//     hSum2[i]->SetLineWidth(2);
//     
//     hData2[i]->SetTitle(("Event " + int2string(i+1)).c_str());
// 
//     hData2[i]->Draw();
//     hPulse2_1[i]->Draw("same");
//     hPulse2_2[i]->Draw("same");
//     hPulse2_3[i]->Draw("same");
//     hPedestals2[i]->Draw("same");
//     hSum2[i]->Draw("same");
//     
//     //leg->Draw();
//     
//     lC0->SaveAs(("Plots/newWindow_EventNo_"+int2string(i+1)+".png").c_str());
//    // lC0->Clear();
//   }
  
  for(int i = 0; i < nHistos; i++) {
      delete hData[i];
      delete hPulse1[i];
      delete hPulse2[i];
      delete hPulse3[i];
      delete hPedestals[i];
    }
    
  for(int i = 0; i < nHistos; i++) {
      delete hData2[i];
      delete hPulse2_1[i];
      delete hPulse2_2[i];
      delete hPulse2_3[i];
      delete hPedestals2[i];
    }
    
  delete hData;
  delete hPulse1;
  delete hPulse2;
  delete hPulse3;
  delete hPedestals;
  
  delete hData2;
  delete hPulse2_1;
  delete hPulse2_2;
  delete hPulse2_3;
  delete hPedestals2;
   
   fclose(fp);

   f->Write();
}
