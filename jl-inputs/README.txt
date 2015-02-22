#How to prepare an input card

#point to text file containing a list of files to run over
in_filelist jl-inputs/Files_Ttbar.txt

#point to output ROOT file
Output_File Output_Ttbar_noPU_all.root

#how many entries to run over; just like cmsDriver/Run, -1 means all entries
Entries -1

#point to directory to put performance plots
Plot_Dir Ttbar_noPU_all

#which part of the hcal:
#0=all, 1=Barrel, 2=Endcap
Region 0

#what pileup condition?
#used for the quantile-based baseline/pedestal correction
#options are 0, 25, 50 (0PU, 25ns/20PU, 50ns/20PU)
Condition 0

#how to deal with pedestal/baseline:
#0=pedestal subtraction only, 
#1=P.S. + baseline subtraction from avg[min[TS,Threshold]],
#2=no pedestal subtraction: avg[TS]: use to estimate noise for strategy 4,
#3=no pedestal subtraction: avg[min[TS,Threshold]],
#4=no pedestal subtraction: quantile method at Q=Quantile
Baseline 4

#which time slew parameterization
#0=Test stand, medium WP
#need to implement more params.
Time_Slew 0

#do we care about negative charge values?
#0=No
#need to implement things where we do care
Neg_Charges 0

#threshold in fC for pedestal/baseline strategies 1,3
#set to 0.0 for other strategies (or at least know it's not used)
Threshold 0.0

#quantile for pedestal/baseline strategy 4
#set to 0.0 for other strategies (or at least know it's not used)
Quantile 0.20
