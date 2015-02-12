//
{
/*{
  
  if (gSystem->Getenv("CMSSW_VERSION")) {
    
    TString rfitpath("/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include");
    TString path = gSystem->GetIncludePath();
    path += "-I. -I$ROOTSYS/src -I";
    path += rfitpath;
    gSystem->SetIncludePath(path.Data());
 
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }
  }*/

  gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  gSystem->AddIncludePath("-I$CMSSW_RELEASE_BASE/src/");
  gInterpreter->AddIncludePath(TString(gSystem->Getenv("CMSSW_BASE"))+"/src/");
  gInterpreter->AddIncludePath(TString(gSystem->Getenv("CMSSW_RELEASE_BASE"))+"/src/");



  TString cmsswbase = getenv("CMSSW_BASE");
  if (cmsswbase.Length() > 0) {
   
       // The CMSSW environment is defined (this is true even for FW Lite)
           // so set up the rest.
               //
                   cout << "Loading FW Lite setup." << endl;
                       gSystem->Load("libFWCoreFWLite.so");
                           AutoLibraryLoader::enable();
                               gSystem->Load("libDataFormatsFWLite.so");
                                   gSystem->Load("libDataFormatsPatCandidates.so");
                     }


 // gSystem->Load("libCMSAnaDataTree.so");
 // gSystem->Load("libCMSAnaUtils.so");
 // gSystem->Load("libCMSAnaJetEnergyCorrections.so");
//  gSystem->Load("$HOME/Delphes-3.0.9/libDelphes.so");
//  gSystem->Load("libFWCoreFWLite.so");
//  AutoLibraryLoader::enable();
}
