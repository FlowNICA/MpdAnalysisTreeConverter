// Standard c++ headers
#include <iostream>
#include <string>
//#include <utility>
#include <set>
#include <map>

// Standard ROOT headers
#include <Rtypes.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TClonesArray.h>
#include <TObject.h>

// FairRoot/MpdRoot headers
#include "FairMCEventHeader.h"
#ifdef _MCSTACK_
#include "FairMCTrack.h"
#endif
#ifdef _MPDMCSTACK_
#include "MpdMCTrack.h"
#endif
#include "MpdEvent.h"
#include "MpdZdcDigi.h"
#include "MpdPid.h"
#include "MpdTrack.h"
#include "MpdKalmanTrack.h"
#include "MpdVertex.h"

// AnalysisTree headers
#include "AnalysisTree/Configuration.hpp"
#include "AnalysisTree/Detector.hpp"

int main(int argc, char **argv)
{
  TString iFileName, oFileName;

  if (argc < 5)
  {
    std::cerr << "./MpdDst2AnalysisTree -i inputfile -o outputfile" << std::endl;
    return 1;
  }
  for (int i = 1; i < argc; i++)
  {
    if (std::string(argv[i]) != "-i" &&
      std::string(argv[i]) != "-o")
    {
      std::cerr << "\n[ERROR]: Unknown parameter " << i << ": " << argv[i] << std::endl;
      return 1;
    }
    else
    {
      if (std::string(argv[i]) == "-i" && i != argc - 1)
      {
        iFileName = argv[++i];
        continue;
      }
      if (std::string(argv[i]) == "-i" && i == argc - 1)
      {
        std::cerr << "\n[ERROR]: Input file name was not specified " << std::endl;
        return 1;
      }
      if (std::string(argv[i]) == "-o" && i != argc - 1)
      {
        oFileName = argv[++i];
        continue;
      }
      if (std::string(argv[i]) == "-o" && i == argc - 1)
      {
        std::cerr << "\n[ERROR]: Output file name was not specified " << std::endl;
        return 1;
      }
    }
  }

  TStopwatch timer;
  timer.Start();

  // Open input dst tree
  TFile *fi = new TFile(iFileName.Data(),"read");
  if (!fi || fi->IsZombie())
  {
    std::cerr << "ERROR: Input file probably is empty. Exit the program." << std::endl;
    return 100;
  }

  TTree  *dstTree = (TTree*) fi->Get("mpdsim");

  // PID related parameters
  const Double_t PIDsigM = 4.0;
  const Double_t PIDsigE = 4.0;
  const Double_t PIDenergy = 11.;
  const Double_t PIDkoeff = 1.;
  const TString PIDgenerator = "URQMD";
  const TString PIDtracking = "CF";
  const TString PIDparticles = "pikapr";
  const Int_t   Num_Of_Modules = 90;

  MpdPid *pid = new MpdPid(PIDsigM, PIDsigE, PIDenergy, PIDkoeff, PIDgenerator, PIDtracking, PIDparticles);

  // Prepare to read input dst
  FairMCEventHeader *MCHeader;
  TClonesArray *MCTracks;
  MpdEvent *MPDEvent;
  TClonesArray *FHCalHits;
  TClonesArray *MpdGlobalTracks;
  MpdZdcDigi *FHCalHit;
  //TClonesArray *mpdKalmanTracks;
  TClonesArray *vertexes;

  MCHeader = nullptr;
  MCTracks = nullptr;
  MPDEvent = nullptr;
  FHCalHits = nullptr;
  //mpdKalmanTracks = nullptr;
  vertexes = nullptr;

  dstTree->SetBranchAddress("MCEventHeader.", &MCHeader);
  dstTree->SetBranchAddress("MCTrack", &MCTracks);
  dstTree->SetBranchAddress("MPDEvent.", &MPDEvent);
  dstTree->SetBranchAddress("ZdcDigi", &FHCalHits);
  //dstTree->SetBranchAddress("TpcKalmanTrack", &mpdKalmanTracks);
  dstTree->SetBranchAddress("Vertex", &vertexes);

  // Set up output dst
  TFile *outFile = new TFile(oFileName.Data(), "RECREATE");
  TTree *outTree = new TTree("mpd_analysistree","AnalysisTree Dst at MPD");

  // Set up AnalysisTree configureation
  //AnalysisTree::Configuration *out_config = new AnalysisTree::Configuration;
  AnalysisTree::ModuleDetector *fhcal_modules = new AnalysisTree::ModuleDetector( Num_Of_Modules );

  // Create branches in the output tree
  std::string str_fhcal_branch = "FHCalModules";
  //AnalysisTree::BranchConfig fhcal_branch(str_fhcal_branch.c_str(), AnalysisTree::DetType::kModule); 
  //out_config->AddBranchConfig(std::move(fhcal_branch));
  outTree->Branch(str_fhcal_branch.c_str(), "AnalysisTree::ModuleDetector", &fhcal_modules);

  // Starting event loop
  TVector3 primaryVertex;
  std::set <Int_t> UsedMCTracks;        // using to remap mc-reco track matching
  std::map <Int_t,Int_t> InitMcNewMcId; // map[old-mc-id] = new-mc-id
  Float_t FHCalSumEnergy[Num_Of_Modules];
  Int_t   FHCalNumOfHits[Num_Of_Modules];
  Int_t n_entries = dstTree->GetEntriesFast();
  Bool_t isGoodPID;
  Short_t charge_mpd;

  for (int iEv = 0; iEv < n_entries; iEv++)
  {
    std::cout << "Event [" << iEv << "/" << n_entries << "]" << std::endl;
    dstTree->GetEntry(iEv);

    UsedMCTracks.clear();
    InitMcNewMcId.clear();
    for (int i=0; i<Num_Of_Modules; i++)
    {
      FHCalSumEnergy[i] = 0.;
      FHCalNumOfHits[i] = 0;
    }

    // Read energy in FHCal modules 
    fhcal_modules->Reserve(Num_Of_Modules);
    for (int imodule=0; imodule<Num_Of_Modules; imodule++)
    {
      auto *module = fhcal_modules->AddChannel();
      module->SetSignal(0.f);
    }
    Int_t number_of_FHCal_hits = FHCalHits->GetEntriesFast();
    for(int ihit=0; ihit<number_of_FHCal_hits; ihit++)
    {
      FHCalHit = (MpdZdcDigi*) FHCalHits->UncheckedAt(ihit);
      Int_t DetId = FHCalHit->GetDetectorID();
      Int_t ModId = FHCalHit->GetModuleID()-1;
      Int_t ModNumber = ModId + (Num_Of_Modules/2) * (DetId-1);

      FHCalSumEnergy[ModNumber] += FHCalHit->GetELoss();
      FHCalNumOfHits[ModNumber]++;
    }
    for (int imodule=0; imodule<Num_Of_Modules; imodule++)
    {
      auto& module = fhcal_modules->GetChannel(imodule);
      module.SetNumber(FHCalNumOfHits[imodule]); // Number of hits that got in the module
      module.SetSignal(FHCalSumEnergy[imodule]); // Total energy from hits in the module
      //if (iEv == 0) std::cout << "Module iD: " << imodule << ": Energy loss recorded in a tree = " << module.GetSignal() << std::endl; 
    }

    outTree->Fill();
  } // End of the event loop
  
  outFile->cd();
  outTree->Print();
  outTree->Write();
  outFile->Close();

  //delete fhcal_modules;

  timer.Stop();
  timer.Print();

  return 0;
}
