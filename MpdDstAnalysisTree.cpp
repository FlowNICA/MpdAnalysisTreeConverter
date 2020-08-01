// Standard c++ headers
#include <iostream>
#include <string>
#include <utility>
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
#include <TMath.h>

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
#include "AnalysisTree/DataHeader.hpp"
#include "AnalysisTree/EventHeader.hpp"
#include "AnalysisTree/Detector.hpp"
#include "AnalysisTree/Matching.hpp"

float get_beamP(float sqrtSnn, float m_target = 0.938, float m_beam = 0.938);
float GetFHCalPhi(int iModule);

int main(int argc, char **argv)
{
  TString iFileName, oFileName;
  std::string system="";
  float sqrtSnn = -1.;
  bool isDataHeaderDefined = false;

  if (argc < 5)
  {
    std::cerr << "./MpdDst2AnalysisTree -i inputfile -o outputfile [OPTIONAL: --sqrtSnn sqrtSnn --system colliding_system_name]" << std::endl;
    return 1;
  }
  for (int i = 1; i < argc; i++)
  {
    if (std::string(argv[i]) != "-i" &&
      std::string(argv[i]) != "-o" &&
      std::string(argv[i]) != "--sqrtSnn" &&
      std::string(argv[i]) != "--system")
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
      if (std::string(argv[i]) == "--system" && i != argc - 1)
      {
        system = std::string(argv[++i]);
        isDataHeaderDefined = true;
        continue;
       }
      if (std::string(argv[i]) == "--system" && i == argc - 1)
      {
        std::cerr << "\n[ERROR]: Output file name was not specified " << std::endl;
        return 1;
      }
      if (std::string(argv[i]) == "--sqrtSnn" && i != argc - 1)
      {
        sqrtSnn = std::stof(std::string(argv[++i]));
        isDataHeaderDefined = true;
        continue;
       }
      if (std::string(argv[i]) == "--sqrtSnn" && i == argc - 1)
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
  TTree *outTree = new TTree("aTree","AnalysisTree Dst at MPD");

  // Set up AnalysisTree data header
  AnalysisTree::DataHeader *dataHeader = new AnalysisTree::DataHeader;
  if (system != "") dataHeader->SetSystem(system);
  if (sqrtSnn > 0) dataHeader->SetBeamMomentum(get_beamP(sqrtSnn));

  // Set up AnalysisTree configureation
  AnalysisTree::Configuration *out_config = new AnalysisTree::Configuration;

  // Set up AnalysisTree configuration
  std::string str_reco_event_branch = "RecoEvent.";
  std::string str_mc_event_branch = "McEvent.";
  std::string str_tpc_tracks_branch = "TpcTracks.";
  std::string str_fhcal_branch = "FHCalModules.";
  std::string str_mc_tracks_branch = "McTracks.";
  std::string str_tpc2mc_tracks_branch = "TpcTracks2McTracks";

  AnalysisTree::BranchConfig reco_event_branch(str_reco_event_branch.c_str(), AnalysisTree::DetType::kEventHeader); 
  AnalysisTree::BranchConfig mc_event_branch(str_mc_event_branch.c_str(), AnalysisTree::DetType::kEventHeader); 
  AnalysisTree::BranchConfig fhcal_branch(str_fhcal_branch.c_str(), AnalysisTree::DetType::kModule); 
  AnalysisTree::BranchConfig tpc_tracks_branch(str_tpc_tracks_branch.c_str(), AnalysisTree::DetType::kTrack); 
  AnalysisTree::BranchConfig mc_tracks_branch(str_mc_tracks_branch.c_str(), AnalysisTree::DetType::kParticle);

  mc_event_branch.AddField<float>("B");
  mc_event_branch.AddField<float>("PhiRp");

  tpc_tracks_branch.AddField<int>("nhits");
  tpc_tracks_branch.AddField<int>("nhits_poss");
  tpc_tracks_branch.AddField<int>("charge");
  tpc_tracks_branch.AddFields<float>({"dca_x","dca_y","dca_z"});
  tpc_tracks_branch.AddField<float>("chi2");
  tpc_tracks_branch.AddField<float>("tof_mass2");
  tpc_tracks_branch.AddField<float>("tof_flag");
  tpc_tracks_branch.AddField<float>("dedx");
  tpc_tracks_branch.AddField<float>("pid_prob_pion");
  tpc_tracks_branch.AddField<float>("pid_prob_kaon");
  tpc_tracks_branch.AddField<float>("pid_prob_proton");

  fhcal_branch.AddField<float>("phi");

  mc_tracks_branch.AddField<int>("mother_id");

  // mc_event's additional field ids
  const int iB = mc_event_branch.GetFieldId("B");
  const int iPhiRp = mc_event_branch.GetFieldId("PhiRp");

  // fhcal_modules' additional field ids
  const int ifhcalphi = fhcal_branch.GetFieldId("phi");

  // tpc_tracks' additional field ids
  const int inhits = tpc_tracks_branch.GetFieldId("nhits");
  const int inhits_poss = tpc_tracks_branch.GetFieldId("nhits_poss");
  const int icharge = tpc_tracks_branch.GetFieldId("charge");
  const int idcax = tpc_tracks_branch.GetFieldId("dca_x");
  const int idcay = tpc_tracks_branch.GetFieldId("dca_y");
  const int idcaz = tpc_tracks_branch.GetFieldId("dca_z");
  const int ichi2 = tpc_tracks_branch.GetFieldId("chi2");
  const int itof_mass2 = tpc_tracks_branch.GetFieldId("tof_mass2");
  const int itof_flag = tpc_tracks_branch.GetFieldId("tof_flag");
  const int idedx = tpc_tracks_branch.GetFieldId("dedx");
  const int ipid_prob_pion = tpc_tracks_branch.GetFieldId("pid_prob_pion");
  const int ipid_prob_kaon = tpc_tracks_branch.GetFieldId("pid_prob_kaon");
  const int ipid_prob_proton = tpc_tracks_branch.GetFieldId("pid_prob_proton");

  // mc_tracks' additional field ids
  const int imother_id = mc_tracks_branch.GetFieldId("mother_id");

  // Initialize AnalysisTree Dst components
  out_config->AddBranchConfig(std::move(reco_event_branch));
  AnalysisTree::EventHeader *reco_event = new AnalysisTree::EventHeader( out_config->GetLastId() );
  out_config->AddBranchConfig(std::move(mc_event_branch));
  AnalysisTree::EventHeader *mc_event = new AnalysisTree::EventHeader( out_config->GetLastId() );
  out_config->AddBranchConfig(std::move(tpc_tracks_branch));
  AnalysisTree::TrackDetector *tpc_tracks = new AnalysisTree::TrackDetector( out_config->GetLastId() ); 
  out_config->AddBranchConfig(std::move(fhcal_branch));
  AnalysisTree::ModuleDetector *fhcal_modules = new AnalysisTree::ModuleDetector( out_config->GetLastId() );
  out_config->AddBranchConfig(std::move(mc_tracks_branch));
  AnalysisTree::Particles *mc_tracks = new AnalysisTree::Particles( out_config->GetLastId() ); 
  AnalysisTree::Matching *tpc2mc_tracks = new AnalysisTree::Matching(out_config->GetBranchConfig(str_tpc_tracks_branch).GetId(), out_config->GetBranchConfig(str_mc_tracks_branch).GetId());
  out_config->AddMatch(str_tpc_tracks_branch, str_mc_tracks_branch, str_tpc2mc_tracks_branch);

  // Create branches in the output tree
  outTree->Branch(str_reco_event_branch.c_str(), "AnalysisTree::EventHeader", &reco_event, 32000, 99);
  outTree->Branch(str_mc_event_branch.c_str(), "AnalysisTree::EventHeader", &mc_event, 32000, 99);
  outTree->Branch(str_tpc_tracks_branch.c_str(), "AnalysisTree::TrackDetector", &tpc_tracks, 256000, 99);
  outTree->Branch(str_fhcal_branch.c_str(), "AnalysisTree::ModuleDetector",  &fhcal_modules, 128000, 99);
  outTree->Branch(str_mc_tracks_branch.c_str(), "AnalysisTree::Particles",  &mc_tracks, 256000, 99);
  outTree->Branch(str_tpc2mc_tracks_branch.c_str(), "AnalysisTree::Matching", &tpc2mc_tracks, 32000, 99);

  // Printout basic configuration info
  out_config->Print();

  // Starting event loop
  TVector3 primaryVertex;
  std::set <Int_t> UsedMCTracks;        // using to remap mc-reco track matching
  std::map <Int_t,Int_t> InitMcNewMcId; // map[old-mc-id] = new-mc-id
  Float_t FHCalSumEnergy[Num_Of_Modules];
  Int_t   FHCalNumOfHits[Num_Of_Modules];
  Int_t n_entries = dstTree->GetEntriesFast();

  bool isGoodPID;
  int charge;

  for (int iEv = 0; iEv < n_entries; iEv++)
  {
    std::cout << "Event [" << iEv << "/" << n_entries << "]" << std::endl;
    dstTree->GetEntry(iEv);

    UsedMCTracks.clear();
    InitMcNewMcId.clear();
    tpc2mc_tracks->Clear();
    for (int i=0; i<Num_Of_Modules; i++)
    {
      FHCalSumEnergy[i] = 0.;
      FHCalNumOfHits[i] = 0;
    }

    tpc_tracks->ClearChannels();
    fhcal_modules->ClearChannels();
    mc_tracks->ClearChannels();

    // Reading Reco Event
    MpdVertex *vertex = (MpdVertex *)vertexes->First();
    vertex->Position(primaryVertex);
    reco_event->SetVertexPosition3(primaryVertex);

    // Read MC Event
    mc_event->Init(out_config->GetBranchConfig(mc_event->GetId()));
    mc_event->SetVertexX(MCHeader->GetX());
    mc_event->SetVertexY(MCHeader->GetY());
    mc_event->SetVertexZ(MCHeader->GetZ());
    mc_event->SetField(float(MCHeader->GetB()), iB);
    mc_event->SetField(float(MCHeader->GetRotZ()), iPhiRp);

    // Read energy in FHCal modules 
    fhcal_modules->Reserve(Num_Of_Modules);
    for (int imodule=0; imodule<Num_Of_Modules; imodule++)
    {
      auto *module = fhcal_modules->AddChannel();
      module->Init(out_config->GetBranchConfig(fhcal_modules->GetId()));
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
      module.SetField(float(GetFHCalPhi(imodule)), ifhcalphi);
    }

    // Reading Reco Tracks
    MpdGlobalTracks = (TClonesArray*)MPDEvent->GetGlobalTracks();
    Int_t Num_of_tpc_tracks = MpdGlobalTracks->GetEntriesFast();
    tpc_tracks->Reserve(Num_of_tpc_tracks);

    for (int itrack=0; itrack<Num_of_tpc_tracks; itrack++)
    {
      MpdTrack* mpdtrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(itrack);
      UsedMCTracks.insert(mpdtrack->GetID());
      auto *track = tpc_tracks->AddChannel();
      track->Init(out_config->GetBranchConfig(tpc_tracks->GetId()));

      track->SetMomentum(mpdtrack->GetPx(),mpdtrack->GetPy(),mpdtrack->GetPz());
      track->SetField(int(mpdtrack->GetNofHits()), inhits);
      track->SetField(int(mpdtrack->GetNofHitsPossTpc()), inhits_poss);
      track->SetField(float(mpdtrack->GetTofMass2()), itof_mass2);
      track->SetField(float(mpdtrack->GetTofFlag()), itof_flag);
      track->SetField(float(mpdtrack->GetdEdXTPC()), idedx);
      track->SetField(float(mpdtrack->GetChi2()), ichi2);
      track->SetField(float(mpdtrack->GetDCAX()), idcax);
      track->SetField(float(mpdtrack->GetDCAY()), idcay);
      track->SetField(float(mpdtrack->GetDCAZ()), idcaz);
      charge = (mpdtrack->GetPt() < 0) ? 1 : -1;
      track->SetField(int(charge), icharge);

      // PID
      if (mpdtrack->GetTofFlag() == 2 || mpdtrack->GetTofFlag() == 6)
      {
        // TOF + TPC
        isGoodPID = pid->FillProbs(TMath::Abs(mpdtrack->GetPt()) * TMath::CosH(mpdtrack->GetEta()),
                                   mpdtrack->GetdEdXTPC(), mpdtrack->GetTofMass2(), charge);
      }
      else
      {
        // TPC only
        isGoodPID = pid->FillProbs(TMath::Abs(mpdtrack->GetPt()) * TMath::CosH(mpdtrack->GetEta()),
                                   mpdtrack->GetdEdXTPC(), charge);
      }

      if (isGoodPID)
      {
        track->SetField(float(pid->GetProbPi()), ipid_prob_pion);
        track->SetField(float(pid->GetProbKa()), ipid_prob_kaon);
        track->SetField(float(pid->GetProbPr()), ipid_prob_proton);
      }
      else
      {
        track->SetField(float(-999.), ipid_prob_pion);
        track->SetField(float(-999.), ipid_prob_kaon);
        track->SetField(float(-999.), ipid_prob_proton);
      }
    } // End of the tpc track loop

    // Read Mc tracks
    Int_t Num_of_mc_tracks = MCTracks->GetEntriesFast();
    mc_tracks->Reserve(Num_of_mc_tracks);

    for (int imctrack=0; imctrack<Num_of_mc_tracks; imctrack++)
    {
#ifdef _MCSTACK_
      FairMCTrack *mctrack = (FairMCTrack*) MCTracks->UncheckedAt(imctrack);
#endif
#ifdef _MPDMCSTACK_
      MpdMCTrack *mctrack = (MpdMCTrack*) MCTracks->UncheckedAt(imctrack);
#endif
      bool isUsed = (UsedMCTracks.count(imctrack));

      // If motherId != 1 and mc track doesn't relate to any reco track - skip
      if (mctrack->GetMotherId() != -1 && !isUsed) continue;

      auto *track = mc_tracks->AddChannel();
      track->Init(out_config->GetBranchConfig(mc_tracks->GetId()));

      // Collect new Mc Ids
      InitMcNewMcId[imctrack] = track->GetId();

      track->SetMomentum(mctrack->GetPx(), mctrack->GetPy(), mctrack->GetPz());
      track->SetPid(int(mctrack->GetPdgCode()));
      track->SetMass(float(mctrack->GetMass()));
      track->SetField(int(mctrack->GetMotherId()), imother_id);
    } // End of the mc track loop

    // reco-mc tracks matching
    for (int itrack=0; itrack<Num_of_tpc_tracks; itrack++)
    {
      MpdTrack* mpdtrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(itrack);
      const int tpc_track_id = itrack;
      const int mc_track_id  = InitMcNewMcId[mpdtrack->GetID()];
      tpc2mc_tracks->AddMatch(tpc_track_id, mc_track_id);
    }

    outTree->Fill();
  } // End of the event loop
  
  outFile->cd();
  outTree->Print();
  outTree->Write();
  if (isDataHeaderDefined) dataHeader->Write("DataHeader");
  out_config->Write("Configuration");
  outFile->Close();

  timer.Stop();
  timer.Print();

  return 0;
}


float get_beamP(float sqrtSnn, float m_target, float m_beam)
{
  return sqrt( pow((pow(sqrtSnn,2) - pow(m_target,2) - pow(m_beam,2))/(2*m_target), 2) - pow(m_beam,2) );
}

float GetFHCalPhi(int iModule)
{
  const int Nmodules = 45;
  int xAxisSwitch = (iModule < Nmodules) ? 1 : -1;
  int module = (iModule < Nmodules) ? iModule : iModule - Nmodules;
  float x, y, phi = -999.;
  if (module >= 0 && module <= 4)
  {
    y = 45.;
    x = (module - 2) * 15.;
    phi = TMath::ATan2(y, x * xAxisSwitch);
  }
  else if ((module >= 5) && (module <= 39))
  {
    y = (3 - (module + 2) / 7) * 15.;
    x = (3 - (module + 2) % 7) * 15.;
    phi = TMath::ATan2(y, x * xAxisSwitch);
  }
  else if ((module >= 40) && (module <= 44))
  {
    y = -45.;
    x = (module - 42) * 15.;
    phi = TMath::ATan2(y, x * xAxisSwitch);
  }
  return phi;
}
