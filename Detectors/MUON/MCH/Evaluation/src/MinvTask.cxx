// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MinvTask.h"
#include "Framework/CallbackService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/InputRecord.h"
#include "Framework/ProcessingContext.h"
#include "MCHEvaluation/ExtendedTrack.h"
#include "MCHTracking/TrackExtrap.h"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <iostream>
#include <Math/GenVector/VectorUtil.h>

using namespace o2::framework;
using namespace std;

namespace o2::mch::eval
{

// ===== CONSTRUCTOR =====
MinvTask::MinvTask(std::shared_ptr<o2::base::GRPGeomRequest> req) : mCcdbRequest(req) {}

// ===== INITIALIZATION OF MAGNETIC FIELD =====
void MinvTask::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj) 
{
  if (mCcdbRequest && base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    if (matcher == framework::ConcreteDataMatcher("GLO", "GRPMAGFIELD", 0)) {
      TrackExtrap::setField();
    }
  }
}

// ===== OUTPUT =====
void MinvTask::init(InitContext& ic)
{
  auto outputFileName = ic.options().get<std::string>("outfile"); 
  o2::base::GRPGeomHelper::instance().setRequest(mCcdbRequest);

  mOutputRootFile = std::make_unique<TFile>(outputFileName.c_str(), "RECREATE"); 

  createHistos();                                                                

  auto stop = [this]() {                                                         
    mOutputRootFile->cd();
    for (auto h : mHistos) {
      h->Write();
    }
    for (auto h2 : mHistos_2) {
      h2->Write();
    }
    mOutputRootFile->Close();
  };
  ic.services().get<CallbackService>().set<CallbackService::Id::Stop>(stop);
}

// ===== INITIALIZATION HISTOS =====
void MinvTask::createHistos()
{
  mHistos.emplace_back(new TH1F("eta", "Pseudo-Rapidity", 200, -5., -2.)); 
  mHistos.emplace_back(new TH1F("pT", "p_{T}", 128, 0., 15.));
  mHistos.emplace_back(new TH1F("M", "Mass", 40, -5., 5.));
  mHistos.emplace_back(new TH1F("Tracks", "Nb of tracks", 10, -1., 9.));
  mHistos.emplace_back(new TH1F("minv1", "Invariant Mass 1", 500, 0., 5.));
  mHistos.emplace_back(new TH1F("minv2", "InvariantMass 2", 500, 0., 5.));
  mHistos.emplace_back(new TH1F("minv3", "Test Track size", 8, -1., 8.));

  Double_t binRange_pT[] = {
    0.,
    1.,
    2.,
    3.,
    4.,
    5.,
    6.,
    8.,
    12.,
    16.,
  };
  mHistos_2.emplace_back(new TH2F("Minv pT", "M_{inv} depending of p_{T}", 9, binRange_pT, 200, 0., 5.));
  mHistos_2.emplace_back(new TH2F("Minv rap", "M_{inv} depending of rapidity", 8, -4.5, -2., 200, 0., 5.));
}

// ===== FILLING HISTOS =====
void MinvTask::fillHistos(gsl::span<const ExtendedTrack> tracks)
{
  TH1* rap = mHistos[0]; 
  TH1* pt = mHistos[1];
  TH1* mass = mHistos[2];
  TH1* tr = mHistos[3];
  // TH1* minv1 = mHistos[4];
  TH1* minv = mHistos[5];
  // TH1* minv3 = mHistos[6];

  TH2* M_pT = mHistos_2[0];
  TH2* M_y = mHistos_2[1];

  for (const auto& t : tracks) { 
    // rap->Fill(t.P().Rapidity());
    rap->Fill(t.P().Eta());
    rap->GetXaxis()->SetTitle("eta");
    rap->GetYaxis()->SetTitle("#");
    pt->Fill(t.P().Pt());
    pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    pt->GetYaxis()->SetTitle("#");
    mass->Fill(t.P().M());
    mass->GetXaxis()->SetTitle("mass (GeV/c^{2})");
    mass->GetYaxis()->SetTitle("#");
    tr->Fill(tracks.size());
    tr->GetYaxis()->SetTitle("#");
  }

  for (auto trk1 = 0; trk1 < tracks.size(); trk1++) { 
    // cout << "BOUCLE J" << endl;
    for (auto trk2 = trk1 + 1; trk2 < tracks.size(); trk2++) { 
      // cout << "BOUCLE J, H : " << j << " et " << h << endl;
      if (tracks[trk1].P().Eta() > -4. && tracks[trk1].P().Eta() < -2.5 && tracks[trk2].P().Eta() > -4. && tracks[trk2].P().Eta() < -2.5) { 
        //Extended tracks
        const auto& test_Inv = ROOT::Math::VectorUtil::InvariantMass(tracks[trk1].P(), tracks[trk2].P());
        const auto& test_P = tracks[trk1].P() + tracks[trk2].P();
        //Lorentz vectors
        ROOT::Math::PxPyPzMVector Lorentz_Reco1(tracks[trk1].P().Px(), tracks[trk1].P().Py(), tracks[trk1].P().Pz(), tracks[trk1].P().M());
        ROOT::Math::PxPyPzMVector Lorentz_Reco2(tracks[trk2].P().Px(), tracks[trk2].P().Py(), tracks[trk2].P().Pz(), tracks[trk2].P().M());
        auto minv_Lorentz = ROOT::Math::VectorUtil::InvariantMass(Lorentz_Reco1, Lorentz_Reco2);
        auto Lorentz_Reco_Fusion = Lorentz_Reco1 + Lorentz_Reco2;
        if (Lorentz_Reco_Fusion.Rapidity() > -4 && Lorentz_Reco_Fusion.Rapidity() < -2.5) {
          minv->Fill(minv_Lorentz, 1);
          M_pT->Fill(Lorentz_Reco_Fusion.Pt(), minv_Lorentz, 1); 
          M_y->Fill(Lorentz_Reco_Fusion.Rapidity(), minv_Lorentz, 1);
        } else {
          continue;
        }
      } else {
        continue;
      }
      minv->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
      M_pT->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      M_pT->GetYaxis()->SetTitle("invariant mass (GeV/c^{2})");
      M_y->GetXaxis()->SetTitle("y");
      M_y->GetYaxis()->SetTitle("invariant mass (GeV/c^{2})");

      //Test
      /*const auto& Mass1 = tracks[j].P().M();
      const auto& Mass2 = tracks[h].P().M();
      const auto& E1 = tracks[j].P().E();
      const auto& E2 = tracks[h].P().E();
      const auto& Ptot = tracks[j].P().Px() * tracks[h].P().Px() + tracks[j].P().Py() * tracks[h].P().Py() + tracks[j].P().Pz() * tracks[h].P().Pz();
      const auto& test_Form = sqrt(Mass1 * Mass1 + Mass2 * Mass2 + 2 * (E1 * E2 - Ptot));  // calcul 2 : formule litterale
      test3->Fill(test_Form);
      test3->GetXaxis()->SetTitle("invariant mass (GeV/c^2)");*/
    }
  }
}

// ===== CREATION OF A EXTENDED TRACK FROM A TRACKMCH =====
std::vector<ExtendedTrack> MinvTask::convert(gsl::span<const TrackMCH> mchTracks,
                                             gsl::span<const Cluster> clusters) const 
{
  std::vector<ExtendedTrack> tracks;
  constexpr double vx{0.0}; 
  constexpr double vy{0.0};
  constexpr double vz{0.0};
  for (const auto& mchTrack : mchTracks) {
    tracks.emplace_back(ExtendedTrack(mchTrack, clusters, vx, vy, vz));
  }
  return tracks;
}

// ===== WARNING lEVEL IN TERMINAL =====
void MinvTask::dump(gsl::span<const ExtendedTrack> tracks) const 
{
  LOGP(warning, "# of tracks = {}", tracks.size());
  for (const auto& t : tracks) {
    LOGP(warning, "Track {}", t.asString());
  }
}

// ===== EXTENDED TRACKS CALL =====
std::vector<ExtendedTrack> MinvTask::getExtendedTracks(const ROFRecord& rof,
                                                       gsl::span<const TrackMCH> tfTracks,
                                                       gsl::span<const Cluster> tfClusters) const
{
  const auto mchTracks = tfTracks.subspan(rof.getFirstIdx(), rof.getNEntries()); 
  return convert(mchTracks, tfClusters);
}

// ===== FINAL =====
void MinvTask::run(ProcessingContext& pc)
{
  if (mCcdbRequest) {
    o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  }
  auto rofs = pc.inputs().get<gsl::span<ROFRecord>>("rofs"); 
  auto tracks = pc.inputs().get<gsl::span<TrackMCH>>("tracks");
  auto clusters = pc.inputs().get<gsl::span<Cluster>>("clusters");
  int compt = 0;
  int comptbis = 0;

  ofstream file_rof, file_tracks;
  file_rof.open("rof.txt");
  file_rof << rofs.size() << endl;
  file_tracks.open("tracks.txt");

  for (auto i = 0; i < rofs.size(); i++) {                     
    auto etracks = getExtendedTracks(rofs[i], tracks, clusters);
    // cout << "DANS LA BOUCLE ROF -> ROF.SIZE() : " << rofs.size() << endl;
    // cout << "DANS LA BOUCLE ROF -> CLUSTER.SIZE() : " << clusters.size() << endl;
    // cout << "DANS LA BOUCLE ROF -> TRACK.SIZE() : " << tracks.size() << endl;
    // cout << "DANS LA BOUCLE ROF -> ETRACK.SIZE() : " << etracks.size() << endl;
    file_tracks << etracks.size() << "\n";
    if (etracks.size() != 0) {
      compt += 1;
    } else if (etracks.size() == 0) {
      comptbis += 1;
    }
    dump(etracks);
    fillHistos(etracks);
  }
  file_tracks.close();
  // cout << "Nb Traks : " << compt << endl;
  // cout << "Nb Traks : " << comptbis << endl;
}
} // namespace o2::mch::eval
