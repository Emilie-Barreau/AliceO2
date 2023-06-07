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

#include "MCHBase/Trackable.h"
#include "Histogrammer.h"
#include "MinvTask.h"
#include "Framework/CallbackService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/InputRecord.h"
#include "Framework/ProcessingContext.h"
#include "MCHEvaluation/ExtendedTrack.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHGlobalMapping/DsIndex.h"
#include "MCHGlobalMapping/ChannelCode.h"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <Math/GenVector/VectorUtil.h>

using namespace o2::framework;
using namespace std;

namespace o2::mch::eval
{

MinvTask::MinvTask(std::shared_ptr<o2::base::GRPGeomRequest> req) : mCcdbRequest(req) {}

// initialize the magnetic field and check if it is on or not
void MinvTask::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (mCcdbRequest && base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    if (matcher == framework::ConcreteDataMatcher("GLO", "GRPMAGFIELD", 0)) {
      TrackExtrap::setField();
    }
  }
}

// function oustide classes, gives a Lorentzvector from an Extended Track
ROOT::Math::PxPyPzMVector getLorentzVector(o2::mch::eval::ExtendedTrack t)
{
  return ROOT::Math::PxPyPzMVector(t.P().Px(), t.P().Py(), t.P().Pz(), t.P().M());
}

// save and create the output
void MinvTask::init(InitContext& ic)
{
  auto outputFileName = ic.options().get<std::string>("outfile");
  o2::base::GRPGeomHelper::instance().setRequest(mCcdbRequest);

  mOutputRootFile = std::make_unique<TFile>(outputFileName.c_str(), "RECREATE");

  auto stop = [this]() {
    mOutputRootFile->cd();
    for (auto h : mHistos) {
      h->Write();
    }
    for (auto h2 : mHistos_2) {
      h2->Write();
    }
    mOutputRootFile->Close();
    // mHistogrammer.save("Histos_reco.root");
    mHistogrammer.save("Histos_reco_cut.root");
  };

  ic.services().get<CallbackService>().set<CallbackService::Id::Stop>(stop);
}

// ===== FILLING HISTOS =====
void MinvTask::fillHistos(gsl::span<const ExtendedTrack> tracks)
{

  // loop for 1D histograms
  for (const auto& t : tracks) {
    auto lv = getLorentzVector(t);
    // mHistogrammer.fillSingleParticleHistos(lv);
  }

  // loop for 2D histograms
  for (auto trk1 = 0; trk1 < tracks.size(); trk1++) {
    auto t1 = tracks[trk1];
    for (auto trk2 = trk1 + 1; trk2 < tracks.size(); trk2++) {
      auto t2 = tracks[trk2];
      auto lv1 = getLorentzVector(t1);
      auto lv2 = getLorentzVector(t2);
      auto lv12 = lv1 + lv2;
      if (tracks[trk1].P().Eta() > -4. && tracks[trk1].P().Eta() < -2.5 && tracks[trk2].P().Eta() > -4. && tracks[trk2].P().Eta() < -2.5) {
        // mHistogrammer.fillDoubleParticleHistos(lv1, lv2);
        if (lv12.Rapidity() > -4 && lv12.Rapidity() < -2.5) {
          mHistogrammer.fillSingleParticleHistos(lv12);
          mHistogrammer.fillDoubleParticleHistos(lv1, lv2);
        } else {
          continue;
        }
      } else {
        continue;
      }
    }
  }
}

// create an Extended Track from a TrackMCH + clusters
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

// warning level in terminal
void MinvTask::dump(gsl::span<const ExtendedTrack> tracks) const
{
  LOGP(warning, "# of tracks = {}", tracks.size());
  for (const auto& t : tracks) {
    LOGP(warning, "Track {}", t.asString());
  }
}

// calling Extended Tracks
std::vector<ExtendedTrack> MinvTask::getExtendedTracks(const ROFRecord& rof,
                                                       gsl::span<const TrackMCH> tfTracks,
                                                       gsl::span<const Cluster> tfClusters) const
{
  const auto mchTracks = tfTracks.subspan(rof.getFirstIdx(), rof.getNEntries());
  // std::cout << "CLUSTERS SIZE TF : " << tfClusters.size() << std::endl;  // tous les clusters
  // std::cout << "TRACKS SIZE TF : " << tfTracks.size() << std::endl;
  // std::cout << "MCHTRACKS SIZE TF : " << mchTracks.size() << std::endl;
  return convert(mchTracks, tfClusters);
}

// loop on ROF
void MinvTask::run(ProcessingContext& pc)
{
  if (mCcdbRequest) {
    o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  }
  auto rofs = pc.inputs().get<gsl::span<ROFRecord>>("rofs");
  auto tracks = pc.inputs().get<gsl::span<TrackMCH>>("tracks");
  auto clusters = pc.inputs().get<gsl::span<Cluster>>("clusters");

  ofstream file;
  file.open("IDENTITY.txt");

  //{1, 1, 1, 1, 1, 1, 2, 2, 2, 2} is trackable
  //{1, 1, 1, 1, 1, 1, 1, 1, 1, 1} is trackable
  //{1, 1, 1, 1, 1, 1, 1, 1, 1, 0} is trackable (idem Ch 7, 8, 9)
  //{1, 0, 1, 0, 1, 0, 1, 0, 1, 0} is trackable (idem inverse)
  //{1, 1, 1, 1, 1, 1, 1, 1, 0, 0} is not trackable (idem St 4)
  //{0, 0, 1, 1, 1, 1, 1, 1, 1, 1} is not trackable (idem St 2, 3 ...)
  //{0, 0, 1, 1, 1, 1, 0, 2, 0, 0} is not trackable
  //{0, 0, 1, 1, 1, 1, 2, 2, 0, 0} is not trackable
  //{1, 1, 1, 1, 1, 1, 0, 2, 0, 0} is not trackable (moreCandidates true and false)
  //{1, 1, 1, 1, 1, 1, 0, 1, 0, 1} is not trackable (moreCandidates = false)
  //{1, 1, 1, 1, 1, 1, 0, 1, 0, 1} is trackable (moreCandidates = true)
  // need at least one hit per station to be trackable
  // for moreCandidates = false : need one hit per chamber for nb 7, 8, 9, 10

  std::vector<Cluster> clust;
  int De_ID;
  int Ch;

  double compt_t1 = 0.;
  double compt_t2 = 0.;
  double compt_clust1 = 0.;
  double compt_clust2 = 0.;

  std::ifstream jsn("/Users/emiliebarreau/alice/TEST_50000evt/geom.json");
  auto transformcreat = o2::mch::geo::transformationFromJSON(jsn);

  std::array<int, 10> ClustperChamber = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // ClusterperChamber

  for (auto i = 0; i < rofs.size(); i++) {
    // auto sppan = clusters.subspan(0, 621500);
    auto etracks = getExtendedTracks(rofs[i], tracks, clusters);
    std::vector<o2::mch::eval::ExtendedTrack> etrackables;
    // std::vector<o2::mch::eval::ExtendedTrack> etracksbis;
    for (auto ext = 0; ext < etracks.size(); ext++) {
      compt_t1 += 1;
      auto etrack = etracks[ext];
      clust = etrack.getClusters();
      ClustperChamber = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      std::vector<Cluster> clst_survivors;
      // clust_rof.insert(clust_rof.end(), clust.begin(), clust.end());
      for (auto clst = 0; clst < clust.size(); clst++) {
        compt_clust1 += 1;
        // file << "================================="
        //  << "\n";
        // file << "Ch " << Ch_ID[clst] << "\n";
        // file << clust[clst].getIdAsString() << "\n";
        // file << clust[clst] << "\n";
        De_ID = clust[clst].getDEId();
        Ch = clust[clst].getChamberId();

        auto T3D = transformcreat(clust[clst].getDEId());
        math_utils::Point3D<double> coord_glob(clust[clst].getX(), clust[clst].getY(), clust[clst].getZ());
        math_utils::Point3D<double> coord_loc(0., 0., 0.);
        T3D.MasterToLocal(coord_glob, coord_loc);
        auto& s = o2::mch::mapping::segmentation(clust[clst].getDEId());
        // int numberofpadsperDE = s.nofPads();
        int PadIndexB = 0;
        int PadIndexNB = 0;
        bool test = s.findPadPairByPosition(coord_loc.X(), coord_loc.Y(), PadIndexB, PadIndexNB); // jusqu'ici ça va
        int compteur_DE = 0;
        if (PadIndexB == -1 || PadIndexNB == -1 || (PadIndexB == -1 && PadIndexNB == -1)) {
          mHistogrammer.DEtest(De_ID);
          compteur_DE += 1;
          file << "DE Id : " << De_ID << "\n";
        }
        ChannelCode c_B;
        ChannelCode c_NB;
        uint16_t SolarIndexB;
        uint16_t SolarIndexNB;
        if (test == true) {
          int DualSampaIdB = s.padDualSampaId(PadIndexB);
          int DualSampaIdNB = s.padDualSampaId(PadIndexNB);
          uint16_t DsIndexB = o2::mch::getDsIndex({clust[clst].getDEId(), DualSampaIdB});
          uint16_t DsIndexNB = o2::mch::getDsIndex({clust[clst].getDEId(), DualSampaIdNB});

          c_B = ChannelCode(De_ID, PadIndexB);
          c_NB = ChannelCode(De_ID, PadIndexNB);
          SolarIndexB = c_B.getSolarIndex();
          SolarIndexNB = c_NB.getSolarIndex();

          file << "Chamber : " << Ch << "\n";
          file << "DE : " << De_ID << "\n";
          file << "SolarIndex INDEX B ET NB : " << SolarIndexB << " & " << SolarIndexNB << "\n";
          file << "SolarIndex ID B ET NB : " << c_B.getSolarId() << " & " << c_NB.getSolarId() << "\n";

          // file << "SOLAR INDEX B ET NB : " << SolarIndexB << " & " << SolarIndexNB << "\n";
        } else {
          continue;
        }

        // file << "SOLAR ID B : " << c_B.getSolarId() << "\n";
        // file << "DE ID : " << De_ID << "\n";
        // file << "DS INDEX B : " << DsIndexB << " & DS INDEX NB : " << DsIndexNB << "\n";
        // file << "SOLAR INDEX B : " << SolarIndexB << " & SOLAR INDEX NB : " << SolarIndexNB << "\n";
        // file << "---------------------------------------------" << "\n";

        // condition to remove a SOLAR
        if (c_B.getSolarId() == 203 || c_B.getSolarId() == 204 || c_NB.getSolarId() == 91 || c_NB.getSolarId() == 92 ||
            c_B.getSolarId() == 259 || c_B.getSolarId() == 260 || c_NB.getSolarId() == 283 || c_NB.getSolarId() == 284) {
          ClustperChamber[Ch] += 0;
        } else {
          ClustperChamber[Ch] += 1;
          clst_survivors.emplace_back(clust[clst]);
          compt_clust2 += 1;
        }

        // condition to remove a DE
        /*if (De_ID == 101 || De_ID == 201) {
          ClustperChamber[Ch] += 0;
        } else {
          ClustperChamber[Ch] += 1;
          clst_survivors.emplace_back(clust[clst]);
          compt_clust2 += 1;
        }*/

        // condition to remove a Chamber
        /*if (Ch == 9) {
          ClustperChamber[Ch] += 0;
        } else {
          ClustperChamber[Ch] += 1;
          clst_survivors.emplace_back(clust[clst]);
          compt_clust2 += 1;
        }*/
      }
      auto etrackbis = ExtendedTrack(clst_survivors, 0., 0., 0.);
      bool trackable = o2::mch::isTrackable(ClustperChamber, {true, true, true, true, true}, false);
      if (trackable == true) {
        // etrackables.emplace_back(etrack); //methode 1
        etrackables.emplace_back(etrackbis);
        compt_t2 += 1;
      }
    }
    dump(etrackables);
    fillHistos(etrackables);
  }
  file << compt_t2 << " are trackables on " << compt_t1 << "\n";
  file << compt_clust2 << " clusters have survived over " << compt_clust1 << "\n";
  auto rapport = compt_t2 * 100 / compt_t1;
  double limit = 20.;
  if (100 - (compt_t2 * 100 / compt_t1) > limit) {
    file << Form("%.2f PERCENT OF TRACKS ARE KEPT", rapport) << "\n";
    file << Form("WARNING : OVER %.f PERCENT OF TRACKS ARE LOST -> CRITICAL EFFICIENCY EXPECTED", limit)
         << "\n";
    LOGP(warning, "KEPT TRACKS : {} %", rapport);
    LOGP(warning, "LOSSES OF TRACK OVER {} %", limit);
  } else {
    file << Form("%.2f PERCENT OF TRACKS ARE KEPT", rapport) << "\n";
    LOGP(warning, "KEPT TRACKS : {} %", rapport);
  }
}
} // namespace o2::mch::eval
