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
    mHistogrammer.save("Histos_reco.root");
  };

  ic.services().get<CallbackService>().set<CallbackService::Id::Stop>(stop);
}

// ===== FILLING HISTOS =====
void MinvTask::fillHistos(gsl::span<const ExtendedTrack> tracks)
{

  // loop for 1D histograms
  for (const auto& t : tracks) {
    auto lv = getLorentzVector(t);
    mHistogrammer.fillSingleParticleHistos(lv);
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
        if (lv12.Rapidity() > -4 && lv12.Rapidity() < -2.5) {
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
  int compt = 0;
  int comptbis = 0;

  std::array<int, 10> ClustperCh;
  std::vector<int> Test;

  ofstream file;
  file.open("IDENTITY.txt");

  for (auto j = 0; j < tracks.size(); j++) {
    ClustperCh.fill(tracks[j].getNClusters());
    for (auto h = 0; h > tracks[j].getNClusters(); h++) {
      Test.push_back(tracks[j].getFirstClusterIdx());
      Test.push_back(tracks[j].getLastClusterIdx());
      // file << "TEST ID CLUSTERS : " << Test[h] << "\n";
    }
  }

  std::vector<Cluster> clust;
  std::vector<int> Clust_ID;
  std::vector<int> DE_Id;
  std::vector<int> Ch_ID;

  // std::array<int, 10> ClustperTrack;
  std::array<int, 10> TestClust = {1, 1, 1, 1, 1, 1, 0, 1, 0, 1};
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

  std::array<int, 10> ClustperTrack = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  std::vector<std::array<int, 10>> Total_Clusters;

  int compt_clust = 0;
  for (auto i = 0; i < rofs.size(); i++) {
    auto etracks = getExtendedTracks(rofs[i], tracks, clusters);
    for (auto ext = 0; ext < etracks.size(); ext++) {
      compt +=1;
      auto etracksbis = etracks[ext];
      clust = etracksbis.getClusters();
      // file << "CLUST SIZE : " << clust.size() << "\n";
      for (auto clst = 0; clst < clust.size(); clst++) {
        ClustperTrack = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        compt_clust += 1;
        Clust_ID.emplace_back(clust[clst].getClusterIndex());
        DE_Id.emplace_back(clust[clst].getDEId());
        Ch_ID.emplace_back(clust[clst].getChamberId());
        file << "=================================" << "\n";
        file << "Ch " << Ch_ID[clst] << "\n";
        file << clust[clst].getIdAsString() << "\n";

        switch (Ch_ID[clst]) {
          case 0:
            ClustperTrack[0] += 1; //ClustperCh
            break;
          case 1:
            ClustperTrack[1] += 1;
            break;
          case 2:
            ClustperTrack[2] += 1;
            break;
          case 3:
            ClustperTrack[3] += 1;
            break;
          case 4:
            ClustperTrack[4] += 1;
            break;
          case 5:
            ClustperTrack[5] += 1;
            break;
          case 6:
            ClustperTrack[6] += 1;
            break;
          case 7:
            ClustperTrack[7] += 1;
            break;
          case 8:
            ClustperTrack[8] += 1;
            break;
          case 9:
            ClustperTrack[9] += 1;
            break;
        }
      }
      Total_Clusters.emplace_back(ClustperCh);
      file << "ISTRACKABLE : " << o2::mch::isTrackable(Total_Clusters[ext], {true, true, true, true, true}, false) << "\n";
    }
    //file << "TOTAL CLUSTERS : " << Total_Clusters.size() << "\n"; //verification : size vector = nb of tracks
    //file << "TOTAL TRACKS : " << compt << "\n";
    dump(etracks);
    fillHistos(etracks);
  }
}
} // namespace o2::mch::eval
