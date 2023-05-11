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
    mHistogrammer.save("Histos_bis.root");
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
         if (lv12.M() > 3.05 && lv12.M() < 3.15) {
            mHistogrammer.fillDoubleParticleHistos(lv1, lv2);
          } else {
            continue;
          }
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

  for (auto i = 0; i < rofs.size(); i++) {
    auto etracks = getExtendedTracks(rofs[i], tracks, clusters);
    if (etracks.size() != 0) {
      compt += 1;
    } else if (etracks.size() == 0) {
      comptbis += 1;
    }
    dump(etracks);
    fillHistos(etracks);
  }
}
} // namespace o2::mch::eval
