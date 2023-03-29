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

using namespace o2::framework;

namespace o2::mch::eval
{

MinvTask::MinvTask(std::shared_ptr<o2::base::GRPGeomRequest> req) : mCcdbRequest(req) {}

void MinvTask::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (mCcdbRequest && base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    if (matcher == framework::ConcreteDataMatcher("GLO", "GRPMAGFIELD", 0)) {
      TrackExtrap::setField();
    }
  }
}

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
    mOutputRootFile->Close();
  };
  ic.services().get<CallbackService>().set<CallbackService::Id::Stop>(stop);
}

void MinvTask::createHistos()
{
        mHistos.emplace_back(new TH1F("y","rapidity",20,-5.,5.));
}

void MinvTask::fillHistos(gsl::span<const ExtendedTrack> tracks) 
{
        TH1* rap = mHistos[0];

        for (const auto& t: tracks) {
                rap->Fill(t.P().Rapidity());
        }
}

std::vector<ExtendedTrack> MinvTask::convert(gsl::span<const TrackMCH> mchTracks,
                                              gsl::span<const Cluster> clusters) const
{
  std::vector<ExtendedTrack> tracks;
  constexpr double vx{0.0};
  constexpr double vy{0.0};
  constexpr double vz{0.0};
  for (const auto& mchTrack : mchTracks) {
    tracks.emplace_back(mchTrack, clusters, vx, vy, vz);
  }
  return tracks;
}

void MinvTask::dump(gsl::span<const ExtendedTrack> tracks) const
{
    LOGP(warning, "# of tracks = {}", tracks.size());
    for (const auto& t : tracks) {
      LOGP(warning, "Track {}", t.asString());
    }
}

std::vector<ExtendedTrack> MinvTask::getExtendedTracks(const ROFRecord& rof,
                                                        gsl::span<const TrackMCH> tfTracks,
                                                        gsl::span<const Cluster> tfClusters) const
{
  const auto mchTracks = tfTracks.subspan(rof.getFirstIdx(), rof.getNEntries());
  return convert(mchTracks, tfClusters);
}

void MinvTask::run(ProcessingContext& pc)
{
  if (mCcdbRequest) {
    o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  }

  auto rofs = pc.inputs().get<gsl::span<ROFRecord>>("rofs");
  auto tracks = pc.inputs().get<gsl::span<TrackMCH>>("tracks");
  auto clusters = pc.inputs().get<gsl::span<Cluster>>("clusters");

  for (auto i = 0; i < rofs.size(); i++) {
    auto etracks = getExtendedTracks(rofs[i], tracks, clusters);
    dump(etracks);
    fillHistos(etracks);
  }
}
} // namespace o2::mch::eval
