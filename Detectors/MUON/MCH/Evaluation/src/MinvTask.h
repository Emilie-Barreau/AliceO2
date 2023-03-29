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

#ifndef O2_MCH_EVALUATION_MINV_TASK_H__
#define O2_MCH_EVALUATION_MINV_TASK_H__

#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/InitContext.h"
#include "MCHEvaluation/ExtendedTrack.h"
#include <TFile.h>
#include <gsl/span>
#include <memory>
#include <vector>

class TH1;

namespace o2::mch::eval
{
class MinvTask
{
 public:
  MinvTask(std::shared_ptr<o2::base::GRPGeomRequest> req);

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj);

  void init(o2::framework::InitContext& ic);

  std::vector<ExtendedTrack> convert(gsl::span<const TrackMCH> mchTracks,
                                   gsl::span<const Cluster> clusters) const;

  void dump(gsl::span<const ExtendedTrack> tracks) const;

  std::vector<ExtendedTrack> getExtendedTracks(const ROFRecord& rof,
                                             gsl::span<const TrackMCH> tfTracks,
                                             gsl::span<const Cluster> tfClusters) const;

  void run(o2::framework::ProcessingContext& pc);

 private:
   void createHistos();
   void fillHistos(gsl::span<const ExtendedTrack> tracks);

 private:
  std::shared_ptr<o2::base::GRPGeomRequest> mCcdbRequest;
  std::unique_ptr<TFile> mOutputRootFile;
  std::vector<TH1*> mHistos;
};
} // namespace o2::mch::eval

#endif
