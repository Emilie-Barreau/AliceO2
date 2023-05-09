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

//==========FABRIQUE L'EXECUTABLE==========

#include "CommonUtils/ConfigurableParam.h"
#include "MinvTask.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/CallbacksPolicy.h"
#include "Framework/CompletionPolicyHelpers.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/ConfigContext.h"
#include "Framework/Logger.h"
#include "Framework/Task.h"
#include "Framework/Variant.h"
#include "Framework/WorkflowSpec.h"
#include <string>

using namespace o2::framework;
using namespace o2::mch;

void customize(std::vector<ConfigParamSpec>& workflowOptions)  //options sur le workflow
{
  std::vector<ConfigParamSpec> options{
    {"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings"}}};  //separation des colonnes ?
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"  //Main

WorkflowSpec defineDataProcessing(ConfigContext const& configcontext)
{
  WorkflowSpec specs;

  o2::conf::ConfigurableParam::updateFromString(configcontext.options().get<std::string>("configKeyValues"));

  Inputs inputs{};  //entrees utiles pour le code : les rof, les clusters et les traces
  inputs.emplace_back("rofs", "MCH", "TRACKROFS", 0, Lifetime::Timeframe);
  inputs.emplace_back("tracks", "MCH", "TRACKS", 0, Lifetime::Timeframe);
  inputs.emplace_back("clusters", "MCH", "TRACKCLUSTERS", 0, Lifetime::Timeframe);

  auto ccdbRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                                false,                             // GRPECS=true
                                                                false,                             // GRPLHCIF
                                                                true,                              // GRPMagField
                                                                false,                             // askMatLUT
                                                                o2::base::GRPGeomRequest::Aligned, // geometry
                                                                inputs);
  specs.emplace_back(DataProcessorSpec{  //configure le device : nom du fichier, outputs, inputs
    "mch-minv",
    inputs,
    Outputs{},
    AlgorithmSpec{adaptFromTask<eval::MinvTask>(ccdbRequest)},
    Options{
      {"outfile", VariantType::String, "minv.root", {"output Root filename"}}
    }});
      

  return specs;
}
