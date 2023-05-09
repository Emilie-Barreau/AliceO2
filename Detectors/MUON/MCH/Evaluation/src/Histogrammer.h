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

#ifndef HISTOGRAMMER_H
#define HISTOGRAMMER_H

#include "Math/Vector4D.h"
#include <vector>

class TH1;
class TH2;

class Histogrammer {

public:

Histogrammer();

void fillSingleParticleHistos(const ROOT::Math::PxPyPzMVector& lor);

void fillDoubleParticleHistos(const ROOT::Math::PxPyPzMVector& lor1, const ROOT::Math::PxPyPzMVector& lor2);

void save(const char* filename);

double GetInvariantMass(const ROOT::Math::PxPyPzMVector& lor1, const ROOT::Math::PxPyPzMVector& lor2);

private:
  //std::unique_ptr<TFile> mRootFile;
  std::vector<TH1*> mHistos;
  std::vector<TH2*> mHistos2;
  double minv;
};

#endif