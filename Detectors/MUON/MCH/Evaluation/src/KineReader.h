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

#ifndef O2_MCH_EVALUATION_KINE_READER_H
#define O2_MCH_EVALUATION_KINE_READER_H

#include "Math/Vector4D.h"

namespace o2::mch::eval
{

// ===== I AM TRYING =====
class KineReader  //struct = all public
{

 public:
  KineReader();                                             // constructor \°w°/

  void createHistos();                               // static function because no object 
  void fillHistos(ROOT::Math::PxPyPzMVector lor);           // function filling histos with Lorentz vectors

 public:
  ROOT::Math::PxPyPzMVector Lorentz_Gen1(const double a, const double b, const double c, const double d);
  ROOT::Math::PxPyPzMVector Lorentz_Gen2(const double a, const double b, const double c, const double d);
  int Nevents;
};

} // namespace o2::mch::eval

#endif