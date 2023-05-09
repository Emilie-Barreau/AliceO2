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

//==========EXTRAPOLE LA TRACE AU VERTEX PRIMAIRE==========

#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "MCHEvaluation/ExtendedTrack.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackFitter.h"
#include "TGeoGlobalMagField.h"
#include <fmt/format.h>

constexpr double muMass = 0.10565800;  //initialise la masse du muon

namespace o2::mch::eval
{
TrackFitter& trackFitter()  //verifie que le champ magnetique est actif ?
{
  static TrackFitter trackFitter;
  static bool isInitialized{false};
  if (!isInitialized) {
    if (!TGeoGlobalMagField::Instance()->GetField()) {
      throw std::runtime_error("Magnetic field should have been initialized before this call to trackFitter()");
    }
    trackFitter.smoothTracks(true);
    TrackExtrap::useExtrapV2();
    isInitialized = true;
  }
  return trackFitter;
}

ExtendedTrack::ExtendedTrack(const TrackMCH& track,
                             gsl::span<const Cluster>& clusters,
                             double x, double y, double z)  //extended depend des clusters choisis, de la trace et du vertex
{
  auto s = clusters.subspan(track.getFirstClusterIdx(),  //recupere l'ID du premier cluster de la trace et le nombre de clusters
                            track.getNClusters());

  mClusters.assign(begin(s), end(s));  //sequence des clusters pour chaque trace
  for (const auto& cluster : mClusters) {
    mTrack.createParamAtCluster(cluster);
  }

  trackFitter().fit(mTrack);  //fit la trace
  extrapToVertex(x, y, z);  //extrapolation au vertex
}

void ExtendedTrack::extrapToVertex(double x, double y, double z)
{  //compute the track parameters at vertex, at DCA and at the end of the absorber

  //extrapolate to vertex
  TrackParam trackParamAtVertex(mTrack.first());  //mTrack.first() = parameters at first cluster
  TrackExtrap::extrapToVertex(trackParamAtVertex, x, y, z, 0., 0.);
  mMomentum4D.SetPx(trackParamAtVertex.px());  //initialise les paramètres du 4-vecteur au vertex ?
  mMomentum4D.SetPy(trackParamAtVertex.py());
  mMomentum4D.SetPz(trackParamAtVertex.pz());
  mMomentum4D.SetM(muMass);

  //extrapolate to DCA (Distance of Closest Approach)
  TrackParam trackParamAtDCA(mTrack.first());
  TrackExtrap::extrapToVertexWithoutBranson(trackParamAtDCA, z);
  double dcaX = trackParamAtDCA.getNonBendingCoor() - x;  //sur plan non courbe
  double dcaY = trackParamAtDCA.getBendingCoor() - y;  //sur plan courbe
  mDCA = std::sqrt(dcaX * dcaX + dcaY * dcaY);

  //extrapolate to the end of the absorber
  TrackParam trackParam(mTrack.first());
  TrackExtrap::extrapToZ(trackParam, -505.);
  double xAbs = trackParam.getNonBendingCoor();
  double yAbs = trackParam.getBendingCoor();
  mRabs = std::sqrt(xAbs * xAbs + yAbs * yAbs);
}

bool ExtendedTrack::operator==(const ExtendedTrack& track) const  //PAS COMPRIS
{
  return areEqual(*this, track, sChi2Max);
}

const TrackParam& ExtendedTrack::param() const
{
  return mTrack.first();
}

const Track& ExtendedTrack::track() const
{
  return mTrack;
}

bool ExtendedTrack::isMatching(const ExtendedTrack& track) const
{
  return areMatching(*this, track, sChi2Max);
}

bool areEqual(const ExtendedTrack& t1, const ExtendedTrack& t2, double chi2Max)  //compare deux traces avec un Chi2
{
  const auto& clusters1 = t1.getClusters();  //recupere les clusters du 1
  const auto& clusters2 = t2.getClusters();  //recupere les clusters du 2
  if (clusters1.size() != clusters2.size()) {  //si la taille est differente ?
    return false;
  }
  for (size_t iCl = 0; iCl != clusters1.size(); ++iCl) {  //sur l'ensemble des clusters
    const auto& cl1 = clusters1[iCl];
    const auto& cl2 = clusters2[iCl];
    if (cl1.getDEId() != cl2.getDEId()) {  //si l'ID des clusters est different ?
      return false;
    }
    double dx = cl1.getX() - cl2.getX();  //delta x et y des clusters
    double dy = cl1.getY() - cl2.getY();
    double chi2 = dx * dx / (cl1.getEx2() + cl2.getEx2()) + dy * dy / (cl1.getEy2() + cl2.getEy2());  //calcul du Chi2
    if (chi2 > chi2Max) {  //ne doit pas dépasser le Chi2Max
      return false;
    }
  }
  return true;
}

bool areMatching(const ExtendedTrack& t1, const ExtendedTrack& t2, double chi2Max)  //verification si 2 traces matchent ?
{
  size_t nMatchClusters(0);
  bool matchCluster[10] = {false, false, false, false, false, false, false, false, false, false};

  const auto& clusters1 = t1.getClusters();  //recuperation des clusters des traces 1 et 2
  const auto& clusters2 = t2.getClusters();

  for (const auto& cl1 : clusters1) {  //sur les clusters 1
    for (const auto& cl2 : clusters2) {  //sur les clusters 2
      if (cl1.getDEId() == cl2.getDEId()) {  //si l'ID des clusters est identique 
        double dx = cl1.getX() - cl2.getX();  //delta x et y des clusters
        double dy = cl1.getY() - cl2.getY();
        double chi2 = dx * dx / (cl1.getEx2() + cl2.getEx2()) + dy * dy / (cl1.getEy2() + cl2.getEy2());  //calcul du Chi2
        if (chi2 <= chi2Max) {  //si inferieur au Chi2Max
          matchCluster[cl1.getChamberId()] = true;
          ++nMatchClusters;
          break;
        }
      }
    }
  }

  return ((matchCluster[0] || matchCluster[1] || matchCluster[2] || matchCluster[3]) &&
          (matchCluster[6] || matchCluster[7] || matchCluster[8] || matchCluster[9]) &&
          (2 * nMatchClusters > clusters1.size() || 2 * nMatchClusters > clusters2.size()));
}

std::ostream& operator<<(std::ostream& out, const ExtendedTrack& track)
{
  out << "{" << track.asString() << "}";
  return out;
}

double ExtendedTrack::getNormalizedChi2() const  //normalisation du Chi2
{
  double chi2 = mTrack.first().getTrackChi2();
  double ndf = 2.0 * mClusters.size() - 6;
  return chi2 / ndf;
}

std::string ExtendedTrack::asString() const
{
  const auto& param = mTrack.first();
  return fmt::format("x = {:7.2f}, y = {:7.2f}, z = {:7.2f}, px = {:7.2f}, py = {:7.2f}, pz = {:7.2f}, sign = {}",
                     param.getNonBendingCoor(), param.getBendingCoor(), param.getZ(),
                     param.px(), param.py(), param.pz(), param.getCharge());
}
} // namespace o2::mch::eval
