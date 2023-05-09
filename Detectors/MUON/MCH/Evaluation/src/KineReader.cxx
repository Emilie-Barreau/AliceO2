#include "KineReader.h"
#include "Histogrammer.h"
#include "Steer/MCKinematicsReader.h"
#include "MinvTask.h"
#include "SimulationDataFormat/MCTrack.h"
#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <vector>
#include <Math/Vector4D.h>
#include <Math/GenVector/VectorUtil.h>

using namespace std;
using namespace o2::steer;
using namespace o2::mch::eval;

KineReader::KineReader() {}   

//function oustide classes, gives a Lorentzvector from a MCtrack
ROOT::Math::PxPyPzMVector getLorentzVector(const o2::MCTrack& t)   
{
  return ROOT::Math::PxPyPzMVector(t.Px(), t.Py(), t.Pz(), t.GetMass());
}

int main()
{
  //calling objects from other classes
  MCKinematicsReader r("sgn", MCKinematicsReader::Mode::kMCKine); 
  Histogrammer histogrammer;

  //text files to know how many muons are primary or not (primary -> J/Psi decay)
  ofstream myfile, myfile2;
  myfile.open("Primary.txt");
  myfile2.open("Second.txt");

  //loop for 1D histograms -> all tracks from all event
  for (int evt = 0; evt < r.getNEvents(0); evt++) {
    for (int trk = 0; trk < r.getTracks(evt).size(); trk++) {
      auto t = r.getTrack(evt, trk);
      auto lv = getLorentzVector(*t);
      if (t->GetPdgCode() == 13 || t->GetPdgCode() == -13) {   
        histogrammer.fillSingleParticleHistos(lv);
      } else if (t->GetPdgCode() == 443) {
        // rap->Fill(r.getTrack(evt, trk)->GetRapidity());
      } else {
        continue;
      }
      if (t->isPrimary() == 1) {
        myfile << r.getTrack(evt, trk)->GetPdgCode() << "\n";
      } else if (r.getTrack(evt, trk)->isSecondary() == 1) {
        myfile2 << r.getTrack(evt, trk)->GetPdgCode() << "\n";
      }
    }
  }
  myfile.close();
  myfile2.close();

  //loop for 2D histograms -> dimuons coming from J/Psi
  for (int Evt = 0; Evt < r.getNEvents(0); Evt++) {
    for (int mu1 = 0; mu1 < r.getTracks(Evt).size(); mu1++) {
      auto t1 = r.getTrack(Evt, mu1);
      for (int mu2 = mu1 + 1; mu2 < r.getTracks(Evt).size(); mu2++) {
        auto t2 = r.getTrack(Evt, mu2);
        if (t1->GetPdgCode() == 13 && t2->GetPdgCode() == -13 || t1->GetPdgCode() == -13 && t2->GetPdgCode() == 13) {
          if (t1->isPrimary() == 1 && t2->isPrimary() == 1) {
            auto lv1 = getLorentzVector(*t1);
            auto lv2 = getLorentzVector(*t2);
            histogrammer.fillDoubleParticleHistos(lv1, lv2);
          } else {
            continue;
          }
        } else {
          continue;
        }
      }
    }
  }

  //saving the histograms
  histogrammer.save("Histos.root");

  return 0;
}
