#include "Histogrammer.h"
#include <iostream>
#include "TH1F.h"

Histogrammer::Histogrammer() {
  mHistos.emplace_back(new TH1F("hphi","hphi",100,0,180));
}
void Histogrammer::fillSingleParticleHisto(const ROOT::Math::PxPyPzMVector& p){
  std::cout << "TEST" << p << std::endl;
  mHistos[0]->Fill(p.Phi());
}

/*
#include "src/MinvTask.h"
#include "src/KineReader.h"
#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <vector>
#include <Math/Vector4D.h>

using namespace std;
using namespace o2::mch::eval;
//using namespace o2::steer; ->reconnais pas ?

class TH1;
class TH2;
class TF1;
class KineReader;
class MinvTask;

int main()
{

  unique_ptr<TFile> myFile(TFile::Open("Reco_gener.root", "RECREATE"));

}
*/