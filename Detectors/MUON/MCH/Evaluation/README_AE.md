# Acceptance-efficiency MCH tracks quick estimation

Goal : From MC simulation, this code can estimate the acceptance-efficiency of the MCH part of the Spectrometer, considering the chambers breakdowns for different levels

Src files :
- ExtendedTrack = MCH tracks extended to primary vertex (no need to modify)
- minv-workflow = build the executable (no need to modify)
- Histogrammer = manage the histograms and root tree (initialization of the range, bin...)
- KineReader = manage the generated particles (cuts to select the primary muons) -> generate the Histos_gen.root
- MinvTask = manage the reconstructed muons from the generated ones (tracks and pair selections) -> generate the Histo_reco.root
           = manage the detections elements, Dual Sampa, Pads... to remove -> each level (DE, DS, Solars) is implemented, you just need to fill the ID 
           = provide some informations about the ratio of lost tracks over reconstructed tracks (open the IDENTITY.txt file)

Macro file :
- acceff = manage the histograms showing the acceptance-efficiency (with and without breakdowns)

TIPS : Keep a clean .root file for the reconstructed muons (without breakdowns) then create news .root files for each modification and specify in the file name (ex : Histo_reco_cut_DE3.root)

============================================================================
HOW TO USE IT :
- Import the root files from your simulation in a repository in $HOME/alice/
- Make the modifications you need
- Compile the code (aliBuild or ninja)
- Run the code **inside your simulation repository**
  KineReader code : `$HOME/alice/sw/BUILD/O2-latest/O2/stage/bin/o2-mch-kine-reader`
  MinvTask code : `o2-mch-tracks-reader-workflow | o2-mch-minv-workflow -b`
- In the acceff.C code, check the name of your fgen, freco and freco_cut (should be the name of your root files from previous codes)
- In root environment : `.x $HOME/alice/O2/Detectors/MUON/MCH/Evaluation/src/acceff.C`
- Change the histogram parameters if needed and root again
**BE CAREFUL ABOUT PACKAGES AND O2 VERSION !!!**

============================================================================