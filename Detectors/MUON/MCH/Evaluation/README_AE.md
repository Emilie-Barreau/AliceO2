# Acceptance-efficiency MCH tracks quick estimation

Goal : From MC simulation, this code can estimate the acceptance-efficiency of the MCH part of the Spectrometer, considering the chambers breakdowns for different levels

Src files :
- ExtendedTrack = MCH tracks extended to primary vertex (no need to modify)
- minv-workflow = build the executable (no need to modify)
- Histogrammer = manage the histograms and root tree (initialization of the range, bin...)
- KineReader = manage the generated particles (cuts to select the primary muons) -> generate the Histos_gen.root
- MinvTask = manage the reconstructed muons from the generated ones (tracks and pair selections) -> generate the Histo_reco.root
           = manage the detections elements, Dual Sampa, Pads... to remove -> each level (DE, DS, Solars) is implemented, you just need to fill the ID 
           = provide some informations about the ratio of lost tracks over reconstructed tracks (open the IDENTITY.txt file -> created from MinvTask, inside the simulation repository)

Macro file :
- acceff = manage the histograms showing the acceptance-efficiency (with and without breakdowns)

TIPS : Keep a clean .root file for the reconstructed muons (without breakdowns) then create news .root files for each modification and specify in the file name (ex : Histo_reco_cut_DE3.root)

**FOR THE COMMAND : check your paths, some files/executables might be in different locations**

============================================================================

HOW TO USE IT :
- Import the root files from your simulation in a repository in $HOME/alice/ (see below for full simulation tutorial)
- Make the modifications you need
- Compile the code (aliBuild or ninja)
- Run the code **inside your simulation repository**
  - KineReader code : `$HOME/alice/sw/BUILD/O2-latest/O2/stage/bin/o2-mch-kine-reader`
  - MinvTask code : `o2-mch-tracks-reader-workflow | o2-mch-minv-workflow -b`
- In the acceff.C code, check the name of your fgen, freco and freco_cut (should be the name of your root files from previous codes)
- In root environment : `.x $HOME/alice/O2/Detectors/MUON/MCH/Evaluation/src/acceff.C`
- Change the histogram parameters if needed and root again
**BE CAREFUL ABOUT PACKAGES AND O2 VERSION !!!**

============================================================================
 
SIMULATION WITH LXPLUS :
- Go to lxplus : `ssh -X youlogin@lxplus.cern.ch`
- Source the environment : `source /cvmfs/alice-nightlies.cern.ch/bin/alienv enter VO_ALICE@O2sim::v20230413-1`
- Create a new repository for your simulation
- Check for the generator file : `O2DPG_ROOT/MC/config/PWGDQ/external/generator/GeneratorParamPromptJpsiToMuonEvtGen_pp13TeV.C` **choose the right one**
- Run the command inside the simulation repository : `o2-sim --timestamp 1669594219618 -j 4 -n 5000 -g external -m HALL MAG DIPO COMP PIPE ABSO SHIL MCH MID -o sgn  --configKeyValues "GeneratorExternal.fileName=$O2DPG_ROOT/MC/config/PWGDQ/external/generator/GeneratorParamPromptJpsiToMuonEvtGen_pp13TeV.C;GeneratorExternal.funcName=GeneratorParamPromptJpsiToMuonEvtGen_pp13TeV()"` (change the parameters if needed)
**Active your token for the CCDB access**
- Digitalization step : `o2-sim-digitizer-workflow -b --sims=sgn`
- Reconstruction step : `o2-mch-reco-workflow -b`
- Exit lxplus, and go back to the O2/latest environment
- Copy your lxplus simulation repository : `scp -r youlogin@lxplus.cern.ch:/afs/cern.ch/user/y/youlogin/simulationrepository ./copyoftherepository`

============================================================================