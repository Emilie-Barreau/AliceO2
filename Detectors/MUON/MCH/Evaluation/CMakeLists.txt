# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

o2_add_library(MCHEvaluation
               SOURCES
                 src/CompareTask.cxx
                 src/CompareTracks.cxx
                 src/Draw.cxx
                 src/MinvTask.cxx
                 src/ExtendedTrack.cxx
                 src/Histos.cxx
                 src/Histogrammer.cxx
               PUBLIC_LINK_LIBRARIES
                  O2::DataFormatsMCH
                  O2::DetectorsRaw
                  O2::DPLUtils
                  O2::Framework
                  O2::MCHTracking
                  O2::MCHGeometryTransformer
                  O2::MCHMappingInterface
                  O2::MCHGlobalMapping
                  fmt::fmt)

o2_add_executable(compare-tracks-worfklow
                  SOURCES
                    src/compare-tracks-workflow.cxx
                  COMPONENT_NAME mch
                  PUBLIC_LINK_LIBRARIES O2::MCHEvaluation)
                  
o2_add_executable(minv-workflow
                  SOURCES
                    src/minv-workflow.cxx
                  COMPONENT_NAME mch
                  PUBLIC_LINK_LIBRARIES
                    O2::MCHEvaluation
                    O2::MCHMappingImpl4
                  )

o2_add_executable(kine-reader
                  SOURCES 
                    src/KineReader.cxx
                  COMPONENT_NAME mch
                  PUBLIC_LINK_LIBRARIES
                    O2::Steer
                    O2::MCHEvaluation
                  )

o2_target_root_dictionary(MCHEvaluation
                          HEADERS
                            include/MCHEvaluation/ExtendedTrack.h
                            include/MCHEvaluation/Draw.h
                          LINKDEF src/MCHEvaluationLinkDef.h)


o2_add_test_root_macro(src/acceff.C
                       #PUBLIC_LINK_LIBRARIES O2::Field
                       LABELS mch eval)
