# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_dpl_workflow(mc-spectra-efficiency
                    SOURCES mcspectraefficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2::AnalysisDataModel O2::AnalysisCore
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(spectra-tof
                    SOURCES spectraTOF.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2::AnalysisDataModel
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(spectra-tpc
                    SOURCES spectraTPC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2::AnalysisDataModel
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(nuclei-spectra
                    SOURCES NucleiSpectraTask.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2::AnalysisDataModel
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(lambdakzeroproducer
                    SOURCES lambdakzeroproducer
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2::AnalysisDataModel O2::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(lambdakzeroconsumer
                    SOURCES lambdakzeroconsumer
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2::AnalysisDataModel O2::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(lambdakzerofinder
                    SOURCES lambdakzerofinder
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2::AnalysisDataModel O2::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(cascadeproducer
                    SOURCES cascadeproducer
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2::AnalysisDataModel O2::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(cascadeconsumer
                    SOURCES cascadeconsumer
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2::AnalysisDataModel O2::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(raacharged
                    SOURCES raacharged.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisDataModel O2::AnalysisCore
                    COMPONENT_NAME Analysis)
