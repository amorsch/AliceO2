// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "FairRunSim.h"
#include "TFluka.h"
#include "SimulationDataFormat/Stack.h"
#include "SimulationDataFormat/StackParam.h"
#include <iostream>
#include "FairLogger.h"
#include "FairModule.h"
#include <DetectorsPassive/Cave.h>
#include "DetectorsBase/MaterialManager.h"
#include "SimSetup/GlobalProcessCutSimParam.h"
#include "Generators/DecayerPythia8.h"

//using declarations here since SetCuts.C and g3Config.C are included within namespace
// these are needed for SetCuts.C inclusion
using o2::GlobalProcessCutSimParam;
using o2::base::ECut;
using o2::base::EProc;
using o2::base::MaterialManager;
// these are used in g3Config.C
using std::cout;
using std::endl;
// these are used in commonConfig.C
using o2::eventgen::DecayerPythia8;
#include <SimSetup/SimSetup.h>

namespace o2
{
namespace flukaconfig
{
#include "../flukaConfig.C"
#include "../SetCuts.h"

void FlukaConfig()
{
  LOG(INFO) << "Setting up FLUKA sim from library code";
  Config();
  SetCuts();
}
} // namespace flukaconfig
} // namespace o2
