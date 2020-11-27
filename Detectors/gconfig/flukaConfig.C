// Configuration macro for FLUKA VirtualMC

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "TFluka.h"
#include "FairRunSim.h"
#include <iostream>
#endif
#include "commonConfig.C"

void Config()
{
  FairRunSim* run = FairRunSim::Instance();
  TString* gModel = run->GetGeoModel();
  TFluka* fluka = new TFluka("C++ Interface to Fluka", 0);
  stackSetup(fluka, run);

  // setup decayer
  decayerSetup(fluka);

  // ******* FLUKA  specific configuration for simulated Runs  *******
}
