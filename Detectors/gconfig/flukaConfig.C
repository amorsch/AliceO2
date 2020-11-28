// Configuration macro for FLUKA VirtualMC

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "TFluka.h"
#include "TSystem.h"
#include "FairRunSim.h"
#include <iostream>
#endif
#include "commonConfig.C"
void linkFlukaFiles();
void Config()
{
  linkFlukaFiles();
  FairRunSim* run = FairRunSim::Instance();
  TString* gModel = run->GetGeoModel();
  TFluka* fluka = new TFluka("C++ Interface to Fluka", 0);
  stackSetup(fluka, run);

  // setup decayer
  decayerSetup(fluka);

  // ******* FLUKA  specific configuration for simulated Runs  *******
}

void linkFlukaFiles()
{
  // Link here some special Fluka files needed
  gSystem->Exec("ln -s $FLUKADATA/neuxsc.bin  .");
  gSystem->Exec("ln -s $FLUKADATA/elasct.bin  .");
  gSystem->Exec("ln -s $FLUKADATA/gxsect.bin  .");
  gSystem->Exec("ln -s $FLUKADATA/nuclear.bin .");
  gSystem->Exec("ln -s $FLUKADATA/sigmapi.bin .");
  gSystem->Exec("ln -s $FLUKADATA/brems_fin.bin .");
  gSystem->Exec("ln -s $FLUKADATA/cohff.bin .");
  gSystem->Exec("ln -s $FLUKADATA/fluodt.dat  .");
  gSystem->Exec("ln -s $FLUKADATA/random.dat  .");
  // Copy the random seed
  gSystem->Exec("cp $FLUKADATA/random.dat old.seed");
  // Give some meaningfull name to the output
  gSystem->Exec("ln -s fluka.out fort.11");
  gSystem->Exec("ln -s fluka.err fort.15");
  gSystem->Exec("ln -fs $ALICE_ROOT/TFluka/macro/FlukaConfig.C Config.C");
  gSystem->Exec("ln -fs $O2_ROOT/share/Detectors/gconfig/data/coreFlukaVmc.inp .");
}
