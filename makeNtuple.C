
#include <iostream>

#include <TFile.h>
#include "TTree.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TLeaf.h"
#include "TBranch.h"
//#include "datamodel/ParticleData.h"
#include <ROOT/RDataFrame.hxx>
#include <cmath>
#include <vector>
#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
//#include "datamodel/MET.h"

using namespace std;

int main()

{

  cout << "RDataFrame: opening root file and getting the events tree" << endl;

  //ROOT::RDataFrame d("events", "/eos/experiment/fcc/ee/generation/DelphesEvents/fcc_v01/p8_ee_ZH_ecm240/events_000038119.root");
  ROOT::RDataFrame d("events", "../FCCSW/ZXRootFiles/ee_ZX_97GeV.root");
  
  auto d2 = d.Define("muons_px","muons.core.p4.px")
    .Define("muons_py","muons.core.p4.py")
    .Define("muons_pz","muons.core.p4.pz")
    .Define("jets_px","jets.core.p4.px")
    .Define("jets_py","jets.core.p4.py")
    .Define("jets_pz","jets.core.p4.pz")
    .Define("jets_M","jets.core.p4.mass");
           

  cout << "RDataFrame: Writing contents to a new tree" << endl;
  d2.Snapshot("events","./ZXRootfiles/FCCee_ZX_97GeV.root",{"muons_px","muons_py","muons_pz","jets_px","jets_py","jets_pz","jets_M"});

  cout << "opening new file and reading its contents" << endl;
  TFile ff("./ZXRootfiles/FCCee_Test_ZX100k.root");
  TTree *tf = (TTree*)ff.Get("events");
  cout<<"New Tree Entries:"<<tf->GetEntries()<<endl;

  /*vector<double> *mpx = new vector<double>();  
  tf->SetBranchAddress("muons_px",&mpx);
  tf->GetEntry(0);
  double k;
  k  = (*mpx)[0];
  cout << mpx->size() << endl;*/
  


  return 0;
}
