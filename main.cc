#include "TChain.h"
#include "looper.h"

int main() {

  looper *l = new looper();
  TChain *chain = new TChain("Events");
  //chain->Add("./TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola-PU20bx25_POSTLS170_V5-v2-cms3-skimss-469438.root");
  //chain->Add("./TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola-PU20bx25_POSTLS170_V5-v2-cms3-skimss-469438_newskim.root");
  //chain->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");

  // chain->Add("./TTbar_SkimSS/merged_ntuple_*_newskim.root");
  // l->ScanChain(chain,"ttbar",0);

  // chain->Add("/hadoop/cms/store/user/namin/CMS3_V07-00-03/QCD_Pt-300to470_Tune4C_13TeV_pythia8_StoreResults-Spring14dr_PU_S14_POSTLS170_V6AN1_miniAOD706p1_814812ec83fce2f620905d2bb30e9100-v1/merged/merged_ntuple_*.root");
  // chain->Add("/hadoop/cms/store/user/namin/CMS3_V07-00-03/QCD_Pt-600to800_Tune4C_13TeV_pythia8_StoreResults-Spring14dr_PU20Bx50_POSTLS170_V7AN1_v1_814812ec83fce2f620905d2bb30e9100-v1/merged/merged_ntuple_*.root");
  // chain->Add("/hadoop/cms/store/user/namin/CMS3_V07-00-03/QCD_Pt-800to1000_Tune4C_13TeV_pythia8_StoreResults-Spring14dr_PU20Bx50_POSTLS170_V7AN1_v1_814812ec83fce2f620905d2bb30e9100-v1/merged/merged_ntuple_*.root");

  chain->Add("./QCD_Pt-300to470/merged_ntuple_*_skimQCD.root");
  chain->Add("./QCD_Pt-470to600/merged_ntuple_*_skimQCD.root");
  chain->Add("./QCD_Pt-600to800/merged_ntuple_*_skimQCD.root");
  chain->Add("./QCD_Pt-800to1000/merged_ntuple_*_skimQCD.root");
  chain->Add("./QCD_Pt-1000to1400/merged_ntuple_*_skimQCD.root");
  chain->Add("./QCD_Pt-1400to1800/merged_ntuple_*_skimQCD.root");
  chain->Add("./QCD_Pt-1800/merged_ntuple_*_skimQCD.root");
  l->ScanChain(chain,"qcd",0);

  // chain->Add("./DYJetsLL/merged_ntuple_1.root");
  // l->ScanChain(chain,"dy",0);

}
