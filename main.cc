#include "TChain.h"
#include "looper.h"

int main() {

  looper *l = new looper();
  TChain *chain = new TChain("Events");
  //chain->Add("./TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola-PU20bx25_POSTLS170_V5-v2-cms3-skimss-469438.root");
  //chain->Add("./TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola-PU20bx25_POSTLS170_V5-v2-cms3-skimss-469438_newskim.root");
  //chain->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");

  // chain->Add("./TTbar_SkimSS/merged_ntuple_*_newskim.root");

  // chain->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-300to470_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-470to600_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-600to800_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-1000to1400_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-1400to1800_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
   chain->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-1800_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");


  l->ScanChain(chain,"test",0);

}
