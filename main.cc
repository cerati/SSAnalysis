#include "TChain.h"
#include "looper.h"

int main() {

  looper *l = new looper();
  TChain *chain = new TChain("Events");
  chain->Add("/tas/cerati/SSAnalysis/SSAnalysis/CMSSW_7_0_6_patch1/src/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola-PU20bx25_POSTLS170_V5-v2-cms3-25k.root");
  //chain->Add("/tas/gzevi/files/miniAOD/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_PU20bx25_POSTLS170_V5-v1_PUiso.root");
  //chain->Add("/tas/gzevi/files/miniAOD/DYJetsToLL_M-50_13TeV-madgraph-pythia8_PAT_CMS3.root");
  //chain->Add("/tas/cerati/SSAnalysis/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_PU20bx25_POSTLS170_V5-v1_BIG.root");
  //chain->Add("/tas/cerati/SSAnalysis/MT2ntuple1and2_2000.root");
  l->ScanChain(chain,"test",0);

}
