#include "TChain.h"
#include "looper.h"

int main() {

  looper *l = new looper();
  TChain *chain = new TChain("Events");
  chain->Add("./TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola-PU20bx25_POSTLS170_V5-v2-cms3-skimss-469438.root");
  l->ScanChain(chain,"test",0);

}
