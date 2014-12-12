#include "TChain.h"
#include "looper.h"

int main() {

  looper *l = new looper();

  bool useSkim = false;

  bool skimAll = false;
  bool runAll  = false;
  bool runLepEff = false;

  //looper::ScanChain( TChain* chain, TString prefix, TString postfix, bool isData, TString whatTest, int nEvents)

  // TChain *chain_ttbar2 = new TChain("Events");
  // chain_ttbar2->Add("./TTJets_skimSS/merged_ntuple_*.root");
  // //chain_ttbar->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
  // l->ScanChain(chain_ttbar2,"ttbar2","",0,"",-1);

    // TChain *chain_TTZJets = new TChain("Events");
    // chain_TTZJets->Add("./TTZJets_skimSS/merged_ntuple_*.root");
    // l->ScanChain(chain_TTZJets,"TTZJets_test","",0,"",-1);

  TChain *chain_T1ttttG1200 = new TChain("Events");
  TChain *chain_T1ttttG1500 = new TChain("Events");
  TChain *chain_T5Full1200 = new TChain("Events");
  TChain *chain_T5Full1500 = new TChain("Events");
  TChain *chain_TTJets = new TChain("Events");
  TChain *chain_TTWJets = new TChain("Events");
  TChain *chain_TTZJets = new TChain("Events");
  TChain *chain_WHZH = new TChain("Events");
  TChain *chain_WW = new TChain("Events");
  TChain *chain_WZJets = new TChain("Events");
  TChain *chain_DY = new TChain("Events");

  if (useSkim) {
    chain_TTJets->Add("./TTJets_skimSS/merged_ntuple_*.root");
    chain_TTWJets->Add("./TTWJets_skimSS/merged_ntuple_*.root");
    chain_TTZJets->Add("./TTZJets_skimSS/merged_ntuple_*.root");
    chain_WW->Add("./WW_skimSS/merged_ntuple_*.root");
    chain_WZJets->Add("./WZJets_skimSS/merged_ntuple_*.root");
    chain_WHZH->Add("./WHZH_skimSS/merged_ntuple_*.root");
    chain_T1ttttG1200->Add("./T1ttttG1200_skimSS/merged_ntuple_*.root");
    chain_T1ttttG1500->Add("./T1ttttG1500_skimSS/merged_ntuple_*.root");
    chain_T5Full1200->Add("./T5Full1200_skimSS/merged_ntuple_*.root");
    chain_T5Full1500->Add("./T5Full1500_skimSS/merged_ntuple_*.root");
    chain_DY->Add("./DYJetsLL/merged_ntuple_1.root");//fixme
  } else {
    chain_T1ttttG1200->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    chain_T1ttttG1500->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    chain_T5Full1200->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/T5Full_T5Full-1200-1000-800-Decay-MGMMatch50/merged/merged_ntuple_*.root");
    chain_T5Full1500->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/T5Full_T5Full-1500-800-100-Decay-MGMMatch50/merged/merged_ntuple_*.root");
    chain_TTJets->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    chain_TTWJets->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/TTWJets_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    chain_TTZJets->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/TTZJets_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    chain_WHZH->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/WH_ZH_HToWW_M-125_13TeV_pythia6/merged/merged_ntuple_*.root");
    chain_WW->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/merged/merged_ntuple_*.root");
    chain_WZJets->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    chain_DY->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/DYJetsToLL_M-50_13TeV-madgraph-pythia8/merged/merged_ntuple_1.root");//fixme
  }

  if (skimAll) {
    l->ScanChain(chain_T1ttttG1200,"T1ttttG1200","",0,"SSskim",-1);
    l->ScanChain(chain_T1ttttG1500,"T1ttttG1500","",0,"SSskim",-1);
    l->ScanChain(chain_T5Full1200,"T5Full1200","",0,"SSskim",-1);
    l->ScanChain(chain_T5Full1500,"T5Full1500","",0,"SSskim",-1);
    l->ScanChain(chain_TTJets,"TTJets","",0,"SSskim",-1);
    l->ScanChain(chain_TTWJets,"TTWJets","",0,"SSskim",-1);
    l->ScanChain(chain_TTZJets,"TTZJets","",0,"SSskim",-1);
    l->ScanChain(chain_WHZH,"WHZH","",0,"SSskim",-1);
    l->ScanChain(chain_WW,"WW","",0,"SSskim",-1);
    l->ScanChain(chain_WZJets,"WZJets","",0,"SSskim",-1);
  }

  //l->ScanChain(chain_T1ttttG1500,"T1ttttG1500","",0,"",-1);
  //l->ScanChain(chain_TTJets,"ttbar","",0,"",-1);
  if (runAll) {
    l->ScanChain(chain_TTJets,"ttbar","",0,"",-1);
    l->ScanChain(chain_TTWJets,"TTWJets","",0,"",-1);
    l->ScanChain(chain_TTZJets,"TTZJets","",0,"",-1);
    l->ScanChain(chain_WHZH,"WHZH","",0,"",-1);
    l->ScanChain(chain_WW,"WW","",0,"",-1);
    l->ScanChain(chain_WZJets,"WZJets","",0,"",-1);
    l->ScanChain(chain_T1ttttG1200,"T1ttttG1200","",0,"",-1);
    l->ScanChain(chain_T1ttttG1500,"T1ttttG1500","",0,"",-1);
    l->ScanChain(chain_T5Full1200,"T5Full1200","",0,"",-1);
    l->ScanChain(chain_T5Full1500,"T5Full1500","",0,"",-1);
  }
  
  if (runLepEff) {
    l->ScanChain(chain_DY,"dy","effic",0,"DYtest",-1);
    l->ScanChain(chain_T1ttttG1200,"T1ttttG1200","effic",0,"DYtest",-1);
    l->ScanChain(chain_T1ttttG1500,"T1ttttG1500","effic",0,"DYtest",-1);
    l->ScanChain(chain_T5Full1200,"T5Full1200","effic",0,"DYtest",-1);
    l->ScanChain(chain_T5Full1500,"T5Full1500","effic",0,"DYtest",-1);
    l->ScanChain(chain_TTWJets,"TTWJets","effic",0,"DYtest",-1);
  }

}
