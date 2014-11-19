#include "TChain.h"
#include "looper.h"

int main() {

  looper *l = new looper();

  bool skimAll = false;
  bool runAll  = true;

  //looper::ScanChain( TChain* chain, TString prefix, TString postfix, bool isData, TString whatTest, int nEvents)

  if (skimAll) {
    TChain *chain_T1ttttG1200skim = new TChain("Events");
    chain_T1ttttG1200skim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    l->ScanChain(chain_T1ttttG1200skim,"T1ttttG1200","",0,"SSskim",-1);
    TChain *chain_T1ttttG1500skim = new TChain("Events");
    chain_T1ttttG1500skim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    l->ScanChain(chain_T1ttttG1500skim,"T1ttttG1500","",0,"SSskim",-1);
    TChain *chain_T5Full1200skim = new TChain("Events");
    chain_T5Full1200skim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/T5Full_T5Full-1200-1000-800-Decay-MGMMatch50/merged/merged_ntuple_*.root");
    l->ScanChain(chain_T5Full1200skim,"T5Full1200","",0,"SSskim",-1);
    TChain *chain_T5Full1500skim = new TChain("Events");
    chain_T5Full1500skim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/T5Full_T5Full-1500-800-100-Decay-MGMMatch50/merged/merged_ntuple_*.root");
    l->ScanChain(chain_T5Full1500skim,"T5Full1500","",0,"SSskim",-1);
    TChain *chain_TTJetsskim = new TChain("Events");
    chain_TTJetsskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    l->ScanChain(chain_TTJetsskim,"TTJets","",0,"SSskim",-1);
    TChain *chain_TTWJetsskim = new TChain("Events");
    chain_TTWJetsskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/TTWJets_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    l->ScanChain(chain_TTWJetsskim,"TTWJets","",0,"SSskim",-1);
    TChain *chain_TTZJetsskim = new TChain("Events");
    chain_TTZJetsskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/TTZJets_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    l->ScanChain(chain_TTZJetsskim,"TTZJets","",0,"SSskim",-1);
    TChain *chain_WHZHskim = new TChain("Events");
    chain_WHZHskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/WH_ZH_HToWW_M-125_13TeV_pythia6/merged/merged_ntuple_*.root");
    l->ScanChain(chain_WHZHskim,"WHZH","",0,"SSskim",-1);
    TChain *chain_WWskim = new TChain("Events");
    chain_WWskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/WWTo2L2Nu_CT10_13TeV-powheg-pythia8-tauola/merged/merged_ntuple_*.root");
    l->ScanChain(chain_WWskim,"WW","",0,"SSskim",-1);
    TChain *chain_WZJetsskim = new TChain("Events");
    chain_WZJetsskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    l->ScanChain(chain_WZJetsskim,"WZJets","",0,"SSskim",-1);
  }

  // TChain *chain_ttbarskim = new TChain("Events");
  // chain_ttbarskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
  // l->ScanChain(chain_ttbarskim,"ttbarskim","",0,"SSskim",-1);

  // TChain *chain_qcdskim = new TChain("Events");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-5to10_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-10to15_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-15to30_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-30to50_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-50to80_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-80to120_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-120to170_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-170to300_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-300to470_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-470to600_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-600to800_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-1000to1400_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-1400to1800_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // chain_qcdskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-1800_Tune4C_13TeV_pythia8/merged/merged_ntuple_*.root");
  // l->ScanChain(chain_qcdskim,"qcdskim","",0,"QCDskim",-1);

  // TChain *chain_qcdmuskim = new TChain("Events");
  // chain_qcdmuskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-15to20_MuEnrichedPt5_TuneZ2star_13TeV_pythia6/merged/merged_ntuple_*.root");
  // chain_qcdmuskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-20to30_MuEnrichedPt5_TuneZ2star_13TeV_pythia6/merged/merged_ntuple_*.root");
  // chain_qcdmuskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-30to50_MuEnrichedPt5_TuneZ2star_13TeV_pythia6/merged/merged_ntuple_*.root");
  // chain_qcdmuskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-50to80_MuEnrichedPt5_TuneZ2star_13TeV_pythia6/merged/merged_ntuple_*.root");
  // chain_qcdmuskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/QCD_Pt-80to120_MuEnrichedPt5_TuneZ2star_13TeV_pythia6/merged/merged_ntuple_*.root");
  // l->ScanChain(chain_qcdmuskim,"qcdmuskim",0);

  if (runAll) {
    TChain *chain_ttbar = new TChain("Events");
    chain_ttbar->Add("./TTJets_skimSS/merged_ntuple_*.root");
    //chain_ttbar->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-04/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
    l->ScanChain(chain_ttbar,"ttbar","",0,"",-1);
    
    TChain *chain_TTWJets = new TChain("Events");
    chain_TTWJets->Add("./TTWJets_skimSS/merged_ntuple_*.root");
    l->ScanChain(chain_TTWJets,"TTWJets","",0,"",-1);
    
    TChain *chain_TTZJets = new TChain("Events");
    chain_TTZJets->Add("./TTZJets_skimSS/merged_ntuple_*.root");
    l->ScanChain(chain_TTZJets,"TTZJets","",0,"",-1);
    
    TChain *chain_WHZHJets = new TChain("Events");
    chain_WHZHJets->Add("./WHZHJets_skimSS/merged_ntuple_*.root");
    l->ScanChain(chain_WHZHJets,"WHZHJets","",0,"",-1);
    
    TChain *chain_WW = new TChain("Events");
    chain_WW->Add("./WW_skimSS/merged_ntuple_*.root");
    l->ScanChain(chain_WW,"WW","",0,"",-1);
    
    TChain *chain_WZJets = new TChain("Events");
    chain_WZJets->Add("./WZJets_skimSS/merged_ntuple_*.root");
    l->ScanChain(chain_WZJets,"WZJets","",0,"",-1);
    
    TChain *chain_T1ttttG1200 = new TChain("Events");
    chain_T1ttttG1200->Add("./T1ttttG1200_skimSS/merged_ntuple_*.root");
    l->ScanChain(chain_T1ttttG1200,"T1ttttG1200","",0,"",-1);
    
    TChain *chain_T1ttttG1500 = new TChain("Events");
    chain_T1ttttG1500->Add("./T1ttttG1500_skimSS/merged_ntuple_*.root");
    l->ScanChain(chain_T1ttttG1500,"T1ttttG1500","",0,"",-1);
    
    TChain *chain_T5Full1200 = new TChain("Events");
    chain_T5Full1200->Add("./T5Full1200_skimSS/merged_ntuple_*.root");
    l->ScanChain(chain_T5Full1200,"T5Full1200","",0,"",-1);
    
    TChain *chain_T5Full1500 = new TChain("Events");
    chain_T5Full1500->Add("./T5Full1500_skimSS/merged_ntuple_*.root");
    l->ScanChain(chain_T5Full1500,"T5Full1500","",0,"",-1);
  }
  
  // TChain *chain_qcd80to120 = new TChain("Events");
  // chain_qcd80to120->Add("./QCD_Pt-80to120/merged_ntuple_*_skimQCD.root");
  // l->ScanChain(chain_qcd80to120,"qcd","pt-80to120",0,"QCDtest",-1);

  // TChain *chain_qcdmupt5_50to80 = new TChain("Events");
  // chain_qcdmupt5_50to80->Add("./QCD_Pt-50to80_MuEnrichedPt5/merged_ntuple_*_skimQCD.root");
  // l->ScanChain(chain_qcdmupt5_50to80,"qcd_mupt5_50to80","",0,"QCDtest",-1);

  // TChain *chain_qcdmupt5 = new TChain("Events");
  // chain_qcdmupt5->Add("./QCD_Pt-15to20_MuEnrichedPt5/merged_ntuple_*_skimQCD.root");
  // chain_qcdmupt5->Add("./QCD_Pt-20to30_MuEnrichedPt5/merged_ntuple_*_skimQCD.root");
  // chain_qcdmupt5->Add("./QCD_Pt-30to50_MuEnrichedPt5/merged_ntuple_*_skimQCD.root");
  // chain_qcdmupt5->Add("./QCD_Pt-50to80_MuEnrichedPt5/merged_ntuple_*_skimQCD.root");
  // chain_qcdmupt5->Add("./QCD_Pt-80to120_MuEnrichedPt5/merged_ntuple_*_skimQCD.root");
  // l->ScanChain(chain_qcdmupt5,"qcd_mupt5","",0,"QCDtest",-1);

  // TChain *chain_qcd = new TChain("Events");
  // // chain_qcd->Add("./QCD_Pt-5to10/merged_ntuple_*_skimQCD.root");
  // // chain_qcd->Add("./QCD_Pt-10to15/merged_ntuple_*_skimQCD.root");
  // chain_qcd->Add("./QCD_Pt-30to50/merged_ntuple_*_skimQCD.root");
  // chain_qcd->Add("./QCD_Pt-50to80/merged_ntuple_*_skimQCD.root");
  // chain_qcd->Add("./QCD_Pt-80to120/merged_ntuple_*_skimQCD.root");
  // chain_qcd->Add("./QCD_Pt-120to170/merged_ntuple_*_skimQCD.root");
  // // chain_qcd->Add("./QCD_Pt-300to470/merged_ntuple_*_skimQCD.root");
  // // chain_qcd->Add("./QCD_Pt-470to600/merged_ntuple_*_skimQCD.root");
  // // chain_qcd->Add("./QCD_Pt-600to800/merged_ntuple_*_skimQCD.root");
  // // chain_qcd->Add("./QCD_Pt-800to1000/merged_ntuple_*_skimQCD.root");
  // // chain_qcd->Add("./QCD_Pt-1000to1400/merged_ntuple_*_skimQCD.root");
  // // chain_qcd->Add("./QCD_Pt-1400to1800/merged_ntuple_*_skimQCD.root");
  // // chain_qcd->Add("./QCD_Pt-1800/merged_ntuple_*_skimQCD.root");
  // l->ScanChain(chain_qcd,"qcd","pt-30to170",0,"QCDtest",-1);

  // TChain *chain_wj = new TChain("Events");
  // chain_wj->Add("./W1234JetsToLNu/merged_ntuple_1.root");
  // l->ScanChain(chain_wj,"wj","fakes",0,"QCDtest",-1);

  // TChain *chain_dy = new TChain("Events");
  // //chain_dy->Add("./DYJetsLL/merged_ntuple_1.root");
  // chain_dy->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/DYJetsToLL_M-50_13TeV-madgraph-pythia8/merged/merged_ntuple_*.root");
  // chain_dy->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
  // chain_dy->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
  // chain_dy->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
  // chain_dy->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
  // //l->ScanChain(chain_dy,"dy","fakes",0,"QCDtest",-1);
  // l->ScanChain(chain_dy,"dy","effic",0,"DYtest",-1);

}
