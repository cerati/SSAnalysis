#include "TChain.h"
#include "looper.h"

int main() {

  looper *l = new looper();


  // TChain *chain_ttbarskim = new TChain("Events");
  // chain_ttbarskim->Add("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/merged/merged_ntuple_*.root");
  // l->ScanChain(chain,"ttbarskim_ttbarskim","",0,"SSskim",-1);

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

  TChain *chain_ttbar = new TChain("Events");
  chain_ttbar->Add("./TTbar_SkimSS/merged_ntuple_*.root");
  l->ScanChain(chain_ttbar,"ttbar","",0,"",-1);

  TChain *chain_qcd80to120 = new TChain("Events");
  chain_qcd80to120->Add("./QCD_Pt-80to120/merged_ntuple_*_skimQCD.root");
  l->ScanChain(chain_qcd80to120,"qcd","pt-80to120",0,"QCDtest",-1);

  TChain *chain_qcdmupt5_50to80 = new TChain("Events");
  chain_qcdmupt5_50to80->Add("./QCD_Pt-50to80_MuEnrichedPt5/merged_ntuple_*_skimQCD.root");
  l->ScanChain(chain_qcdmupt5_50to80,"qcd_mupt5_50to80","",0,"QCDtest",-1);

  TChain *chain_qcdmupt5 = new TChain("Events");
  chain_qcdmupt5->Add("./QCD_Pt-15to20_MuEnrichedPt5/merged_ntuple_*_skimQCD.root");
  chain_qcdmupt5->Add("./QCD_Pt-20to30_MuEnrichedPt5/merged_ntuple_*_skimQCD.root");
  chain_qcdmupt5->Add("./QCD_Pt-30to50_MuEnrichedPt5/merged_ntuple_*_skimQCD.root");
  chain_qcdmupt5->Add("./QCD_Pt-50to80_MuEnrichedPt5/merged_ntuple_*_skimQCD.root");
  chain_qcdmupt5->Add("./QCD_Pt-80to120_MuEnrichedPt5/merged_ntuple_*_skimQCD.root");
  l->ScanChain(chain_qcdmupt5,"qcd_mupt5","",0,"QCDtest",-1);

  TChain *chain_qcd = new TChain("Events");
  // chain_qcd->Add("./QCD_Pt-5to10/merged_ntuple_*_skimQCD.root");
  // chain_qcd->Add("./QCD_Pt-10to15/merged_ntuple_*_skimQCD.root");
  chain_qcd->Add("./QCD_Pt-30to50/merged_ntuple_*_skimQCD.root");
  chain_qcd->Add("./QCD_Pt-50to80/merged_ntuple_*_skimQCD.root");
  chain_qcd->Add("./QCD_Pt-80to120/merged_ntuple_*_skimQCD.root");
  chain_qcd->Add("./QCD_Pt-120to170/merged_ntuple_*_skimQCD.root");
  // chain_qcd->Add("./QCD_Pt-300to470/merged_ntuple_*_skimQCD.root");
  // chain_qcd->Add("./QCD_Pt-470to600/merged_ntuple_*_skimQCD.root");
  // chain_qcd->Add("./QCD_Pt-600to800/merged_ntuple_*_skimQCD.root");
  // chain_qcd->Add("./QCD_Pt-800to1000/merged_ntuple_*_skimQCD.root");
  // chain_qcd->Add("./QCD_Pt-1000to1400/merged_ntuple_*_skimQCD.root");
  // chain_qcd->Add("./QCD_Pt-1400to1800/merged_ntuple_*_skimQCD.root");
  // chain_qcd->Add("./QCD_Pt-1800/merged_ntuple_*_skimQCD.root");
  l->ScanChain(chain_qcd,"qcd","pt-30to170",0,"QCDtest",-1);

  TChain *chain_wj = new TChain("Events");
  chain_wj->Add("./W1234JetsToLNu/merged_ntuple_1.root");
  l->ScanChain(chain_wj,"wj","fakes",0,"QCDtest",-1);

  TChain *chain_dy = new TChain("Events");
  chain_dy->Add("./DYJetsLL/merged_ntuple_1.root");
  l->ScanChain(chain_dy,"dy","fakes",0,"QCDtest",-1);
  l->ScanChain(chain_dy,"dy","effic",0,"DYtest",-1);

}
