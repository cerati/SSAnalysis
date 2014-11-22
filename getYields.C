{

  TFile *ttbar_file = TFile::Open("ttbar_histos.root");
  cout << "ttbar BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR28=" << hyp_highpt_sr->GetBinContent(29) << endl;

  TFile *TTWJets_file = TFile::Open("TTWJets_histos.root");
  cout << "TTWJets BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR28=" << hyp_highpt_sr->GetBinContent(29) << endl;

  TFile *TTZJets_file = TFile::Open("TTZJets_histos.root");
  cout << "TTZJets BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR28=" << hyp_highpt_sr->GetBinContent(29) << endl;

  TFile *WZJets_file = TFile::Open("WZJets_histos.root");
  cout << "WZJets BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR28=" << hyp_highpt_sr->GetBinContent(29) << endl;

  TFile *WW_file = TFile::Open("WW_histos.root");
  cout << "WW BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR28=" << hyp_highpt_sr->GetBinContent(29) << endl;

  TFile *T1ttttG1200_file = TFile::Open("T1ttttG1200_histos.root");
  cout << "T1ttttG1200 BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR28=" << hyp_highpt_sr->GetBinContent(29) << endl;

  TFile *T1ttttG1500_file = TFile::Open("T1ttttG1500_histos.root");
  cout << "T1ttttG1500 BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR28=" << hyp_highpt_sr->GetBinContent(29) << endl;

  TFile *T5Full1200_file = TFile::Open("T5Full1200_histos.root");
  cout << "T5Full1200 BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR28=" << hyp_highpt_sr->GetBinContent(29) << endl;

  TFile *T5Full1500_file = TFile::Open("T5Full1500_histos.root");
  cout << "T5Full100 BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR28=" << hyp_highpt_sr->GetBinContent(29) << endl;

}
