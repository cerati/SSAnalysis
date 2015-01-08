{

  //TString dir = "results_PHYS14";
  TString dir = ".";

  TFile *ttbar_file = TFile::Open(dir+"/ttbar_histos.root");
  cout << "ttbar BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21)+hyp_highpt_sr->GetBinContent(31)
       << " BR0E=" << hyp_highpt_br->GetBinContent(1) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR30=" << hyp_highpt_sr->GetBinContent(31) << endl;

  TFile *TTWJets_file = TFile::Open(dir+"/TTWJets_histos.root");
  cout << "TTWJets BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21)+hyp_highpt_sr->GetBinContent(31)
       << " BR0E=" << hyp_highpt_br->GetBinContent(1) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR30=" << hyp_highpt_sr->GetBinContent(31) << endl;

  TFile *TTZJets_file = TFile::Open(dir+"/TTZJets_histos.root");
  cout << "TTZJets BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21)+hyp_highpt_sr->GetBinContent(31) 
       << " BR0E=" << hyp_highpt_br->GetBinContent(1) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR30=" << hyp_highpt_sr->GetBinContent(31) << endl;

  TFile *WZJets_file = TFile::Open(dir+"/WZJets_histos.root");
  cout << "WZJets BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21)+hyp_highpt_sr->GetBinContent(31) 
       << " BR0E=" << hyp_highpt_br->GetBinContent(1) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR30=" << hyp_highpt_sr->GetBinContent(31) << endl;

  // TFile *WW_file = TFile::Open(dir+"/WW_histos.root");
  // cout << "WW BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21)+hyp_highpt_sr->GetBinContent(31) 
  //<< " BR0E=" << hyp_highpt_br->GetBinContent(1) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR30=" << hyp_highpt_sr->GetBinContent(31) << endl;

  // TFile *WHZH_file = TFile::Open(dir+"/WHZH_histos.root");
  // cout << "WHZH BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21)+hyp_highpt_sr->GetBinContent(31) 
  //<< " BR0E=" << hyp_highpt_br->GetBinContent(1) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR30=" << hyp_highpt_sr->GetBinContent(31) << endl;

  TFile *T1ttttG1200_file = TFile::Open(dir+"/T1ttttG1200_histos.root");
  cout << "T1ttttG1200 BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21)+hyp_highpt_sr->GetBinContent(31) 
       << " BR0E=" << hyp_highpt_br->GetBinContent(1) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR30=" << hyp_highpt_sr->GetBinContent(31) << endl;

  TFile *T1ttttG1500_file = TFile::Open(dir+"/T1ttttG1500_histos.root");
  cout << "T1ttttG1500 BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21)+hyp_highpt_sr->GetBinContent(31) 
       << " BR0E=" << hyp_highpt_br->GetBinContent(1) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR30=" << hyp_highpt_sr->GetBinContent(31) << endl;

  // TFile *T5Full1200_file = TFile::Open(dir+"/T5Full1200_histos.root");
  // cout << "T5Full1200 BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21)+hyp_highpt_sr->GetBinContent(31) 
  //<< " BR0E=" << hyp_highpt_br->GetBinContent(1) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR30=" << hyp_highpt_sr->GetBinContent(31) << endl;

  // TFile *T5Full1500_file = TFile::Open(dir+"/T5Full1500_histos.root");
  // cout << "T5Full100 BR0I=" << hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21)+hyp_highpt_sr->GetBinContent(31) 
  //<< " BR0E=" << hyp_highpt_br->GetBinContent(1) << " BR10=" << hyp_highpt_br->GetBinContent(11) << " BR20=" << hyp_highpt_br->GetBinContent(21) << " SR30=" << hyp_highpt_sr->GetBinContent(31) << endl;


  // TString dataset = "TTZJets";
  // TFile *_file = TFile::Open(dir+"/"+dataset+"_histos.root");
  // cout << "Dataset = " << dataset << endl
  //      << Form("BR0I=%.2f",hyp_highpt_br->GetBinContent(1)+hyp_highpt_br->GetBinContent(11)+hyp_highpt_br->GetBinContent(21))  << endl
  //      << Form("SR01=%.2f",hyp_highpt_sr->GetBinContent(2)) << endl
  //      << Form("SR02=%.2f",hyp_highpt_sr->GetBinContent(3)) << endl
  //      << Form("SR03=%.2f",hyp_highpt_sr->GetBinContent(4)) << endl
  //      << Form("SR04=%.2f",hyp_highpt_sr->GetBinContent(5)) << endl
  //      << Form("SR05=%.2f",hyp_highpt_sr->GetBinContent(6)) << endl
  //      << Form("SR06=%.2f",hyp_highpt_sr->GetBinContent(7)) << endl
  //      << Form("SR07=%.2f",hyp_highpt_sr->GetBinContent(8)) << endl
  //      << Form("SR08=%.2f",hyp_highpt_sr->GetBinContent(9)) << endl
  //      << Form("BR10=%.2f",hyp_highpt_br->GetBinContent(11))  << endl
  //      << Form("SR11=%.2f",hyp_highpt_sr->GetBinContent(12)) << endl
  //      << Form("SR12=%.2f",hyp_highpt_sr->GetBinContent(13)) << endl
  //      << Form("SR13=%.2f",hyp_highpt_sr->GetBinContent(14)) << endl
  //      << Form("SR14=%.2f",hyp_highpt_sr->GetBinContent(15)) << endl
  //      << Form("SR15=%.2f",hyp_highpt_sr->GetBinContent(16)) << endl
  //      << Form("SR16=%.2f",hyp_highpt_sr->GetBinContent(17)) << endl
  //      << Form("SR17=%.2f",hyp_highpt_sr->GetBinContent(18)) << endl
  //      << Form("SR18=%.2f",hyp_highpt_sr->GetBinContent(19)) << endl
  //      << Form("BR20=%.2f",hyp_highpt_br->GetBinContent(21))  << endl
  //      << Form("SR21=%.2f",hyp_highpt_sr->GetBinContent(22)) << endl
  //      << Form("SR22=%.2f",hyp_highpt_sr->GetBinContent(23)) << endl
  //      << Form("SR23=%.2f",hyp_highpt_sr->GetBinContent(24)) << endl
  //      << Form("SR24=%.2f",hyp_highpt_sr->GetBinContent(25)) << endl
  //      << Form("SR25=%.2f",hyp_highpt_sr->GetBinContent(26)) << endl
  //      << Form("SR26=%.2f",hyp_highpt_sr->GetBinContent(27)) << endl
  //      << Form("SR27=%.2f",hyp_highpt_sr->GetBinContent(28)) << endl
  //      << Form("SR28=%.2f",hyp_highpt_sr->GetBinContent(29)) << endl;
  //      << Form("BR30=%.2f",hyp_highpt_br->GetBinContent(31))  << endl
  //      << Form("SR31=%.2f",hyp_highpt_sr->GetBinContent(32)) << endl
  //      << Form("SR32=%.2f",hyp_highpt_sr->GetBinContent(33)) << endl
  //      << Form("SR33=%.2f",hyp_highpt_sr->GetBinContent(34)) << endl
  //      << Form("SR34=%.2f",hyp_highpt_sr->GetBinContent(35)) << endl
  //      << Form("SR35=%.2f",hyp_highpt_sr->GetBinContent(36)) << endl
  //      << Form("SR36=%.2f",hyp_highpt_sr->GetBinContent(37)) << endl
  //      << Form("SR37=%.2f",hyp_highpt_sr->GetBinContent(38)) << endl
  //      << Form("SR38=%.2f",hyp_highpt_sr->GetBinContent(39)) << endl;


}
