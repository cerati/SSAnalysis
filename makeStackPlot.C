{

  gROOT.Reset();
  TString plot = "hyp_highpt_sr";

  TCanvas c1;
  c1.SetLogy();

  TFile *ttbar_file = TFile::Open("ttbar_histos.root");
  TH1F* ttbar_h = (TH1F*) ttbar_file->Get(plot);
  ttbar_h->SetFillColor(kBlue);

  TFile *TTW_file = TFile::Open("TTWJets_histos.root");
  TH1F* TTW_h = (TH1F*) TTW_file->Get(plot);
  TTW_h->SetFillColor(kYellow);

  TFile *TTZ_file = TFile::Open("TTZJets_histos.root");
  TH1F* TTZ_h = (TH1F*) TTZ_file->Get(plot);
  TTZ_h->SetFillColor(kYellow+1);

  TFile *WZ_file = TFile::Open("WZJets_histos.root");
  TH1F* WZ_h = (TH1F*) WZ_file->Get(plot);
  WZ_h->SetFillColor(kYellow+2);

  TFile *WW_file = TFile::Open("WW_histos.root");
  TH1F* WW_h = (TH1F*) WW_file->Get(plot);
  if (WW_h) WW_h->SetFillColor(kYellow-1);

  THStack hs;
  hs.SetTitle(plot);
  hs.Add(ttbar_h);
  hs.Add(TTW_h);
  hs.Add(TTZ_h);
  hs.Add(WZ_h);
  if (WW_h) hs.Add(WW_h);
  hs.Draw("HIST");

  TFile *T1ttttG1500_file = TFile::Open("T1ttttG1500_histos.root");
  TH1F* T1ttttG1500_h = (TH1F*) T1ttttG1500_file->Get(plot);
  T1ttttG1500_h->SetLineColor(kRed);
  T1ttttG1500_h->SetLineStyle(1);
  T1ttttG1500_h->SetLineWidth(2);
  T1ttttG1500_h->Draw("SAME,HIST");

  TFile *T1ttttG1200_file = TFile::Open("T1ttttG1200_histos.root");
  TH1F* T1ttttG1200_h = (TH1F*) T1ttttG1200_file->Get(plot);
  T1ttttG1200_h->SetLineColor(kRed);
  T1ttttG1200_h->SetLineStyle(2);
  T1ttttG1200_h->SetLineWidth(2);
  T1ttttG1200_h->Draw("SAME,HIST");

  TFile *T5Full1500_file = TFile::Open("T5Full1500_histos.root");
  TH1F* T5Full1500_h = (TH1F*) T5Full1500_file->Get(plot);
  T5Full1500_h->SetLineColor(kMagenta);
  T5Full1500_h->SetLineStyle(1);
  T5Full1500_h->SetLineWidth(2);
  T5Full1500_h->Draw("SAME,HIST");

  TFile *T5Full1200_file = TFile::Open("T5Full1200_histos.root");
  TH1F* T5Full1200_h = (TH1F*) T5Full1200_file->Get(plot);
  T5Full1200_h->SetLineColor(kMagenta);
  T5Full1200_h->SetLineStyle(2);
  T5Full1200_h->SetLineWidth(2);
  T5Full1200_h->Draw("SAME,HIST");


}
