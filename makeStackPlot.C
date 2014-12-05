{

  gROOT.Reset();
  gStyle.SetOptStat(0);
  TString plot = "";
  plot = "hyp_ss_mass_mindRjj";
  plot = "hyp_mtmin";
  plot = "hyp_highpt_genbjets_pt40";
  plot = "hyp_highpt_genbjets_pt1";
  plot = "hyp_highpt_gennobjets_csv";
  plot = "hyp_highpthtmet_nbtag";
  plot = "hyp_ss_foFromWtrail_mu_relIso03";
  plot = "hyp_ss_foFromWlead_mu_relIso03";
  plot = "hyp_highpt_mtmin";
  plot = "hyp_highpthtmet_mtmin";
  plot = "hyp_highpt_genbjets_pt20";
  plot = "hyp_highpthtmetmt_nbtag";
  plot = "hyp_highpthtmetmt_pttrai";
  plot = "hyp_highpthtmetmt_sr";
  plot = "hyp_highpt_sr";
  bool norm = false;
  bool save = false;
  bool logy = false;

  TCanvas c1;
  if (logy) c1.SetLogy();

  TFile *ttbar_file = TFile::Open("ttbar_histos.root");
  TH1F* ttbar_h = (ttbar_file->GetListOfKeys()->Contains(plot) ? (TH1F*) ttbar_file->Get(plot) : 0);
  if (ttbar_h) {
    if (norm) ttbar_h->SetLineColor(kBlue);
    else ttbar_h->SetFillColor(kBlue);
  }

  TFile *TTW_file = TFile::Open("TTWJets_histos.root");
  TH1F* TTW_h = (TTW_file->GetListOfKeys()->Contains(plot) ? (TH1F*) TTW_file->Get(plot) : 0);
  if (TTW_h) {
    if (norm) TTW_h->SetLineColor(kOrange+1);
    else TTW_h->SetFillColor(kOrange+1);
  }

  TFile *TTZ_file = TFile::Open("TTZJets_histos.root");
  TH1F* TTZ_h = (TTZ_file->GetListOfKeys()->Contains(plot) ? (TH1F*) TTZ_file->Get(plot) : 0);
  if (TTZ_h) {
    if (norm) TTZ_h->SetLineColor(kYellow+1);
    else TTZ_h->SetFillColor(kYellow+1);
  }

  TFile *WZ_file = TFile::Open("WZJets_histos.root");
  TH1F* WZ_h = (WZ_file->GetListOfKeys()->Contains(plot) ? (TH1F*) WZ_file->Get(plot) : 0);
  if (WZ_h) {
    if (norm) WZ_h->SetLineColor(kGreen+2);
    else WZ_h->SetFillColor(kGreen+2);
  }

  TFile *WW_file = TFile::Open("WW_histos.root");
  TH1F* WW_h = (WW_file->GetListOfKeys()->Contains(plot) ? (TH1F*) WW_file->Get(plot) : 0);
  if (WW_h) {
    if (norm) WW_h->SetLineColor(kGreen-1);
    else WW_h->SetFillColor(kGreen-1);
  }

  if (ttbar_h && norm && ttbar_h->GetEntries()<=20) ttbar_h=0;
  if (TTW_h && norm && TTW_h->GetEntries()<=20) TTW_h=0;
  if (TTZ_h && norm && TTZ_h->GetEntries()<=20) TTZ_h=0;
  if (WZ_h  && norm && WZ_h->GetEntries()<=20) WZ_h=0;
  if (WW_h  && norm && WW_h->GetEntries()<=20) WW_h=0;

  float maxy = -1.;
  TH1* h0 = 0;
  if (norm) {
    if (ttbar_h) ttbar_h->Scale(1./ttbar_h->Integral());
    if (TTW_h)   TTW_h->Scale(1./TTW_h->Integral());
    if (TTZ_h)   TTZ_h->Scale(1./TTZ_h->Integral());
    if (WZ_h)    WZ_h->Scale(1./WZ_h->Integral());
    if (WW_h)    WW_h->Scale(1./WW_h->Integral());

    if (ttbar_h) ttbar_h->SetLineWidth(2);
    if (TTW_h)   TTW_h->SetLineWidth(2);
    if (TTZ_h)   TTZ_h->SetLineWidth(2);
    if (WZ_h)    WZ_h->SetLineWidth(2);
    if (WW_h)    WW_h->SetLineWidth(2);

    if (ttbar_h)    h0 = ttbar_h;
    else if (TTW_h) h0 = TTW_h;
    else if (TTZ_h) h0 = TTZ_h;
    else if (WZ_h)  h0 = WZ_h;
    else if (WW_h)  h0 = WW_h;

    if (ttbar_h && ttbar_h->GetMaximum()>maxy) maxy = ttbar_h->GetMaximum();
    if (TTW_h && TTW_h->GetMaximum()>maxy) maxy = TTW_h->GetMaximum();
    if (TTZ_h && TTZ_h->GetMaximum()>maxy) maxy = TTZ_h->GetMaximum();
    if (WZ_h  && WZ_h->GetMaximum()>maxy) maxy = WZ_h->GetMaximum();
    if (WW_h  && WW_h->GetMaximum()>maxy) maxy = WW_h->GetMaximum();

    if (h0) h0->Draw("HIST");

    if (ttbar_h) TTW_h->Draw("SAME,HIST");
    if (TTW_h)   TTW_h->Draw("SAME,HIST");
    if (TTZ_h)   TTZ_h->Draw("SAME,HIST");
    if (WZ_h)    WZ_h->Draw("SAME,HIST");
    if (WW_h)    WW_h->Draw("SAME,HIST");

  } else {
    THStack hs;
    hs.SetTitle(plot);
    if (ttbar_h) hs.Add(ttbar_h);
    if (TTW_h) hs.Add(TTW_h);
    if (TTZ_h) hs.Add(TTZ_h);
    if (WZ_h) hs.Add(WZ_h);
    if (WW_h) hs.Add(WW_h);
    hs.Draw("HIST");
  }

  TFile *T1ttttG1500_file = TFile::Open("T1ttttG1500_histos.root");
  TH1F* T1ttttG1500_h = (TH1F*) T1ttttG1500_file->Get(plot);
  T1ttttG1500_h->SetLineColor(kRed);
  T1ttttG1500_h->SetLineStyle(1);
  T1ttttG1500_h->SetLineWidth(2);
  if (norm) {
    T1ttttG1500_h->Scale(1./T1ttttG1500_h->Integral());
    if (h0==0) {
      h0=T1ttttG1500_h;
      h0->Draw("HIST");
    }
  }
  T1ttttG1500_h->Draw("SAME,HIST");
  if (T1ttttG1500_h && T1ttttG1500_h->GetMaximum()>maxy) maxy = T1ttttG1500_h->GetMaximum();

  TFile *T1ttttG1200_file = TFile::Open("T1ttttG1200_histos.root");
  TH1F* T1ttttG1200_h = (TH1F*) T1ttttG1200_file->Get(plot);
  T1ttttG1200_h->SetLineColor(kRed);
  T1ttttG1200_h->SetLineStyle(7);
  T1ttttG1200_h->SetLineWidth(2);
  if (norm) T1ttttG1200_h->Scale(1./T1ttttG1200_h->Integral());
  T1ttttG1200_h->Draw("SAME,HIST");
  if (T1ttttG1200_h && T1ttttG1200_h->GetMaximum()>maxy) maxy = T1ttttG1200_h->GetMaximum();

  TFile *T5Full1500_file = TFile::Open("T5Full1500_histos.root");
  TH1F* T5Full1500_h = (TH1F*) T5Full1500_file->Get(plot);
  if (T5Full1500_h && T5Full1500_h->GetEntries()>20) {
    T5Full1500_h->SetLineColor(kBlack);
    T5Full1500_h->SetLineStyle(1);
    T5Full1500_h->SetLineWidth(2);
    if (norm) T5Full1500_h->Scale(1./T5Full1500_h->Integral());
    T5Full1500_h->Draw("SAME,HIST");
    if (T5Full1500_h && T5Full1500_h->GetMaximum()>maxy) maxy = T5Full1500_h->GetMaximum();
  }

  TFile *T5Full1200_file = TFile::Open("T5Full1200_histos.root");
  TH1F* T5Full1200_h = (TH1F*) T5Full1200_file->Get(plot);
  if (T5Full1200_h && T5Full1200_h->GetEntries()>20) {
    T5Full1200_h->SetLineColor(kBlack);
    T5Full1200_h->SetLineStyle(7);
    T5Full1200_h->SetLineWidth(2);
    T5Full1200_h->Draw("SAME,HIST");
    if (norm) T5Full1200_h->Scale(1./T5Full1200_h->Integral());
    T5Full1200_h->Draw("SAME,HIST");
    if (T5Full1200_h && T5Full1200_h->GetMaximum()>maxy) maxy = T5Full1200_h->GetMaximum();
  }

  if (norm) {
    if (logy) h0->GetYaxis()->SetRangeUser(0.001,1.1*maxy);
    else h0->GetYaxis()->SetRangeUser(0.,1.1*maxy);
  }

  TLegend* leg  = new TLegend(0.73,0.54,0.89,0.89);
  leg->SetNColumns(1);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  if (ttbar_h) leg->AddEntry(ttbar_h,"ttbar",(norm ? "L" : "F"));
  if (TTW_h)   leg->AddEntry(TTW_h,"TTW",(norm ? "L" : "F"));
  if (TTZ_h)   leg->AddEntry(TTZ_h,"TTZ",(norm ? "L" : "F"));
  if (WZ_h)    leg->AddEntry(WZ_h,"WZ",(norm ? "L" : "F"));
  if (WW_h)    leg->AddEntry(WW_h,"WW",(norm ? "L" : "F"));
  if (T1ttttG1500_h) leg->AddEntry(T1ttttG1500_h,"T1_1500","L");
  if (T1ttttG1200_h) leg->AddEntry(T1ttttG1200_h,"T1_1200","L");
  if (T5Full1500_h)  leg->AddEntry(T5Full1500_h,"T5_1500","L");
  if (T5Full1200_h)  leg->AddEntry(T5Full1200_h,"T5_1200","L");
  leg->Draw();

  if (save) c1.SaveAs(plot+".png");

}
