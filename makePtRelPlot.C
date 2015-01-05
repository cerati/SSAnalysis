{

  gROOT.Reset();
  gStyle.SetOptStat(0);

  TString dir = "results_PHYS14";

  TString file1 = dir+"/T1ttttG1500_histos.root";
  TString file2 = dir+"/T1ttttG1200_histos.root";
  TString file3 = dir+"/ttbar_histos.root";

  TString plot1 = "hyp_ss_foFromWlead_mu_ptrel";
  TString plot2 = "hyp_ss_foFromWlead_mu_ptrel";
  TString plot3 = "hyp_ss_foFakelead_mu_ptrel";
 
  bool norm = true;
  bool save = true;
  bool logy = false;

  TCanvas c1;c1.SetGridy();c1.SetGridx();

  TFile *_file1 = TFile::Open(file1);
  TH1F* h1 = (TH1F*) _file1->Get(plot1);
  h1->SetLineColor(kRed);
  h1->SetLineStyle(1);
  h1->SetLineWidth(2);
  if (norm) h1->Scale(1./h1->Integral());
  h1->GetYaxis()->SetRangeUser(0,1.0);
  h1->GetXaxis()->SetTitle("pTrel [GeV]");
  h1->GetYaxis()->SetTitle("fraction");
  h1->Draw("HIST");

  TFile *_file2 = TFile::Open(file2);
  TH1F* h2 = (TH1F*) _file2->Get(plot2);
  h2->SetLineColor(kRed);
  h2->SetLineStyle(2);
  h2->SetLineWidth(2);
  if (norm) h2->Scale(1./h2->Integral());
  h2->Draw("HIST,SAME");

  TFile *_file3 = TFile::Open(file3);
  TH1F* h3 = (TH1F*) _file3->Get(plot3);
  h3->SetLineColor(kBlue);
  h3->SetLineStyle(1);
  h3->SetLineWidth(3);
  if (norm) h3->Scale(1./h3->Integral());
  h3->Draw("HIST,SAME");

  TLegend* leg  = new TLegend(0.65,0.54,0.89,0.89);
  leg->SetNColumns(1);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->AddEntry(h3,"ttbar fakes",(norm ? "L" : "F"));
  leg->AddEntry(h1,"T1_1500 fromW","L");
  leg->AddEntry(h2,"T1_1200 fromW","L");
  leg->Draw();

  if (save) c1.SaveAs(dir+"/"+plot1+".png");

}
