{

  gROOT.Reset();
  gStyle.SetOptStat(0);

  TString file1 = "T1ttttG1200_histos.root";
  TString file2 = "T1ttttG1200_histos.root";

  TString plot1 = "hyp_highpt_gennobjets_csv";
  TString plot2 = "hyp_highpt_genbjets_csv";

  bool norm = true;
  bool save = false;
  bool logy = false;

  TCanvas c1;c1.SetGridy();c1.SetGridx();

  TFile *_file1 = TFile::Open(file1);
  TH1F* h1 = (TH1F*) _file1->Get(plot1);
  h1->SetLineColor(kRed);
  h1->SetLineStyle(2);
  h1->SetLineWidth(2);
  if (norm) h1->Scale(1./h1->Integral());
  TH1F* h1_integral = h1->Clone("integral");
  for (unsigned int bin=1;bin<=h1->GetNbinsX();++bin) {
    h1_integral->SetBinContent(bin,h1->Integral(bin,h1->GetNbinsX()+1));
  }
  h1_integral->Draw("HIST");

  TFile *_file2 = TFile::Open(file2);
  TH1F* h2 = (TH1F*) _file2->Get(plot2);
  h2->SetLineColor(kRed);
  h2->SetLineStyle(2);
  h2->SetLineWidth(2);
  if (norm) h2->Scale(1./h2->Integral());
  TH1F* h2_integral = h2->Clone("integral");
  for (unsigned int bin=1;bin<=h2->GetNbinsX();++bin) {
    h2_integral->SetBinContent(bin,h2->Integral(bin,h2->GetNbinsX()+1));
  }
  h2_integral->Draw("HIST,SAME");



}
