{

  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.2f");

  TString sample = "qcd";

  TFile* f = TFile::Open(sample+"_histos.root");

  TH2F* mud = (TH2F*) f->Get("fr_mu_den");
  TH2F* mun = (TH2F*) f->Get("fr_mu_num");
  TH2F* muf = (TH2F*) mud->Clone("fr_mu");
  muf->Reset();
  muf->SetTitle("muon fake rate");
  muf->GetXaxis()->SetTitle("p_{T} [GeV]");
  muf->GetYaxis()->SetTitle("|#eta|");
  for (int i=1;i<=muf->GetXaxis()->GetNbins();++i) {
    for (int j=1;j<=muf->GetYaxis()->GetNbins();++j) {
      if (mud->GetBinContent(i,j)>0) { 
	float fr = mun->GetBinContent(i,j)/mud->GetBinContent(i,j);
	muf->SetBinContent(i,j,fr);
	muf->SetBinError(i,j,sqrt( fr*(1-fr)/mud->GetBinContent(i,j)) );
      }
    }
  }
  TCanvas c1;
  muf->Draw("texte,colz");
  c1.RedrawAxis();

  TH2F* eld = (TH2F*) f->Get("fr_el_den");
  TH2F* eln = (TH2F*) f->Get("fr_el_num");
  TH2F* elf = (TH2F*) eld->Clone("fr_el");
  elf->Reset();
  elf->SetTitle("electron fake rate");
  elf->GetXaxis()->SetTitle("p_{T} [GeV]");
  elf->GetYaxis()->SetTitle("|#eta|");
  for (int i=1;i<=elf->GetXaxis()->GetNbins();++i) {
    for (int j=1;j<=elf->GetYaxis()->GetNbins();++j) {
      if (eld->GetBinContent(i,j)>0) { 
	float fr = eln->GetBinContent(i,j)/eld->GetBinContent(i,j);
	elf->SetBinContent(i,j,fr);
	elf->SetBinError(i,j,sqrt( fr*(1-fr)/eld->GetBinContent(i,j)) );
      }
    }
  }
  TCanvas c2;
  elf->Draw("texte,colz");
  c2.RedrawAxis();

  c1.SaveAs(sample+"_mu_fr.png");
  c2.SaveAs(sample+"_el_fr.png");

}
