{

  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.2f");

  TString sample = "dy";

  TFile* f = TFile::Open(sample+"_histos_effic.root");

  TH2F* mud = (TH2F*) f->Get("ef_mu_den");
  TH2F* mun = (TH2F*) f->Get("ef_mu_num");
  TH2F* muf = (TH2F*) mud->Clone("ef_mu");
  muf->Reset();
  muf->SetTitle("muon efficiency");
  muf->GetXaxis()->SetTitle("p_{T} [GeV]");
  muf->GetYaxis()->SetTitle("|#eta|");
  for (int i=1;i<=muf->GetXaxis()->GetNbins();++i) {
    for (int j=1;j<=muf->GetYaxis()->GetNbins();++j) {
      if (mud->GetBinContent(i,j)>0) { 
	float ef = mun->GetBinContent(i,j)/mud->GetBinContent(i,j);
	muf->SetBinContent(i,j,ef);
	muf->SetBinError(i,j,sqrt( ef*(1-ef)/mud->GetBinContent(i,j)) );
      }
    }
  }
  TCanvas c1;
  muf->GetZaxis()->SetRangeUser(0.,1.);
  muf->Draw("texte,colz");
  c1.RedrawAxis();

  TH2F* eld = (TH2F*) f->Get("ef_el_den");
  TH2F* eln = (TH2F*) f->Get("ef_el_num");
  TH2F* elf = (TH2F*) eld->Clone("ef_el");
  elf->Reset();
  elf->SetTitle("electron efficiency");
  elf->GetXaxis()->SetTitle("p_{T} [GeV]");
  elf->GetYaxis()->SetTitle("|#eta|");
  for (int i=1;i<=elf->GetXaxis()->GetNbins();++i) {
    for (int j=1;j<=elf->GetYaxis()->GetNbins();++j) {
      if (eld->GetBinContent(i,j)>0) { 
	float ef = eln->GetBinContent(i,j)/eld->GetBinContent(i,j);
	elf->SetBinContent(i,j,ef);
	elf->SetBinError(i,j,sqrt( ef*(1-ef)/eld->GetBinContent(i,j)) );
      }
    }
  }
  TCanvas c2;
  elf->GetZaxis()->SetRangeUser(0.,1.);
  elf->Draw("texte,colz");
  c2.RedrawAxis();

  c1.SaveAs(sample+"_mu_ef.png");
  c2.SaveAs(sample+"_el_ef.png");

}
