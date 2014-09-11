#ifndef looper_h
#define looper_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>

#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"


class TChain;

class looper
{
 public:
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector; 

  looper() {
    p4_ = new LorentzVector();
  };
  ~looper() {
    delete babyFile_;
    delete babyTree_;
  };
  
  int ScanChain ( TChain* chain, TString prefix, TString postfix, bool isData, TString whatTest = "", int nEvents = -1);
  
  void MakeBabyNtuple (const char *);
  void InitBabyNtuple ();
  void FillBabyNtuple (){babyTree_->Fill();}
  void CloseBabyNtuple () { babyFile_->cd();babyTree_->Write();babyFile_->Close();}
  
  void CreateOutputFile(const char * name){
    outf = TFile::Open(name,"RECREATE");
  }
  void SaveHistos(){outf->cd();outf->Write();outf->Close();}
  void fillUnderOverFlow(TH1F *h1, float value, float weight = 1);
  
  void printEvent(  ostream& ostr = std::cout );
  
  float deltaPhi( float phi1 , float phi2 ) {
    float dphi = fabs( phi1 - phi2 );
    if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
    return dphi;
  }

  float deltaR( LorentzVector lv1, LorentzVector lv2 ) {
    return sqrt( pow(deltaPhi(lv1.phi(),lv2.phi()),2) + pow(lv1.eta()-lv2.eta(),2) );
  }

  float mt(float pt1, float pt2, float dphi){
    return 2*sqrt(pt1*pt2)*fabs(sin(dphi/2));
  }

  template<class T, class U> void makeFillHisto1D(const char* name,const char* title,
						  int nbins,U minx,U maxx,U value, 
						  float weight = 1.) {
    T* h = (T*) outf->Get(name);
    if (!h) {
      outf->cd();
      h = new T(name, title, nbins, minx, maxx);
      h->Sumw2();
    }
    h->Fill(std::max(minx,std::min(value,U(h->GetBinCenter(nbins)))),weight);
  }
  template<class T, class U> void makeFillHisto2D(const char* name,const char* title,
						  int nbinsx,U minx,U maxx,U valuex,
						  int nbinsy,U miny,U maxy,U valuey, 
						  float weight = 1.) {
    T* h = (T*) outf->Get(name);
    if (!h) {
      outf->cd();
      h = new T(name, title, nbinsx, minx, maxx, nbinsy, miny, maxy);
      h->Sumw2();
    }
    h->Fill(std::max(minx,std::min(valuex,U(h->GetXaxis()->GetBinCenter(nbinsx)))),std::max(miny,std::min(valuey,U(h->GetYaxis()->GetBinCenter(nbinsy)))),weight);
    float newvaluex = std::max(minx,std::min(valuex,U(h->GetXaxis()->GetBinCenter(nbinsx))));
    float newvaluey = std::max(miny,std::min(valuey,U(h->GetYaxis()->GetBinCenter(nbinsy))));
    if (TString(name)=="evt_mt_vs_pt" && valuey>80. && newvaluey<80.) std::cout << valuex << " " << newvaluex << " " << valuey << " " << newvaluey << std::endl;
  }
  
 private:
  
  bool makebaby;
  bool makehist;
  bool maketext;
  bool makeSSskim;
  bool makeQCDskim;
  bool makeQCDtest;
  bool makeDYtest;
  
  //ntuple, file
  TFile *babyFile_;
  TTree *babyTree_;
  TFile *outf;
  
  //baby vars
  int run_,ls_,evt_;
  float weight_;

  LorentzVector* p4_;

  //histos
  TH1F* h_dummy;
  
  ofstream ofile;
};

#endif
