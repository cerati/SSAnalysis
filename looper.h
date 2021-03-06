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
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

  looper() {
    p4_ = new LorentzVector();
  };
  ~looper() {
    delete babyFile_;
    delete babyTree_;
  };
  
  int ScanChain (TChain*, const char*, bool isData, int nEvents = -1);
  
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

  template<class T, class U> void makeFillHisto1D(const char* name,const char* title,int nbins,U minx,U maxx,U value) {
    T* h = (T*) outf->Get(name);
    if (!h) {
      outf->cd();
      h = new T(name, title, nbins, minx, maxx);
    }
    h->Fill(value);
  }
  
 private:
  
  bool makebaby;
  bool makehist;
  bool maketext;
  
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
