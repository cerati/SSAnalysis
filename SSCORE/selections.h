#ifndef SELECTIONS_H
#define SELECTIONS_H
#include "CMS2.h"
#include "TString.h"

const static float ptCutHigh = 25.;
const static float ptCutLow = 10.;

enum AnalysisBit { HighPt = 0, LowPt = 1, VeryLowPt = 2 };

//fixme: put WF and FSR in different categories
enum LeptonCategories { Prompt = 0, PromptWS = 1, PromptWF = 2, PromptFSR = 2, 
			FakeLightTrue = 3, FakeC = 4, FakeB = 5, FakeLightFake = 6, FakeHiPtGamma = 7, 
			FakeUnknown = 8, FakeLowPtGamma = 9, All9999 = 10,
			Other = 11, End = 12};

float muRelIso03(unsigned int);
float eleRelIso03(unsigned int);
float muRelIso03EA(unsigned int);
float eleRelIso03EA(unsigned int);
float muRelIso03DB(unsigned int);
float muRelIso04DB(unsigned int);
float eleRelIso03DB(unsigned int);

float muRelIsoTest(unsigned int, float dr, float deltaZCut=0.1);
float muRelIsoTestDB(unsigned int, float dr, float deltaZCut=0.1);

float elRelIsoTest(unsigned int, float dr, float deltaZCut=0.1);
float elRelIsoTestDB(unsigned int, float dr, float deltaZCut=0.1);

struct Lep {
  Lep(int pdgid, int idxx):pdgid_(pdgid),idx_(idxx){}
  int charge() {return -1*pdgid_/abs(pdgid_);}
  int pdgId() {return pdgid_;}
  int idx() {return idx_;}
  float pt() {return abs(pdgid_)==11 ? cms2.els_p4().at(idx_).pt() : cms2.mus_p4().at(idx_).pt();}
  float eta() {return abs(pdgid_)==11 ? cms2.els_p4().at(idx_).eta() : cms2.mus_p4().at(idx_).eta();}
  LorentzVector p4() {return abs(pdgid_)==11 ? cms2.els_p4().at(idx_) : cms2.mus_p4().at(idx_);}
  int mc3_id() {return abs(pdgid_)==11 ? cms2.els_mc3_id().at(idx_) : cms2.mus_mc3_id().at(idx_);}
  int mc3idx() {return abs(pdgid_)==11 ? cms2.els_mc3idx().at(idx_) : cms2.mus_mc3idx().at(idx_);}
  int mc3_motherid() {return abs(pdgid_)==11 ? cms2.els_mc3_motherid().at(idx_) : cms2.mus_mc3_motherid().at(idx_);}
  int mc3_motheridx() {return abs(pdgid_)==11 ? cms2.els_mc3_motheridx().at(idx_) : cms2.mus_mc3_motheridx().at(idx_);}
  int mc_id() { return abs(pdgid_)==11 ? cms2.els_mc_id().at(idx_) : cms2.mus_mc_id().at(idx_);}
  int mcidx() { return abs(pdgid_)==11 ? cms2.els_mcidx().at(idx_) : cms2.mus_mcidx().at(idx_);}
  int mc_motherid() {return abs(pdgid_)==11 ? cms2.els_mc_motherid().at(idx_) : cms2.mus_mc_motherid().at(idx_);}
  LorentzVector mc_p4() { return abs(pdgid_)==11 ? cms2.els_mc_p4().at(idx_) : cms2.mus_mc_p4().at(idx_);}
  float relIso03() { return abs(pdgid_)==11 ? eleRelIso03(idx_) : muRelIso03(idx_);}
  float dxyPV() { return abs(pdgid_)==11 ? cms2.els_dxyPV().at(idx_) : cms2.mus_dxyPV().at(idx_);}
  float dzPV() { return abs(pdgid_)==11 ? cms2.els_dzPV().at(idx_) : cms2.mus_dzPV().at(idx_);}
private:
  int pdgid_, idx_;
};

struct DilepHyp {
  DilepHyp(Lep lepone_, Lep leptwo_):
    leadlep(lepone_),trailep(leptwo_) {
    if (lepone_.pt()<leptwo_.pt()) {
      leadlep = leptwo_;
      trailep = lepone_;
    }
  }
  int charge() {return leadlep.charge()+trailep.charge();}
  LorentzVector p4() {return leadlep.p4()+trailep.p4();}
  Lep leadLep() {return leadlep;}
  Lep traiLep() {return trailep;}
private:
  Lep leadlep, trailep;
};

struct Jet {
  Jet(int idxx):idx_(idxx){}
  LorentzVector p4() {return cms2.pfjets_p4()[idx_]/**cms2.pfjets_corL1FastL2L3()[idx_]*/;}//fixme
  float pt() {return p4().pt();}
  float eta() {return p4().eta();}
  float phi() {return p4().phi();}
  float csv() {return cms2.pfjets_pfCombinedSecondaryVertexBJetTag()[idx_];}
  float csvivf() {return cms2.pfjets_combinedInclusiveSecondaryVertexV2BJetTag()[idx_];}
  bool isBtag() {return csvivf()>0.814;}
  int   mc3_id() {return cms2.pfjets_mc3_id()[idx_];}
  LorentzVector genjet_p4() {return cms2.pfjets_mc_p4()[idx_];}
  LorentzVector genps_p4() {return cms2.pfjets_mc_gp_p4()[idx_];}
  int idx() {return idx_;}
private:
  int idx_;
};

bool ptsort (int i,int j);
bool lepsort (Lep i,Lep j);
bool jetptsort (Jet i,Jet j);

float computePtRel(Lep& lep, vector<Jet> lepjets);

inline float deltaPhi( float phi1 , float phi2 ) {
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

inline float deltaR( LorentzVector lv1, LorentzVector lv2 ) {
  return sqrt( pow(deltaPhi(lv1.phi(),lv2.phi()),2) + pow(lv1.eta()-lv2.eta(),2) );
}

inline float mt(float pt1, float pt2, float dphi){
  return 2*sqrt(pt1*pt2)*fabs(sin(dphi/2));
}


//functions for veto, FO, Tight selection
bool isGoodLepton(int id, int idx);
bool isGoodLeptonNoIso(int id, int idx);
bool isDenominatorLepton(int id, int idx);
bool isGoodVetoElectron(unsigned int);
bool isGoodVetoMuon(unsigned int);
bool isFakableElectron(unsigned int);
bool isFakableMuon(unsigned int);
bool isGoodElectron(unsigned int);
bool isGoodMuon(unsigned int);

bool isGoodVetoElectronNoIso(unsigned int);
bool isGoodVetoMuonNoIso(unsigned int);
bool isFakableElectronNoIso(unsigned int);
bool isFakableMuonNoIso(unsigned int);
bool isGoodElectronNoIso(unsigned int);
bool isGoodMuonNoIso(unsigned int);

bool isGoodVertex(size_t ivtx);
int firstGoodVertex();

bool isLoosePFJet(unsigned int pfJetIdx);
bool isMediumPFJet(unsigned int pfJetIdx);
bool isTightPFJet(unsigned int pfJetIdx);
bool isVetoElectron(unsigned int);
bool isElectronFO(unsigned int);
int  isElectronFO_debug(unsigned int);
bool threeChargeAgree(unsigned int);
bool isMediumElectron(unsigned int);
bool isLooseMuon(unsigned int);
bool isTightMuon(unsigned int);
bool isMuonFO(unsigned int);

TString triggerName(TString);
bool passHLTTriggerPattern(const char*);

bool idIsCharm(int id);
bool idIsBeauty(int id);
bool isFromWZ(Lep lep);
bool isFromW(Lep lep);
bool isFromZ(Lep lep);
bool isFromB(Lep lep);
bool isFromC(Lep lep);
bool isFromLight(Lep lep);
bool isFromLightFake(Lep lep);

unsigned int analysisCategory(Lep lep1, Lep lep2);
void passesBaselineCuts(int njets, int nbtag, float met, float ht, unsigned int& analysisBitMask);
int baselineRegion(int nbtag);
void passesSignalRegionCuts(float ht, float met, unsigned int& analysisBitMask);
int signalRegion(int njets, int nbtag, float met, float ht);

bool makesExtraZ(int idx);
bool makesExtraGammaStar(int idx);
bool hypsFromFirstGoodVertex(size_t hypIdx, float dz_cut = 1.0);

struct metStruct{
  metStruct() : met(-999.), metphi(-999.), metx(-999.), mety(-999.), sumet(-999.)  {}
  float met;
  float metphi;
  float metx;
  float mety;
  float sumet;
};
metStruct trackerMET(float deltaZCut = 0.2 /*, const std::vector<LorentzVector>* jets = 0*/);

#endif
