#ifndef SELECTIONS_H
#define SELECTIONS_H
#include "CMS2.h"
#include "TString.h"

enum AnalysisBit { HighPt = 0, LowPt = 1, VeryLowPt = 2 };

float muRelIso03(unsigned int);
float muRelIso04(unsigned int);
float eleRelIso03(unsigned int);

struct Lep {
  Lep(int pdgid, int idx):pdgid_(pdgid),idx_(idx){}
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
Jet(int idx):idx_(idx){}
  LorentzVector p4() {return cms2.pfjets_p4()[idx_]*cms2.pfjets_corL1FastL2L3()[idx_];}
  float pt() {return p4().pt();}
  float eta() {return p4().eta();}
  float csv() {return cms2.pfjets_combinedSecondaryVertexBJetTag()[idx_];}
private:
  int idx_;
};

bool isFakableElectron(unsigned int);
bool isFakableMuon(unsigned int);
bool isGoodElectron(unsigned int);
bool isGoodMuon(unsigned int);
bool isGoodVetoElectron(unsigned int);
bool isGoodVetoMuon(unsigned int);

bool isGoodVertex(size_t ivtx);
int firstGoodVertex();

bool isLoosePFJet(unsigned int pfJetIdx);
bool isMediumPFJet(unsigned int pfJetIdx);
bool isTightPFJet(unsigned int pfJetIdx);
bool isVetoElectron(unsigned int);
bool isLooseElectron(unsigned int);
bool isMediumElectron(unsigned int);
bool isElectronFO(unsigned int);
bool isTightElectron(unsigned int);
bool isLooseMuon(unsigned int);
bool isTightMuon(unsigned int);
bool isMuonFO(unsigned int);
bool threeChargeAgree(unsigned int);
int eleTightID(unsigned int);
int muTightID(unsigned int);
TString triggerName(TString);
bool passHLTTriggerPattern(const char*);

bool idIsCharm(int id);
bool idIsBeauty(int id);
bool isFromW(Lep lep);
bool isFromB(Lep lep);
bool isFromC(Lep lep);
bool isFromLight(Lep lep);
bool isFromLightFake(Lep lep);

unsigned int analysisCategory(Lep lep1, Lep lep2);
void passesBaselineCuts(int njets, int nbtag, float met, float ht, unsigned int& analysisBitMask);
int baselineRegion(int nbtag);
void passesSignalRegionCuts(float ht, float met, unsigned int& analysisBitMask);
int signalRegion(int njets, int nbtag, float met, float ht);

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
