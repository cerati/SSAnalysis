#ifndef SELECTIONS_H
#define SELECTIONS_H
#include "CMS2.h"
#include "TString.h"

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

bool isFakableElectron(unsigned int);
bool isFakableMuon(unsigned int);
bool isGoodElectron(unsigned int);
bool isGoodMuon(unsigned int);

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
float muRelIso03(unsigned int);
float muRelIso04(unsigned int);
float eleRelIso03(unsigned int);
int eleTightID(unsigned int);
int muTightID(unsigned int);
TString triggerName(TString);
bool passHLTTriggerPattern(const char*);

bool isFromW(Lep lep);
bool idIsCharm(int id);
bool idIsBeauty(int id);

#endif
