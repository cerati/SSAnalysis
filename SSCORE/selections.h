#ifndef SELECTIONS_H
#define SELECTIONS_H
#include "../CORE/CMS3.h"
#include "../CORE/SSSelections.h"
#include "TString.h"

//Enums
enum hyp_type_t { EE, MM, EM, UNASSIGNED }; 
//fixme: put WF and FSR in different categories
enum LeptonCategories { Prompt = 0, PromptWS = 1, PromptWF = 2, PromptFSR = 2, 
			FakeLightTrue = 3, FakeC = 4, FakeB = 5, FakeLightFake = 6, FakeHiPtGamma = 7, 
			FakeUnknown = 8, FakeLowPtGamma = 9, All9999 = 10,
			Other = 11, End = 12};

//Structs
struct hyp_result_t { int best_hyp; int hyp_class; };
struct particle_t { int id; LorentzVector p4; int idx; };

struct Lep {
  Lep(int pdgid, int idxx):pdgid_(pdgid),idx_(idxx){}
  int charge() {return -1*pdgid_/abs(pdgid_);}
  int pdgId() {return pdgid_;}
  int idx() {return idx_;}
  float pt() {return abs(pdgid_)==11 ? cms3.els_p4().at(idx_).pt() : cms3.mus_p4().at(idx_).pt();}
  float eta() {return abs(pdgid_)==11 ? cms3.els_p4().at(idx_).eta() : cms3.mus_p4().at(idx_).eta();}
  LorentzVector p4() {return abs(pdgid_)==11 ? cms3.els_p4().at(idx_) : cms3.mus_p4().at(idx_);}
  float relIso03() { return abs(pdgid_)==11 ? eleRelIso03(idx_, SS) : muRelIso03(idx_, SS);}
  float dxyPV() { return abs(pdgid_)==11 ? cms3.els_dxyPV().at(idx_) : cms3.mus_dxyPV().at(idx_);}
  float dzPV() { return abs(pdgid_)==11 ? cms3.els_dzPV().at(idx_) : cms3.mus_dzPV().at(idx_);}
  float d0Err() { return abs(pdgid_)==11 ? cms3.els_d0Err().at(idx_) : cms3.mus_d0Err().at(idx_);}
  float ip3d() { return abs(pdgid_)==11 ? cms3.els_ip3d().at(idx_) : cms3.mus_ip3d().at(idx_);}
  float ip3dErr() { return abs(pdgid_)==11 ? cms3.els_ip3derr().at(idx_) : cms3.mus_ip3derr().at(idx_);}
  int mc3_id() {return abs(pdgid_)==11 ? cms3.els_mc3_id().at(idx_) : cms3.mus_mc3_id().at(idx_);}
  int mc3idx() {return abs(pdgid_)==11 ? cms3.els_mc3idx().at(idx_) : cms3.mus_mc3idx().at(idx_);}
  int mc3_motherid() {return abs(pdgid_)==11 ? cms3.els_mc3_motherid().at(idx_) : cms3.mus_mc3_motherid().at(idx_);}
  int mc3_motheridx() {return abs(pdgid_)==11 ? cms3.els_mc3_motheridx().at(idx_) : cms3.mus_mc3_motheridx().at(idx_);}
  int mc_id() { return abs(pdgid_)==11 ? cms3.els_mc_id().at(idx_) : cms3.mus_mc_id().at(idx_);}
  int mcidx() { return abs(pdgid_)==11 ? cms3.els_mcidx().at(idx_) : cms3.mus_mcidx().at(idx_);}
  int mc_motherid() {return abs(pdgid_)==11 ? cms3.els_mc_motherid().at(idx_) : cms3.mus_mc_motherid().at(idx_);}
  LorentzVector mc_p4() { return abs(pdgid_)==11 ? cms3.els_mc_p4().at(idx_) : cms3.mus_mc_p4().at(idx_);}
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
  LorentzVector p4() {return cms3.pfjets_p4()[idx_]/**cms3.pfjets_corL1FastL2L3()[idx_]*/;}//fixme
  float pt() {return p4().pt();}
  float eta() {return p4().eta();}
  float phi() {return p4().phi();}
  float csv() {return cms3.pfjets_pfCombinedSecondaryVertexBJetTag()[idx_];}
  float csvivf() {return cms3.pfjets_combinedInclusiveSecondaryVertexV2BJetTag()[idx_];}
  bool isBtag() {return csvivf()>0.814;}
  int   mc3_id() {return cms3.pfjets_mc3_id()[idx_];}
  LorentzVector genjet_p4() {return cms3.pfjets_mc_p4()[idx_];}
  LorentzVector genps_p4() {return cms3.pfjets_mc_gp_p4()[idx_];}
  int idx() {return idx_;}
private:
  int idx_;
};


float muRelIsoTest(unsigned int, float dr, float deltaZCut=0.1);
float muRelIsoTestDB(unsigned int, float dr, float deltaZCut=0.1);

float elRelIsoTest(unsigned int, float dr, float deltaZCut=0.1);
float elRelIsoTestDB(unsigned int, float dr, float deltaZCut=0.1);

bool ptsort (int i,int j);
bool lepsort (Lep i,Lep j);
bool jetptsort (Jet i,Jet j);

std::vector<Lep> getBestSSLeps(std::vector<Lep> leps);

float computePtRel(Lep& lep, vector<Jet> lepjets);

float computeLD(DilepHyp hyp, vector<Jet> alljets, float met, float minmt);

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

bool isGoodVertex(size_t ivtx);
int firstGoodVertex();

int  isElectronFO_debug(unsigned int);

TString triggerName(TString);
bool passHLTTriggerPattern(const char*);

unsigned int analysisCategory(Lep lep1, Lep lep2);

struct metStruct{
  metStruct() : met(-999.), metphi(-999.), metx(-999.), mety(-999.), sumet(-999.)  {}
  float met;
  float metphi;
  float metx;
  float mety;
  float sumet;
};
metStruct trackerMET(float deltaZCut = 0.2 /*, const std::vector<LorentzVector>* jets = 0*/);

hyp_result_t chooseBestHyp(bool verbose=false);
int isGoodHyp(int iHyp, int analType = 2, bool verbose=false);
std::pair<particle_t, int> getThirdLepton(int hyp);
std::vector<particle_t> getGenPair(bool verbose=false);

template <typename T> int sgn(T val){
    return (T(0) < val) - (val < T(0));
}
double calculateMt(const LorentzVector p4, double met, double met_phi);
int lepMotherID(Lep lep);
float DileptonTriggerScaleFactor(int hyp_type, anal_type_t anal_type, LorentzVector trailing_p4);
float TagAndProbeScaleFactor(int id, float pt, float eta);
float DileptonTagAndProbeScaleFactor(const int hyp_idx);

#endif
