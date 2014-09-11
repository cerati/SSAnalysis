#include "selections.h"

using namespace tas;

bool isGoodVertex(size_t ivtx) {
  if (cms2.vtxs_isFake()[ivtx]) return false;
  if (cms2.vtxs_ndof()[ivtx] <= 4.) return false;
  if (cms2.vtxs_position()[ivtx].Rho() > 2.0) return false;
  if (fabs(cms2.vtxs_position()[ivtx].Z()) > 24.0) return false;
  return true;
}

int firstGoodVertex () {
    for (unsigned int vidx = 0; vidx < cms2.vtxs_position().size(); vidx++) {
        if (isGoodVertex(vidx))
            return vidx;
    }
    return -1;
}


bool isLoosePFJet(unsigned int pfJetIdx){

    float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / pfjets_p4()[pfJetIdx].energy();
    float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / pfjets_p4()[pfJetIdx].energy();
    float pfjet_cef_  = pfjets_chargedEmE()[pfJetIdx] / pfjets_p4()[pfJetIdx].energy();
    float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / pfjets_p4()[pfJetIdx].energy();
    int   pfjet_cm_   = pfjets_chargedMultiplicity()[pfJetIdx];
    float pfjet_eta   = fabs(pfjets_p4()[pfJetIdx].eta());

    if (pfjets_pfcandIndicies().size() < 2)
        return false;
    if (pfjet_nef_ >= 0.99)
        return false;
    if (pfjet_nhf_ >= 0.99)
        return false;

    if (pfjet_eta < 2.4){
      if (pfjet_cm_ < 1)
          return false;
      if (pfjet_chf_ < 1e-6)
          return false;
      if (pfjet_cef_ >= 0.99)
          return false;
    }

    return true;
}

bool isMediumPFJet(unsigned int pfJetIdx){


    float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / pfjets_p4()[pfJetIdx].energy();
    float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / pfjets_p4()[pfJetIdx].energy();

    if (pfjet_nef_ >= 0.95)
        return false;
    if (pfjet_nhf_ >= 0.95)
        return false;

    if (!isLoosePFJet(pfJetIdx)) return false;

    return true;
}

bool isTightPFJet(unsigned int pfJetIdx){


    float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / pfjets_p4()[pfJetIdx].energy();
    float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / pfjets_p4()[pfJetIdx].energy();

    if (pfjet_nef_ >= 0.90)
        return false;
    if (pfjet_nhf_ >= 0.90)
        return false;

    if (!isLoosePFJet(pfJetIdx)) return false;

    return true;
}


//2012 Electron IDs
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification#Electron_ID_Working_Points
bool isVetoElectron(unsigned int elIdx){

  return false;

}

bool isLooseElectron(unsigned int elIdx){

  return false;

}

bool isElectronFO(unsigned int elIdx){//fixme
  //same as medium but removed dxy cut and iso<0.6
  if(fabs(els_etaSC().at(elIdx)) <= 1.479){
    //csa14 cuts 50ns
    if(fabs(els_dEtaIn().at(elIdx)) >= 0.015) return false; 
    if(fabs(els_dPhiIn().at(elIdx)) >= 0.051) return false; 
    if(els_sigmaIEtaIEta().at(elIdx) >= 0.01) return false; 
    if(els_hOverE().at(elIdx) >= 0.10) return false; 
    //if(fabs(els_dxyPV().at(elIdx)) >= 0.012) return false; //is this wrt the correct PV?
    if(fabs(els_dzPV().at(elIdx)) >= 0.030) return false; //is this wrt the correct PV?
    if( fabs( (1.0/els_ecalEnergy().at(elIdx)) - (els_eOverPIn().at(elIdx)/els_ecalEnergy().at(elIdx)) ) >= 0.053) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.6) return false; //fixme add dBeta PU correction
    if( els_conv_vtx_flag().at(elIdx) ) return false;
    if( els_lost_pixelhits().at(elIdx) > 1) return false;
    //csa14 cuts 50ns
    return true;

  } else if((fabs(els_etaSC().at(elIdx)) > 1.479) && (fabs(els_etaSC().at(elIdx)) < 2.5)){
    //csa14 cuts 50ns
    if(fabs(els_dEtaIn().at(elIdx)) >= 0.023) return false; 
    if(fabs(els_dPhiIn().at(elIdx)) >= 0.056) return false; 
    if(els_sigmaIEtaIEta().at(elIdx) >= 0.030) return false; 
    if(els_hOverE().at(elIdx) >= 0.099) return false; 
    //if(fabs(els_dxyPV().at(elIdx)) >= 0.068) return false; //is this wrt the correct PV?
    if(fabs(els_dzPV().at(elIdx)) >= 0.78) return false; //is this wrt the correct PV?
    if( fabs( (1.0/els_ecalEnergy().at(elIdx)) - (els_eOverPIn().at(elIdx)/els_ecalEnergy().at(elIdx)) ) >= 0.11) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.6 ) return false; //fixme add dBeta PU correction
    if( els_conv_vtx_flag().at(elIdx) ) return false;
    if( els_lost_pixelhits().at(elIdx) > 1) return false;
    //csa14 cuts 50ns
    return true;

  } else return false;

}

bool isMediumElectron(unsigned int elIdx){

  if(fabs(els_etaSC().at(elIdx)) <= 1.479){
    //csa14 cuts 50ns
    if(fabs(els_dEtaIn().at(elIdx)) >= 0.015) return false; 
    if(fabs(els_dPhiIn().at(elIdx)) >= 0.051) return false; 
    if(els_sigmaIEtaIEta().at(elIdx) >= 0.01) return false; 
    if(els_hOverE().at(elIdx) >= 0.10) return false; 
    if(fabs(els_dxyPV().at(elIdx)) >= 0.012) return false; //is this wrt the correct PV?
    if(fabs(els_dzPV().at(elIdx)) >= 0.030) return false; //is this wrt the correct PV?
    if( fabs( (1.0/els_ecalEnergy().at(elIdx)) - (els_eOverPIn().at(elIdx)/els_ecalEnergy().at(elIdx)) ) >= 0.053) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.14) return false; //fixme add dBeta PU correction
    if( els_conv_vtx_flag().at(elIdx) ) return false;
    if( els_lost_pixelhits().at(elIdx) > 1) return false;
    //csa14 cuts 50ns
    return true;

  } else if((fabs(els_etaSC().at(elIdx)) > 1.479) && (fabs(els_etaSC().at(elIdx)) < 2.5)){
    //csa14 cuts 50ns
    if(fabs(els_dEtaIn().at(elIdx)) >= 0.023) return false; 
    if(fabs(els_dPhiIn().at(elIdx)) >= 0.056) return false; 
    if(els_sigmaIEtaIEta().at(elIdx) >= 0.030) return false; 
    if(els_hOverE().at(elIdx) >= 0.099) return false; 
    if(fabs(els_dxyPV().at(elIdx)) >= 0.068) return false; //is this wrt the correct PV?
    if(fabs(els_dzPV().at(elIdx)) >= 0.78) return false; //is this wrt the correct PV?
    if( fabs( (1.0/els_ecalEnergy().at(elIdx)) - (els_eOverPIn().at(elIdx)/els_ecalEnergy().at(elIdx)) ) >= 0.11) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.15 ) return false; //fixme add dBeta PU correction
    if( els_conv_vtx_flag().at(elIdx) ) return false;
    if( els_lost_pixelhits().at(elIdx) > 1) return false;
    //csa14 cuts 50ns
    return true;

  } else return false;

}

bool isTightElectron(unsigned int elIdx){

  if(fabs(els_etaSC().at(elIdx)) <= 1.479){
    //csa14 cuts 50ns
    if(fabs(els_dEtaIn().at(elIdx)) >= 0.012) return false; 
    if(fabs(els_dPhiIn().at(elIdx)) >= 0.024) return false; 
    if(els_sigmaIEtaIEta().at(elIdx) >= 0.01) return false; 
    if(els_hOverE().at(elIdx) >= 0.074) return false; 
    if(fabs(els_dxyPV().at(elIdx)) >= 0.0091) return false; //is this wrt the correct PV?
    if(fabs(els_dzPV().at(elIdx)) >= 0.017) return false; //is this wrt the correct PV?
    if( fabs( (1.0/els_ecalEnergy().at(elIdx)) - (els_eOverPIn().at(elIdx)/els_ecalEnergy().at(elIdx)) ) >= 0.026) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.10) return false; //fixme add dBeta PU correction
    if( els_conv_vtx_flag().at(elIdx) ) return false;
    if( els_lost_pixelhits().at(elIdx) > 1) return false;
    //csa14 cuts 50ns
    return true;

  } else if((fabs(els_etaSC().at(elIdx)) > 1.479) && (fabs(els_etaSC().at(elIdx)) < 2.5)){
    //csa14 cuts 50ns
    if(fabs(els_dEtaIn().at(elIdx)) >= 0.019) return false; 
    if(fabs(els_dPhiIn().at(elIdx)) >= 0.043) return false; 
    if(els_sigmaIEtaIEta().at(elIdx) >= 0.029) return false; 
    if(els_hOverE().at(elIdx) >= 0.08) return false; 
    if(fabs(els_dxyPV().at(elIdx)) >= 0.037) return false; //is this wrt the correct PV?
    if(fabs(els_dzPV().at(elIdx)) >= 0.065) return false; //is this wrt the correct PV?
    if( fabs( (1.0/els_ecalEnergy().at(elIdx)) - (els_eOverPIn().at(elIdx)/els_ecalEnergy().at(elIdx)) ) >= 0.076) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.14 ) return false; //fixme add dBeta PU correction
    if( els_conv_vtx_flag().at(elIdx) ) return false;
    if( els_lost_pixelhits().at(elIdx) > 1) return false;
    //csa14 cuts 50ns
    return true;

  } else return false;

}


bool isLooseMuon(unsigned int muIdx){

  if(!mus_pid_PFMuon().at(muIdx)) return false; 
   
  bool isGlobal  = true;
  bool isTracker = true;
  if (((mus_type().at(muIdx)) & (1<<1)) == 0) isGlobal  = false;
  if (((mus_type().at(muIdx)) & (1<<2)) == 0) isTracker = false;

  if(!(isGlobal || isTracker)) return false;
  
  return true;

}

bool isMuonFO(unsigned int muIdx){
  if (!isLooseMuon(muIdx)) return false;
  if (mus_gfit_chi2().at(muIdx)/mus_gfit_ndof().at(muIdx) >= 10)      return false; 
  if (mus_gfit_validSTAHits().at(muIdx) == 0)                         return false; 
  if (mus_numberOfMatchedStations().at(muIdx) < 2)                    return false;
  if (mus_validPixelHits().at(muIdx) == 0)                            return false;
  if (mus_nlayers().at(muIdx) < 6)                                    return false;
  if (mus_dxyPV().at(muIdx) > 0.2)                                    return false;
  if (mus_dzPV().at(muIdx) > 0.2)                                     return false;
  return true;
}

bool isTightMuon(unsigned int muIdx){
  if (!isMuonFO(muIdx)) return false;
  //fixme not applying MIP requirement in calo
  if (mus_dxyPV().at(muIdx) > 0.01)                                  return false;
  return true;
}


bool threeChargeAgree(unsigned int elIdx){

  if(els_charge().at(elIdx) != els_trk_charge().at(elIdx)) return false;
  if(els_charge().at(elIdx) != els_sccharge().at(elIdx))   return false;
  return true;

}

float muRelIso03(unsigned int muIdx){

  float chiso     = mus_isoR03_pf_ChargedHadronPt().at(muIdx);
  float nhiso     = mus_isoR03_pf_NeutralHadronEt().at(muIdx);
  float emiso     = mus_isoR03_pf_PhotonEt().at(muIdx);
  float deltaBeta = mus_isoR03_pf_PUPt().at(muIdx);
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(mus_p4().at(muIdx).pt());

}

float muRelIso04(unsigned int muIdx){

  float chiso     = mus_isoR04_pf_ChargedHadronPt().at(muIdx);
  float nhiso     = mus_isoR04_pf_NeutralHadronEt().at(muIdx);
  float emiso     = mus_isoR04_pf_PhotonEt().at(muIdx);
  float deltaBeta = mus_isoR04_pf_PUPt().at(muIdx);
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(mus_p4().at(muIdx).pt());

}
float eleRelIso03(unsigned int elIdx){

  float chiso     = els_pfChargedHadronIso().at(elIdx);
  float nhiso     = els_pfNeutralHadronIso().at(elIdx);
  float emiso     = els_pfPhotonIso().at(elIdx);
  float deltaBeta = els_pfPUIso().at(elIdx);
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(els_p4().at(elIdx).pt());
}
int eleTightID(unsigned int elIdx){
  if (isTightElectron(elIdx))  return 3;
  if (isMediumElectron(elIdx)) return 2;
  if (isLooseElectron(elIdx))  return 1;
  if (isVetoElectron(elIdx))   return 0;
  return -1;
}
int muTightID(unsigned int muIdx){
  if (isTightMuon(muIdx))  return 1;
  if (isLooseMuon(muIdx))   return 0;
  return -1;
}
//-------------------------------------------------------
// get exact trigger name corresponding to given pattern
//-------------------------------------------------------
TString triggerName(TString triggerPattern){

  bool    foundTrigger  = false;
  TString exact_hltname = "";

  for( unsigned int itrig = 0 ; itrig < hlt_trigNames().size() ; ++itrig ){
    if( TString( hlt_trigNames().at(itrig) ).Contains( triggerPattern ) ){
      foundTrigger  = true;
      exact_hltname = hlt_trigNames().at(itrig);
      break;
    }
  }

  if( !foundTrigger) return "TRIGGER_NOT_FOUND";

  return exact_hltname;

}
//---------------------------------------------
// Check if trigger passes
//---------------------------------------------
bool passHLTTriggerPattern(const char* arg){

  TString HLTTriggerPattern(arg);
  TString HLTTrigger = triggerName( HLTTriggerPattern );

  if( HLTTrigger.Contains("TRIGGER_NOT_FOUND")){
    return false;
  }
  return passHLTTrigger( HLTTrigger );
}


bool idIsCharm(int id) {
  id = abs(id);
  if (
      id == 4       ||
      id == 411     ||
      id == 421     ||
      id == 10411   ||
      id == 10421   ||
      id == 413     ||
      id == 423     ||
      id == 10413   ||
      id == 10423   ||
      id == 20413   ||
      id == 20423   ||
      id == 415     ||
      id == 425     ||
      id == 431     ||
      id == 10431   ||
      id == 433     ||
      id == 10433   ||
      id == 20433   ||
      id == 435     ||
      id == 441     ||
      id == 10441   ||
      id == 100441  ||
      id == 443     ||
      id == 10443   ||
      id == 20443   ||
      id == 100443  ||
      id == 30443   ||
      id == 9000443 ||
      id == 9010443 ||
      id == 9020443 ||
      id == 445     ||
      id == 9000445 ||
      id == 4122    ||
      id == 4222    ||
      id == 4212    ||
      id == 4112    ||
      id == 4224    ||
      id == 4214    ||
      id == 4114    ||
      id == 4232    ||
      id == 4132    ||
      id == 4322    ||
      id == 4312    ||
      id == 4324    ||
      id == 4314    ||
      id == 4332    ||
      id == 4334    ||
      id == 4412    ||
      id == 4422    ||
      id == 4414    ||
      id == 4424    ||
      id == 4432    ||
      id == 4434    ||
      id == 4444
      ) {
    return true;
  }
  else return false;
}

bool idIsBeauty(int id) {
  id = abs(id);
  if (
      id == 5       ||
      id == 511     ||
      id == 521     ||
      id == 10511   ||
      id == 10521   ||
      id == 513     ||
      id == 523     ||
      id == 10513   ||
      id == 10523   ||
      id == 20513   ||
      id == 20523   ||
      id == 515     ||
      id == 525     ||
      id == 531     ||
      id == 10531   ||
      id == 533     ||
      id == 10533   ||
      id == 20533   ||
      id == 535     ||
      id == 541     ||
      id == 10541   ||
      id == 543     ||
      id == 10543   ||
      id == 20543   ||
      id == 545     ||
      id == 551     ||
      id == 10551   ||
      id == 100551  ||
      id == 110551  ||
      id == 200551  ||
      id == 210551  ||
      id == 553     ||
      id == 10553   ||
      id == 20553   ||
      id == 30553   ||
      id == 100553  ||
      id == 110553  ||
      id == 120553  ||
      id == 130553  ||
      id == 200553  ||
      id == 210553  ||
      id == 220553  ||
      id == 300553  ||
      id == 9000553 ||
      id == 9010553 ||
      id == 555     ||
      id == 10555   ||
      id == 20555   ||
      id == 100555  ||
      id == 110555  ||
      id == 120555  ||
      id == 200555  ||
      id == 557     ||
      id == 100557  ||
      id == 5122    || 
      id == 5112    ||
      id == 5212    ||
      id == 5222    ||
      id == 5114    ||
      id == 5214    ||
      id == 5224    ||
      id == 5132    ||
      id == 5232    ||
      id == 5312    ||
      id == 5322    ||
      id == 5314    ||
      id == 5324    ||
      id == 5332    ||
      id == 5334    ||
      id == 5142    ||
      id == 5242    ||
      id == 5412    ||
      id == 5422    ||
      id == 5414    ||
      id == 5424    ||
      id == 5342    ||
      id == 5432    ||
      id == 5434    ||
      id == 5442    ||
      id == 5444    ||
      id == 5512    ||
      id == 5522    ||
      id == 5514    ||
      id == 5524    ||
      id == 5532    ||
      id == 5534    ||
      id == 5542    ||
      id == 5544    ||
      id == 5554 
      ) {
    return true;
  }
  else return false;
}


bool isFakableElectron(unsigned int elidx){
  if (els_p4().at(elidx).pt()<10.) return false;
  if (isElectronFO(elidx)==0) return false;
  if (threeChargeAgree(elidx)==0) return false;
  return true;
}

bool isGoodElectron(unsigned int elidx){
  if (isFakableElectron(elidx)==0) return false;
  if (isMediumElectron(elidx)==0) return false;
  if(fabs(els_dxyPV().at(elidx)) >= 0.01) return false;
  return true;
}

bool isFakableMuon(unsigned int muidx){
  if (mus_p4().at(muidx).pt()<5.) return false;
  if (isMuonFO(muidx)==0) return false;
  if (muRelIso03(muidx)>1.0 ) return false;
  return true;
}
bool isGoodMuon(unsigned int muidx){
  if (isFakableMuon(muidx)==0) return false;
  if (isTightMuon(muidx)==0) return false;
  if (muRelIso03(muidx)>0.1 ) return false;
  return true;
}

bool isFromW(Lep lep) {
  //status 1 leptons with mother W or with mother tau whose mother is W
  if ( (abs(lep.mc_id())==11 || abs(lep.mc_id())==13) && (abs(lep.mc_motherid())==24 || (abs(lep.mc_motherid())==15 && abs(genps_id_mother()[lep.mc3_motheridx()])==24) ) ) return true;
  //status 1 photons whose mother is a lepton, matching status 3 leptons whose mother is W
  if ( abs(lep.mc_id())==22 && abs(lep.mc_motherid())==abs(lep.pdgId()) && abs(lep.mc3_id())==abs(lep.pdgId()) && abs(lep.mc3_motherid())==24 ) return true;
  //everything else
  return false;
}
bool isFromB(Lep lep) {
  //true lepton from b quark
  return (abs(lep.mc_id())==11 || abs(lep.mc_id())==13) && idIsBeauty(lep.mc_motherid());
}
bool isFromC(Lep lep) {
  //true lepton from c quark
  return (abs(lep.mc_id())==11 || abs(lep.mc_id())==13) && idIsCharm(lep.mc_motherid());
}
bool isFromLight(Lep lep) {
  //true lepton from light quark
  return (abs(lep.mc_id())==11 || abs(lep.mc_id())==13) && ( (abs(lep.mc_motherid())>200 && abs(lep.mc_motherid())<400) || (abs(lep.mc_motherid())>0 && abs(lep.mc_motherid())<4) );
}
bool isFromLightFake(Lep lep) {
  //fake lepton from light quark: match light hadron or in any case does not match a lepton and mother is light hadron
  return (abs(lep.mc_id())>200 && abs(lep.mc_id())<400) || 
    ( abs(lep.mc_id())!=11 && abs(lep.mc_id())!=13 && ( (abs(lep.mc_motherid())>200 && abs(lep.mc_motherid())<400) || (abs(lep.mc_motherid())>0 && abs(lep.mc_motherid())<4) ) );
}

unsigned int analysisCategory(Lep lep1, Lep lep2) {
  unsigned int result = 0;
  if (lep1.pt()>20 && lep2.pt()>20) {
    result |= 1<<HighPt;
    result |= 1<<LowPt;
    result |= 1<<VeryLowPt;
  } else if (lep1.pt()>10 && lep2.pt()>10) {
    result |= 1<<LowPt;
    result |= 1<<VeryLowPt;
  } else if (lep1.pt()>5 && lep2.pt()>5) {
    //electrons must be at leat 10 GeV
    if ( (abs(lep1.pdgId())==13 || lep1.pt()>10) && (abs(lep2.pdgId())==13 || lep2.pt()>10) ) result |= 1<<VeryLowPt;
  }
  return result;
}

void passesBaselineCuts(int njets, int nbtag, float met, float ht, unsigned int& analysisBitMask) {
  if (analysisBitMask & 1<<HighPt) {
    if (!(ht>80 && njets>=2 && (met>30 || ht>500)))  analysisBitMask &= ~(1<<HighPt);
  } else {
    if (!(ht>250 && njets>=2 && (met>30 || ht>500))) {
      analysisBitMask &= ~(1<<LowPt);
      analysisBitMask &= ~(1<<VeryLowPt);
    }
  }
}

int baselineRegion(int nbtag) {
  if (nbtag==0) return 0;
  else if (nbtag==1) return 10;
  else return 20;
}

void passesSignalRegionCuts(float ht, unsigned int& analysisBitMask) {
  if (analysisBitMask & 1<<HighPt) {
    if (ht<200)  analysisBitMask &= ~(1<<HighPt);
  } else {
    if (ht<250) {
      analysisBitMask &= ~(1<<LowPt);
      analysisBitMask &= ~(1<<VeryLowPt);
    }
  }
}

//this assumes that the event has already passed all selections (including min ht)
int signalRegion(int njets, int nbtag, float met, float ht) {
  int result = 1;
  if (nbtag==1) result+=10;
  else if (nbtag>=2) result+=20;
  if (met>120) result+=4;
  if (njets>=4) result+=2;
  if (ht>400) result+=1;
  return result;
}
