#include "selections.h"
#include "CMS2.h"

using namespace tas;

bool isLoosePFJet(unsigned int pfJetIdx){

    float pfjet_chf_  = cms2.pfjets_chargedHadronE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    float pfjet_nhf_  = cms2.pfjets_neutralHadronE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    float pfjet_cef_  = cms2.pfjets_chargedEmE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    float pfjet_nef_  = cms2.pfjets_neutralEmE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    int   pfjet_cm_   = cms2.pfjets_chargedMultiplicity()[pfJetIdx];
    float pfjet_eta   = fabs(cms2.pfjets_p4()[pfJetIdx].eta());

    if (cms2.pfjets_pfcandIndicies().size() < 2)
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


    float pfjet_nhf_  = cms2.pfjets_neutralHadronE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    float pfjet_nef_  = cms2.pfjets_neutralEmE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();

    if (pfjet_nef_ >= 0.95)
        return false;
    if (pfjet_nhf_ >= 0.95)
        return false;

    if (!isLoosePFJet(pfJetIdx)) return false;

    return true;
}

bool isTightPFJet(unsigned int pfJetIdx){


    float pfjet_nhf_  = cms2.pfjets_neutralHadronE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();
    float pfjet_nef_  = cms2.pfjets_neutralEmE()[pfJetIdx] / cms2.pfjets_p4()[pfJetIdx].energy();

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

  if(fabs(cms2.els_etaSC().at(elIdx)) <= 1.479){
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.007) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.8) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.01) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.15) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.04) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.2) return false; //is this wrt the correct PV?
    if(eleRelIso03(elIdx) >= 0.15) return false; 
    return true;

  } else if((fabs(cms2.els_etaSC().at(elIdx)) > 1.479) && (fabs(cms2.els_etaSC().at(elIdx)) < 2.5)){
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.01) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.7) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.03) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.04) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.2) return false; //is this wrt the correct PV?
    if( eleRelIso03(elIdx) >= 0.15) return false; 
    return true;

  } else return false;

}

bool isLooseElectron(unsigned int elIdx){

  if(fabs(cms2.els_etaSC().at(elIdx)) <= 1.479){
    /*
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.007) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.15) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.01) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.12) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.02) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.2) return false; //is this wrt the correct PV?
    if( fabs( (1.0/cms2.els_ecalEnergy().at(elIdx)) - (cms2.els_eOverPIn().at(elIdx)/cms2.els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.15) return false; 
    if( cms2.els_conv_vtx_flag().at(elIdx) ) return false;
    if( cms2.els_lost_pixelhits().at(elIdx) > 1) return false;
    */
    //csa14 cuts
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.012) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.15) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.01) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.12) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.02) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.2) return false; //is this wrt the correct PV?
    if( fabs( (1.0/cms2.els_ecalEnergy().at(elIdx)) - (cms2.els_eOverPIn().at(elIdx)/cms2.els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.23) return false; 
    if( cms2.els_conv_vtx_flag().at(elIdx) ) return false;
    if( cms2.els_lost_pixelhits().at(elIdx) > 1) return false;
    //csa14 cuts
    return true;

  } else if((fabs(cms2.els_etaSC().at(elIdx)) > 1.479) && (fabs(cms2.els_etaSC().at(elIdx)) < 2.5)){
    /*
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.009) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.10) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.03) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.10) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.02) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.2) return false; //is this wrt the correct PV?
    if( fabs( (1.0/cms2.els_ecalEnergy().at(elIdx)) - (cms2.els_eOverPIn().at(elIdx)/cms2.els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.15 || (eleRelIso03(elIdx) >= 0.1 && cms2.els_p4().at(elIdx).pt() < 20) ) return false; 
    if( cms2.els_conv_vtx_flag().at(elIdx) ) return false;
    if( cms2.els_lost_pixelhits().at(elIdx) > 1) return false;
    */
    //csa14 cuts
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.014) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.10) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.033) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.12) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.02) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.2) return false; //is this wrt the correct PV?
    if( fabs( (1.0/cms2.els_ecalEnergy().at(elIdx)) - (cms2.els_eOverPIn().at(elIdx)/cms2.els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.25 /*|| (eleRelIso03(elIdx) >= 0.1 && cms2.els_p4().at(elIdx).pt() < 20)*/ ) return false; 
    if( cms2.els_conv_vtx_flag().at(elIdx) ) return false;
    if( cms2.els_lost_pixelhits().at(elIdx) > 1) return false;
    //csa14 cuts
    return true;

  } else return false;

}

bool isMediumElectron(unsigned int elIdx){

  if(fabs(cms2.els_etaSC().at(elIdx)) <= 1.479){
    /*
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.004) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.06) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.01) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.12) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.02) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.1) return false; //is this wrt the correct PV?
    if( fabs( (1.0/cms2.els_ecalEnergy().at(elIdx)) - (cms2.els_eOverPIn().at(elIdx)/cms2.els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.15) return false; 
    if( cms2.els_conv_vtx_flag().at(elIdx) ) return false;
    if( cms2.els_lost_pixelhits().at(elIdx) > 1) return false;
    */
    //csa14 cuts
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.009) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.06) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.01) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.12) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.02) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.1) return false; //is this wrt the correct PV?
    if( fabs( (1.0/cms2.els_ecalEnergy().at(elIdx)) - (cms2.els_eOverPIn().at(elIdx)/cms2.els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.23) return false; 
    if( cms2.els_conv_vtx_flag().at(elIdx) ) return false;
    if( cms2.els_lost_pixelhits().at(elIdx) > 1) return false;
    //csa14 cuts
    return true;

  } else if((fabs(cms2.els_etaSC().at(elIdx)) > 1.479) && (fabs(cms2.els_etaSC().at(elIdx)) < 2.5)){
    /*
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.007) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.03) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.03) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.10) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.02) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.1) return false; //is this wrt the correct PV?
    if( fabs( (1.0/cms2.els_ecalEnergy().at(elIdx)) - (cms2.els_eOverPIn().at(elIdx)/cms2.els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.15 || (eleRelIso03(elIdx) >= 0.1 && cms2.els_p4().at(elIdx).pt() < 20) ) return false; 
    if( cms2.els_conv_vtx_flag().at(elIdx) ) return false;
    if( cms2.els_lost_pixelhits().at(elIdx) > 1) return false;
    */
    //csa14 cuts
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.012) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.03) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.031) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.12) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.02) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.1) return false; //is this wrt the correct PV?
    if( fabs( (1.0/cms2.els_ecalEnergy().at(elIdx)) - (cms2.els_eOverPIn().at(elIdx)/cms2.els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.21 /*|| (eleRelIso03(elIdx) >= 0.1 && cms2.els_p4().at(elIdx).pt() < 20)*/ ) return false; 
    if( cms2.els_conv_vtx_flag().at(elIdx) ) return false;
    if( cms2.els_lost_pixelhits().at(elIdx) > 1) return false;
    //csa14 cuts
    return true;

  } else return false;

}

bool isTightElectron(unsigned int elIdx){

  if(fabs(cms2.els_etaSC().at(elIdx)) <= 1.479){
    /*
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.004) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.03) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.01) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.12) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.02) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.1) return false; //is this wrt the correct PV?
    if( fabs( (1.0/cms2.els_ecalEnergy().at(elIdx)) - (cms2.els_eOverPIn().at(elIdx)/cms2.els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.1) return false; 
    if( cms2.els_conv_vtx_flag().at(elIdx) ) return false;
    if( cms2.els_lost_pixelhits().at(elIdx) > 0) return false;
    */
    //csa14 cuts
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.009) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.03) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.01) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.12) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.02) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.1) return false; //is this wrt the correct PV?
    if( fabs( (1.0/cms2.els_ecalEnergy().at(elIdx)) - (cms2.els_eOverPIn().at(elIdx)/cms2.els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.18) return false; 
    if( cms2.els_conv_vtx_flag().at(elIdx) ) return false;
    if( cms2.els_lost_pixelhits().at(elIdx) > 0) return false;
    //csa14 cuts
    return true;

  } else if((fabs(cms2.els_etaSC().at(elIdx)) > 1.479) && (fabs(cms2.els_etaSC().at(elIdx)) < 2.5)){
    /*
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.005) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.02) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.03) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.10) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.02) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.1) return false; //is this wrt the correct PV?
    if( fabs( (1.0/cms2.els_ecalEnergy().at(elIdx)) - (cms2.els_eOverPIn().at(elIdx)/cms2.els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.1 || (eleRelIso03(elIdx) >= 0.07 && cms2.els_p4().at(elIdx).pt() < 20) ) return false; 
    if( cms2.els_conv_vtx_flag().at(elIdx) ) return false;
    if( cms2.els_lost_pixelhits().at(elIdx) > 0) return false;
    */
    //csa14 cuts
    if(fabs(cms2.els_dEtaIn().at(elIdx)) >= 0.010) return false; 
    if(fabs(cms2.els_dPhiIn().at(elIdx)) >= 0.02) return false; 
    if(cms2.els_sigmaIEtaIEta().at(elIdx) >= 0.031) return false; 
    if(cms2.els_hOverE().at(elIdx) >= 0.12) return false; 
    if(fabs(cms2.els_dxyPV().at(elIdx)) >= 0.02) return false; //is this wrt the correct PV?
    if(fabs(cms2.els_dzPV().at(elIdx)) >= 0.1) return false; //is this wrt the correct PV?
    if( fabs( (1.0/cms2.els_ecalEnergy().at(elIdx)) - (cms2.els_eOverPIn().at(elIdx)/cms2.els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( eleRelIso03(elIdx) >= 0.16 /*|| (eleRelIso03(elIdx) >= 0.07 && cms2.els_p4().at(elIdx).pt() < 20)*/ ) return false; 
    if( cms2.els_conv_vtx_flag().at(elIdx) ) return false;
    if( cms2.els_lost_pixelhits().at(elIdx) > 0) return false;
    //csa14 cuts
    return true;

  } else return false;

}


bool isLooseMuon(unsigned int muIdx){

  if(!cms2.mus_pid_PFMuon().at(muIdx)) return false; 
   
  bool isGlobal  = true;
  bool isTracker = true;
  if (((cms2.mus_type().at(muIdx)) & (1<<1)) == 0) isGlobal  = false;
  if (((cms2.mus_type().at(muIdx)) & (1<<2)) == 0) isTracker = false;

  if(!(isGlobal || isTracker)) return false;
  
  return true;

}

bool isTightMuon(unsigned int muIdx){
  if (!isLooseMuon(muIdx)) return false;
  if (cms2.mus_gfit_chi2().at(muIdx)/cms2.mus_gfit_ndof().at(muIdx) >= 10) return false; 
  if (cms2.mus_gfit_validSTAHits().at(muIdx) == 0)                         return false; 
  if (cms2.mus_numberOfMatchedStations().at(muIdx) < 2)                    return false;
  if (cms2.mus_dxyPV().at(muIdx) > 0.2)                                    return false;
  if (cms2.mus_dzPV().at(muIdx) > 0.5)                                     return false;
  if (cms2.mus_validPixelHits().at(muIdx) == 0)                            return false;
  if (cms2.mus_nlayers().at(muIdx) < 6)                                    return false;
  return true;

}


bool threeChargeAgree(unsigned int elIdx){

  if(cms2.els_charge().at(elIdx) != cms2.els_trk_charge().at(elIdx)) return false;
  if(cms2.els_charge().at(elIdx) != cms2.els_sccharge().at(elIdx))   return false;
  return true;

}

float muRelIso03(unsigned int muIdx){

  float chiso     = cms2.mus_isoR03_pf_ChargedHadronPt().at(muIdx);
  float nhiso     = cms2.mus_isoR03_pf_NeutralHadronEt().at(muIdx);
  float emiso     = cms2.mus_isoR03_pf_PhotonEt().at(muIdx);
  float deltaBeta = cms2.mus_isoR03_pf_PUPt().at(muIdx);
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(cms2.mus_p4().at(muIdx).pt());

}

float muRelIso04(unsigned int muIdx){

  float chiso     = cms2.mus_isoR04_pf_ChargedHadronPt().at(muIdx);
  float nhiso     = cms2.mus_isoR04_pf_NeutralHadronEt().at(muIdx);
  float emiso     = cms2.mus_isoR04_pf_PhotonEt().at(muIdx);
  float deltaBeta = cms2.mus_isoR04_pf_PUPt().at(muIdx);
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(cms2.mus_p4().at(muIdx).pt());

}
float eleRelIso03(unsigned int elIdx){

  float chiso     = cms2.els_pfChargedHadronIso().at(elIdx);
  float nhiso     = cms2.els_pfNeutralHadronIso().at(elIdx);
  float emiso     = cms2.els_pfPhotonIso().at(elIdx);
  float deltaBeta = cms2.els_pfPUIso().at(elIdx);
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(cms2.els_p4().at(elIdx).pt());
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
