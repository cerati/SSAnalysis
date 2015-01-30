#include "selections.h"
#include "../CORE/ElectronSelections.h"
#include "../CORE/MuonSelections.h"

using namespace tas;

bool ptsort (int i,int j) { return (genps_p4()[i].pt()>genps_p4()[j].pt()); }

bool lepsort (Lep i,Lep j) { 
  if ( abs(i.pdgId())==abs(j.pdgId()) ) return ( i.pt()>j.pt() ); //sort by pt if same flavor
  else return ( abs(i.pdgId())>abs(j.pdgId()) ); //prefer muons over electrons, but check that mu have pt>25//fixme, need to sync // && i.pt()>ptCutHigh
}

bool jetptsort (Jet i,Jet j) { return (i.pt()>j.pt()); }

std::vector<Lep> getBestSSLeps(std::vector<Lep> leps) {
  vector<Lep> hypleps;
  //select best hyp in case >3 good leps
  if (leps.size()>2) {
    vector<Lep> lepsp, lepsn;
    for (unsigned int gl=0;gl<leps.size();++gl) {
      if (leps[gl].pdgId()>0) lepsn.push_back(leps[gl]);
      else lepsp.push_back(leps[gl]);
    }
    //sort leps by muon and then by pt
    std::sort(lepsp.begin(),lepsp.end(),lepsort);
    std::sort(lepsn.begin(),lepsn.end(),lepsort);
    //take first SS hyp
    if (lepsp.size()<2) {
      hypleps.push_back(lepsn[0]);
      hypleps.push_back(lepsn[1]);
    } else if (lepsn.size()<2) {
      hypleps.push_back(lepsp[0]);
      hypleps.push_back(lepsp[1]);
    } else {
      if ( abs(lepsn[0].pdgId()+lepsn[1].pdgId())>abs(lepsp[0].pdgId()+lepsp[1].pdgId()) ) {
	hypleps.push_back(lepsn[0]);
	hypleps.push_back(lepsn[1]);
      } else if ( abs(lepsn[0].pdgId()+lepsn[1].pdgId())<abs(lepsp[0].pdgId()+lepsp[1].pdgId()) ) {
	hypleps.push_back(lepsp[0]);
	hypleps.push_back(lepsp[1]);
      } else if ( (lepsn[0].pt()+lepsn[1].pt())>(lepsp[0].pt()+lepsp[1].pt()) ) {
	hypleps.push_back(lepsn[0]);
	hypleps.push_back(lepsn[1]);
      } else {
	hypleps.push_back(lepsp[0]);
	hypleps.push_back(lepsp[1]);
      }
    }
  } else if (leps.size()==2) {
    hypleps.push_back(leps[0]);
    hypleps.push_back(leps[1]);	
  }      
  return hypleps;
}

float computeLD(DilepHyp hyp, vector<Jet> alljets, float met, float minmt) {
  //fixme: should variables be truncated?
  int njets25 = 0;
  float ht25 = 0;
  float htratio25 = 0;
  for (unsigned int j=0;j<alljets.size();j++) {
    float jetpt = alljets[j].pt();
    if (jetpt<25) continue;
    njets25++;
    ht25+=jetpt;
    if (fabs(alljets[j].eta())<1.2) htratio25+=jetpt;
  }
  htratio25/=ht25;
  float maxlepeta = std::max(fabs(hyp.leadLep().eta()),fabs(hyp.traiLep().eta()));
  if (hyp.leadLep().pt()>ptCutHigh) {
    return 0.147*met/100. + 0.178*ht25/1000. + 0.045*minmt/100. + 0.036*njets25 - 0.105*maxlepeta + 0.196*htratio25;
  } else {
    return 0.099*met/100. + 0.80*ht25/1000. + 0.004*njets25 - 0.046*maxlepeta + 0.094*htratio25 - 0.5;
  }
}


bool isGoodVertex(size_t ivtx) {
  if (vtxs_isFake()[ivtx]) return false;
  if (vtxs_ndof()[ivtx] <= 4.) return false;
  if (vtxs_position()[ivtx].Rho() > 2.0) return false;
  if (fabs(vtxs_position()[ivtx].Z()) > 24.0) return false;
  return true;
}

int firstGoodVertex () {
    for (unsigned int vidx = 0; vidx < vtxs_position().size(); vidx++) {
        if (isGoodVertex(vidx))
            return vidx;
    }
    return -1;
}

bool isElectronFO(unsigned int elIdx){//fixme
  //same as medium but removed dxy cut and iso<0.6
  if(fabs(els_etaSC().at(elIdx)) <= 1.479){
    if(fabs(els_dEtaIn().at(elIdx)) >= 0.004) return false;
    if(fabs(els_dPhiIn().at(elIdx)) >= 0.06) return false; 
    if(els_sigmaIEtaIEta_full5x5().at(elIdx) >= 0.01) return false; 
    if(els_hOverE().at(elIdx) >= 0.12) return false; 
    if( fabs( (1.0/els_ecalEnergy().at(elIdx)) - (els_eOverPIn().at(elIdx)/els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( els_conv_vtx_flag().at(elIdx) ) return false;
    if( els_exp_innerlayers().at(elIdx) > 0) return false;
    return true;
  } else if((fabs(els_etaSC().at(elIdx)) > 1.479) && (fabs(els_etaSC().at(elIdx)) < 2.5)){
    if(fabs(els_dEtaIn().at(elIdx)) >= 0.007) return false;
    if(fabs(els_dPhiIn().at(elIdx)) >= 0.03) return false; 
    if(els_sigmaIEtaIEta_full5x5().at(elIdx) >= 0.03) return false; 
    if(els_hOverE().at(elIdx) >= 0.1) return false; 
    if( fabs( (1.0/els_ecalEnergy().at(elIdx)) - (els_eOverPIn().at(elIdx)/els_ecalEnergy().at(elIdx)) ) >= 0.05) return false; // |1/E - 1/p|
    if( els_conv_vtx_flag().at(elIdx) ) return false;
    if( els_exp_innerlayers().at(elIdx) > 0) return false;
    return true;
  } else return false;
}

int isElectronFO_debug(unsigned int elIdx){//fixme
  if(fabs(els_etaSC().at(elIdx)) <= 1.479){
    if(fabs(els_dEtaIn().at(elIdx)) >= 0.004) return 1;
    if(fabs(els_dPhiIn().at(elIdx)) >= 0.06) return 2; 
    if(els_sigmaIEtaIEta_full5x5().at(elIdx) >= 0.01) return 3; 
    if(els_hOverE().at(elIdx) >= 0.12) return 4; 
    if(fabs(els_dxyPV().at(elIdx)) >= 0.05) return 5; //is this wrt the correct PV?
    if(fabs(els_dzPV().at(elIdx)) >= 0.1) return 5; //is this wrt the correct PV?
    if( fabs( (1.0/els_ecalEnergy().at(elIdx)) - (els_eOverPIn().at(elIdx)/els_ecalEnergy().at(elIdx)) ) >= 0.05) return 6; // |1/E - 1/p|
    if( eleRelIso03(elIdx, SS) >= 0.5 ) return 7; 
    if( els_conv_vtx_flag().at(elIdx) ) return 8;
    if( els_exp_innerlayers().at(elIdx) > 0) return 9;
    return 0;
  } else if((fabs(els_etaSC().at(elIdx)) > 1.479) && (fabs(els_etaSC().at(elIdx)) < 2.5)){
    if(fabs(els_dEtaIn().at(elIdx)) >= 0.007) return 1;
    if(fabs(els_dPhiIn().at(elIdx)) >= 0.03) return 2; 
    if(els_sigmaIEtaIEta_full5x5().at(elIdx) >= 0.03) return 3; 
    if(els_hOverE().at(elIdx) >= 0.1) return 4; 
    if(fabs(els_dxyPV().at(elIdx)) >= 0.05) return 5; //is this wrt the correct PV?
    if(fabs(els_dzPV().at(elIdx)) >= 0.1) return 5; //is this wrt the correct PV?
    if( fabs( (1.0/els_ecalEnergy().at(elIdx)) - (els_eOverPIn().at(elIdx)/els_ecalEnergy().at(elIdx)) ) >= 0.05) return 6; // |1/E - 1/p|
    if( eleRelIso03(elIdx, SS) >= 0.5 ) return 7; 
    if( els_conv_vtx_flag().at(elIdx) ) return 8;
    if( els_exp_innerlayers().at(elIdx) > 0) return 9;
    return 0;
  } else return -1;
}

float muRelIsoTest(unsigned int muIdx, float dr, float deltaZCut){
  float chiso     = 0.;
  float nhiso     = 0.;
  float emiso     = 0.;
  for (unsigned int i=0; i<pfcands_particleId().size(); ++i){
    if ( fabs(deltaR(pfcands_p4().at(i),mus_p4().at(muIdx)))>dr ) continue;  
    if ( fabs(pfcands_particleId().at(i))==211 && fabs(pfcands_dz().at(i)) < deltaZCut ) chiso+=pfcands_p4().at(i).pt();
    if ( fabs(pfcands_particleId().at(i))==130 ) nhiso+=pfcands_p4().at(i).pt();
    if ( fabs(pfcands_particleId().at(i))==22  ) emiso+=pfcands_p4().at(i).pt();
  }
  float absiso = chiso + std::max(float(0.0), nhiso + emiso);
  return absiso/(mus_p4().at(muIdx).pt());
}
float muRelIsoTestDB(unsigned int muIdx, float dr, float deltaZCut){
  float chiso     = 0.;
  float nhiso     = 0.;
  float emiso     = 0.;
  float deltaBeta = 0.;
  for (unsigned int i=0; i<pfcands_particleId().size(); ++i){
    if ( fabs(deltaR(pfcands_p4().at(i),mus_p4().at(muIdx)))>dr ) continue;  
    if ( fabs(pfcands_particleId().at(i))==211 ) {
      if (fabs(pfcands_dz().at(i)) < deltaZCut) chiso+=pfcands_p4().at(i).pt();
      else deltaBeta+=pfcands_p4().at(i).pt();
    }
    if ( fabs(pfcands_particleId().at(i))==130 ) nhiso+=pfcands_p4().at(i).pt();
    if ( fabs(pfcands_particleId().at(i))==22  ) emiso+=pfcands_p4().at(i).pt();
  }
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(mus_p4().at(muIdx).pt());
}

float elRelIsoTest(unsigned int elIdx, float dr, float deltaZCut){
  float chiso     = 0.;
  float nhiso     = 0.;
  float emiso     = 0.;
  for (unsigned int i=0; i<pfcands_particleId().size(); ++i){
    if ( fabs(deltaR(pfcands_p4().at(i),els_p4().at(elIdx)))>dr ) continue;  
    if ( fabs(pfcands_particleId().at(i))==211 && fabs(pfcands_dz().at(i)) < deltaZCut ) chiso+=pfcands_p4().at(i).pt();
    if ( fabs(pfcands_particleId().at(i))==130 ) nhiso+=pfcands_p4().at(i).pt();
    if ( fabs(pfcands_particleId().at(i))==22  ) emiso+=pfcands_p4().at(i).pt();
  }
  float absiso = chiso + std::max(float(0.0), nhiso + emiso);
  return absiso/(els_p4().at(elIdx).pt());
}
float elRelIsoTestDB(unsigned int elIdx, float dr, float deltaZCut){
  float chiso     = 0.;
  float nhiso     = 0.;
  float emiso     = 0.;
  float deltaBeta = 0.;
  for (unsigned int i=0; i<pfcands_particleId().size(); ++i){
    if ( fabs(deltaR(pfcands_p4().at(i),els_p4().at(elIdx)))>dr ) continue;  
    if ( fabs(pfcands_particleId().at(i))==211 ) {
      if (fabs(pfcands_dz().at(i)) < deltaZCut) chiso+=pfcands_p4().at(i).pt();
      else deltaBeta+=pfcands_p4().at(i).pt();
    }
    if ( fabs(pfcands_particleId().at(i))==130 ) nhiso+=pfcands_p4().at(i).pt();
    if ( fabs(pfcands_particleId().at(i))==22  ) emiso+=pfcands_p4().at(i).pt();
  }
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(els_p4().at(elIdx).pt());
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

unsigned int analysisCategory(Lep lep1, Lep lep2) {
  unsigned int result = 0;
  if (lep1.pt()>ptCutHigh && lep2.pt()>ptCutHigh) {
    result |= 1<<HighHigh;
  } else if (lep1.pt()>ptCutHigh && lep2.pt()>ptCutLow) {
    result |= 1<<HighLow;
  } else if (lep1.pt()>ptCutLow && lep2.pt()>ptCutLow) {
    result |= 1<<HighLow;
  }
  return result;
}

float computePtRel(Lep& lep, vector<Jet> lepjets) {

  if (abs(lep.pdgId())==13 && isGoodMuonNoIso(lep.idx())==0) return 0.;
  if (abs(lep.pdgId())==11 && isGoodElectronNoIso(lep.idx())==0) return 0.;
  if (lep.relIso03()<0.1) return 0.;//ok, this is inverted here
  int lepjetidx = -1;
  float mindr = 0.7;
  for (unsigned int j=0;j<lepjets.size();++j) {
    float dr = deltaR(lepjets[j].p4(),lep.p4());
    if (dr<mindr) {
      mindr = dr;
      lepjetidx = j;
    }
  } 
  if (lepjetidx>=0) {
    float sinA = fabs(lep.p4().x()*lepjets[lepjetidx].p4().y()-lep.p4().y()*lepjets[lepjetidx].p4().x())/(sqrt(lep.p4().x()*lep.p4().x()+lep.p4().y()*lep.p4().y())*sqrt(lepjets[lepjetidx].p4().x()*lepjets[lepjetidx].p4().x()+lepjets[lepjetidx].p4().y()*lepjets[lepjetidx].p4().y()));//fixme fabs? 
    return lep.pt()*sinA;
  } else return 0.;

}


//#include "Math/VectorUtil.h"
metStruct trackerMET(float deltaZCut/*, const std::vector<LorentzVector>* jets */) {

  if ( vtxs_sumpt().empty() ) return metStruct();
  double pX(0), pY(0);
  
  for (unsigned int i=0; i<pfcands_particleId().size(); ++i){
    if ( pfcands_charge().at(i)==0 ) continue;
    // if ( jets ){
    //   bool matched = false;
    //   for ( std::vector<LorentzVector>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet )
    // 	if ( fabs(ROOT::Math::VectorUtil::DeltaR(pfcands_p4().at(i),*jet))<0.5 ) matched=true;
    //   if (matched) continue;
    // }
    
    if ( fabs(pfcands_dz().at(i)) > deltaZCut) continue;
    
    pX -= pfcands_p4().at(i).px();
    pY -= pfcands_p4().at(i).py();
  }
  
  // if (jets){
  //   for ( std::vector<LorentzVector>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet ){
  //     pX -= jet->px();
  //     pY -= jet->py();
  //   }
  // }
  metStruct met;
  met.met     = sqrt(pX * pX + pY * pY);
  met.metphi  = atan2(pY, pX);
  met.metx = pX;
  met.mety = pY;
  return met;
}
