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

double calculateMt(const LorentzVector p4, double met, double met_phi){
  float phi1 = p4.Phi();
  float phi2 = met_phi;
  float Et1  = p4.Et();
  float Et2  = met;
  return sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
}

int lepMotherID(Lep lep){
  if (tas::evt_isRealData()) return 1;
  else if (isFromZ(lep.pdgId(),lep.idx()) || isFromW(lep.pdgId(),lep.idx())){
    if (sgn(lep.pdgId()) == sgn(lep.mc_id())) return 1;
    else return 2;
  }
  else if (isFromB(lep.pdgId(),lep.idx())) return -1;
  else if (isFromC(lep.pdgId(),lep.idx())) return -2;
  return 0;
}

float DileptonTriggerScaleFactor(int hyp_type, anal_type_t anal_type, LorentzVector trailing_p4){

  //Local variables
  const bool is_barrel = (fabs(trailing_p4.eta()) < 1.0);
  const float lep_pt   = trailing_p4.pt();
  
  // values supplied by ETH (LHCP value)
  if (anal_type == HighHigh){
      if (hyp_type == 3) return (lep_pt < 30 ? 0.92 : 0.96); 
      else if (hyp_type == 0) return (is_barrel   ? 0.90 : 0.81); 
      else return 0.93;                        
  }
  if (anal_type == HighLow){
      if (hyp_type == 0) return 0.93;                      
      else if (hyp_type == 3) return (is_barrel ? 0.94 : 0.90); 
      else return 0.93;                      
  }
  if (anal_type == LowLow){
      // not measured yet --> using vlow_pt
      if (hyp_type == 0) return 0.93;                      
      else if (hyp_type == 3) return (is_barrel ? 0.94 : 0.90); 
      else return 0.93;                     
      return 0.0;
  }

  return 0.0;
}

float TagAndProbeScaleFactor(int id, float pt, float eta){
  // haven't measured very-low pt
  if (pt < 10) return 1.0;
  const double aeta = fabs(eta);

  if (abs(id)==11){
    if (1.4442 < aeta && aeta < 1.566) return 0;
    if (10 < pt && pt <= 15){
      if (0.00   < aeta && aeta <= 0.80  ) return 0.834;
      if (0.80   < aeta && aeta <= 1.4442) return 0.973;
      if (1.566  < aeta && aeta <= 2.00  ) return 0.954;
      if (2.00   < aeta)                   return 1.119;
    }
    else if (15 < pt && pt <= 20){
      if (0.00   < aeta && aeta <= 0.80  ) return 0.918;
      if (0.80   < aeta && aeta <= 1.4442) return 0.906;
      if (1.566  < aeta && aeta <= 2.00  ) return 0.909;
      if (2.00   < aeta)                   return 0.944;
    }
    else if (20 < pt && pt <= 30){
      if (0.00   < aeta && aeta <= 0.80  ) return 0.954;
      if (0.80   < aeta && aeta <= 1.4442) return 0.923;
      if (1.566  < aeta && aeta <= 2.00  ) return 0.921;
      if (2.00   < aeta)                   return 0.993;
    }
    else if (30 < pt && pt <= 40){
      if (0.00   < aeta && aeta <= 0.80  ) return 0.960;
      if (0.80   < aeta && aeta <= 1.4442) return 0.935;
      if (1.566  < aeta && aeta <= 2.00  ) return 0.924;
      if (2.00   < aeta)                   return 0.959;
    }
    else if (40 < pt && pt <= 50){
      if (0.00   < aeta && aeta <= 0.80  ) return 0.972;
      if (0.80   < aeta && aeta <= 1.4442) return 0.955;
      if (1.566  < aeta && aeta <= 2.00  ) return 0.950;
      if (2.00   < aeta)                   return 0.968;
    }
    else if (pt > 50){
      if (0.00   < aeta && aeta <= 0.80  ) return 0.969;
      if (0.80   < aeta && aeta <= 1.4442) return 0.956;
      if (1.566  < aeta && aeta <= 2.00  ) return 0.995;
      if (2.00   < aeta)                   return 0.969;
    }
    return -999999.0;
  }//electrons

  if (abs(id)==13){
    if (10 < pt && pt <= 15){
      if (0.00  < aeta && aeta <= 1.20) return 0.973;
      if (1.20  < aeta)                 return 0.954;
    }
    else if (15 < pt && pt <= 20){
      if (0.00  < aeta && aeta <= 1.20) return 0.957;
      if (1.20  < aeta)                 return 0.971;
    }
    else if (20 < pt && pt <= 30){
      if (0.00  < aeta && aeta <= 1.20) return 0.964;
      if (1.20  < aeta)                 return 0.981;
    }
    else if (30 < pt && pt <= 40){
      if (0.00  < aeta && aeta <= 1.20) return 0.971;
      if (1.20  < aeta)                 return 0.978;
    }
    else if (40 < pt && pt <= 50){
      if (0.00  < aeta && aeta <= 1.20) return 0.978;
      if (1.20  < aeta)                 return 0.984;
    }
    else if (pt > 50){
      if (0.00  < aeta && aeta <= 1.20) return 0.974;
      if (1.20  < aeta)                 return 0.977;
    }
    return -999999.0;
  }//muons

  // if we get here, return bogus value
  return -999999.0;
}

float DileptonTagAndProbeScaleFactor(const int hyp_idx){
  float eta_1 = (abs(tas::hyp_ll_id().at(hyp_idx)) == 13) ? tas::hyp_ll_p4().at(hyp_idx).eta() : tas::els_etaSC().at(tas::hyp_ll_index().at(hyp_idx));
  float eta_2 = (abs(tas::hyp_lt_id().at(hyp_idx)) == 13) ? tas::hyp_lt_p4().at(hyp_idx).eta() : tas::els_etaSC().at(tas::hyp_lt_index().at(hyp_idx));
  const float sf1 = TagAndProbeScaleFactor(tas::hyp_ll_id().at(hyp_idx), tas::hyp_ll_p4().at(hyp_idx).pt(), eta_1);
  const float sf2 = TagAndProbeScaleFactor(tas::hyp_lt_id().at(hyp_idx), tas::hyp_lt_p4().at(hyp_idx).pt(), eta_2);
  return sf1*sf2;
}

int isGoodHyp(int iHyp, int analType, bool verbose){

  //Bunch o' variables
  //bool isData = tas::evt_isRealData();
  float pt_ll = tas::hyp_ll_p4().at(iHyp).pt(); 
  float pt_lt = tas::hyp_lt_p4().at(iHyp).pt(); 
  float eta_ll = tas::hyp_ll_p4().at(iHyp).eta();
  float eta_lt = tas::hyp_lt_p4().at(iHyp).eta();
  int idx_ll = tas::hyp_ll_index().at(iHyp);
  int idx_lt = tas::hyp_lt_index().at(iHyp);
  int id_ll = tas::hyp_ll_id().at(iHyp);
  int id_lt = tas::hyp_lt_id().at(iHyp);
  bool isss = false;
  if (sgn(id_ll) == sgn(id_lt)) isss = true;  
  bool passed_id_numer_ll = isGoodLepton(id_ll, idx_ll);
  bool passed_id_numer_lt = isGoodLepton(id_lt, idx_lt);
  bool passed_id_denom_ll = isDenominatorLepton(id_ll, idx_ll);
  bool passed_id_denom_lt = isDenominatorLepton(id_lt, idx_lt);
  bool extraZ = makesExtraZ(iHyp);
  bool extraGammaStar = makesExtraGammaStar(iHyp);

  //Verbose info:
  if (verbose && pt_ll > 10 && pt_lt > 10){
    cout << "hyp " << iHyp << "leptons: " << id_ll << " " << pt_ll << " " << id_lt << " " << pt_lt << endl;
    cout << "   isss: " << isss << endl;
    cout << "   extraZ: " << extraZ << endl;
    cout << "   extraG: " << extraGammaStar << endl;
    cout << "   invt mass: " << (tas::hyp_ll_p4().at(iHyp) + tas::hyp_lt_p4().at(iHyp)).M() << endl;
    cout << "   passes eta: " << (fabs(eta_ll) < 2.4 && fabs(eta_lt) < 2.4) << " etas are " << eta_ll << " and " << eta_lt << endl;
    cout << "   passes hypsFromFirstGoodVertex: " << hypsFromFirstGoodVertex(iHyp) << endl;
    //cout << "   passes triggers: " << passesTriggerVeryLowPt(tas::hyp_type().at(iHyp)) << endl; 
    cout << "   lepton with pT " << pt_ll << " passes id: " << passed_id_numer_ll << endl;
    cout << "   lepton with pT " << pt_lt << " passes id: " << passed_id_numer_lt << endl;
    if (abs(id_ll) == 11) cout << "   lepton with pT " << pt_ll << " passes 3chg: " << threeChargeAgree(idx_ll) << endl;
    if (abs(id_lt) == 11) cout << "   lepton with pT " << pt_lt << " passes 3chg: " << threeChargeAgree(idx_lt) << endl;
  }

  //Cuts:
  if (analType == 0 && pt_ll < 20) return 0;
  else if (analType == 1 && pt_ll < 10) return 0;
  else if (analType == 2 && abs(id_ll) == 11 && pt_ll < 7) return 0;
  else if (analType == 2 && abs(id_ll) == 13 && pt_ll < 5) return 0;
  if (analType == 0 && pt_lt < 20) return 0;
  else if (analType == 1 && pt_lt < 10) return 0;
  else if (analType == 2 && abs(id_lt) == 11 && pt_lt < 7) return 0;
  else if (analType == 2 && abs(id_lt) == 13 && pt_lt < 5) return 0;
  if (fabs(eta_ll) > 2.4) return 0;
  if (fabs(eta_lt) > 2.4) return 0;
  if (extraZ) return 0;
  if (extraGammaStar) return 0;
  if ((tas::hyp_ll_p4().at(iHyp) + tas::hyp_lt_p4().at(iHyp)).M() < 8) return 0; 
  if (!hypsFromFirstGoodVertex(iHyp)) return 0;
  //if (isData == true && passesTriggerVeryLowPt(tas::hyp_type().at(iHyp)) == 0) return 0;

  //Results
  if (passed_id_numer_ll == 0 && passed_id_denom_ll == 0) return 0;      // 0 if ll fails
  if (passed_id_numer_lt == 0 && passed_id_denom_lt == 0) return 0; // 0 if lt fails
  else if (passed_id_numer_lt == 1 && passed_id_numer_ll == 1 && isss == 1) return 3;  // 3 if both high pass, SS
  else if (passed_id_numer_lt == 1 && passed_id_numer_ll == 1 && isss == 0) return 4;  // 4 if both high pass, OS
  else if (passed_id_numer_lt == 0 && passed_id_numer_ll == 0 && passed_id_denom_lt == 1 && passed_id_denom_ll == 1 && isss == true) return 1; // 1 if both low pass
  else if (isss == true) return 2; //2 for lowpass/highpass
  else return 0; //non-highpass OS
}

hyp_result_t chooseBestHyp(bool verbose){

  //List of good hyps
  vector <int> good_hyps_ss; 
  vector <int> good_hyps_sf; 
  vector <int> good_hyps_df; 
  vector <int> good_hyps_os; 
  for (unsigned int i = 0; i < tas::hyp_type().size(); i++){
    int good_hyp_result = isGoodHyp(i, 2,verbose);
    if (good_hyp_result == 3) good_hyps_ss.push_back(i); 
    if (good_hyp_result == 2) good_hyps_sf.push_back(i); 
    else if (good_hyp_result == 1) good_hyps_df.push_back(i); 
    else if (good_hyp_result == 4) good_hyps_os.push_back(i); 
  }

  //hyp_class_ to track SS(3), SF(2), DF(1), OS(4), or none(0)
  int hyp_class_;

  //Load good hyps in, SS then SF then DF then OS
  vector <int> good_hyps;
  if (good_hyps_ss.size() != 0){
     good_hyps = good_hyps_ss;
     hyp_class_ = 3;
  }
  else if (good_hyps_sf.size() != 0){
    good_hyps = good_hyps_sf;
    hyp_class_ = 2;
  }
  else if (good_hyps_df.size() != 0){
     good_hyps = good_hyps_df;
     hyp_class_ = 1;
  }
  else if (good_hyps_os.size() != 0){
    good_hyps = good_hyps_os;
    hyp_class_ = 4;
  }
  else hyp_class_ = 0; 

  //If no hyps or one hyps, know what to do
  int best_hyp_ = -1;
  if (good_hyps.size() == 1) best_hyp_ = good_hyps.at(0);

  //Otherwise, pick ones with more muons, then highest pT  
  if (good_hyps.size() > 1){
    best_hyp_ = good_hyps.at(0); 
    for (unsigned int i = 1; i < good_hyps.size(); i++){
      int hyp = good_hyps.at(i);
      if (tas::hyp_type().at(hyp) < tas::hyp_type().at(best_hyp_)) best_hyp_ = hyp;
      else if (tas::hyp_type().at(hyp) == tas::hyp_type().at(best_hyp_) && (tas::hyp_ll_p4().at(hyp)+tas::hyp_lt_p4().at(hyp)).pt() > (tas::hyp_ll_p4().at(best_hyp_) + tas::hyp_lt_p4().at(best_hyp_)).pt()) best_hyp_ = hyp;
    }
  }

  if (best_hyp_ < 0){
    hyp_result_t null = { -1, -1 };
    return null;
  }

  hyp_result_t temp;
  temp.best_hyp = best_hyp_;
  temp.hyp_class = hyp_class_; 
  return temp;
}

vector <particle_t> getGenPair(bool verbose){

  vector <particle_t> gen_particles;

  //First get all gen leptons 
  for (unsigned int gidx = 0; gidx < tas::genps_p4().size(); gidx++){
    if (tas::genps_status().at(gidx) != 1) continue;
    int id = tas::genps_id().at(gidx);
    if (abs(id) != 11 && abs(id) != 13 && abs(id) != 15) continue;
    float eta = tas::genps_p4().at(gidx).eta();
    if (fabs(eta) > 2.4) continue;
    int did = -1;
    int didx = -1;
    if (abs(id) == 15){  
      float dpt = -1;
      for (unsigned int didx_ = 0; didx_ < tas::genps_lepdaughter_id().at(gidx).size(); didx_++){
        int did_ = tas::genps_lepdaughter_id().at(gidx).at(didx_);
        if (abs(did_) != 11 && abs(did_) != 13) continue;
        if (fabs(tas::genps_lepdaughter_p4().at(gidx).at(didx_).eta()) > 2.4) continue;
        float dpt_ = tas::genps_lepdaughter_p4().at(gidx).at(didx_).pt();
        if (dpt_ <= dpt) continue;
        dpt = dpt_;
        didx = didx_;
        did = did_;
      }
    }
    particle_t temp;
    temp.id = abs(id) == 15 ? did : id;
    if (abs(id) == 15 && didx < 0) continue;
    temp.idx = abs(id) == 15 ? tas::genps_lepdaughter_idx().at(gidx).at(didx) : gidx;
    temp.p4 = abs(id) == 15 ? tas::genps_lepdaughter_p4().at(gidx).at(didx) : tas::genps_p4().at(gidx);
    if (temp.id != -1) gen_particles.push_back(temp);
  }
 
  if (gen_particles.size() < 2) return gen_particles;

  //Now loop over gen hyps and pick the best
  int gen_hyp_class = 5; 
  int type = 5; 
  particle_t lep1; 
  particle_t lep2;
  lep1.id = -1;
  lep2.id = -1;
  for (unsigned int idx1 = 0; idx1 < gen_particles.size(); idx1++){
    if (verbose) cout << "gen lep " << idx1 << ": " << gen_particles.at(idx1).id << " " << gen_particles.at(idx1).p4.pt() << endl;
    for (unsigned int idx2 = idx1+1; idx2 < gen_particles.size(); idx2++){
      int id1 = gen_particles.at(idx1).id;
      int id2 = gen_particles.at(idx2).id;
      int gen_hyp_class_ = 5;
      int type_ = 5;
      LorentzVector p41 = gen_particles.at(idx1).p4;
      LorentzVector p42 = gen_particles.at(idx2).p4;
      if (id1*id2 < 0) gen_hyp_class_ = 4; 
      if (id1*id2 > 0) gen_hyp_class_ = 3; 
      if (min(p41.pt(), p42.pt()) < 5) continue;
      if (max(p41.pt(), p42.pt()) < 5) continue;
      if (abs(id1) == 11 && abs(id2) == 11) type_ = 3; 
      if (abs(id1) == 13 && abs(id2) == 13) type_ = 0; 
      if ((abs(id1) == 13 && abs(id2) == 11) || (abs(id1) == 11 && abs(id2) == 13)) type_ = 1; 
      if (gen_hyp_class_ < gen_hyp_class){
        gen_hyp_class = gen_hyp_class_;
        type = type_;
        lep1 = gen_particles.at(idx1);
        lep2 = gen_particles.at(idx2);
      }
      else if (gen_hyp_class_ == gen_hyp_class && type_ < type){
        gen_hyp_class = gen_hyp_class_;
        type = type_;
        lep1 = gen_particles.at(idx1);
        lep2 = gen_particles.at(idx2);
      } 
      else if (gen_hyp_class_ == gen_hyp_class && type_ == type && p41.pt()+p42.pt() > lep1.p4.pt()+lep2.p4.pt()){
        gen_hyp_class = gen_hyp_class_;
        type = type_;
        lep1 = gen_particles.at(idx1);
        lep2 = gen_particles.at(idx2);
      }
    }
  }
  
  //Return leptons you found
  vector <particle_t> selected_gen_leptons;
  if (lep1.p4.pt() > lep2.p4.pt()){
    selected_gen_leptons.push_back(lep1);
    selected_gen_leptons.push_back(lep2);
  }
  else{ 
    selected_gen_leptons.push_back(lep2);
    selected_gen_leptons.push_back(lep1);
  }

  return selected_gen_leptons;

}

pair<particle_t, int> getThirdLepton(int hyp){

  //Selected Lepton Information
  int ll_id = tas::hyp_ll_id().at(hyp);
  int lt_id = tas::hyp_lt_id().at(hyp);
  unsigned int ll_idx = tas::hyp_ll_index().at(hyp);
  unsigned int lt_idx = tas::hyp_lt_index().at(hyp);

  //Store best lepton
  int lep3_id_ = -1;
  int lep3_idx_ = -1;
  int quality = 0;
  LorentzVector lep3_p4_; 

  //Electron Loop 
  for (unsigned int i = 0; i < tas::els_p4().size(); i++){

    //Remove electrons already selected
    if (abs(ll_id) == 11 && ll_idx == i) continue; 
    if (abs(lt_id) == 11 && lt_idx == i) continue; 

    //Remove electrons that fail kinematically
    if (tas::els_p4().at(i).pt() < 20) continue;
    if (fabs(tas::els_p4().at(i).eta()) > 2.4) continue;

    //Remove electrons that fail loosest ID, determine tighter IDs
    int quality_ = 0;
    if (!isGoodVetoElectron(i)) continue;
    if (isFakableElectron(i)) quality_ = 1;
    if (isGoodElectron(i)) quality_ = 2;

    //Choose the highest-quality, highest-pT electron 
    if (quality_ > quality || (quality_ == quality && tas::els_p4().at(i).pt() > lep3_p4_.pt())){
       quality = quality_;
       lep3_p4_ = tas::els_p4().at(i); 
       lep3_id_ = -11*tas::els_charge().at(i);
    } 
  }
  
  //Muon Loop
  for (unsigned int i = 0; i < tas::mus_p4().size(); i++){

    //Remove electrons already selected
    if (abs(ll_id) == 13 && ll_idx == i) continue; 
    if (abs(lt_id) == 13 && lt_idx == i) continue; 
   
    //Remove electrons that fail kinematically
    if (tas::mus_p4().at(i).pt() < 20) continue;
    if (fabs(tas::mus_p4().at(i).eta()) > 2.4) continue;

    //Remove muons that fail ID
    int quality_ = 0; 
    if (!isGoodVetoMuon(i)) continue;
    if (isFakableMuon(i)) quality_ = 1;
    if (isGoodMuon(i)) quality_ = 2;

    //Choose the highest-quality, highest-pT electron 
    if (quality_ > quality || (quality_ == quality && tas::mus_p4().at(i).pt() > lep3_p4_.pt())){
       quality = quality_;
       lep3_p4_ = tas::mus_p4().at(i); 
       lep3_id_ = -11*tas::mus_charge().at(i);
    } 

  }//Muon loop

  particle_t result;
  result.id = lep3_id_;
  result.p4 = lep3_p4_;
  result.idx = lep3_idx_;

  return pair<particle_t, int>(result, quality);

}
