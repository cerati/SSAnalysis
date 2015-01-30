#include "helper_babymaker.h"

using namespace tas;

//Main functions
void babyMaker::MakeBabyNtuple(const char* output_name){

  //Create Baby
  //TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  //rootdir->cd();
  BabyFile = new TFile(Form("%s/%s", path.Data(), output_name), "RECREATE");
  BabyFile->cd();
  BabyTree = new TTree("t", "SS2015 Baby Ntuple");

  //Define Branches
  BabyTree->Branch("met", &met);
  BabyTree->Branch("metPhi", &metPhi);
  BabyTree->Branch("event", &event);
  BabyTree->Branch("lumi", &lumi);
  BabyTree->Branch("run", &run);
  BabyTree->Branch("filt_csc", &filt_csc);
  BabyTree->Branch("filt_hbhe", &filt_hbhe);
  BabyTree->Branch("filt_hcallaser", &filt_hcallaser);
  BabyTree->Branch("filt_ecaltp", &filt_ecaltp);
  BabyTree->Branch("filt_trkfail", &filt_trkfail);
  BabyTree->Branch("filt_eebadsc", &filt_eebadsc);
  BabyTree->Branch("is_real_data", &is_real_data);
  BabyTree->Branch("scale1fb", &scale1fb);
  BabyTree->Branch("xsec", &xsec);
  BabyTree->Branch("kfactor", &kfactor);
  BabyTree->Branch("gen_met", &gen_met);
  BabyTree->Branch("gen_met_phi", &gen_met_phi);
  BabyTree->Branch("njets", &njets);
  BabyTree->Branch("hyp_class", &hyp_class);
  BabyTree->Branch("lep1_p4", &lep1_p4);
  BabyTree->Branch("lep2_p4", &lep2_p4);
  BabyTree->Branch("ht", &ht);
  BabyTree->Branch("lep1_motherID", &lep1_motherID);
  BabyTree->Branch("lep2_motherID", &lep2_motherID);
  BabyTree->Branch("lep1_mc_id", &lep1_mc_id);
  BabyTree->Branch("lep2_mc_id", &lep2_mc_id);
  BabyTree->Branch("lep1_id", &lep1_id);
  BabyTree->Branch("lep2_id", &lep2_id);
  BabyTree->Branch("lep1_idx", &lep1_idx);
  BabyTree->Branch("lep2_idx", &lep2_idx);
  BabyTree->Branch("jets", &jets);
  BabyTree->Branch("btags_disc", &btags_disc);
  BabyTree->Branch("jets_disc", &jets_disc);
  BabyTree->Branch("btags", &btags);
  BabyTree->Branch("nbtags", &nbtags);
  BabyTree->Branch("sf_dilepTrig_hpt", &sf_dilepTrig_hpt);
  BabyTree->Branch("sf_dilepTrig_lpt", &sf_dilepTrig_lpt);
  BabyTree->Branch("sf_dilepTrig_vlpt", &sf_dilepTrig_vlpt);
  BabyTree->Branch("hyp_type", &hyp_type);
  BabyTree->Branch("sf_dilep_eff", &sf_dilep_eff);
  BabyTree->Branch("mt", &mt);
  BabyTree->Branch("mt_l2", &mt_l2);
  BabyTree->Branch("mt2", &mt2);
  BabyTree->Branch("mGluino", &mGluino);
  BabyTree->Branch("mLSP", &mLSP);
  BabyTree->Branch("mSbottom", &mSbottom);
  BabyTree->Branch("mChargino", &mChargino);
  BabyTree->Branch("lep1_id_gen", &lep1_id_gen);
  BabyTree->Branch("lep2_id_gen", &lep2_id_gen);
  BabyTree->Branch("lep1_p4_gen", &lep1_p4_gen);
  BabyTree->Branch("lep2_p4_gen", &lep2_p4_gen);
  BabyTree->Branch("lep3_id", &lep3_id);
  BabyTree->Branch("lep3_idx", &lep3_idx);
  BabyTree->Branch("lep3_p4", &lep3_p4);
  BabyTree->Branch("lep3_quality", &lep3_quality);
  BabyTree->Branch("lep1_iso", &lep1_iso);
  BabyTree->Branch("lep2_iso", &lep2_iso);
  BabyTree->Branch("dilep_p4", &dilep_p4);
  BabyTree->Branch("genps_p4", &genps_p4);
  BabyTree->Branch("genps_id", &genps_id);
  BabyTree->Branch("genps_id_mother", &genps_id_mother);
  BabyTree->Branch("genps_status", &genps_status);
  BabyTree->Branch("genps_id_grandma", &genps_id_grandma);
  BabyTree->Branch("lep1_passes_id", &lep1_passes_id);
  BabyTree->Branch("lep2_passes_id", &lep2_passes_id);
  BabyTree->Branch("lep1_dxyPV", &lep1_dxyPV);
  BabyTree->Branch("lep2_dxyPV", &lep2_dxyPV);
  BabyTree->Branch("lep1_dZ", &lep1_dZ);
  BabyTree->Branch("lep2_dZ", &lep2_dZ);
  BabyTree->Branch("lep1_d0_err", &lep1_d0_err);
  BabyTree->Branch("lep2_d0_err", &lep2_d0_err);
  BabyTree->Branch("lep1_ip3d", &lep1_ip3d);
  BabyTree->Branch("lep2_ip3d", &lep2_ip3d);
  BabyTree->Branch("lep1_ip3d_err", &lep1_ip3d_err);
  BabyTree->Branch("lep2_ip3d_err", &lep2_ip3d_err);
  BabyTree->Branch("nGoodElectrons7", &nGoodElectrons7);
  BabyTree->Branch("nGoodElectrons10", &nGoodElectrons10);
  BabyTree->Branch("nGoodElectrons25", &nGoodElectrons25);
  BabyTree->Branch("nGoodMuons5", &nGoodMuons5);
  BabyTree->Branch("nGoodMuons10", &nGoodMuons10);
  BabyTree->Branch("nGoodMuons25", &nGoodMuons25);
  BabyTree->Branch("filename", &filename);

  //Print warning!
  cout << "Careful!! Path is " << path << endl;

}

void babyMaker::InitBabyNtuple(){

    //Gen variables
    met = -1;
    metPhi = -1;
    event = -1;
    lumi = -1;
    run = -1;
    filt_csc = 0;
    filt_hbhe = 0;
    filt_hcallaser = 0;
    filt_ecaltp = 0;
    filt_trkfail = 0;
    filt_eebadsc = 0;
    is_real_data = 0;
    scale1fb = -1;
    xsec = -1;
    kfactor = -1;
    gen_met = -1;
    gen_met_phi = -1;
    njets = -1;
    hyp_class = -1;
    ht = -1;
    lep1_motherID = 0;
    lep2_motherID = 0;
    lep1_mc_id = -1;
    lep2_mc_id = -1; 
    lep1_id = -1;
    lep2_id = -1;
    lep1_idx = -1;
    lep2_idx = -1;
    jets.clear();
    btags_disc.clear();
    jets_disc.clear();
    btags.clear();
    nbtags = -1;
    sf_dilepTrig_hpt = -1;
    sf_dilepTrig_lpt = -1;
    sf_dilepTrig_vlpt = -1;
    hyp_type = -1;
    sf_dilep_eff = -1;
    mt = -1;
    mt_l2 = -1;
    mt2 = -1;
    mGluino = -1;
    mLSP = -1;
    mSbottom = -1;
    mChargino = -1;
    lep1_id_gen = -1;
    lep2_id_gen = -1;
    lep3_id = -1;
    lep3_idx = -1;
    lep3_quality = -1;
    lep1_iso = -1;
    lep2_iso = -1;
    genps_p4.clear();
    genps_id.clear();
    genps_id_mother.clear();
    genps_status.clear();
    genps_id_grandma.clear();
    lep1_passes_id = false;
    lep2_passes_id = false;
    lep1_dxyPV = -999998;
    lep2_dxyPV = -999998;
    lep1_dZ = -999998;
    lep2_dZ = -999998;
    lep1_d0_err = -999998;
    lep2_d0_err = -999998;
    lep1_ip3d = -999998;
    lep2_ip3d = -999998;
    lep1_ip3d_err = -999998;
    lep2_ip3d_err = -999998;
    nGoodElectrons7 = 0;
    nGoodElectrons10 = 0;
    nGoodElectrons25 = 0;
    nGoodMuons5 = 0;
    nGoodMuons10 = 0;
    nGoodMuons25 = 0;
    filename = "";
} 

template <typename T> int sgn(T val){
    return (T(0) < val) - (val < T(0));
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

float leptonDxy(const int id, const int idx){
  if (abs(id)==13) return tas::mus_dxyPV().at(idx); 
  else if (abs(id)==11) return tas::els_dxyPV().at(idx);   
  return -999999.0;
}

val_err_t leptonIp3d(const int id, const int idx){
  val_err_t ip3d;
  if (abs(id) == 11){
    ip3d.value = tas::els_ip3d().at(idx); 
    ip3d.error = tas::els_ip3derr().at(idx); 
  }
  else if (abs(id) == 13){
    ip3d.value = tas::mus_ip3d().at(idx); 
    ip3d.error = tas::mus_ip3derr().at(idx); 
  }
  else{
    ip3d.value = -999998;
    ip3d.error = -999998;
  }
  return ip3d;
}

int babyMaker::isGoodHyp(int iHyp, int analType){

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

hyp_result_t babyMaker::chooseBestHyp(){

  //List of good hyps
  vector <int> good_hyps_ss; 
  vector <int> good_hyps_sf; 
  vector <int> good_hyps_df; 
  vector <int> good_hyps_os; 
  for (unsigned int i = 0; i < tas::hyp_type().size(); i++){
    int good_hyp_result = isGoodHyp(i, 2);
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

float DileptonTriggerScaleFactor(int hyp_type, anal_type_t anal_type, LorentzVector trailing_p4){

  //Local variables
  const bool is_barrel = (fabs(trailing_p4.eta()) < 1.0);
  const float lep_pt   = trailing_p4.pt();
  
  // values supplied by ETH (LHCP value)
  if (anal_type == HIGH_PT){
      if (hyp_type == 3) return (lep_pt < 30 ? 0.92 : 0.96); 
      else if (hyp_type == 0) return (is_barrel   ? 0.90 : 0.81); 
      else return 0.93;                        
  }
  if (anal_type == LOW_PT){
      if (hyp_type == 0) return 0.93;                      
      else if (hyp_type == 3) return (is_barrel ? 0.94 : 0.90); 
      else return 0.93;                      
  }
  if (anal_type == VLOW_PT){
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


vector <particle_t> babyMaker::getGenPair(){

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

float leptonD0err(const int id, const int idx){
    if (abs(id)==13) return tas::mus_d0Err().at(idx); 
    else if (abs(id)==11) return tas::els_d0Err().at(idx);   
    return -999999.0;
}

float leptonDZ(const int id, const int idx){
    if (abs(id) == 11) return tas::els_dzPV().at(idx);
    if (abs(id) == 13) return tas::mus_dzPV().at(idx);
    return -999999.0;
}

pair<particle_t, int> babyMaker::getThirdLepton(int hyp){

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

bool leptonIsFromWag(int idx, int id){

  LorentzVector lep = abs(id) == 11 ? tas::els_p4().at(idx) : tas::mus_p4().at(idx); 
  float minDR = 9999;
  int genPart = -1;

  for (unsigned int i = 0; i < tas::genps_id().size(); i++){
    if (tas::genps_status().at(i) != 23) continue;
    if (abs(tas::genps_id().at(i)) != 11 && abs(tas::genps_id().at(i)) != 13 && abs(tas::genps_id().at(i) != 15)) continue;
    float minDR_ = ROOT::Math::VectorUtil::DeltaR(lep, tas::genps_p4().at(i)); 
    if (minDR_ < minDR){
      minDR = minDR_;
      genPart = i;
    }
  }

  if (genPart < 0) return false;
  if (minDR > 0.4) return false;
  if (abs(tas::genps_id_mother().at(genPart)) != 24) return false;
  return true;

}

//Main function
int babyMaker::ProcessBaby(){

  //Initialize variables
  InitBabyNtuple();
  
  //Local variables
  bool isData = tas::evt_isRealData();
  
  //Debug mode
  if (verbose && tas::evt_event() != evt_cut) return -1;
  //if (verbose) cout << "file name is " << file->GetName() << endl;

  //Preliminary stuff
  if (tas::hyp_type().size() < 1) return -1;
  if (tas::mus_dxyPV().size() != tas::mus_dzPV().size()) return -1;

  //Fill Easy Variables
  //filename = currentFile->GetTitle();//fixme
  met = evt_pfmet();
  metPhi = evt_pfmetPhi();
  event = tas::evt_event();
  lumi = tas::evt_lumiBlock();
  run = tas::evt_run();
  is_real_data = tas::evt_isRealData();
  xsec = tas::evt_xsec_incl();
  kfactor = tas::evt_kfactor();
  gen_met = tas::gen_met();
  gen_met_phi = tas::gen_metPhi();

  //Fill data vs. mc variables
  filt_csc = is_real_data ? tas::evt_cscTightHaloId() : 1;
  filt_hbhe = is_real_data ? (tas::evt_hbheFilter() && tas::hcalnoise_isolatedNoiseSumE() < 50.0 && tas::hcalnoise_isolatedNoiseSumEt() < 25.0 && tas::hcalnoise_numIsolatedNoiseChannels() < 10) : 1;
  filt_hcallaser = is_real_data ? tas::filt_hcalLaser() : 1;
  filt_ecaltp = is_real_data ? tas::filt_ecalTP() : 1;
  filt_trkfail = is_real_data ? tas::filt_trackingFailure() : 1;
  filt_eebadsc = is_real_data ? tas::filt_eeBadSc() : 1;
  scale1fb = is_real_data ? 1 : tas::evt_scale1fb();
  
  //Fill signal variables
  /*
    mGluino = signal == T1TTTT ? tas::sparm_values().at(0): -1;
    mLSP = signal == T1TTTT ? tas::sparm_values().at(1) : -1;
    mSbottom = signal == T4TW ? tas::sparm_values().at(0) : -1;
    mChargino = signal == T4TW ? tas::sparm_values().at(1) : -1;
    if (signal == T1TTTT && (mGluino != 1500 || mLSP != 100)) continue;
    if (signal == T4TW && (mSbottom != 500 || mChargino != 200)) continue;
    if (signal == T1TTTT) scale1fb = 1000*gluino_cs(mGluino).value/20000;
    if (signal == T4TW) scale1fb = 1000*stop_sbottom_cs(mSbottom).value/47500;
    if (signal != NONE && verbose) cout << scale1fb << endl;
  */
  
  //Fill lepton variables
  hyp_result_t best_hyp_info = chooseBestHyp();
  hyp_class = best_hyp_info.hyp_class;
  int best_hyp = best_hyp_info.best_hyp;
  if (verbose) cout << "chose hyp: " << best_hyp << " of class" << hyp_class << endl;
  if (hyp_class == 0 || hyp_class == -1) return -1;
  lep1_p4 = (tas::hyp_ll_p4().at(best_hyp).pt() > tas::hyp_lt_p4().at(best_hyp).pt()) ? tas::hyp_ll_p4().at(best_hyp) : tas::hyp_lt_p4().at(best_hyp);
  lep2_p4 = (tas::hyp_ll_p4().at(best_hyp).pt() <= tas::hyp_lt_p4().at(best_hyp).pt()) ? tas::hyp_ll_p4().at(best_hyp) : tas::hyp_lt_p4().at(best_hyp);
  lep1_id = (tas::hyp_ll_p4().at(best_hyp).pt() > tas::hyp_lt_p4().at(best_hyp).pt()) ? tas::hyp_ll_id().at(best_hyp) : tas::hyp_lt_id().at(best_hyp);
  lep2_id = (tas::hyp_ll_p4().at(best_hyp).pt() <= tas::hyp_lt_p4().at(best_hyp).pt()) ? tas::hyp_ll_id().at(best_hyp) : tas::hyp_lt_id().at(best_hyp);
  lep1_idx = (tas::hyp_ll_p4().at(best_hyp).pt() > tas::hyp_lt_p4().at(best_hyp).pt()) ? tas::hyp_ll_index().at(best_hyp) : tas::hyp_lt_index().at(best_hyp);
  lep2_idx = (tas::hyp_ll_p4().at(best_hyp).pt() <= tas::hyp_lt_p4().at(best_hyp).pt()) ? tas::hyp_ll_index().at(best_hyp) : tas::hyp_lt_index().at(best_hyp);
  lep1_dxyPV = leptonDxy(lep1_id, lep1_idx);
  lep2_dxyPV = leptonDxy(lep2_id, lep2_idx);
  lep1_dZ = leptonDZ(lep1_id, lep1_idx);
  lep2_dZ = leptonDZ(lep2_id, lep2_idx);
  lep1_d0_err = leptonD0err(lep1_id, lep1_idx);
  lep2_d0_err = leptonD0err(lep2_id, lep2_idx);
  lep1_ip3d = leptonIp3d(lep1_id, lep1_idx).value;
  lep2_ip3d = leptonIp3d(lep2_id, lep2_idx).value;
  lep1_ip3d_err = leptonIp3d(lep1_id, lep1_idx).error;
  lep2_ip3d_err = leptonIp3d(lep2_id, lep2_idx).error;
  Lep lep1 = Lep(lep1_id, lep1_idx);
  Lep lep2 = Lep(lep2_id, lep2_idx);
  lep1_motherID = lepMotherID(lep1);
  lep2_motherID = lepMotherID(lep2);
  lep1_mc_id = abs(lep1_id) == 11 ? tas::els_mc_id().at(lep1_idx) : tas::mus_mc_id().at(lep1_idx);
  lep2_mc_id = abs(lep2_id) == 11 ? tas::els_mc_id().at(lep2_idx) : tas::mus_mc_id().at(lep2_idx);
  hyp_type = tas::hyp_type().at(best_hyp);
  pair <particle_t, int> thirdLepton = getThirdLepton(best_hyp);
  lep3_id = thirdLepton.first.id;
  lep3_idx = thirdLepton.first.idx;
  lep3_p4 = thirdLepton.first.p4;
  lep3_quality = thirdLepton.second;
  lep1_iso = isIsolatedLepton(lep1_id, lep1_idx);
  lep2_iso = isIsolatedLepton(lep2_id, lep2_idx);
  dilep_p4 = lep1_p4 + lep2_p4; 
  lep1_passes_id = isGoodLepton(lep1_id, lep1_idx);
  lep2_passes_id = isGoodLepton(lep2_id, lep2_idx);
  
  //Fill generated lepton variables, ignoring reco (matching to reco done above)
  vector <particle_t> genPair = getGenPair();
  if (genPair.size() == 2){
    lep1_id_gen = genPair.at(0).id;
    lep2_id_gen = genPair.at(1).id;
    lep1_p4_gen = genPair.at(0).p4;
    lep2_p4_gen = genPair.at(1).p4;
  }
  
  //Fill all generated particles
  if (!isData){
    genps_p4 = tas::genps_p4();
    genps_id = tas::genps_id();
    genps_id_mother = tas::genps_id_mother();
    genps_status = tas::genps_status(); 
    genps_id_grandma = tas::genps_id_simplegrandma(); 
  }
  
  //Determine and save jet and b-tag variables
  ht = 0;
  for (unsigned int i = 0; i < tas::pfjets_p4().size(); i++){
    LorentzVector jet = tas::pfjets_p4().at(i);
    
    //Kinematic jet cuts
    if (jet.pt() < 40) continue;
    if (fabs(jet.eta()) > 2.4) continue;
    
    //Verbose
    if (verbose) cout << "Possible jet with pT: " << jet.pt() << endl;
    
    //Require loose jet ID
    if (!isLoosePFJet(i)) continue;
    
    //Jet cleaning -- electrons
    bool jetIsLep = false;
    for (unsigned int eidx = 0; eidx < tas::els_p4().size(); eidx++){
      LorentzVector electron = tas::els_p4().at(eidx);
      if (electron.pt() < 7) continue;
      if (!isGoodVetoElectron(eidx)) continue;
      if (ROOT::Math::VectorUtil::DeltaR(jet, electron) > 0.4) continue;
      jetIsLep = true;
    }
    if (jetIsLep == true) continue;
    
    //Jet cleaning -- muons
    for (unsigned int muidx = 0; muidx < tas::mus_p4().size(); muidx++){
      LorentzVector muon = tas::mus_p4().at(muidx);
      if (muon.pt() < 5) continue;
      if (!isGoodVetoMuon(muidx)) continue;
      if (ROOT::Math::VectorUtil::DeltaR(jet, muon) > 0.4) continue;
      jetIsLep = true;
    }
    if (jetIsLep == true) continue;
    
    //Save jets that make it this far
    jets.push_back(jet);
    ht += jet.pt();
    float disc = tas::pfjets_combinedInclusiveSecondaryVertexV2BJetTag().at(i);
    jets_disc.push_back(disc);
    
    //Btag discriminator
    if (disc < 0.814) continue;
    
    //Save btags that make it this far
    btags_disc.push_back(disc);
    btags.push_back(tas::pfjets_p4().at(i));
  }
  njets = jets.size();
  nbtags = btags.size();
  
  //Verbose for jets
  if (verbose){
    cout << "njets: " << njets << endl;
    cout << "nbtags: " <<  nbtags << endl;
    for (unsigned int i = 0; i < jets.size(); i++) cout << i << " " << jets[i].pt() << " " << jets[i].eta() << endl;
  } 
  
  //nGood Leptons
  for (unsigned int eidx = 0; eidx < tas::els_p4().size(); eidx++){
    if (!isGoodVetoElectron(eidx)) continue;
    if (tas::els_p4().at(eidx).pt() < 7) continue;
    nGoodElectrons7++;
    if (verbose) cout << "good elec: " << tas::els_p4().at(eidx).pt() << endl;
    if (tas::els_p4().at(eidx).pt() < 10) continue;
    nGoodElectrons10++;
    if (tas::els_p4().at(eidx).pt() < 25) continue;
    nGoodElectrons25++;
  }
  for (unsigned int muidx = 0; muidx < tas::mus_p4().size(); muidx++){
    if (!isGoodVetoMuon(muidx)) continue;
    if (tas::mus_p4().at(muidx).pt() < 5) continue;
    nGoodMuons5++;
    if (verbose) cout << "good muon: " << tas::mus_p4().at(muidx).pt() << endl;
    if (tas::mus_p4().at(muidx).pt() < 10) continue;
    nGoodMuons10++;
    if (tas::mus_p4().at(muidx).pt() < 25) continue;
    nGoodMuons25++;
  }
  
  //Scale Factors
  sf_dilepTrig_hpt = DileptonTriggerScaleFactor(hyp_type, HIGH_PT, lep2_p4);
  sf_dilepTrig_lpt = DileptonTriggerScaleFactor(hyp_type, LOW_PT, lep2_p4); 
  sf_dilepTrig_vlpt = DileptonTriggerScaleFactor(hyp_type, VLOW_PT, lep2_p4);
  sf_dilep_eff = DileptonTagAndProbeScaleFactor(best_hyp);
  if (verbose) cout << "sf_dilep_eff: " << sf_dilep_eff << endl;
  
  //MT variables
  mt = calculateMt(lep1_p4, met, metPhi);
  mt2 = MT2(met, metPhi, lep1_p4, lep2_p4);
  mt_l2 = calculateMt(lep2_p4, met, metPhi);
  
  //Ht and MET
  if (verbose) cout << "ht: " << ht << endl;
  if (verbose) cout << "met: " << met << endl;
  
  //Fill Baby
  BabyTree->Fill();

  return 0;  

}
