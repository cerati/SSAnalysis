// -*- C++ -*-
#ifndef CMS2_H
#define CMS2_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"
#include <vector> 
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std; 
class CMS2 {
private: 
protected: 
	unsigned int index;
	TBits hlt_bits_;
	TBranch *hlt_bits_branch;
	bool hlt_bits_isLoaded;
	vector<TString> evt_CMS2tag_;
	TBranch *evt_CMS2tag_branch;
	bool evt_CMS2tag_isLoaded;
	vector<TString> evt_dataset_;
	TBranch *evt_dataset_branch;
	bool evt_dataset_isLoaded;
	vector<TString> hlt_trigNames_;
	TBranch *hlt_trigNames_branch;
	bool hlt_trigNames_isLoaded;
	bool hcalnoise_HasBadRBXTS4TS5_;
	TBranch *hcalnoise_HasBadRBXTS4TS5_branch;
	bool hcalnoise_HasBadRBXTS4TS5_isLoaded;
	bool evt_cscTightHaloId_;
	TBranch *evt_cscTightHaloId_branch;
	bool evt_cscTightHaloId_isLoaded;
	bool evt_hbheFilter_;
	TBranch *evt_hbheFilter_branch;
	bool evt_hbheFilter_isLoaded;
	bool filt_cscBeamHalo_;
	TBranch *filt_cscBeamHalo_branch;
	bool filt_cscBeamHalo_isLoaded;
	bool filt_ecalLaser_;
	TBranch *filt_ecalLaser_branch;
	bool filt_ecalLaser_isLoaded;
	bool filt_ecalTP_;
	TBranch *filt_ecalTP_branch;
	bool filt_ecalTP_isLoaded;
	bool filt_eeBadSc_;
	TBranch *filt_eeBadSc_branch;
	bool filt_eeBadSc_isLoaded;
	bool filt_goodVertices_;
	TBranch *filt_goodVertices_branch;
	bool filt_goodVertices_isLoaded;
	bool filt_hbheNoise_;
	TBranch *filt_hbheNoise_branch;
	bool filt_hbheNoise_isLoaded;
	bool filt_hcalLaser_;
	TBranch *filt_hcalLaser_branch;
	bool filt_hcalLaser_isLoaded;
	bool filt_metfilter_;
	TBranch *filt_metfilter_branch;
	bool filt_metfilter_isLoaded;
	bool filt_trackingFailure_;
	TBranch *filt_trackingFailure_branch;
	bool filt_trackingFailure_isLoaded;
	bool filt_trkPOGFilters_;
	TBranch *filt_trkPOGFilters_branch;
	bool filt_trkPOGFilters_isLoaded;
	bool filt_trkPOG_logErrorTooManyClusters_;
	TBranch *filt_trkPOG_logErrorTooManyClusters_branch;
	bool filt_trkPOG_logErrorTooManyClusters_isLoaded;
	bool filt_trkPOG_manystripclus53X_;
	TBranch *filt_trkPOG_manystripclus53X_branch;
	bool filt_trkPOG_manystripclus53X_isLoaded;
	bool filt_trkPOG_toomanystripclus53X_;
	TBranch *filt_trkPOG_toomanystripclus53X_branch;
	bool filt_trkPOG_toomanystripclus53X_isLoaded;
	vector<bool> els_conv_vtx_flag_;
	TBranch *els_conv_vtx_flag_branch;
	bool els_conv_vtx_flag_isLoaded;
	vector<bool> els_isGsfCtfScPixChargeConsistent_;
	TBranch *els_isGsfCtfScPixChargeConsistent_branch;
	bool els_isGsfCtfScPixChargeConsistent_isLoaded;
	vector<bool> els_passingMvaPreselection_;
	TBranch *els_passingMvaPreselection_branch;
	bool els_passingMvaPreselection_isLoaded;
	vector<bool> els_passingPflowPreselection_;
	TBranch *els_passingPflowPreselection_branch;
	bool els_passingPflowPreselection_isLoaded;
	vector<bool> mus_isRPCMuon_;
	TBranch *mus_isRPCMuon_branch;
	bool mus_isRPCMuon_isLoaded;
	vector<bool> mus_tightMatch_;
	TBranch *mus_tightMatch_branch;
	bool mus_tightMatch_isLoaded;
	vector<bool> mus_updatedSta_;
	TBranch *mus_updatedSta_branch;
	bool mus_updatedSta_isLoaded;
	vector<bool> photons_haspixelSeed_;
	TBranch *photons_haspixelSeed_branch;
	bool photons_haspixelSeed_isLoaded;
	float evt_bs_Xwidth_;
	TBranch *evt_bs_Xwidth_branch;
	bool evt_bs_Xwidth_isLoaded;
	float evt_bs_XwidthErr_;
	TBranch *evt_bs_XwidthErr_branch;
	bool evt_bs_XwidthErr_isLoaded;
	float evt_bs_Ywidth_;
	TBranch *evt_bs_Ywidth_branch;
	bool evt_bs_Ywidth_isLoaded;
	float evt_bs_YwidthErr_;
	TBranch *evt_bs_YwidthErr_branch;
	bool evt_bs_YwidthErr_isLoaded;
	float evt_bs_dxdz_;
	TBranch *evt_bs_dxdz_branch;
	bool evt_bs_dxdz_isLoaded;
	float evt_bs_dxdzErr_;
	TBranch *evt_bs_dxdzErr_branch;
	bool evt_bs_dxdzErr_isLoaded;
	float evt_bs_dydz_;
	TBranch *evt_bs_dydz_branch;
	bool evt_bs_dydz_isLoaded;
	float evt_bs_dydzErr_;
	TBranch *evt_bs_dydzErr_branch;
	bool evt_bs_dydzErr_isLoaded;
	float evt_bs_sigmaZ_;
	TBranch *evt_bs_sigmaZ_branch;
	bool evt_bs_sigmaZ_isLoaded;
	float evt_bs_sigmaZErr_;
	TBranch *evt_bs_sigmaZErr_branch;
	bool evt_bs_sigmaZErr_isLoaded;
	float evt_bs_xErr_;
	TBranch *evt_bs_xErr_branch;
	bool evt_bs_xErr_isLoaded;
	float evt_bs_yErr_;
	TBranch *evt_bs_yErr_branch;
	bool evt_bs_yErr_isLoaded;
	float evt_bs_zErr_;
	TBranch *evt_bs_zErr_branch;
	bool evt_bs_zErr_isLoaded;
	float evt_bField_;
	TBranch *evt_bField_branch;
	bool evt_bField_isLoaded;
	float evt_fixgrid_all_rho_;
	TBranch *evt_fixgrid_all_rho_branch;
	bool evt_fixgrid_all_rho_isLoaded;
	float evt_fixgridfastjet_allcalo_rho_;
	TBranch *evt_fixgridfastjet_allcalo_rho_branch;
	bool evt_fixgridfastjet_allcalo_rho_isLoaded;
	float evt_fixgridfastjet_all_rho_;
	TBranch *evt_fixgridfastjet_all_rho_branch;
	bool evt_fixgridfastjet_all_rho_isLoaded;
	float evt_fixgridfastjet_centralcalo_rho_;
	TBranch *evt_fixgridfastjet_centralcalo_rho_branch;
	bool evt_fixgridfastjet_centralcalo_rho_isLoaded;
	float evt_fixgridfastjet_centralchargedpileup_rho_;
	TBranch *evt_fixgridfastjet_centralchargedpileup_rho_branch;
	bool evt_fixgridfastjet_centralchargedpileup_rho_isLoaded;
	float evt_fixgridfastjet_centralneutral_rho_;
	TBranch *evt_fixgridfastjet_centralneutral_rho_branch;
	bool evt_fixgridfastjet_centralneutral_rho_isLoaded;
	float evt_kfactor_;
	TBranch *evt_kfactor_branch;
	bool evt_kfactor_isLoaded;
	float evt_scale1fb_;
	TBranch *evt_scale1fb_branch;
	bool evt_scale1fb_isLoaded;
	float evt_xsec_excl_;
	TBranch *evt_xsec_excl_branch;
	bool evt_xsec_excl_isLoaded;
	float evt_xsec_incl_;
	TBranch *evt_xsec_incl_branch;
	bool evt_xsec_incl_isLoaded;
	float gen_met_;
	TBranch *gen_met_branch;
	bool gen_met_isLoaded;
	float gen_metPhi_;
	TBranch *gen_metPhi_branch;
	bool gen_metPhi_isLoaded;
	float genps_alphaQCD_;
	TBranch *genps_alphaQCD_branch;
	bool genps_alphaQCD_isLoaded;
	float genps_pthat_;
	TBranch *genps_pthat_branch;
	bool genps_pthat_isLoaded;
	float genps_qScale_;
	TBranch *genps_qScale_branch;
	bool genps_qScale_isLoaded;
	float genps_weight_;
	TBranch *genps_weight_branch;
	bool genps_weight_isLoaded;
	float gen_sumEt_;
	TBranch *gen_sumEt_branch;
	bool gen_sumEt_isLoaded;
	float hcalnoise_GetRecHitEnergy_;
	TBranch *hcalnoise_GetRecHitEnergy_branch;
	bool hcalnoise_GetRecHitEnergy_isLoaded;
	float hcalnoise_GetRecHitEnergy15_;
	TBranch *hcalnoise_GetRecHitEnergy15_branch;
	bool hcalnoise_GetRecHitEnergy15_isLoaded;
	float hcalnoise_GetTotalCalibCharge_;
	TBranch *hcalnoise_GetTotalCalibCharge_branch;
	bool hcalnoise_GetTotalCalibCharge_isLoaded;
	float hcalnoise_TS4TS5NoiseSumE_;
	TBranch *hcalnoise_TS4TS5NoiseSumE_branch;
	bool hcalnoise_TS4TS5NoiseSumE_isLoaded;
	float hcalnoise_TS4TS5NoiseSumEt_;
	TBranch *hcalnoise_TS4TS5NoiseSumEt_branch;
	bool hcalnoise_TS4TS5NoiseSumEt_isLoaded;
	float hcalnoise_eventChargeFraction_;
	TBranch *hcalnoise_eventChargeFraction_branch;
	bool hcalnoise_eventChargeFraction_isLoaded;
	float hcalnoise_eventEMEnergy_;
	TBranch *hcalnoise_eventEMEnergy_branch;
	bool hcalnoise_eventEMEnergy_isLoaded;
	float hcalnoise_eventEMFraction_;
	TBranch *hcalnoise_eventEMFraction_branch;
	bool hcalnoise_eventEMFraction_isLoaded;
	float hcalnoise_eventHadEnergy_;
	TBranch *hcalnoise_eventHadEnergy_branch;
	bool hcalnoise_eventHadEnergy_isLoaded;
	float hcalnoise_eventTrackEnergy_;
	TBranch *hcalnoise_eventTrackEnergy_branch;
	bool hcalnoise_eventTrackEnergy_isLoaded;
	float hcalnoise_flatNoiseSumE_;
	TBranch *hcalnoise_flatNoiseSumE_branch;
	bool hcalnoise_flatNoiseSumE_isLoaded;
	float hcalnoise_flatNoiseSumEt_;
	TBranch *hcalnoise_flatNoiseSumEt_branch;
	bool hcalnoise_flatNoiseSumEt_isLoaded;
	float hcalnoise_isolatedNoiseSumE_;
	TBranch *hcalnoise_isolatedNoiseSumE_branch;
	bool hcalnoise_isolatedNoiseSumE_isLoaded;
	float hcalnoise_isolatedNoiseSumEt_;
	TBranch *hcalnoise_isolatedNoiseSumEt_branch;
	bool hcalnoise_isolatedNoiseSumEt_isLoaded;
	float hcalnoise_max10GeVHitTime_;
	TBranch *hcalnoise_max10GeVHitTime_branch;
	bool hcalnoise_max10GeVHitTime_isLoaded;
	float hcalnoise_max25GeVHitTime_;
	TBranch *hcalnoise_max25GeVHitTime_branch;
	bool hcalnoise_max25GeVHitTime_isLoaded;
	float hcalnoise_maxE10TS_;
	TBranch *hcalnoise_maxE10TS_branch;
	bool hcalnoise_maxE10TS_isLoaded;
	float hcalnoise_maxE2Over10TS_;
	TBranch *hcalnoise_maxE2Over10TS_branch;
	bool hcalnoise_maxE2Over10TS_isLoaded;
	float hcalnoise_maxE2TS_;
	TBranch *hcalnoise_maxE2TS_branch;
	bool hcalnoise_maxE2TS_isLoaded;
	float hcalnoise_min10GeVHitTime_;
	TBranch *hcalnoise_min10GeVHitTime_branch;
	bool hcalnoise_min10GeVHitTime_isLoaded;
	float hcalnoise_min25GeVHitTime_;
	TBranch *hcalnoise_min25GeVHitTime_branch;
	bool hcalnoise_min25GeVHitTime_isLoaded;
	float hcalnoise_minE10TS_;
	TBranch *hcalnoise_minE10TS_branch;
	bool hcalnoise_minE10TS_isLoaded;
	float hcalnoise_minE2Over10TS_;
	TBranch *hcalnoise_minE2Over10TS_branch;
	bool hcalnoise_minE2Over10TS_isLoaded;
	float hcalnoise_minE2TS_;
	TBranch *hcalnoise_minE2TS_branch;
	bool hcalnoise_minE2TS_isLoaded;
	float hcalnoise_minHPDEMF_;
	TBranch *hcalnoise_minHPDEMF_branch;
	bool hcalnoise_minHPDEMF_isLoaded;
	float hcalnoise_minRBXEMF_;
	TBranch *hcalnoise_minRBXEMF_branch;
	bool hcalnoise_minRBXEMF_isLoaded;
	float hcalnoise_rms10GeVHitTime_;
	TBranch *hcalnoise_rms10GeVHitTime_branch;
	bool hcalnoise_rms10GeVHitTime_isLoaded;
	float hcalnoise_rms25GeVHitTime_;
	TBranch *hcalnoise_rms25GeVHitTime_branch;
	bool hcalnoise_rms25GeVHitTime_isLoaded;
	float hcalnoise_spikeNoiseSumE_;
	TBranch *hcalnoise_spikeNoiseSumE_branch;
	bool hcalnoise_spikeNoiseSumE_isLoaded;
	float hcalnoise_spikeNoiseSumEt_;
	TBranch *hcalnoise_spikeNoiseSumEt_branch;
	bool hcalnoise_spikeNoiseSumEt_isLoaded;
	float hcalnoise_triangleNoiseSumE_;
	TBranch *hcalnoise_triangleNoiseSumE_branch;
	bool hcalnoise_triangleNoiseSumE_isLoaded;
	float hcalnoise_triangleNoiseSumEt_;
	TBranch *hcalnoise_triangleNoiseSumEt_branch;
	bool hcalnoise_triangleNoiseSumEt_isLoaded;
	float pdfinfo_pdf1_;
	TBranch *pdfinfo_pdf1_branch;
	bool pdfinfo_pdf1_isLoaded;
	float pdfinfo_pdf2_;
	TBranch *pdfinfo_pdf2_branch;
	bool pdfinfo_pdf2_isLoaded;
	float pdfinfo_scale_;
	TBranch *pdfinfo_scale_branch;
	bool pdfinfo_scale_isLoaded;
	float pdfinfo_x1_;
	TBranch *pdfinfo_x1_branch;
	bool pdfinfo_x1_isLoaded;
	float pdfinfo_x2_;
	TBranch *pdfinfo_x2_branch;
	bool pdfinfo_x2_isLoaded;
	float evt_fixgrid_rho_all_;
	TBranch *evt_fixgrid_rho_all_branch;
	bool evt_fixgrid_rho_all_isLoaded;
	float evt_fixgrid_rho_ctr_;
	TBranch *evt_fixgrid_rho_ctr_branch;
	bool evt_fixgrid_rho_ctr_isLoaded;
	float evt_fixgrid_rho_fwd_;
	TBranch *evt_fixgrid_rho_fwd_branch;
	bool evt_fixgrid_rho_fwd_isLoaded;
	float evt_pfmet_;
	TBranch *evt_pfmet_branch;
	bool evt_pfmet_isLoaded;
	float evt_pfmetPhi_;
	TBranch *evt_pfmetPhi_branch;
	bool evt_pfmetPhi_isLoaded;
	float evt_pfmetSig_;
	TBranch *evt_pfmetSig_branch;
	bool evt_pfmetSig_isLoaded;
	float evt_pfmetSignificance_;
	TBranch *evt_pfmetSignificance_branch;
	bool evt_pfmetSignificance_isLoaded;
	float evt_pfsumet_;
	TBranch *evt_pfsumet_branch;
	bool evt_pfsumet_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  evt_bsp4_;
	TBranch *evt_bsp4_branch;
	bool evt_bsp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_mc_motherp4_;
	TBranch *els_mc_motherp4_branch;
	bool els_mc_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_mc_p4_;
	TBranch *els_mc_p4_branch;
	bool els_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_mc_motherp4_;
	TBranch *mus_mc_motherp4_branch;
	bool mus_mc_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_mc_p4_;
	TBranch *mus_mc_p4_branch;
	bool mus_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfjets_mc_gp_p4_;
	TBranch *pfjets_mc_gp_p4_branch;
	bool pfjets_mc_gp_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfjets_mc_motherp4_;
	TBranch *pfjets_mc_motherp4_branch;
	bool pfjets_mc_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfjets_mc_p4_;
	TBranch *pfjets_mc_p4_branch;
	bool pfjets_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > photons_mc_motherp4_;
	TBranch *photons_mc_motherp4_branch;
	bool photons_mc_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > photons_mc_p4_;
	TBranch *photons_mc_p4_branch;
	bool photons_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_conv_pos_p4_;
	TBranch *els_conv_pos_p4_branch;
	bool els_conv_pos_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_inner_position_;
	TBranch *els_inner_position_branch;
	bool els_inner_position_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_mc_patMatch_p4_;
	TBranch *els_mc_patMatch_p4_branch;
	bool els_mc_patMatch_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_outer_position_;
	TBranch *els_outer_position_branch;
	bool els_outer_position_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_p4_;
	TBranch *els_p4_branch;
	bool els_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_p4In_;
	TBranch *els_p4In_branch;
	bool els_p4In_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_p4Out_;
	TBranch *els_p4Out_branch;
	bool els_p4Out_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_trk_p4_;
	TBranch *els_trk_p4_branch;
	bool els_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_trk_vertex_p4_;
	TBranch *els_trk_vertex_p4_branch;
	bool els_trk_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_vertex_p4_;
	TBranch *els_vertex_p4_branch;
	bool els_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > genjets_p4NoMuNoNu_;
	TBranch *genjets_p4NoMuNoNu_branch;
	bool genjets_p4NoMuNoNu_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > genps_p4_;
	TBranch *genps_p4_branch;
	bool genps_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > genps_prod_vtx_;
	TBranch *genps_prod_vtx_branch;
	bool genps_prod_vtx_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_gfit_p4_;
	TBranch *mus_gfit_p4_branch;
	bool mus_gfit_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_gfit_vertex_p4_;
	TBranch *mus_gfit_vertex_p4_branch;
	bool mus_gfit_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_mc_patMatch_p4_;
	TBranch *mus_mc_patMatch_p4_branch;
	bool mus_mc_patMatch_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_p4_;
	TBranch *mus_p4_branch;
	bool mus_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_pfp4_;
	TBranch *mus_pfp4_branch;
	bool mus_pfp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_sta_p4_;
	TBranch *mus_sta_p4_branch;
	bool mus_sta_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_sta_vertex_p4_;
	TBranch *mus_sta_vertex_p4_branch;
	bool mus_sta_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_trk_p4_;
	TBranch *mus_trk_p4_branch;
	bool mus_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_vertex_p4_;
	TBranch *mus_vertex_p4_branch;
	bool mus_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfcands_p4_;
	TBranch *pfcands_p4_branch;
	bool pfcands_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfjets_p4_;
	TBranch *pfjets_p4_branch;
	bool pfjets_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > taus_pf_lead_chargecand_p4_;
	TBranch *taus_pf_lead_chargecand_p4_branch;
	bool taus_pf_lead_chargecand_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > taus_pf_lead_neutrcand_p4_;
	TBranch *taus_pf_lead_neutrcand_p4_branch;
	bool taus_pf_lead_neutrcand_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > taus_pf_p4_;
	TBranch *taus_pf_p4_branch;
	bool taus_pf_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > photons_p4_;
	TBranch *photons_p4_branch;
	bool photons_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > vtxs_position_;
	TBranch *vtxs_position_branch;
	bool vtxs_position_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > genps_lepdaughter_p4_;
	TBranch *genps_lepdaughter_p4_branch;
	bool genps_lepdaughter_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > hlt_trigObjs_p4_;
	TBranch *hlt_trigObjs_p4_branch;
	bool hlt_trigObjs_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > taus_pf_isocands_p4_;
	TBranch *taus_pf_isocands_p4_branch;
	bool taus_pf_isocands_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > taus_pf_signalcands_p4_;
	TBranch *taus_pf_signalcands_p4_branch;
	bool taus_pf_signalcands_p4_isLoaded;
	vector<float> evt_bs_covMatrix_;
	TBranch *evt_bs_covMatrix_branch;
	bool evt_bs_covMatrix_isLoaded;
	vector<float> els_mc3dr_;
	TBranch *els_mc3dr_branch;
	bool els_mc3dr_isLoaded;
	vector<float> els_mcdr_;
	TBranch *els_mcdr_branch;
	bool els_mcdr_isLoaded;
	vector<float> mus_mc3dr_;
	TBranch *mus_mc3dr_branch;
	bool mus_mc3dr_isLoaded;
	vector<float> mus_mcdr_;
	TBranch *mus_mcdr_branch;
	bool mus_mcdr_isLoaded;
	vector<float> pfjets_mc3dr_;
	TBranch *pfjets_mc3dr_branch;
	bool pfjets_mc3dr_isLoaded;
	vector<float> pfjets_mcdr_;
	TBranch *pfjets_mcdr_branch;
	bool pfjets_mcdr_isLoaded;
	vector<float> pfjets_mc_emEnergy_;
	TBranch *pfjets_mc_emEnergy_branch;
	bool pfjets_mc_emEnergy_isLoaded;
	vector<float> pfjets_mc_gpdr_;
	TBranch *pfjets_mc_gpdr_branch;
	bool pfjets_mc_gpdr_isLoaded;
	vector<float> pfjets_mc_hadEnergy_;
	TBranch *pfjets_mc_hadEnergy_branch;
	bool pfjets_mc_hadEnergy_isLoaded;
	vector<float> pfjets_mc_invEnergy_;
	TBranch *pfjets_mc_invEnergy_branch;
	bool pfjets_mc_invEnergy_isLoaded;
	vector<float> pfjets_mc_otherEnergy_;
	TBranch *pfjets_mc_otherEnergy_branch;
	bool pfjets_mc_otherEnergy_isLoaded;
	vector<float> photons_mc3dr_;
	TBranch *photons_mc3dr_branch;
	bool photons_mc3dr_isLoaded;
	vector<float> photons_mcdr_;
	TBranch *photons_mcdr_branch;
	bool photons_mcdr_isLoaded;
	vector<float> els_chi2_;
	TBranch *els_chi2_branch;
	bool els_chi2_isLoaded;
	vector<float> els_ckf_chi2_;
	TBranch *els_ckf_chi2_branch;
	bool els_ckf_chi2_isLoaded;
	vector<float> els_ckf_ndof_;
	TBranch *els_ckf_ndof_branch;
	bool els_ckf_ndof_isLoaded;
	vector<float> els_conv_dcot_;
	TBranch *els_conv_dcot_branch;
	bool els_conv_dcot_isLoaded;
	vector<float> els_conv_dist_;
	TBranch *els_conv_dist_branch;
	bool els_conv_dist_isLoaded;
	vector<float> els_conv_old_dcot_;
	TBranch *els_conv_old_dcot_branch;
	bool els_conv_old_dcot_isLoaded;
	vector<float> els_conv_old_dist_;
	TBranch *els_conv_old_dist_branch;
	bool els_conv_old_dist_isLoaded;
	vector<float> els_conv_old_radius_;
	TBranch *els_conv_old_radius_branch;
	bool els_conv_old_radius_isLoaded;
	vector<float> els_conv_radius_;
	TBranch *els_conv_radius_branch;
	bool els_conv_radius_isLoaded;
	vector<float> els_d0_;
	TBranch *els_d0_branch;
	bool els_d0_isLoaded;
	vector<float> els_d0Err_;
	TBranch *els_d0Err_branch;
	bool els_d0Err_isLoaded;
	vector<float> els_d0corr_;
	TBranch *els_d0corr_branch;
	bool els_d0corr_isLoaded;
	vector<float> els_d0corrPhi_;
	TBranch *els_d0corrPhi_branch;
	bool els_d0corrPhi_isLoaded;
	vector<float> els_d0phiCov_;
	TBranch *els_d0phiCov_branch;
	bool els_d0phiCov_isLoaded;
	vector<float> els_dEtaIn_;
	TBranch *els_dEtaIn_branch;
	bool els_dEtaIn_isLoaded;
	vector<float> els_dEtaOut_;
	TBranch *els_dEtaOut_branch;
	bool els_dEtaOut_isLoaded;
	vector<float> els_dPhiIn_;
	TBranch *els_dPhiIn_branch;
	bool els_dPhiIn_isLoaded;
	vector<float> els_dPhiInPhiOut_;
	TBranch *els_dPhiInPhiOut_branch;
	bool els_dPhiInPhiOut_isLoaded;
	vector<float> els_dPhiOut_;
	TBranch *els_dPhiOut_branch;
	bool els_dPhiOut_isLoaded;
	vector<float> els_deltaEtaEleClusterTrackAtCalo_;
	TBranch *els_deltaEtaEleClusterTrackAtCalo_branch;
	bool els_deltaEtaEleClusterTrackAtCalo_isLoaded;
	vector<float> els_deltaPhiEleClusterTrackAtCalo_;
	TBranch *els_deltaPhiEleClusterTrackAtCalo_branch;
	bool els_deltaPhiEleClusterTrackAtCalo_isLoaded;
	vector<float> els_dxyPV_;
	TBranch *els_dxyPV_branch;
	bool els_dxyPV_isLoaded;
	vector<float> els_dzPV_;
	TBranch *els_dzPV_branch;
	bool els_dzPV_isLoaded;
	vector<float> els_e1x5_;
	TBranch *els_e1x5_branch;
	bool els_e1x5_isLoaded;
	vector<float> els_e1x5_full5x5_;
	TBranch *els_e1x5_full5x5_branch;
	bool els_e1x5_full5x5_isLoaded;
	vector<float> els_e2x5Max_;
	TBranch *els_e2x5Max_branch;
	bool els_e2x5Max_isLoaded;
	vector<float> els_e2x5Max_full5x5_;
	TBranch *els_e2x5Max_full5x5_branch;
	bool els_e2x5Max_full5x5_isLoaded;
	vector<float> els_e3x3_;
	TBranch *els_e3x3_branch;
	bool els_e3x3_isLoaded;
	vector<float> els_e5x5_;
	TBranch *els_e5x5_branch;
	bool els_e5x5_isLoaded;
	vector<float> els_e5x5_full5x5_;
	TBranch *els_e5x5_full5x5_branch;
	bool els_e5x5_full5x5_isLoaded;
	vector<float> els_eMax_;
	TBranch *els_eMax_branch;
	bool els_eMax_isLoaded;
	vector<float> els_eOverPIn_;
	TBranch *els_eOverPIn_branch;
	bool els_eOverPIn_isLoaded;
	vector<float> els_eOverPOut_;
	TBranch *els_eOverPOut_branch;
	bool els_eOverPOut_isLoaded;
	vector<float> els_eSC_;
	TBranch *els_eSC_branch;
	bool els_eSC_isLoaded;
	vector<float> els_eSCPresh_;
	TBranch *els_eSCPresh_branch;
	bool els_eSCPresh_isLoaded;
	vector<float> els_eSCRaw_;
	TBranch *els_eSCRaw_branch;
	bool els_eSCRaw_isLoaded;
	vector<float> els_eSeed_;
	TBranch *els_eSeed_branch;
	bool els_eSeed_isLoaded;
	vector<float> els_eSeedOverPIn_;
	TBranch *els_eSeedOverPIn_branch;
	bool els_eSeedOverPIn_isLoaded;
	vector<float> els_eSeedOverPOut_;
	TBranch *els_eSeedOverPOut_branch;
	bool els_eSeedOverPOut_isLoaded;
	vector<float> els_ecalEnergy_;
	TBranch *els_ecalEnergy_branch;
	bool els_ecalEnergy_isLoaded;
	vector<float> els_ecalEnergyError_;
	TBranch *els_ecalEnergyError_branch;
	bool els_ecalEnergyError_isLoaded;
	vector<float> els_ecalIso_;
	TBranch *els_ecalIso_branch;
	bool els_ecalIso_isLoaded;
	vector<float> els_ecalIso04_;
	TBranch *els_ecalIso04_branch;
	bool els_ecalIso04_isLoaded;
	vector<float> els_etaErr_;
	TBranch *els_etaErr_branch;
	bool els_etaErr_isLoaded;
	vector<float> els_etaSC_;
	TBranch *els_etaSC_branch;
	bool els_etaSC_isLoaded;
	vector<float> els_etaSCwidth_;
	TBranch *els_etaSCwidth_branch;
	bool els_etaSCwidth_isLoaded;
	vector<float> els_fbrem_;
	TBranch *els_fbrem_branch;
	bool els_fbrem_isLoaded;
	vector<float> els_hOverE_;
	TBranch *els_hOverE_branch;
	bool els_hOverE_isLoaded;
	vector<float> els_hOverEBC_;
	TBranch *els_hOverEBC_branch;
	bool els_hOverEBC_isLoaded;
	vector<float> els_hcalDepth1OverEcal_;
	TBranch *els_hcalDepth1OverEcal_branch;
	bool els_hcalDepth1OverEcal_isLoaded;
	vector<float> els_hcalDepth1TowerSumEt_;
	TBranch *els_hcalDepth1TowerSumEt_branch;
	bool els_hcalDepth1TowerSumEt_isLoaded;
	vector<float> els_hcalDepth1TowerSumEt04_;
	TBranch *els_hcalDepth1TowerSumEt04_branch;
	bool els_hcalDepth1TowerSumEt04_isLoaded;
	vector<float> els_hcalDepth2OverEcal_;
	TBranch *els_hcalDepth2OverEcal_branch;
	bool els_hcalDepth2OverEcal_isLoaded;
	vector<float> els_hcalDepth2TowerSumEt_;
	TBranch *els_hcalDepth2TowerSumEt_branch;
	bool els_hcalDepth2TowerSumEt_isLoaded;
	vector<float> els_hcalDepth2TowerSumEt04_;
	TBranch *els_hcalDepth2TowerSumEt04_branch;
	bool els_hcalDepth2TowerSumEt04_isLoaded;
	vector<float> els_hcalIso_;
	TBranch *els_hcalIso_branch;
	bool els_hcalIso_isLoaded;
	vector<float> els_hcalIso04_;
	TBranch *els_hcalIso04_branch;
	bool els_hcalIso04_isLoaded;
	vector<float> els_ip3d_;
	TBranch *els_ip3d_branch;
	bool els_ip3d_isLoaded;
	vector<float> els_ip3derr_;
	TBranch *els_ip3derr_branch;
	bool els_ip3derr_isLoaded;
	vector<float> els_iso03_pf_;
	TBranch *els_iso03_pf_branch;
	bool els_iso03_pf_isLoaded;
	vector<float> els_iso03_pf2012_ch_;
	TBranch *els_iso03_pf2012_ch_branch;
	bool els_iso03_pf2012_ch_isLoaded;
	vector<float> els_iso03_pf2012_em_;
	TBranch *els_iso03_pf2012_em_branch;
	bool els_iso03_pf2012_em_isLoaded;
	vector<float> els_iso03_pf2012ext_ch_;
	TBranch *els_iso03_pf2012ext_ch_branch;
	bool els_iso03_pf2012ext_ch_isLoaded;
	vector<float> els_iso03_pf2012ext_em_;
	TBranch *els_iso03_pf2012ext_em_branch;
	bool els_iso03_pf2012ext_em_isLoaded;
	vector<float> els_iso03_pf2012ext_nh_;
	TBranch *els_iso03_pf2012ext_nh_branch;
	bool els_iso03_pf2012ext_nh_isLoaded;
	vector<float> els_iso03_pf2012_nh_;
	TBranch *els_iso03_pf2012_nh_branch;
	bool els_iso03_pf2012_nh_isLoaded;
	vector<float> els_iso03_pf_ch_;
	TBranch *els_iso03_pf_ch_branch;
	bool els_iso03_pf_ch_isLoaded;
	vector<float> els_iso03_pf_gamma05_;
	TBranch *els_iso03_pf_gamma05_branch;
	bool els_iso03_pf_gamma05_isLoaded;
	vector<float> els_iso03_pf_nhad05_;
	TBranch *els_iso03_pf_nhad05_branch;
	bool els_iso03_pf_nhad05_isLoaded;
	vector<float> els_iso04_pf_;
	TBranch *els_iso04_pf_branch;
	bool els_iso04_pf_isLoaded;
	vector<float> els_iso04_pf2012_ch_;
	TBranch *els_iso04_pf2012_ch_branch;
	bool els_iso04_pf2012_ch_isLoaded;
	vector<float> els_iso04_pf2012_em_;
	TBranch *els_iso04_pf2012_em_branch;
	bool els_iso04_pf2012_em_isLoaded;
	vector<float> els_iso04_pf2012ext_ch_;
	TBranch *els_iso04_pf2012ext_ch_branch;
	bool els_iso04_pf2012ext_ch_isLoaded;
	vector<float> els_iso04_pf2012ext_em_;
	TBranch *els_iso04_pf2012ext_em_branch;
	bool els_iso04_pf2012ext_em_isLoaded;
	vector<float> els_iso04_pf2012ext_nh_;
	TBranch *els_iso04_pf2012ext_nh_branch;
	bool els_iso04_pf2012ext_nh_isLoaded;
	vector<float> els_iso04_pf2012_nh_;
	TBranch *els_iso04_pf2012_nh_branch;
	bool els_iso04_pf2012_nh_isLoaded;
	vector<float> els_iso04_pf_ch_;
	TBranch *els_iso04_pf_ch_branch;
	bool els_iso04_pf_ch_isLoaded;
	vector<float> els_iso04_pf_gamma05_;
	TBranch *els_iso04_pf_gamma05_branch;
	bool els_iso04_pf_gamma05_isLoaded;
	vector<float> els_iso04_pf_nhad05_;
	TBranch *els_iso04_pf_nhad05_branch;
	bool els_iso04_pf_nhad05_isLoaded;
	vector<float> els_layer1_charge_;
	TBranch *els_layer1_charge_branch;
	bool els_layer1_charge_isLoaded;
	vector<float> els_lh_;
	TBranch *els_lh_branch;
	bool els_lh_isLoaded;
	vector<float> els_mc_patMatch_dr_;
	TBranch *els_mc_patMatch_dr_branch;
	bool els_mc_patMatch_dr_isLoaded;
	vector<float> els_mva_;
	TBranch *els_mva_branch;
	bool els_mva_isLoaded;
	vector<float> els_ndof_;
	TBranch *els_ndof_branch;
	bool els_ndof_isLoaded;
	vector<float> els_pfChargedHadronIso_;
	TBranch *els_pfChargedHadronIso_branch;
	bool els_pfChargedHadronIso_isLoaded;
	vector<float> els_pfNeutralHadronIso_;
	TBranch *els_pfNeutralHadronIso_branch;
	bool els_pfNeutralHadronIso_isLoaded;
	vector<float> els_pfPUIso_;
	TBranch *els_pfPUIso_branch;
	bool els_pfPUIso_isLoaded;
	vector<float> els_pfPhotonIso_;
	TBranch *els_pfPhotonIso_branch;
	bool els_pfPhotonIso_isLoaded;
	vector<float> els_phiErr_;
	TBranch *els_phiErr_branch;
	bool els_phiErr_isLoaded;
	vector<float> els_phiSC_;
	TBranch *els_phiSC_branch;
	bool els_phiSC_isLoaded;
	vector<float> els_phiSCwidth_;
	TBranch *els_phiSCwidth_branch;
	bool els_phiSCwidth_isLoaded;
	vector<float> els_ptErr_;
	TBranch *els_ptErr_branch;
	bool els_ptErr_isLoaded;
	vector<float> els_ptErrGsf_;
	TBranch *els_ptErrGsf_branch;
	bool els_ptErrGsf_isLoaded;
	vector<float> els_r9_;
	TBranch *els_r9_branch;
	bool els_r9_isLoaded;
	vector<float> els_r9_full5x5_;
	TBranch *els_r9_full5x5_branch;
	bool els_r9_full5x5_isLoaded;
	vector<float> els_sigmaEtaEta_;
	TBranch *els_sigmaEtaEta_branch;
	bool els_sigmaEtaEta_isLoaded;
	vector<float> els_sigmaEtaEta_full5x5_;
	TBranch *els_sigmaEtaEta_full5x5_branch;
	bool els_sigmaEtaEta_full5x5_isLoaded;
	vector<float> els_sigmaIEtaIEta_;
	TBranch *els_sigmaIEtaIEta_branch;
	bool els_sigmaIEtaIEta_isLoaded;
	vector<float> els_sigmaIEtaIEtaSC_;
	TBranch *els_sigmaIEtaIEtaSC_branch;
	bool els_sigmaIEtaIEtaSC_isLoaded;
	vector<float> els_sigmaIEtaIEta_full5x5_;
	TBranch *els_sigmaIEtaIEta_full5x5_branch;
	bool els_sigmaIEtaIEta_full5x5_isLoaded;
	vector<float> els_sigmaIEtaIPhi_;
	TBranch *els_sigmaIEtaIPhi_branch;
	bool els_sigmaIEtaIPhi_isLoaded;
	vector<float> els_sigmaIPhiIPhi_;
	TBranch *els_sigmaIPhiIPhi_branch;
	bool els_sigmaIPhiIPhi_isLoaded;
	vector<float> els_sigmaIPhiIPhiSC_;
	TBranch *els_sigmaIPhiIPhiSC_branch;
	bool els_sigmaIPhiIPhiSC_isLoaded;
	vector<float> els_sigmaIPhiIPhi_full5x5_;
	TBranch *els_sigmaIPhiIPhi_full5x5_branch;
	bool els_sigmaIPhiIPhi_full5x5_isLoaded;
	vector<float> els_sigmaIphiIphi_;
	TBranch *els_sigmaIphiIphi_branch;
	bool els_sigmaIphiIphi_isLoaded;
	vector<float> els_sigmaPhiPhi_;
	TBranch *els_sigmaPhiPhi_branch;
	bool els_sigmaPhiPhi_isLoaded;
	vector<float> els_tkIso_;
	TBranch *els_tkIso_branch;
	bool els_tkIso_isLoaded;
	vector<float> els_tkIso04_;
	TBranch *els_tkIso04_branch;
	bool els_tkIso04_isLoaded;
	vector<float> els_trackMomentumError_;
	TBranch *els_trackMomentumError_branch;
	bool els_trackMomentumError_isLoaded;
	vector<float> els_trkdr_;
	TBranch *els_trkdr_branch;
	bool els_trkdr_isLoaded;
	vector<float> els_trkshFrac_;
	TBranch *els_trkshFrac_branch;
	bool els_trkshFrac_isLoaded;
	vector<float> els_z0_;
	TBranch *els_z0_branch;
	bool els_z0_isLoaded;
	vector<float> els_z0Err_;
	TBranch *els_z0Err_branch;
	bool els_z0Err_isLoaded;
	vector<float> els_z0corr_;
	TBranch *els_z0corr_branch;
	bool els_z0corr_isLoaded;
	vector<float> mus_backToBackCompat_;
	TBranch *mus_backToBackCompat_branch;
	bool mus_backToBackCompat_isLoaded;
	vector<float> mus_best_dxyPV_;
	TBranch *mus_best_dxyPV_branch;
	bool mus_best_dxyPV_isLoaded;
	vector<float> mus_best_dzPV_;
	TBranch *mus_best_dzPV_branch;
	bool mus_best_dzPV_isLoaded;
	vector<float> mus_caloCompatibility_;
	TBranch *mus_caloCompatibility_branch;
	bool mus_caloCompatibility_isLoaded;
	vector<float> mus_chi2_;
	TBranch *mus_chi2_branch;
	bool mus_chi2_isLoaded;
	vector<float> mus_chi2LocalMomentum_;
	TBranch *mus_chi2LocalMomentum_branch;
	bool mus_chi2LocalMomentum_isLoaded;
	vector<float> mus_chi2LocalPosition_;
	TBranch *mus_chi2LocalPosition_branch;
	bool mus_chi2LocalPosition_isLoaded;
	vector<float> mus_cosmicCompat_;
	TBranch *mus_cosmicCompat_branch;
	bool mus_cosmicCompat_isLoaded;
	vector<float> mus_d0_;
	TBranch *mus_d0_branch;
	bool mus_d0_isLoaded;
	vector<float> mus_d0Err_;
	TBranch *mus_d0Err_branch;
	bool mus_d0Err_isLoaded;
	vector<float> mus_d0corr_;
	TBranch *mus_d0corr_branch;
	bool mus_d0corr_isLoaded;
	vector<float> mus_d0corrPhi_;
	TBranch *mus_d0corrPhi_branch;
	bool mus_d0corrPhi_isLoaded;
	vector<float> mus_d0phiCov_;
	TBranch *mus_d0phiCov_branch;
	bool mus_d0phiCov_isLoaded;
	vector<float> mus_dxyPV_;
	TBranch *mus_dxyPV_branch;
	bool mus_dxyPV_isLoaded;
	vector<float> mus_dzPV_;
	TBranch *mus_dzPV_branch;
	bool mus_dzPV_isLoaded;
	vector<float> mus_ecal_time_;
	TBranch *mus_ecal_time_branch;
	bool mus_ecal_time_isLoaded;
	vector<float> mus_ecal_timeError_;
	TBranch *mus_ecal_timeError_branch;
	bool mus_ecal_timeError_isLoaded;
	vector<float> mus_e_em_;
	TBranch *mus_e_em_branch;
	bool mus_e_em_isLoaded;
	vector<float> mus_e_emS9_;
	TBranch *mus_e_emS9_branch;
	bool mus_e_emS9_isLoaded;
	vector<float> mus_e_had_;
	TBranch *mus_e_had_branch;
	bool mus_e_had_isLoaded;
	vector<float> mus_e_hadS9_;
	TBranch *mus_e_hadS9_branch;
	bool mus_e_hadS9_isLoaded;
	vector<float> mus_e_ho_;
	TBranch *mus_e_ho_branch;
	bool mus_e_ho_isLoaded;
	vector<float> mus_e_hoS9_;
	TBranch *mus_e_hoS9_branch;
	bool mus_e_hoS9_isLoaded;
	vector<float> mus_emMax_;
	TBranch *mus_emMax_branch;
	bool mus_emMax_isLoaded;
	vector<float> mus_emS25_;
	TBranch *mus_emS25_branch;
	bool mus_emS25_isLoaded;
	vector<float> mus_etaErr_;
	TBranch *mus_etaErr_branch;
	bool mus_etaErr_isLoaded;
	vector<float> mus_gfit_chi2_;
	TBranch *mus_gfit_chi2_branch;
	bool mus_gfit_chi2_isLoaded;
	vector<float> mus_gfit_d0_;
	TBranch *mus_gfit_d0_branch;
	bool mus_gfit_d0_isLoaded;
	vector<float> mus_gfit_d0Err_;
	TBranch *mus_gfit_d0Err_branch;
	bool mus_gfit_d0Err_isLoaded;
	vector<float> mus_gfit_d0corr_;
	TBranch *mus_gfit_d0corr_branch;
	bool mus_gfit_d0corr_isLoaded;
	vector<float> mus_gfit_d0corrPhi_;
	TBranch *mus_gfit_d0corrPhi_branch;
	bool mus_gfit_d0corrPhi_isLoaded;
	vector<float> mus_gfit_d0phiCov_;
	TBranch *mus_gfit_d0phiCov_branch;
	bool mus_gfit_d0phiCov_isLoaded;
	vector<float> mus_gfit_ndof_;
	TBranch *mus_gfit_ndof_branch;
	bool mus_gfit_ndof_isLoaded;
	vector<float> mus_gfit_qoverp_;
	TBranch *mus_gfit_qoverp_branch;
	bool mus_gfit_qoverp_isLoaded;
	vector<float> mus_gfit_qoverpError_;
	TBranch *mus_gfit_qoverpError_branch;
	bool mus_gfit_qoverpError_isLoaded;
	vector<float> mus_gfit_z0_;
	TBranch *mus_gfit_z0_branch;
	bool mus_gfit_z0_isLoaded;
	vector<float> mus_gfit_z0Err_;
	TBranch *mus_gfit_z0Err_branch;
	bool mus_gfit_z0Err_isLoaded;
	vector<float> mus_gfit_z0corr_;
	TBranch *mus_gfit_z0corr_branch;
	bool mus_gfit_z0corr_isLoaded;
	vector<float> mus_glbKink_;
	TBranch *mus_glbKink_branch;
	bool mus_glbKink_isLoaded;
	vector<float> mus_glbTrackProbability_;
	TBranch *mus_glbTrackProbability_branch;
	bool mus_glbTrackProbability_isLoaded;
	vector<float> mus_globalDeltaEtaPhi_;
	TBranch *mus_globalDeltaEtaPhi_branch;
	bool mus_globalDeltaEtaPhi_isLoaded;
	vector<float> mus_hadMax_;
	TBranch *mus_hadMax_branch;
	bool mus_hadMax_isLoaded;
	vector<float> mus_hcal_time_;
	TBranch *mus_hcal_time_branch;
	bool mus_hcal_time_isLoaded;
	vector<float> mus_hcal_timeError_;
	TBranch *mus_hcal_timeError_branch;
	bool mus_hcal_timeError_isLoaded;
	vector<float> mus_ip3d_;
	TBranch *mus_ip3d_branch;
	bool mus_ip3d_isLoaded;
	vector<float> mus_ip3derr_;
	TBranch *mus_ip3derr_branch;
	bool mus_ip3derr_isLoaded;
	vector<float> mus_iso03_emEt_;
	TBranch *mus_iso03_emEt_branch;
	bool mus_iso03_emEt_isLoaded;
	vector<float> mus_iso03_hadEt_;
	TBranch *mus_iso03_hadEt_branch;
	bool mus_iso03_hadEt_isLoaded;
	vector<float> mus_iso03_hoEt_;
	TBranch *mus_iso03_hoEt_branch;
	bool mus_iso03_hoEt_isLoaded;
	vector<float> mus_iso03_sumPt_;
	TBranch *mus_iso03_sumPt_branch;
	bool mus_iso03_sumPt_isLoaded;
	vector<float> mus_iso05_emEt_;
	TBranch *mus_iso05_emEt_branch;
	bool mus_iso05_emEt_isLoaded;
	vector<float> mus_iso05_hadEt_;
	TBranch *mus_iso05_hadEt_branch;
	bool mus_iso05_hadEt_isLoaded;
	vector<float> mus_iso05_hoEt_;
	TBranch *mus_iso05_hoEt_branch;
	bool mus_iso05_hoEt_isLoaded;
	vector<float> mus_iso05_sumPt_;
	TBranch *mus_iso05_sumPt_branch;
	bool mus_iso05_sumPt_isLoaded;
	vector<float> mus_isoMeanDRR03_pf_ChargedHadronPt_;
	TBranch *mus_isoMeanDRR03_pf_ChargedHadronPt_branch;
	bool mus_isoMeanDRR03_pf_ChargedHadronPt_isLoaded;
	vector<float> mus_isoMeanDRR03_pf_ChargedParticlePt_;
	TBranch *mus_isoMeanDRR03_pf_ChargedParticlePt_branch;
	bool mus_isoMeanDRR03_pf_ChargedParticlePt_isLoaded;
	vector<float> mus_isoMeanDRR03_pf_NeutralHadronEt_;
	TBranch *mus_isoMeanDRR03_pf_NeutralHadronEt_branch;
	bool mus_isoMeanDRR03_pf_NeutralHadronEt_isLoaded;
	vector<float> mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_;
	TBranch *mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_branch;
	bool mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_isLoaded;
	vector<float> mus_isoMeanDRR03_pf_PUPt_;
	TBranch *mus_isoMeanDRR03_pf_PUPt_branch;
	bool mus_isoMeanDRR03_pf_PUPt_isLoaded;
	vector<float> mus_isoMeanDRR03_pf_PhotonEt_;
	TBranch *mus_isoMeanDRR03_pf_PhotonEt_branch;
	bool mus_isoMeanDRR03_pf_PhotonEt_isLoaded;
	vector<float> mus_isoMeanDRR03_pf_PhotonEtHighThreshold_;
	TBranch *mus_isoMeanDRR03_pf_PhotonEtHighThreshold_branch;
	bool mus_isoMeanDRR03_pf_PhotonEtHighThreshold_isLoaded;
	vector<float> mus_isoMeanDRR04_pf_ChargedHadronPt_;
	TBranch *mus_isoMeanDRR04_pf_ChargedHadronPt_branch;
	bool mus_isoMeanDRR04_pf_ChargedHadronPt_isLoaded;
	vector<float> mus_isoMeanDRR04_pf_ChargedParticlePt_;
	TBranch *mus_isoMeanDRR04_pf_ChargedParticlePt_branch;
	bool mus_isoMeanDRR04_pf_ChargedParticlePt_isLoaded;
	vector<float> mus_isoMeanDRR04_pf_NeutralHadronEt_;
	TBranch *mus_isoMeanDRR04_pf_NeutralHadronEt_branch;
	bool mus_isoMeanDRR04_pf_NeutralHadronEt_isLoaded;
	vector<float> mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_;
	TBranch *mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_branch;
	bool mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_isLoaded;
	vector<float> mus_isoMeanDRR04_pf_PUPt_;
	TBranch *mus_isoMeanDRR04_pf_PUPt_branch;
	bool mus_isoMeanDRR04_pf_PUPt_isLoaded;
	vector<float> mus_isoMeanDRR04_pf_PhotonEt_;
	TBranch *mus_isoMeanDRR04_pf_PhotonEt_branch;
	bool mus_isoMeanDRR04_pf_PhotonEt_isLoaded;
	vector<float> mus_isoMeanDRR04_pf_PhotonEtHighThreshold_;
	TBranch *mus_isoMeanDRR04_pf_PhotonEtHighThreshold_branch;
	bool mus_isoMeanDRR04_pf_PhotonEtHighThreshold_isLoaded;
	vector<float> mus_isoR03_pf_ChargedHadronPt_;
	TBranch *mus_isoR03_pf_ChargedHadronPt_branch;
	bool mus_isoR03_pf_ChargedHadronPt_isLoaded;
	vector<float> mus_isoR03_pf_ChargedParticlePt_;
	TBranch *mus_isoR03_pf_ChargedParticlePt_branch;
	bool mus_isoR03_pf_ChargedParticlePt_isLoaded;
	vector<float> mus_isoR03_pf_NeutralHadronEt_;
	TBranch *mus_isoR03_pf_NeutralHadronEt_branch;
	bool mus_isoR03_pf_NeutralHadronEt_isLoaded;
	vector<float> mus_isoR03_pf_NeutralHadronEtHighThreshold_;
	TBranch *mus_isoR03_pf_NeutralHadronEtHighThreshold_branch;
	bool mus_isoR03_pf_NeutralHadronEtHighThreshold_isLoaded;
	vector<float> mus_isoR03_pf_PUPt_;
	TBranch *mus_isoR03_pf_PUPt_branch;
	bool mus_isoR03_pf_PUPt_isLoaded;
	vector<float> mus_isoR03_pf_PhotonEt_;
	TBranch *mus_isoR03_pf_PhotonEt_branch;
	bool mus_isoR03_pf_PhotonEt_isLoaded;
	vector<float> mus_isoR03_pf_PhotonEtHighThreshold_;
	TBranch *mus_isoR03_pf_PhotonEtHighThreshold_branch;
	bool mus_isoR03_pf_PhotonEtHighThreshold_isLoaded;
	vector<float> mus_isoR04_pf_ChargedHadronPt_;
	TBranch *mus_isoR04_pf_ChargedHadronPt_branch;
	bool mus_isoR04_pf_ChargedHadronPt_isLoaded;
	vector<float> mus_isoR04_pf_ChargedParticlePt_;
	TBranch *mus_isoR04_pf_ChargedParticlePt_branch;
	bool mus_isoR04_pf_ChargedParticlePt_isLoaded;
	vector<float> mus_isoR04_pf_NeutralHadronEt_;
	TBranch *mus_isoR04_pf_NeutralHadronEt_branch;
	bool mus_isoR04_pf_NeutralHadronEt_isLoaded;
	vector<float> mus_isoR04_pf_NeutralHadronEtHighThreshold_;
	TBranch *mus_isoR04_pf_NeutralHadronEtHighThreshold_branch;
	bool mus_isoR04_pf_NeutralHadronEtHighThreshold_isLoaded;
	vector<float> mus_isoR04_pf_PUPt_;
	TBranch *mus_isoR04_pf_PUPt_branch;
	bool mus_isoR04_pf_PUPt_isLoaded;
	vector<float> mus_isoR04_pf_PhotonEt_;
	TBranch *mus_isoR04_pf_PhotonEt_branch;
	bool mus_isoR04_pf_PhotonEt_isLoaded;
	vector<float> mus_isoR04_pf_PhotonEtHighThreshold_;
	TBranch *mus_isoR04_pf_PhotonEtHighThreshold_branch;
	bool mus_isoR04_pf_PhotonEtHighThreshold_isLoaded;
	vector<float> mus_isoSumDRR03_pf_ChargedHadronPt_;
	TBranch *mus_isoSumDRR03_pf_ChargedHadronPt_branch;
	bool mus_isoSumDRR03_pf_ChargedHadronPt_isLoaded;
	vector<float> mus_isoSumDRR03_pf_ChargedParticlePt_;
	TBranch *mus_isoSumDRR03_pf_ChargedParticlePt_branch;
	bool mus_isoSumDRR03_pf_ChargedParticlePt_isLoaded;
	vector<float> mus_isoSumDRR03_pf_NeutralHadronEt_;
	TBranch *mus_isoSumDRR03_pf_NeutralHadronEt_branch;
	bool mus_isoSumDRR03_pf_NeutralHadronEt_isLoaded;
	vector<float> mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_;
	TBranch *mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_branch;
	bool mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_isLoaded;
	vector<float> mus_isoSumDRR03_pf_PUPt_;
	TBranch *mus_isoSumDRR03_pf_PUPt_branch;
	bool mus_isoSumDRR03_pf_PUPt_isLoaded;
	vector<float> mus_isoSumDRR03_pf_PhotonEt_;
	TBranch *mus_isoSumDRR03_pf_PhotonEt_branch;
	bool mus_isoSumDRR03_pf_PhotonEt_isLoaded;
	vector<float> mus_isoSumDRR03_pf_PhotonEtHighThreshold_;
	TBranch *mus_isoSumDRR03_pf_PhotonEtHighThreshold_branch;
	bool mus_isoSumDRR03_pf_PhotonEtHighThreshold_isLoaded;
	vector<float> mus_isoSumDRR04_pf_ChargedHadronPt_;
	TBranch *mus_isoSumDRR04_pf_ChargedHadronPt_branch;
	bool mus_isoSumDRR04_pf_ChargedHadronPt_isLoaded;
	vector<float> mus_isoSumDRR04_pf_ChargedParticlePt_;
	TBranch *mus_isoSumDRR04_pf_ChargedParticlePt_branch;
	bool mus_isoSumDRR04_pf_ChargedParticlePt_isLoaded;
	vector<float> mus_isoSumDRR04_pf_NeutralHadronEt_;
	TBranch *mus_isoSumDRR04_pf_NeutralHadronEt_branch;
	bool mus_isoSumDRR04_pf_NeutralHadronEt_isLoaded;
	vector<float> mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_;
	TBranch *mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_branch;
	bool mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_isLoaded;
	vector<float> mus_isoSumDRR04_pf_PUPt_;
	TBranch *mus_isoSumDRR04_pf_PUPt_branch;
	bool mus_isoSumDRR04_pf_PUPt_isLoaded;
	vector<float> mus_isoSumDRR04_pf_PhotonEt_;
	TBranch *mus_isoSumDRR04_pf_PhotonEt_branch;
	bool mus_isoSumDRR04_pf_PhotonEt_isLoaded;
	vector<float> mus_isoSumDRR04_pf_PhotonEtHighThreshold_;
	TBranch *mus_isoSumDRR04_pf_PhotonEtHighThreshold_branch;
	bool mus_isoSumDRR04_pf_PhotonEtHighThreshold_isLoaded;
	vector<float> mus_iso_ecalvetoDep_;
	TBranch *mus_iso_ecalvetoDep_branch;
	bool mus_iso_ecalvetoDep_isLoaded;
	vector<float> mus_iso_hcalvetoDep_;
	TBranch *mus_iso_hcalvetoDep_branch;
	bool mus_iso_hcalvetoDep_isLoaded;
	vector<float> mus_iso_hovetoDep_;
	TBranch *mus_iso_hovetoDep_branch;
	bool mus_iso_hovetoDep_isLoaded;
	vector<float> mus_iso_trckvetoDep_;
	TBranch *mus_iso_trckvetoDep_branch;
	bool mus_iso_trckvetoDep_isLoaded;
	vector<float> mus_localDistance_;
	TBranch *mus_localDistance_branch;
	bool mus_localDistance_isLoaded;
	vector<float> mus_mc_patMatch_dr_;
	TBranch *mus_mc_patMatch_dr_branch;
	bool mus_mc_patMatch_dr_isLoaded;
	vector<float> mus_ndof_;
	TBranch *mus_ndof_branch;
	bool mus_ndof_isLoaded;
	vector<float> mus_overlapCompat_;
	TBranch *mus_overlapCompat_branch;
	bool mus_overlapCompat_isLoaded;
	vector<float> mus_phiErr_;
	TBranch *mus_phiErr_branch;
	bool mus_phiErr_isLoaded;
	vector<float> mus_ptErr_;
	TBranch *mus_ptErr_branch;
	bool mus_ptErr_isLoaded;
	vector<float> mus_qoverp_;
	TBranch *mus_qoverp_branch;
	bool mus_qoverp_isLoaded;
	vector<float> mus_qoverpError_;
	TBranch *mus_qoverpError_branch;
	bool mus_qoverpError_isLoaded;
	vector<float> mus_segmCompatibility_;
	TBranch *mus_segmCompatibility_branch;
	bool mus_segmCompatibility_isLoaded;
	vector<float> mus_staRelChi2_;
	TBranch *mus_staRelChi2_branch;
	bool mus_staRelChi2_isLoaded;
	vector<float> mus_sta_chi2_;
	TBranch *mus_sta_chi2_branch;
	bool mus_sta_chi2_isLoaded;
	vector<float> mus_sta_d0_;
	TBranch *mus_sta_d0_branch;
	bool mus_sta_d0_isLoaded;
	vector<float> mus_sta_d0Err_;
	TBranch *mus_sta_d0Err_branch;
	bool mus_sta_d0Err_isLoaded;
	vector<float> mus_sta_d0corr_;
	TBranch *mus_sta_d0corr_branch;
	bool mus_sta_d0corr_isLoaded;
	vector<float> mus_sta_d0corrPhi_;
	TBranch *mus_sta_d0corrPhi_branch;
	bool mus_sta_d0corrPhi_isLoaded;
	vector<float> mus_sta_d0phiCov_;
	TBranch *mus_sta_d0phiCov_branch;
	bool mus_sta_d0phiCov_isLoaded;
	vector<float> mus_sta_ndof_;
	TBranch *mus_sta_ndof_branch;
	bool mus_sta_ndof_isLoaded;
	vector<float> mus_sta_qoverp_;
	TBranch *mus_sta_qoverp_branch;
	bool mus_sta_qoverp_isLoaded;
	vector<float> mus_sta_qoverpError_;
	TBranch *mus_sta_qoverpError_branch;
	bool mus_sta_qoverpError_isLoaded;
	vector<float> mus_sta_z0_;
	TBranch *mus_sta_z0_branch;
	bool mus_sta_z0_isLoaded;
	vector<float> mus_sta_z0Err_;
	TBranch *mus_sta_z0Err_branch;
	bool mus_sta_z0Err_isLoaded;
	vector<float> mus_sta_z0corr_;
	TBranch *mus_sta_z0corr_branch;
	bool mus_sta_z0corr_isLoaded;
	vector<float> mus_timeAtIpInOut_;
	TBranch *mus_timeAtIpInOut_branch;
	bool mus_timeAtIpInOut_isLoaded;
	vector<float> mus_timeAtIpInOutErr_;
	TBranch *mus_timeAtIpInOutErr_branch;
	bool mus_timeAtIpInOutErr_isLoaded;
	vector<float> mus_timeAtIpOutIn_;
	TBranch *mus_timeAtIpOutIn_branch;
	bool mus_timeAtIpOutIn_isLoaded;
	vector<float> mus_timeAtIpOutInErr_;
	TBranch *mus_timeAtIpOutInErr_branch;
	bool mus_timeAtIpOutInErr_isLoaded;
	vector<float> mus_timeCompat_;
	TBranch *mus_timeCompat_branch;
	bool mus_timeCompat_isLoaded;
	vector<float> mus_tower_;
	TBranch *mus_tower_branch;
	bool mus_tower_isLoaded;
	vector<float> mus_towerS9_;
	TBranch *mus_towerS9_branch;
	bool mus_towerS9_isLoaded;
	vector<float> mus_trkKink_;
	TBranch *mus_trkKink_branch;
	bool mus_trkKink_isLoaded;
	vector<float> mus_trkRelChi2_;
	TBranch *mus_trkRelChi2_branch;
	bool mus_trkRelChi2_isLoaded;
	vector<float> mus_vertexCompat_;
	TBranch *mus_vertexCompat_branch;
	bool mus_vertexCompat_isLoaded;
	vector<float> mus_vertexphi_;
	TBranch *mus_vertexphi_branch;
	bool mus_vertexphi_isLoaded;
	vector<float> mus_z0_;
	TBranch *mus_z0_branch;
	bool mus_z0_isLoaded;
	vector<float> mus_z0Err_;
	TBranch *mus_z0Err_branch;
	bool mus_z0Err_isLoaded;
	vector<float> mus_z0corr_;
	TBranch *mus_z0corr_branch;
	bool mus_z0corr_isLoaded;
	vector<float> pfcands_dz_;
	TBranch *pfcands_dz_branch;
	bool pfcands_dz_isLoaded;
	vector<float> pfjets_area_;
	TBranch *pfjets_area_branch;
	bool pfjets_area_isLoaded;
	vector<float> pfjets_chargedEmE_;
	TBranch *pfjets_chargedEmE_branch;
	bool pfjets_chargedEmE_isLoaded;
	vector<float> pfjets_chargedHadronE_;
	TBranch *pfjets_chargedHadronE_branch;
	bool pfjets_chargedHadronE_isLoaded;
	vector<float> pfjets_combinedSecondaryVertexBJetTag_;
	TBranch *pfjets_combinedSecondaryVertexBJetTag_branch;
	bool pfjets_combinedSecondaryVertexBJetTag_isLoaded;
	vector<float> pfjets_cor_;
	TBranch *pfjets_cor_branch;
	bool pfjets_cor_isLoaded;
	vector<float> pfjets_corL1Fast_;
	TBranch *pfjets_corL1Fast_branch;
	bool pfjets_corL1Fast_isLoaded;
	vector<float> pfjets_corL1FastL2L3_;
	TBranch *pfjets_corL1FastL2L3_branch;
	bool pfjets_corL1FastL2L3_isLoaded;
	vector<float> pfjets_corL1FastL2L3residual_;
	TBranch *pfjets_corL1FastL2L3residual_branch;
	bool pfjets_corL1FastL2L3residual_isLoaded;
	vector<float> pfjets_corL1L2L3_;
	TBranch *pfjets_corL1L2L3_branch;
	bool pfjets_corL1L2L3_isLoaded;
	vector<float> pfjets_electronE_;
	TBranch *pfjets_electronE_branch;
	bool pfjets_electronE_isLoaded;
	vector<float> pfjets_hfEmE_;
	TBranch *pfjets_hfEmE_branch;
	bool pfjets_hfEmE_isLoaded;
	vector<float> pfjets_hfHadronE_;
	TBranch *pfjets_hfHadronE_branch;
	bool pfjets_hfHadronE_isLoaded;
	vector<float> pfjets_jetBProbabilityBJetTag_;
	TBranch *pfjets_jetBProbabilityBJetTag_branch;
	bool pfjets_jetBProbabilityBJetTag_isLoaded;
	vector<float> pfjets_jetProbabilityBJetTag_;
	TBranch *pfjets_jetProbabilityBJetTag_branch;
	bool pfjets_jetProbabilityBJetTag_isLoaded;
	vector<float> pfjets_muonE_;
	TBranch *pfjets_muonE_branch;
	bool pfjets_muonE_isLoaded;
	vector<float> pfjets_neutralEmE_;
	TBranch *pfjets_neutralEmE_branch;
	bool pfjets_neutralEmE_isLoaded;
	vector<float> pfjets_neutralHadronE_;
	TBranch *pfjets_neutralHadronE_branch;
	bool pfjets_neutralHadronE_isLoaded;
	vector<float> pfjets_photonE_;
	TBranch *pfjets_photonE_branch;
	bool pfjets_photonE_isLoaded;
	vector<float> pfjets_pileupJetId_;
	TBranch *pfjets_pileupJetId_branch;
	bool pfjets_pileupJetId_isLoaded;
	vector<float> pfjets_simpleSecondaryVertexHighEffBJetTag_;
	TBranch *pfjets_simpleSecondaryVertexHighEffBJetTag_branch;
	bool pfjets_simpleSecondaryVertexHighEffBJetTag_isLoaded;
	vector<float> pfjets_simpleSecondaryVertexHighPurBJetTags_;
	TBranch *pfjets_simpleSecondaryVertexHighPurBJetTags_branch;
	bool pfjets_simpleSecondaryVertexHighPurBJetTags_isLoaded;
	vector<float> pfjets_trackCountingHighEffBJetTag_;
	TBranch *pfjets_trackCountingHighEffBJetTag_branch;
	bool pfjets_trackCountingHighEffBJetTag_isLoaded;
	vector<float> pfjets_trackCountingHighPurBJetTag_;
	TBranch *pfjets_trackCountingHighPurBJetTag_branch;
	bool pfjets_trackCountingHighPurBJetTag_isLoaded;
	vector<float> taus_pf_againstElectronDeadECAL_;
	TBranch *taus_pf_againstElectronDeadECAL_branch;
	bool taus_pf_againstElectronDeadECAL_isLoaded;
	vector<float> taus_pf_againstElectronLoose_;
	TBranch *taus_pf_againstElectronLoose_branch;
	bool taus_pf_againstElectronLoose_isLoaded;
	vector<float> taus_pf_againstElectronMedium_;
	TBranch *taus_pf_againstElectronMedium_branch;
	bool taus_pf_againstElectronMedium_isLoaded;
	vector<float> taus_pf_againstElectronTight_;
	TBranch *taus_pf_againstElectronTight_branch;
	bool taus_pf_againstElectronTight_isLoaded;
	vector<float> taus_pf_againstMuonLoose_;
	TBranch *taus_pf_againstMuonLoose_branch;
	bool taus_pf_againstMuonLoose_isLoaded;
	vector<float> taus_pf_againstMuonLoose2_;
	TBranch *taus_pf_againstMuonLoose2_branch;
	bool taus_pf_againstMuonLoose2_isLoaded;
	vector<float> taus_pf_againstMuonLoose3_;
	TBranch *taus_pf_againstMuonLoose3_branch;
	bool taus_pf_againstMuonLoose3_isLoaded;
	vector<float> taus_pf_againstMuonMedium_;
	TBranch *taus_pf_againstMuonMedium_branch;
	bool taus_pf_againstMuonMedium_isLoaded;
	vector<float> taus_pf_againstMuonMedium2_;
	TBranch *taus_pf_againstMuonMedium2_branch;
	bool taus_pf_againstMuonMedium2_isLoaded;
	vector<float> taus_pf_againstMuonTight_;
	TBranch *taus_pf_againstMuonTight_branch;
	bool taus_pf_againstMuonTight_isLoaded;
	vector<float> taus_pf_againstMuonTight2_;
	TBranch *taus_pf_againstMuonTight2_branch;
	bool taus_pf_againstMuonTight2_isLoaded;
	vector<float> taus_pf_againstMuonTight3_;
	TBranch *taus_pf_againstMuonTight3_branch;
	bool taus_pf_againstMuonTight3_isLoaded;
	vector<float> taus_pf_byCombinedIsolationDeltaBetaCorrRaw_;
	TBranch *taus_pf_byCombinedIsolationDeltaBetaCorrRaw_branch;
	bool taus_pf_byCombinedIsolationDeltaBetaCorrRaw_isLoaded;
	vector<float> taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_;
	TBranch *taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_branch;
	bool taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_isLoaded;
	vector<float> taus_pf_byDecayModeFinding_;
	TBranch *taus_pf_byDecayModeFinding_branch;
	bool taus_pf_byDecayModeFinding_isLoaded;
	vector<float> taus_pf_byLooseCombinedIsolationDeltaBetaCorr_;
	TBranch *taus_pf_byLooseCombinedIsolationDeltaBetaCorr_branch;
	bool taus_pf_byLooseCombinedIsolationDeltaBetaCorr_isLoaded;
	vector<float> taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_;
	TBranch *taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_branch;
	bool taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_isLoaded;
	vector<float> taus_pf_byMediumCombinedIsolationDeltaBetaCorr_;
	TBranch *taus_pf_byMediumCombinedIsolationDeltaBetaCorr_branch;
	bool taus_pf_byMediumCombinedIsolationDeltaBetaCorr_isLoaded;
	vector<float> taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_;
	TBranch *taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_branch;
	bool taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_isLoaded;
	vector<float> taus_pf_byTightCombinedIsolationDeltaBetaCorr_;
	TBranch *taus_pf_byTightCombinedIsolationDeltaBetaCorr_branch;
	bool taus_pf_byTightCombinedIsolationDeltaBetaCorr_isLoaded;
	vector<float> taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_;
	TBranch *taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_branch;
	bool taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_isLoaded;
	vector<float> taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_;
	TBranch *taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_branch;
	bool taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_isLoaded;
	vector<float> photons_e1x5_;
	TBranch *photons_e1x5_branch;
	bool photons_e1x5_isLoaded;
	vector<float> photons_e2x5Max_;
	TBranch *photons_e2x5Max_branch;
	bool photons_e2x5Max_isLoaded;
	vector<float> photons_e3x3_;
	TBranch *photons_e3x3_branch;
	bool photons_e3x3_isLoaded;
	vector<float> photons_e5x5_;
	TBranch *photons_e5x5_branch;
	bool photons_e5x5_isLoaded;
	vector<float> photons_ecalIso03_;
	TBranch *photons_ecalIso03_branch;
	bool photons_ecalIso03_isLoaded;
	vector<float> photons_ecalIso04_;
	TBranch *photons_ecalIso04_branch;
	bool photons_ecalIso04_isLoaded;
	vector<float> photons_full3x3_e3x3_;
	TBranch *photons_full3x3_e3x3_branch;
	bool photons_full3x3_e3x3_isLoaded;
	vector<float> photons_full5x5_e1x5_;
	TBranch *photons_full5x5_e1x5_branch;
	bool photons_full5x5_e1x5_isLoaded;
	vector<float> photons_full5x5_e2x5Max_;
	TBranch *photons_full5x5_e2x5Max_branch;
	bool photons_full5x5_e2x5Max_isLoaded;
	vector<float> photons_full5x5_e5x5_;
	TBranch *photons_full5x5_e5x5_branch;
	bool photons_full5x5_e5x5_isLoaded;
	vector<float> photons_full5x5_hOverE_;
	TBranch *photons_full5x5_hOverE_branch;
	bool photons_full5x5_hOverE_isLoaded;
	vector<float> photons_full5x5_sigmaEtaEta_;
	TBranch *photons_full5x5_sigmaEtaEta_branch;
	bool photons_full5x5_sigmaEtaEta_isLoaded;
	vector<float> photons_full5x5_sigmaIEtaIEta_;
	TBranch *photons_full5x5_sigmaIEtaIEta_branch;
	bool photons_full5x5_sigmaIEtaIEta_isLoaded;
	vector<float> photons_hOverE_;
	TBranch *photons_hOverE_branch;
	bool photons_hOverE_isLoaded;
	vector<float> photons_hcalDepth1TowerSumEtBcConeDR03_;
	TBranch *photons_hcalDepth1TowerSumEtBcConeDR03_branch;
	bool photons_hcalDepth1TowerSumEtBcConeDR03_isLoaded;
	vector<float> photons_hcalDepth1TowerSumEtBcConeDR04_;
	TBranch *photons_hcalDepth1TowerSumEtBcConeDR04_branch;
	bool photons_hcalDepth1TowerSumEtBcConeDR04_isLoaded;
	vector<float> photons_hcalDepth2TowerSumEtBcConeDR03_;
	TBranch *photons_hcalDepth2TowerSumEtBcConeDR03_branch;
	bool photons_hcalDepth2TowerSumEtBcConeDR03_isLoaded;
	vector<float> photons_hcalDepth2TowerSumEtBcConeDR04_;
	TBranch *photons_hcalDepth2TowerSumEtBcConeDR04_branch;
	bool photons_hcalDepth2TowerSumEtBcConeDR04_isLoaded;
	vector<float> photons_hcalIso03_;
	TBranch *photons_hcalIso03_branch;
	bool photons_hcalIso03_isLoaded;
	vector<float> photons_hcalIso04_;
	TBranch *photons_hcalIso04_branch;
	bool photons_hcalIso04_isLoaded;
	vector<float> photons_hcalTowerSumEtBcConeDR03_;
	TBranch *photons_hcalTowerSumEtBcConeDR03_branch;
	bool photons_hcalTowerSumEtBcConeDR03_isLoaded;
	vector<float> photons_hcalTowerSumEtBcConeDR04_;
	TBranch *photons_hcalTowerSumEtBcConeDR04_branch;
	bool photons_hcalTowerSumEtBcConeDR04_isLoaded;
	vector<float> photons_ntkIsoHollow03_;
	TBranch *photons_ntkIsoHollow03_branch;
	bool photons_ntkIsoHollow03_isLoaded;
	vector<float> photons_ntkIsoHollow04_;
	TBranch *photons_ntkIsoHollow04_branch;
	bool photons_ntkIsoHollow04_isLoaded;
	vector<float> photons_ntkIsoSolid03_;
	TBranch *photons_ntkIsoSolid03_branch;
	bool photons_ntkIsoSolid03_isLoaded;
	vector<float> photons_ntkIsoSolid04_;
	TBranch *photons_ntkIsoSolid04_branch;
	bool photons_ntkIsoSolid04_isLoaded;
	vector<float> photons_sigmaEtaEta_;
	TBranch *photons_sigmaEtaEta_branch;
	bool photons_sigmaEtaEta_isLoaded;
	vector<float> photons_sigmaIEtaIEta_;
	TBranch *photons_sigmaIEtaIEta_branch;
	bool photons_sigmaIEtaIEta_isLoaded;
	vector<float> photons_tkIsoHollow03_;
	TBranch *photons_tkIsoHollow03_branch;
	bool photons_tkIsoHollow03_isLoaded;
	vector<float> photons_tkIsoHollow04_;
	TBranch *photons_tkIsoHollow04_branch;
	bool photons_tkIsoHollow04_isLoaded;
	vector<float> photons_tkIsoSolid03_;
	TBranch *photons_tkIsoSolid03_branch;
	bool photons_tkIsoSolid03_isLoaded;
	vector<float> photons_tkIsoSolid04_;
	TBranch *photons_tkIsoSolid04_branch;
	bool photons_tkIsoSolid04_isLoaded;
	vector<float> puInfo_trueNumInteractions_;
	TBranch *puInfo_trueNumInteractions_branch;
	bool puInfo_trueNumInteractions_isLoaded;
	vector<float> convs_chi2_;
	TBranch *convs_chi2_branch;
	bool convs_chi2_isLoaded;
	vector<float> convs_dl_;
	TBranch *convs_dl_branch;
	bool convs_dl_isLoaded;
	vector<float> convs_ndof_;
	TBranch *convs_ndof_branch;
	bool convs_ndof_isLoaded;
	vector<float> vtxs_chi2_;
	TBranch *vtxs_chi2_branch;
	bool vtxs_chi2_isLoaded;
	vector<float> vtxs_ndof_;
	TBranch *vtxs_ndof_branch;
	bool vtxs_ndof_isLoaded;
	vector<float> vtxs_sumpt_;
	TBranch *vtxs_sumpt_branch;
	bool vtxs_sumpt_isLoaded;
	vector<float> vtxs_xError_;
	TBranch *vtxs_xError_branch;
	bool vtxs_xError_isLoaded;
	vector<float> vtxs_yError_;
	TBranch *vtxs_yError_branch;
	bool vtxs_yError_isLoaded;
	vector<float> vtxs_zError_;
	TBranch *vtxs_zError_branch;
	bool vtxs_zError_isLoaded;
	vector<vector<float> > els_convs_dcot_;
	TBranch *els_convs_dcot_branch;
	bool els_convs_dcot_isLoaded;
	vector<vector<float> > els_convs_dist_;
	TBranch *els_convs_dist_branch;
	bool els_convs_dist_isLoaded;
	vector<vector<float> > els_convs_radius_;
	TBranch *els_convs_radius_branch;
	bool els_convs_radius_isLoaded;
	vector<vector<float> > mus_stationShowerDeltaR_;
	TBranch *mus_stationShowerDeltaR_branch;
	bool mus_stationShowerDeltaR_isLoaded;
	vector<vector<float> > mus_stationShowerSizeT_;
	TBranch *mus_stationShowerSizeT_branch;
	bool mus_stationShowerSizeT_isLoaded;
	vector<vector<float> > puInfo_instLumi_;
	TBranch *puInfo_instLumi_branch;
	bool puInfo_instLumi_isLoaded;
	vector<vector<float> > vtxs_covMatrix_;
	TBranch *vtxs_covMatrix_branch;
	bool vtxs_covMatrix_isLoaded;
	int evt_bsType_;
	TBranch *evt_bsType_branch;
	bool evt_bsType_isLoaded;
	int evt_bunchCrossing_;
	TBranch *evt_bunchCrossing_branch;
	bool evt_bunchCrossing_isLoaded;
	int evt_experimentType_;
	TBranch *evt_experimentType_branch;
	bool evt_experimentType_isLoaded;
	int evt_isRealData_;
	TBranch *evt_isRealData_branch;
	bool evt_isRealData_isLoaded;
	int evt_orbitNumber_;
	TBranch *evt_orbitNumber_branch;
	bool evt_orbitNumber_isLoaded;
	int evt_storeNumber_;
	TBranch *evt_storeNumber_branch;
	bool evt_storeNumber_isLoaded;
	int hcalnoise_GetRecHitCount_;
	TBranch *hcalnoise_GetRecHitCount_branch;
	bool hcalnoise_GetRecHitCount_isLoaded;
	int hcalnoise_GetRecHitCount15_;
	TBranch *hcalnoise_GetRecHitCount15_branch;
	bool hcalnoise_GetRecHitCount15_isLoaded;
	int hcalnoise_maxHPDHits_;
	TBranch *hcalnoise_maxHPDHits_branch;
	bool hcalnoise_maxHPDHits_isLoaded;
	int hcalnoise_maxHPDNoOtherHits_;
	TBranch *hcalnoise_maxHPDNoOtherHits_branch;
	bool hcalnoise_maxHPDNoOtherHits_isLoaded;
	int hcalnoise_maxRBXHits_;
	TBranch *hcalnoise_maxRBXHits_branch;
	bool hcalnoise_maxRBXHits_isLoaded;
	int hcalnoise_maxZeros_;
	TBranch *hcalnoise_maxZeros_branch;
	bool hcalnoise_maxZeros_isLoaded;
	int hcalnoise_noiseFilterStatus_;
	TBranch *hcalnoise_noiseFilterStatus_branch;
	bool hcalnoise_noiseFilterStatus_isLoaded;
	int hcalnoise_noiseType_;
	TBranch *hcalnoise_noiseType_branch;
	bool hcalnoise_noiseType_isLoaded;
	int hcalnoise_num10GeVHits_;
	TBranch *hcalnoise_num10GeVHits_branch;
	bool hcalnoise_num10GeVHits_isLoaded;
	int hcalnoise_num25GeVHits_;
	TBranch *hcalnoise_num25GeVHits_branch;
	bool hcalnoise_num25GeVHits_isLoaded;
	int hcalnoise_numFlatNoiseChannels_;
	TBranch *hcalnoise_numFlatNoiseChannels_branch;
	bool hcalnoise_numFlatNoiseChannels_isLoaded;
	int hcalnoise_numIsolatedNoiseChannels_;
	TBranch *hcalnoise_numIsolatedNoiseChannels_branch;
	bool hcalnoise_numIsolatedNoiseChannels_isLoaded;
	int hcalnoise_numProblematicRBXs_;
	TBranch *hcalnoise_numProblematicRBXs_branch;
	bool hcalnoise_numProblematicRBXs_isLoaded;
	int hcalnoise_numSpikeNoiseChannels_;
	TBranch *hcalnoise_numSpikeNoiseChannels_branch;
	bool hcalnoise_numSpikeNoiseChannels_isLoaded;
	int hcalnoise_numTS4TS5NoiseChannels_;
	TBranch *hcalnoise_numTS4TS5NoiseChannels_branch;
	bool hcalnoise_numTS4TS5NoiseChannels_isLoaded;
	int hcalnoise_numTriangleNoiseChannels_;
	TBranch *hcalnoise_numTriangleNoiseChannels_branch;
	bool hcalnoise_numTriangleNoiseChannels_isLoaded;
	int hcalnoise_passHighLevelNoiseFilter_;
	TBranch *hcalnoise_passHighLevelNoiseFilter_branch;
	bool hcalnoise_passHighLevelNoiseFilter_isLoaded;
	int hcalnoise_passLooseNoiseFilter_;
	TBranch *hcalnoise_passLooseNoiseFilter_branch;
	bool hcalnoise_passLooseNoiseFilter_isLoaded;
	int hcalnoise_passTightNoiseFilter_;
	TBranch *hcalnoise_passTightNoiseFilter_branch;
	bool hcalnoise_passTightNoiseFilter_isLoaded;
	int pdfinfo_id1_;
	TBranch *pdfinfo_id1_branch;
	bool pdfinfo_id1_isLoaded;
	int pdfinfo_id2_;
	TBranch *pdfinfo_id2_branch;
	bool pdfinfo_id2_isLoaded;
	vector<int> els_mc3_id_;
	TBranch *els_mc3_id_branch;
	bool els_mc3_id_isLoaded;
	vector<int> els_mc3idx_;
	TBranch *els_mc3idx_branch;
	bool els_mc3idx_isLoaded;
	vector<int> els_mc3_motherid_;
	TBranch *els_mc3_motherid_branch;
	bool els_mc3_motherid_isLoaded;
	vector<int> els_mc3_motheridx_;
	TBranch *els_mc3_motheridx_branch;
	bool els_mc3_motheridx_isLoaded;
	vector<int> els_mc_id_;
	TBranch *els_mc_id_branch;
	bool els_mc_id_isLoaded;
	vector<int> els_mcidx_;
	TBranch *els_mcidx_branch;
	bool els_mcidx_isLoaded;
	vector<int> els_mc_motherid_;
	TBranch *els_mc_motherid_branch;
	bool els_mc_motherid_isLoaded;
	vector<int> mus_mc3_id_;
	TBranch *mus_mc3_id_branch;
	bool mus_mc3_id_isLoaded;
	vector<int> mus_mc3idx_;
	TBranch *mus_mc3idx_branch;
	bool mus_mc3idx_isLoaded;
	vector<int> mus_mc3_motherid_;
	TBranch *mus_mc3_motherid_branch;
	bool mus_mc3_motherid_isLoaded;
	vector<int> mus_mc3_motheridx_;
	TBranch *mus_mc3_motheridx_branch;
	bool mus_mc3_motheridx_isLoaded;
	vector<int> mus_mc_id_;
	TBranch *mus_mc_id_branch;
	bool mus_mc_id_isLoaded;
	vector<int> mus_mcidx_;
	TBranch *mus_mcidx_branch;
	bool mus_mcidx_isLoaded;
	vector<int> mus_mc_motherid_;
	TBranch *mus_mc_motherid_branch;
	bool mus_mc_motherid_isLoaded;
	vector<int> pfjets_mc3_id_;
	TBranch *pfjets_mc3_id_branch;
	bool pfjets_mc3_id_isLoaded;
	vector<int> pfjets_mc3idx_;
	TBranch *pfjets_mc3idx_branch;
	bool pfjets_mc3idx_isLoaded;
	vector<int> pfjets_mc_gpidx_;
	TBranch *pfjets_mc_gpidx_branch;
	bool pfjets_mc_gpidx_isLoaded;
	vector<int> pfjets_mc_id_;
	TBranch *pfjets_mc_id_branch;
	bool pfjets_mc_id_isLoaded;
	vector<int> pfjets_mcidx_;
	TBranch *pfjets_mcidx_branch;
	bool pfjets_mcidx_isLoaded;
	vector<int> pfjets_mc_motherid_;
	TBranch *pfjets_mc_motherid_branch;
	bool pfjets_mc_motherid_isLoaded;
	vector<int> photons_mc3_id_;
	TBranch *photons_mc3_id_branch;
	bool photons_mc3_id_isLoaded;
	vector<int> photons_mc3idx_;
	TBranch *photons_mc3idx_branch;
	bool photons_mc3idx_isLoaded;
	vector<int> photons_mc3_motherid_;
	TBranch *photons_mc3_motherid_branch;
	bool photons_mc3_motherid_isLoaded;
	vector<int> photons_mc3_motheridx_;
	TBranch *photons_mc3_motheridx_branch;
	bool photons_mc3_motheridx_isLoaded;
	vector<int> photons_mc_id_;
	TBranch *photons_mc_id_branch;
	bool photons_mc_id_isLoaded;
	vector<int> photons_mcidx_;
	TBranch *photons_mcidx_branch;
	bool photons_mcidx_isLoaded;
	vector<int> photons_mc_motherid_;
	TBranch *photons_mc_motherid_branch;
	bool photons_mc_motherid_isLoaded;
	vector<int> els_category_;
	TBranch *els_category_branch;
	bool els_category_isLoaded;
	vector<int> els_charge_;
	TBranch *els_charge_branch;
	bool els_charge_isLoaded;
	vector<int> els_ckf_charge_;
	TBranch *els_ckf_charge_branch;
	bool els_ckf_charge_isLoaded;
	vector<int> els_ckf_laywithmeas_;
	TBranch *els_ckf_laywithmeas_branch;
	bool els_ckf_laywithmeas_isLoaded;
	vector<int> els_class_;
	TBranch *els_class_branch;
	bool els_class_isLoaded;
	vector<int> els_conv_delMissHits_;
	TBranch *els_conv_delMissHits_branch;
	bool els_conv_delMissHits_isLoaded;
	vector<int> els_conv_flag_;
	TBranch *els_conv_flag_branch;
	bool els_conv_flag_isLoaded;
	vector<int> els_conv_gsftkidx_;
	TBranch *els_conv_gsftkidx_branch;
	bool els_conv_gsftkidx_isLoaded;
	vector<int> els_conv_old_delMissHits_;
	TBranch *els_conv_old_delMissHits_branch;
	bool els_conv_old_delMissHits_isLoaded;
	vector<int> els_conv_old_flag_;
	TBranch *els_conv_old_flag_branch;
	bool els_conv_old_flag_isLoaded;
	vector<int> els_conv_old_gsftkidx_;
	TBranch *els_conv_old_gsftkidx_branch;
	bool els_conv_old_gsftkidx_isLoaded;
	vector<int> els_conv_old_tkidx_;
	TBranch *els_conv_old_tkidx_branch;
	bool els_conv_old_tkidx_isLoaded;
	vector<int> els_conv_tkidx_;
	TBranch *els_conv_tkidx_branch;
	bool els_conv_tkidx_isLoaded;
	vector<int> els_exp_innerlayers_;
	TBranch *els_exp_innerlayers_branch;
	bool els_exp_innerlayers_isLoaded;
	vector<int> els_exp_outerlayers_;
	TBranch *els_exp_outerlayers_branch;
	bool els_exp_outerlayers_isLoaded;
	vector<int> els_fiduciality_;
	TBranch *els_fiduciality_branch;
	bool els_fiduciality_isLoaded;
	vector<int> els_gsftrkidx_;
	TBranch *els_gsftrkidx_branch;
	bool els_gsftrkidx_isLoaded;
	vector<int> els_layer1_det_;
	TBranch *els_layer1_det_branch;
	bool els_layer1_det_isLoaded;
	vector<int> els_layer1_layer_;
	TBranch *els_layer1_layer_branch;
	bool els_layer1_layer_isLoaded;
	vector<int> els_layer1_sizerphi_;
	TBranch *els_layer1_sizerphi_branch;
	bool els_layer1_sizerphi_isLoaded;
	vector<int> els_layer1_sizerz_;
	TBranch *els_layer1_sizerz_branch;
	bool els_layer1_sizerz_isLoaded;
	vector<int> els_lostHits_;
	TBranch *els_lostHits_branch;
	bool els_lostHits_isLoaded;
	vector<int> els_lost_pixelhits_;
	TBranch *els_lost_pixelhits_branch;
	bool els_lost_pixelhits_isLoaded;
	vector<int> els_mc_patMatch_id_;
	TBranch *els_mc_patMatch_id_branch;
	bool els_mc_patMatch_id_isLoaded;
	vector<int> els_nSeed_;
	TBranch *els_nSeed_branch;
	bool els_nSeed_isLoaded;
	vector<int> els_nlayers_;
	TBranch *els_nlayers_branch;
	bool els_nlayers_isLoaded;
	vector<int> els_nlayers3D_;
	TBranch *els_nlayers3D_branch;
	bool els_nlayers3D_isLoaded;
	vector<int> els_nlayersLost_;
	TBranch *els_nlayersLost_branch;
	bool els_nlayersLost_isLoaded;
	vector<int> els_sccharge_;
	TBranch *els_sccharge_branch;
	bool els_sccharge_isLoaded;
	vector<int> els_scindex_;
	TBranch *els_scindex_branch;
	bool els_scindex_isLoaded;
	vector<int> els_trk_charge_;
	TBranch *els_trk_charge_branch;
	bool els_trk_charge_isLoaded;
	vector<int> els_trkidx_;
	TBranch *els_trkidx_branch;
	bool els_trkidx_isLoaded;
	vector<int> els_type_;
	TBranch *els_type_branch;
	bool els_type_isLoaded;
	vector<int> els_validHits_;
	TBranch *els_validHits_branch;
	bool els_validHits_isLoaded;
	vector<int> els_valid_pixelhits_;
	TBranch *els_valid_pixelhits_branch;
	bool els_valid_pixelhits_isLoaded;
	vector<int> genps_id_;
	TBranch *genps_id_branch;
	bool genps_id_isLoaded;
	vector<int> genps_id_mother_;
	TBranch *genps_id_mother_branch;
	bool genps_id_mother_isLoaded;
	vector<int> genps_status_;
	TBranch *genps_status_branch;
	bool genps_status_isLoaded;
	vector<int> mus_algo_;
	TBranch *mus_algo_branch;
	bool mus_algo_isLoaded;
	vector<int> mus_charge_;
	TBranch *mus_charge_branch;
	bool mus_charge_isLoaded;
	vector<int> mus_ecal_rawId_;
	TBranch *mus_ecal_rawId_branch;
	bool mus_ecal_rawId_isLoaded;
	vector<int> mus_exp_innerlayers_;
	TBranch *mus_exp_innerlayers_branch;
	bool mus_exp_innerlayers_isLoaded;
	vector<int> mus_exp_outerlayers_;
	TBranch *mus_exp_outerlayers_branch;
	bool mus_exp_outerlayers_isLoaded;
	vector<int> mus_gfit_algo_;
	TBranch *mus_gfit_algo_branch;
	bool mus_gfit_algo_isLoaded;
	vector<int> mus_gfit_etaErr_;
	TBranch *mus_gfit_etaErr_branch;
	bool mus_gfit_etaErr_isLoaded;
	vector<int> mus_gfit_exp_innerlayer_;
	TBranch *mus_gfit_exp_innerlayer_branch;
	bool mus_gfit_exp_innerlayer_isLoaded;
	vector<int> mus_gfit_exp_outerlayer_;
	TBranch *mus_gfit_exp_outerlayer_branch;
	bool mus_gfit_exp_outerlayer_isLoaded;
	vector<int> mus_gfit_lostHits_;
	TBranch *mus_gfit_lostHits_branch;
	bool mus_gfit_lostHits_isLoaded;
	vector<int> mus_gfit_lostPixelHits_;
	TBranch *mus_gfit_lostPixelHits_branch;
	bool mus_gfit_lostPixelHits_isLoaded;
	vector<int> mus_gfit_nlayers_;
	TBranch *mus_gfit_nlayers_branch;
	bool mus_gfit_nlayers_isLoaded;
	vector<int> mus_gfit_nlayers3D_;
	TBranch *mus_gfit_nlayers3D_branch;
	bool mus_gfit_nlayers3D_isLoaded;
	vector<int> mus_gfit_nlayersLost_;
	TBranch *mus_gfit_nlayersLost_branch;
	bool mus_gfit_nlayersLost_isLoaded;
	vector<int> mus_gfit_phiErr_;
	TBranch *mus_gfit_phiErr_branch;
	bool mus_gfit_phiErr_isLoaded;
	vector<int> mus_gfit_ptErr_;
	TBranch *mus_gfit_ptErr_branch;
	bool mus_gfit_ptErr_isLoaded;
	vector<int> mus_gfit_qualityMask_;
	TBranch *mus_gfit_qualityMask_branch;
	bool mus_gfit_qualityMask_isLoaded;
	vector<int> mus_gfit_trk_charge_;
	TBranch *mus_gfit_trk_charge_branch;
	bool mus_gfit_trk_charge_isLoaded;
	vector<int> mus_gfit_validHits_;
	TBranch *mus_gfit_validHits_branch;
	bool mus_gfit_validHits_isLoaded;
	vector<int> mus_gfit_validPixelHits_;
	TBranch *mus_gfit_validPixelHits_branch;
	bool mus_gfit_validPixelHits_isLoaded;
	vector<int> mus_gfit_validSTAHits_;
	TBranch *mus_gfit_validSTAHits_branch;
	bool mus_gfit_validSTAHits_isLoaded;
	vector<int> mus_gfit_validSiHits_;
	TBranch *mus_gfit_validSiHits_branch;
	bool mus_gfit_validSiHits_isLoaded;
	vector<int> mus_goodmask_;
	TBranch *mus_goodmask_branch;
	bool mus_goodmask_isLoaded;
	vector<int> mus_hcal_rawId_;
	TBranch *mus_hcal_rawId_branch;
	bool mus_hcal_rawId_isLoaded;
	vector<int> mus_iso03_njets_;
	TBranch *mus_iso03_njets_branch;
	bool mus_iso03_njets_isLoaded;
	vector<int> mus_iso03_ntrk_;
	TBranch *mus_iso03_ntrk_branch;
	bool mus_iso03_ntrk_isLoaded;
	vector<int> mus_iso05_ntrk_;
	TBranch *mus_iso05_ntrk_branch;
	bool mus_iso05_ntrk_isLoaded;
	vector<int> mus_lostHits_;
	TBranch *mus_lostHits_branch;
	bool mus_lostHits_isLoaded;
	vector<int> mus_lostPixelHits_;
	TBranch *mus_lostPixelHits_branch;
	bool mus_lostPixelHits_isLoaded;
	vector<int> mus_mc_patMatch_id_;
	TBranch *mus_mc_patMatch_id_branch;
	bool mus_mc_patMatch_id_isLoaded;
	vector<int> mus_muonBestTrackType_;
	TBranch *mus_muonBestTrackType_branch;
	bool mus_muonBestTrackType_isLoaded;
	vector<int> mus_nOverlaps_;
	TBranch *mus_nOverlaps_branch;
	bool mus_nOverlaps_isLoaded;
	vector<int> mus_nlayers_;
	TBranch *mus_nlayers_branch;
	bool mus_nlayers_isLoaded;
	vector<int> mus_nlayers3D_;
	TBranch *mus_nlayers3D_branch;
	bool mus_nlayers3D_isLoaded;
	vector<int> mus_nlayersLost_;
	TBranch *mus_nlayersLost_branch;
	bool mus_nlayersLost_isLoaded;
	vector<int> mus_nmatches_;
	TBranch *mus_nmatches_branch;
	bool mus_nmatches_isLoaded;
	vector<int> mus_numberOfMatchedStations_;
	TBranch *mus_numberOfMatchedStations_branch;
	bool mus_numberOfMatchedStations_isLoaded;
	vector<int> mus_overlap0_;
	TBranch *mus_overlap0_branch;
	bool mus_overlap0_isLoaded;
	vector<int> mus_overlap1_;
	TBranch *mus_overlap1_branch;
	bool mus_overlap1_isLoaded;
	vector<int> mus_pfcharge_;
	TBranch *mus_pfcharge_branch;
	bool mus_pfcharge_isLoaded;
	vector<int> mus_pfparticleId_;
	TBranch *mus_pfparticleId_branch;
	bool mus_pfparticleId_isLoaded;
	vector<int> mus_pid_PFMuon_;
	TBranch *mus_pid_PFMuon_branch;
	bool mus_pid_PFMuon_isLoaded;
	vector<int> mus_pid_TM2DCompatibilityLoose_;
	TBranch *mus_pid_TM2DCompatibilityLoose_branch;
	bool mus_pid_TM2DCompatibilityLoose_isLoaded;
	vector<int> mus_pid_TM2DCompatibilityTight_;
	TBranch *mus_pid_TM2DCompatibilityTight_branch;
	bool mus_pid_TM2DCompatibilityTight_isLoaded;
	vector<int> mus_pid_TMLastStationLoose_;
	TBranch *mus_pid_TMLastStationLoose_branch;
	bool mus_pid_TMLastStationLoose_isLoaded;
	vector<int> mus_pid_TMLastStationTight_;
	TBranch *mus_pid_TMLastStationTight_branch;
	bool mus_pid_TMLastStationTight_isLoaded;
	vector<int> mus_qualityMask_;
	TBranch *mus_qualityMask_branch;
	bool mus_qualityMask_isLoaded;
	vector<int> mus_sta_algo_;
	TBranch *mus_sta_algo_branch;
	bool mus_sta_algo_isLoaded;
	vector<int> mus_sta_etaErr_;
	TBranch *mus_sta_etaErr_branch;
	bool mus_sta_etaErr_isLoaded;
	vector<int> mus_sta_exp_innerlayer_;
	TBranch *mus_sta_exp_innerlayer_branch;
	bool mus_sta_exp_innerlayer_isLoaded;
	vector<int> mus_sta_exp_outerlayer_;
	TBranch *mus_sta_exp_outerlayer_branch;
	bool mus_sta_exp_outerlayer_isLoaded;
	vector<int> mus_sta_lostHits_;
	TBranch *mus_sta_lostHits_branch;
	bool mus_sta_lostHits_isLoaded;
	vector<int> mus_sta_lostPixelHits_;
	TBranch *mus_sta_lostPixelHits_branch;
	bool mus_sta_lostPixelHits_isLoaded;
	vector<int> mus_sta_nlayers_;
	TBranch *mus_sta_nlayers_branch;
	bool mus_sta_nlayers_isLoaded;
	vector<int> mus_sta_nlayers3D_;
	TBranch *mus_sta_nlayers3D_branch;
	bool mus_sta_nlayers3D_isLoaded;
	vector<int> mus_sta_nlayersLost_;
	TBranch *mus_sta_nlayersLost_branch;
	bool mus_sta_nlayersLost_isLoaded;
	vector<int> mus_sta_phiErr_;
	TBranch *mus_sta_phiErr_branch;
	bool mus_sta_phiErr_isLoaded;
	vector<int> mus_sta_ptErr_;
	TBranch *mus_sta_ptErr_branch;
	bool mus_sta_ptErr_isLoaded;
	vector<int> mus_sta_qualityMask_;
	TBranch *mus_sta_qualityMask_branch;
	bool mus_sta_qualityMask_isLoaded;
	vector<int> mus_sta_trk_charge_;
	TBranch *mus_sta_trk_charge_branch;
	bool mus_sta_trk_charge_isLoaded;
	vector<int> mus_sta_validHits_;
	TBranch *mus_sta_validHits_branch;
	bool mus_sta_validHits_isLoaded;
	vector<int> mus_sta_validPixelHits_;
	TBranch *mus_sta_validPixelHits_branch;
	bool mus_sta_validPixelHits_isLoaded;
	vector<int> mus_timeDirection_;
	TBranch *mus_timeDirection_branch;
	bool mus_timeDirection_isLoaded;
	vector<int> mus_timeNumStationsUsed_;
	TBranch *mus_timeNumStationsUsed_branch;
	bool mus_timeNumStationsUsed_isLoaded;
	vector<int> mus_trk_charge_;
	TBranch *mus_trk_charge_branch;
	bool mus_trk_charge_isLoaded;
	vector<int> mus_trkidx_;
	TBranch *mus_trkidx_branch;
	bool mus_trkidx_isLoaded;
	vector<int> mus_type_;
	TBranch *mus_type_branch;
	bool mus_type_isLoaded;
	vector<int> mus_validHits_;
	TBranch *mus_validHits_branch;
	bool mus_validHits_isLoaded;
	vector<int> mus_validPixelHits_;
	TBranch *mus_validPixelHits_branch;
	bool mus_validPixelHits_isLoaded;
	vector<int> pfcands_charge_;
	TBranch *pfcands_charge_branch;
	bool pfcands_charge_isLoaded;
	vector<int> pfcands_particleId_;
	TBranch *pfcands_particleId_branch;
	bool pfcands_particleId_isLoaded;
	vector<int> pfjets_chargedHadronMultiplicity_;
	TBranch *pfjets_chargedHadronMultiplicity_branch;
	bool pfjets_chargedHadronMultiplicity_isLoaded;
	vector<int> pfjets_chargedMultiplicity_;
	TBranch *pfjets_chargedMultiplicity_branch;
	bool pfjets_chargedMultiplicity_isLoaded;
	vector<int> pfjets_electronMultiplicity_;
	TBranch *pfjets_electronMultiplicity_branch;
	bool pfjets_electronMultiplicity_isLoaded;
	vector<int> pfjets_hfEmMultiplicity_;
	TBranch *pfjets_hfEmMultiplicity_branch;
	bool pfjets_hfEmMultiplicity_isLoaded;
	vector<int> pfjets_hfHadronMultiplicity_;
	TBranch *pfjets_hfHadronMultiplicity_branch;
	bool pfjets_hfHadronMultiplicity_isLoaded;
	vector<int> pfjets_muonMultiplicity_;
	TBranch *pfjets_muonMultiplicity_branch;
	bool pfjets_muonMultiplicity_isLoaded;
	vector<int> pfjets_neutralHadronMultiplicity_;
	TBranch *pfjets_neutralHadronMultiplicity_branch;
	bool pfjets_neutralHadronMultiplicity_isLoaded;
	vector<int> pfjets_neutralMultiplicity_;
	TBranch *pfjets_neutralMultiplicity_branch;
	bool pfjets_neutralMultiplicity_isLoaded;
	vector<int> pfjets_partonFlavour_;
	TBranch *pfjets_partonFlavour_branch;
	bool pfjets_partonFlavour_isLoaded;
	vector<int> pfjets_photonMultiplicity_;
	TBranch *pfjets_photonMultiplicity_branch;
	bool pfjets_photonMultiplicity_isLoaded;
	vector<int> taus_pf_charge_;
	TBranch *taus_pf_charge_branch;
	bool taus_pf_charge_isLoaded;
	vector<int> photons_fiduciality_;
	TBranch *photons_fiduciality_branch;
	bool photons_fiduciality_isLoaded;
	vector<int> puInfo_bunchCrossing_;
	TBranch *puInfo_bunchCrossing_branch;
	bool puInfo_bunchCrossing_isLoaded;
	vector<int> puInfo_nPUvertices_;
	TBranch *puInfo_nPUvertices_branch;
	bool puInfo_nPUvertices_isLoaded;
	vector<int> convs_algo_;
	TBranch *convs_algo_branch;
	bool convs_algo_isLoaded;
	vector<int> convs_isConverted_;
	TBranch *convs_isConverted_branch;
	bool convs_isConverted_isLoaded;
	vector<int> convs_quality_;
	TBranch *convs_quality_branch;
	bool convs_quality_isLoaded;
	vector<int> vtxs_isFake_;
	TBranch *vtxs_isFake_branch;
	bool vtxs_isFake_isLoaded;
	vector<int> vtxs_isValid_;
	TBranch *vtxs_isValid_branch;
	bool vtxs_isValid_isLoaded;
	vector<int> vtxs_tracksSize_;
	TBranch *vtxs_tracksSize_branch;
	bool vtxs_tracksSize_isLoaded;
	vector<vector<int> > els_convs_delMissHits_;
	TBranch *els_convs_delMissHits_branch;
	bool els_convs_delMissHits_isLoaded;
	vector<vector<int> > els_convs_flag_;
	TBranch *els_convs_flag_branch;
	bool els_convs_flag_isLoaded;
	vector<vector<int> > els_convs_gsftkidx_;
	TBranch *els_convs_gsftkidx_branch;
	bool els_convs_gsftkidx_isLoaded;
	vector<vector<int> > els_convs_tkidx_;
	TBranch *els_convs_tkidx_branch;
	bool els_convs_tkidx_isLoaded;
	vector<vector<int> > els_PFCand_idx_;
	TBranch *els_PFCand_idx_branch;
	bool els_PFCand_idx_isLoaded;
	vector<vector<int> > genps_lepdaughter_id_;
	TBranch *genps_lepdaughter_id_branch;
	bool genps_lepdaughter_id_isLoaded;
	vector<vector<int> > genps_lepdaughter_idx_;
	TBranch *genps_lepdaughter_idx_branch;
	bool genps_lepdaughter_idx_isLoaded;
	vector<vector<int> > hlt_trigObjs_id_;
	TBranch *hlt_trigObjs_id_branch;
	bool hlt_trigObjs_id_isLoaded;
	vector<vector<int> > mus_nStationCorrelatedHits_;
	TBranch *mus_nStationCorrelatedHits_branch;
	bool mus_nStationCorrelatedHits_isLoaded;
	vector<vector<int> > mus_nStationHits_;
	TBranch *mus_nStationHits_branch;
	bool mus_nStationHits_isLoaded;
	vector<vector<int> > pfjets_pfcandIndicies_;
	TBranch *pfjets_pfcandIndicies_branch;
	bool pfjets_pfcandIndicies_isLoaded;
	vector<vector<int> > convs_nHitsBeforeVtx_;
	TBranch *convs_nHitsBeforeVtx_branch;
	bool convs_nHitsBeforeVtx_isLoaded;
	vector<vector<int> > convs_tkalgo_;
	TBranch *convs_tkalgo_branch;
	bool convs_tkalgo_isLoaded;
	vector<vector<int> > convs_tkidx_;
	TBranch *convs_tkidx_branch;
	bool convs_tkidx_isLoaded;
	unsigned int els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_;
	TBranch *els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch;
	bool els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_isLoaded;
	unsigned int els_HLT_Ele17_Ele8_LeadingLeg_version_;
	TBranch *els_HLT_Ele17_Ele8_LeadingLeg_version_branch;
	bool els_HLT_Ele17_Ele8_LeadingLeg_version_isLoaded;
	unsigned int els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_;
	TBranch *els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch;
	bool els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_isLoaded;
	unsigned int els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_;
	TBranch *els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch;
	bool els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_isLoaded;
	unsigned int els_HLT_Ele17_Ele8_TrailingLeg_version_;
	TBranch *els_HLT_Ele17_Ele8_TrailingLeg_version_branch;
	bool els_HLT_Ele17_Ele8_TrailingLeg_version_isLoaded;
	unsigned int els_HLT_Ele17_Ele8_version_;
	TBranch *els_HLT_Ele17_Ele8_version_branch;
	bool els_HLT_Ele17_Ele8_version_isLoaded;
	unsigned int els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_;
	TBranch *els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch;
	bool els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_isLoaded;
	unsigned int els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_;
	TBranch *els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch;
	bool els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_isLoaded;
	unsigned int els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_;
	TBranch *els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch;
	bool els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_isLoaded;
	unsigned int els_HLT_Ele27_WP80_version_;
	TBranch *els_HLT_Ele27_WP80_version_branch;
	bool els_HLT_Ele27_WP80_version_isLoaded;
	unsigned int els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_;
	TBranch *els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch;
	bool els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_isLoaded;
	unsigned int els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_;
	TBranch *els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch;
	bool els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_isLoaded;
	unsigned int els_HLT_Mu17_Ele8_TrailingLeg_version_;
	TBranch *els_HLT_Mu17_Ele8_TrailingLeg_version_branch;
	bool els_HLT_Mu17_Ele8_TrailingLeg_version_isLoaded;
	unsigned int els_HLT_Mu17_Ele8_version_;
	TBranch *els_HLT_Mu17_Ele8_version_branch;
	bool els_HLT_Mu17_Ele8_version_isLoaded;
	unsigned int els_HLT_Mu8_Ele17_version_;
	TBranch *els_HLT_Mu8_Ele17_version_branch;
	bool els_HLT_Mu8_Ele17_version_isLoaded;
	unsigned int evt_nels_;
	TBranch *evt_nels_branch;
	bool evt_nels_isLoaded;
	unsigned int evt_detectorStatus_;
	TBranch *evt_detectorStatus_branch;
	bool evt_detectorStatus_isLoaded;
	unsigned int evt_event_;
	TBranch *evt_event_branch;
	bool evt_event_isLoaded;
	unsigned int evt_lumiBlock_;
	TBranch *evt_lumiBlock_branch;
	bool evt_lumiBlock_isLoaded;
	unsigned int evt_run_;
	TBranch *evt_run_branch;
	bool evt_run_isLoaded;
	unsigned int evt_ngenjetsNoMuNoNu_;
	TBranch *evt_ngenjetsNoMuNoNu_branch;
	bool evt_ngenjetsNoMuNoNu_isLoaded;
	unsigned int genps_signalProcessID_;
	TBranch *genps_signalProcessID_branch;
	bool genps_signalProcessID_isLoaded;
	unsigned int mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_;
	TBranch *mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch;
	bool mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_isLoaded;
	unsigned int mus_HLT_IsoMu24_eta2p1_version_;
	TBranch *mus_HLT_IsoMu24_eta2p1_version_branch;
	bool mus_HLT_IsoMu24_eta2p1_version_isLoaded;
	unsigned int mus_HLT_Mu17_Ele8_LeadingLeg_version_;
	TBranch *mus_HLT_Mu17_Ele8_LeadingLeg_version_branch;
	bool mus_HLT_Mu17_Ele8_LeadingLeg_version_isLoaded;
	unsigned int mus_HLT_Mu17_Ele8_version_;
	TBranch *mus_HLT_Mu17_Ele8_version_branch;
	bool mus_HLT_Mu17_Ele8_version_isLoaded;
	unsigned int mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_;
	TBranch *mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch;
	bool mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_isLoaded;
	unsigned int mus_HLT_Mu17_Mu8_LeadingLeg_version_;
	TBranch *mus_HLT_Mu17_Mu8_LeadingLeg_version_branch;
	bool mus_HLT_Mu17_Mu8_LeadingLeg_version_isLoaded;
	unsigned int mus_HLT_Mu17_Mu8_TrailingLeg_version_;
	TBranch *mus_HLT_Mu17_Mu8_TrailingLeg_version_branch;
	bool mus_HLT_Mu17_Mu8_TrailingLeg_version_isLoaded;
	unsigned int mus_HLT_Mu17_Mu8_version_;
	TBranch *mus_HLT_Mu17_Mu8_version_branch;
	bool mus_HLT_Mu17_Mu8_version_isLoaded;
	unsigned int mus_HLT_Mu17_TkMu8_LeadingLeg_version_;
	TBranch *mus_HLT_Mu17_TkMu8_LeadingLeg_version_branch;
	bool mus_HLT_Mu17_TkMu8_LeadingLeg_version_isLoaded;
	unsigned int mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_;
	TBranch *mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch;
	bool mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_isLoaded;
	unsigned int mus_HLT_Mu17_TkMu8_TrailingLeg_version_;
	TBranch *mus_HLT_Mu17_TkMu8_TrailingLeg_version_branch;
	bool mus_HLT_Mu17_TkMu8_TrailingLeg_version_isLoaded;
	unsigned int mus_HLT_Mu17_TkMu8_version_;
	TBranch *mus_HLT_Mu17_TkMu8_version_branch;
	bool mus_HLT_Mu17_TkMu8_version_isLoaded;
	unsigned int mus_HLT_Mu8_Ele17_TrailingLeg_version_;
	TBranch *mus_HLT_Mu8_Ele17_TrailingLeg_version_branch;
	bool mus_HLT_Mu8_Ele17_TrailingLeg_version_isLoaded;
	unsigned int mus_HLT_Mu8_Ele17_version_;
	TBranch *mus_HLT_Mu8_Ele17_version_branch;
	bool mus_HLT_Mu8_Ele17_version_isLoaded;
	unsigned int evt_nphotons_;
	TBranch *evt_nphotons_branch;
	bool evt_nphotons_isLoaded;
	unsigned int evt_nvtxs_;
	TBranch *evt_nvtxs_branch;
	bool evt_nvtxs_isLoaded;
	vector<unsigned int> els_HLT_Ele17_Ele8_;
	TBranch *els_HLT_Ele17_Ele8_branch;
	bool els_HLT_Ele17_Ele8_isLoaded;
	vector<unsigned int> els_HLT_Ele17_Ele8_L1sL1DoubleEG137_;
	TBranch *els_HLT_Ele17_Ele8_L1sL1DoubleEG137_branch;
	bool els_HLT_Ele17_Ele8_L1sL1DoubleEG137_isLoaded;
	vector<unsigned int> els_HLT_Ele17_Ele8_LeadingLeg_;
	TBranch *els_HLT_Ele17_Ele8_LeadingLeg_branch;
	bool els_HLT_Ele17_Ele8_LeadingLeg_isLoaded;
	vector<unsigned int> els_HLT_Ele17_Ele8_Mass50_LeadingLeg_;
	TBranch *els_HLT_Ele17_Ele8_Mass50_LeadingLeg_branch;
	bool els_HLT_Ele17_Ele8_Mass50_LeadingLeg_isLoaded;
	vector<unsigned int> els_HLT_Ele17_Ele8_Mass50_TrailingLeg_;
	TBranch *els_HLT_Ele17_Ele8_Mass50_TrailingLeg_branch;
	bool els_HLT_Ele17_Ele8_Mass50_TrailingLeg_isLoaded;
	vector<unsigned int> els_HLT_Ele17_Ele8_TrailingLeg_;
	TBranch *els_HLT_Ele17_Ele8_TrailingLeg_branch;
	bool els_HLT_Ele17_Ele8_TrailingLeg_isLoaded;
	vector<unsigned int> els_HLT_Ele20_SC4_Mass50_LeadingLeg_;
	TBranch *els_HLT_Ele20_SC4_Mass50_LeadingLeg_branch;
	bool els_HLT_Ele20_SC4_Mass50_LeadingLeg_isLoaded;
	vector<unsigned int> els_HLT_Ele20_SC4_Mass50_TrailingLeg_;
	TBranch *els_HLT_Ele20_SC4_Mass50_TrailingLeg_branch;
	bool els_HLT_Ele20_SC4_Mass50_TrailingLeg_isLoaded;
	vector<unsigned int> els_HLT_Ele27_WP80_;
	TBranch *els_HLT_Ele27_WP80_branch;
	bool els_HLT_Ele27_WP80_isLoaded;
	vector<unsigned int> els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_;
	TBranch *els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_branch;
	bool els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_isLoaded;
	vector<unsigned int> els_HLT_Ele32_SC17_Mass50_LeadingLeg_;
	TBranch *els_HLT_Ele32_SC17_Mass50_LeadingLeg_branch;
	bool els_HLT_Ele32_SC17_Mass50_LeadingLeg_isLoaded;
	vector<unsigned int> els_HLT_Ele32_SC17_Mass50_TrailingLeg_;
	TBranch *els_HLT_Ele32_SC17_Mass50_TrailingLeg_branch;
	bool els_HLT_Ele32_SC17_Mass50_TrailingLeg_isLoaded;
	vector<unsigned int> els_HLT_Mu17_Ele8_;
	TBranch *els_HLT_Mu17_Ele8_branch;
	bool els_HLT_Mu17_Ele8_isLoaded;
	vector<unsigned int> els_HLT_Mu17_Ele8_TrailingLeg_;
	TBranch *els_HLT_Mu17_Ele8_TrailingLeg_branch;
	bool els_HLT_Mu17_Ele8_TrailingLeg_isLoaded;
	vector<unsigned int> els_HLT_Mu8_Ele17_;
	TBranch *els_HLT_Mu8_Ele17_branch;
	bool els_HLT_Mu8_Ele17_isLoaded;
	vector<unsigned int> els_id2012ext_loose_;
	TBranch *els_id2012ext_loose_branch;
	bool els_id2012ext_loose_isLoaded;
	vector<unsigned int> els_id2012ext_medium_;
	TBranch *els_id2012ext_medium_branch;
	bool els_id2012ext_medium_isLoaded;
	vector<unsigned int> els_id2012ext_tight_;
	TBranch *els_id2012ext_tight_branch;
	bool els_id2012ext_tight_isLoaded;
	vector<unsigned int> els_id2012ext_veto_;
	TBranch *els_id2012ext_veto_branch;
	bool els_id2012ext_veto_isLoaded;
	vector<unsigned int> els_id2012_loose_;
	TBranch *els_id2012_loose_branch;
	bool els_id2012_loose_isLoaded;
	vector<unsigned int> els_id2012_medium_;
	TBranch *els_id2012_medium_branch;
	bool els_id2012_medium_isLoaded;
	vector<unsigned int> els_id2012_tight_;
	TBranch *els_id2012_tight_branch;
	bool els_id2012_tight_isLoaded;
	vector<unsigned int> els_id2012_veto_;
	TBranch *els_id2012_veto_branch;
	bool els_id2012_veto_isLoaded;
	vector<unsigned int> hlt_prescales_;
	TBranch *hlt_prescales_branch;
	bool hlt_prescales_isLoaded;
	vector<unsigned int> mus_HLT_IsoMu24_eta2p1_;
	TBranch *mus_HLT_IsoMu24_eta2p1_branch;
	bool mus_HLT_IsoMu24_eta2p1_isLoaded;
	vector<unsigned int> mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_;
	TBranch *mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_branch;
	bool mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_isLoaded;
	vector<unsigned int> mus_HLT_Mu17_Ele8_;
	TBranch *mus_HLT_Mu17_Ele8_branch;
	bool mus_HLT_Mu17_Ele8_isLoaded;
	vector<unsigned int> mus_HLT_Mu17_Ele8_LeadingLeg_;
	TBranch *mus_HLT_Mu17_Ele8_LeadingLeg_branch;
	bool mus_HLT_Mu17_Ele8_LeadingLeg_isLoaded;
	vector<unsigned int> mus_HLT_Mu17_Mu8_;
	TBranch *mus_HLT_Mu17_Mu8_branch;
	bool mus_HLT_Mu17_Mu8_isLoaded;
	vector<unsigned int> mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_;
	TBranch *mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_branch;
	bool mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_isLoaded;
	vector<unsigned int> mus_HLT_Mu17_Mu8_LeadingLeg_;
	TBranch *mus_HLT_Mu17_Mu8_LeadingLeg_branch;
	bool mus_HLT_Mu17_Mu8_LeadingLeg_isLoaded;
	vector<unsigned int> mus_HLT_Mu17_Mu8_TrailingLeg_;
	TBranch *mus_HLT_Mu17_Mu8_TrailingLeg_branch;
	bool mus_HLT_Mu17_Mu8_TrailingLeg_isLoaded;
	vector<unsigned int> mus_HLT_Mu17_TkMu8_;
	TBranch *mus_HLT_Mu17_TkMu8_branch;
	bool mus_HLT_Mu17_TkMu8_isLoaded;
	vector<unsigned int> mus_HLT_Mu17_TkMu8_LeadingLeg_;
	TBranch *mus_HLT_Mu17_TkMu8_LeadingLeg_branch;
	bool mus_HLT_Mu17_TkMu8_LeadingLeg_isLoaded;
	vector<unsigned int> mus_HLT_Mu17_TkMu8_TrailingLeg_;
	TBranch *mus_HLT_Mu17_TkMu8_TrailingLeg_branch;
	bool mus_HLT_Mu17_TkMu8_TrailingLeg_isLoaded;
	vector<unsigned int> mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_;
	TBranch *mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_branch;
	bool mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_isLoaded;
	vector<unsigned int> mus_HLT_Mu8_Ele17_;
	TBranch *mus_HLT_Mu8_Ele17_branch;
	bool mus_HLT_Mu8_Ele17_isLoaded;
	vector<unsigned int> mus_HLT_Mu8_Ele17_TrailingLeg_;
	TBranch *mus_HLT_Mu8_Ele17_TrailingLeg_branch;
	bool mus_HLT_Mu8_Ele17_TrailingLeg_isLoaded;
public: 
void Init(TTree *tree) {
	hlt_bits_branch = 0;
	if (tree->GetAlias("hlt_bits") != 0) {
		hlt_bits_branch = tree->GetBranch(tree->GetAlias("hlt_bits"));
		hlt_bits_branch->SetAddress(&hlt_bits_);
	}
	evt_bsp4_branch = 0;
	if (tree->GetAlias("evt_bsp4") != 0) {
		evt_bsp4_branch = tree->GetBranch(tree->GetAlias("evt_bsp4"));
		evt_bsp4_branch->SetAddress(&evt_bsp4_);
	}
	els_mc_motherp4_branch = 0;
	if (tree->GetAlias("els_mc_motherp4") != 0) {
		els_mc_motherp4_branch = tree->GetBranch(tree->GetAlias("els_mc_motherp4"));
		els_mc_motherp4_branch->SetAddress(&els_mc_motherp4_);
	}
	els_mc_p4_branch = 0;
	if (tree->GetAlias("els_mc_p4") != 0) {
		els_mc_p4_branch = tree->GetBranch(tree->GetAlias("els_mc_p4"));
		els_mc_p4_branch->SetAddress(&els_mc_p4_);
	}
	mus_mc_motherp4_branch = 0;
	if (tree->GetAlias("mus_mc_motherp4") != 0) {
		mus_mc_motherp4_branch = tree->GetBranch(tree->GetAlias("mus_mc_motherp4"));
		mus_mc_motherp4_branch->SetAddress(&mus_mc_motherp4_);
	}
	mus_mc_p4_branch = 0;
	if (tree->GetAlias("mus_mc_p4") != 0) {
		mus_mc_p4_branch = tree->GetBranch(tree->GetAlias("mus_mc_p4"));
		mus_mc_p4_branch->SetAddress(&mus_mc_p4_);
	}
	pfjets_mc_gp_p4_branch = 0;
	if (tree->GetAlias("pfjets_mc_gp_p4") != 0) {
		pfjets_mc_gp_p4_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_gp_p4"));
		pfjets_mc_gp_p4_branch->SetAddress(&pfjets_mc_gp_p4_);
	}
	pfjets_mc_motherp4_branch = 0;
	if (tree->GetAlias("pfjets_mc_motherp4") != 0) {
		pfjets_mc_motherp4_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_motherp4"));
		pfjets_mc_motherp4_branch->SetAddress(&pfjets_mc_motherp4_);
	}
	pfjets_mc_p4_branch = 0;
	if (tree->GetAlias("pfjets_mc_p4") != 0) {
		pfjets_mc_p4_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_p4"));
		pfjets_mc_p4_branch->SetAddress(&pfjets_mc_p4_);
	}
	photons_mc_motherp4_branch = 0;
	if (tree->GetAlias("photons_mc_motherp4") != 0) {
		photons_mc_motherp4_branch = tree->GetBranch(tree->GetAlias("photons_mc_motherp4"));
		photons_mc_motherp4_branch->SetAddress(&photons_mc_motherp4_);
	}
	photons_mc_p4_branch = 0;
	if (tree->GetAlias("photons_mc_p4") != 0) {
		photons_mc_p4_branch = tree->GetBranch(tree->GetAlias("photons_mc_p4"));
		photons_mc_p4_branch->SetAddress(&photons_mc_p4_);
	}
	els_conv_pos_p4_branch = 0;
	if (tree->GetAlias("els_conv_pos_p4") != 0) {
		els_conv_pos_p4_branch = tree->GetBranch(tree->GetAlias("els_conv_pos_p4"));
		els_conv_pos_p4_branch->SetAddress(&els_conv_pos_p4_);
	}
	els_inner_position_branch = 0;
	if (tree->GetAlias("els_inner_position") != 0) {
		els_inner_position_branch = tree->GetBranch(tree->GetAlias("els_inner_position"));
		els_inner_position_branch->SetAddress(&els_inner_position_);
	}
	els_mc_patMatch_p4_branch = 0;
	if (tree->GetAlias("els_mc_patMatch_p4") != 0) {
		els_mc_patMatch_p4_branch = tree->GetBranch(tree->GetAlias("els_mc_patMatch_p4"));
		els_mc_patMatch_p4_branch->SetAddress(&els_mc_patMatch_p4_);
	}
	els_outer_position_branch = 0;
	if (tree->GetAlias("els_outer_position") != 0) {
		els_outer_position_branch = tree->GetBranch(tree->GetAlias("els_outer_position"));
		els_outer_position_branch->SetAddress(&els_outer_position_);
	}
	els_p4_branch = 0;
	if (tree->GetAlias("els_p4") != 0) {
		els_p4_branch = tree->GetBranch(tree->GetAlias("els_p4"));
		els_p4_branch->SetAddress(&els_p4_);
	}
	els_p4In_branch = 0;
	if (tree->GetAlias("els_p4In") != 0) {
		els_p4In_branch = tree->GetBranch(tree->GetAlias("els_p4In"));
		els_p4In_branch->SetAddress(&els_p4In_);
	}
	els_p4Out_branch = 0;
	if (tree->GetAlias("els_p4Out") != 0) {
		els_p4Out_branch = tree->GetBranch(tree->GetAlias("els_p4Out"));
		els_p4Out_branch->SetAddress(&els_p4Out_);
	}
	els_trk_p4_branch = 0;
	if (tree->GetAlias("els_trk_p4") != 0) {
		els_trk_p4_branch = tree->GetBranch(tree->GetAlias("els_trk_p4"));
		els_trk_p4_branch->SetAddress(&els_trk_p4_);
	}
	els_trk_vertex_p4_branch = 0;
	if (tree->GetAlias("els_trk_vertex_p4") != 0) {
		els_trk_vertex_p4_branch = tree->GetBranch(tree->GetAlias("els_trk_vertex_p4"));
		els_trk_vertex_p4_branch->SetAddress(&els_trk_vertex_p4_);
	}
	els_vertex_p4_branch = 0;
	if (tree->GetAlias("els_vertex_p4") != 0) {
		els_vertex_p4_branch = tree->GetBranch(tree->GetAlias("els_vertex_p4"));
		els_vertex_p4_branch->SetAddress(&els_vertex_p4_);
	}
	genjets_p4NoMuNoNu_branch = 0;
	if (tree->GetAlias("genjets_p4NoMuNoNu") != 0) {
		genjets_p4NoMuNoNu_branch = tree->GetBranch(tree->GetAlias("genjets_p4NoMuNoNu"));
		genjets_p4NoMuNoNu_branch->SetAddress(&genjets_p4NoMuNoNu_);
	}
	genps_p4_branch = 0;
	if (tree->GetAlias("genps_p4") != 0) {
		genps_p4_branch = tree->GetBranch(tree->GetAlias("genps_p4"));
		genps_p4_branch->SetAddress(&genps_p4_);
	}
	genps_prod_vtx_branch = 0;
	if (tree->GetAlias("genps_prod_vtx") != 0) {
		genps_prod_vtx_branch = tree->GetBranch(tree->GetAlias("genps_prod_vtx"));
		genps_prod_vtx_branch->SetAddress(&genps_prod_vtx_);
	}
	mus_gfit_p4_branch = 0;
	if (tree->GetAlias("mus_gfit_p4") != 0) {
		mus_gfit_p4_branch = tree->GetBranch(tree->GetAlias("mus_gfit_p4"));
		mus_gfit_p4_branch->SetAddress(&mus_gfit_p4_);
	}
	mus_gfit_vertex_p4_branch = 0;
	if (tree->GetAlias("mus_gfit_vertex_p4") != 0) {
		mus_gfit_vertex_p4_branch = tree->GetBranch(tree->GetAlias("mus_gfit_vertex_p4"));
		mus_gfit_vertex_p4_branch->SetAddress(&mus_gfit_vertex_p4_);
	}
	mus_mc_patMatch_p4_branch = 0;
	if (tree->GetAlias("mus_mc_patMatch_p4") != 0) {
		mus_mc_patMatch_p4_branch = tree->GetBranch(tree->GetAlias("mus_mc_patMatch_p4"));
		mus_mc_patMatch_p4_branch->SetAddress(&mus_mc_patMatch_p4_);
	}
	mus_p4_branch = 0;
	if (tree->GetAlias("mus_p4") != 0) {
		mus_p4_branch = tree->GetBranch(tree->GetAlias("mus_p4"));
		mus_p4_branch->SetAddress(&mus_p4_);
	}
	mus_pfp4_branch = 0;
	if (tree->GetAlias("mus_pfp4") != 0) {
		mus_pfp4_branch = tree->GetBranch(tree->GetAlias("mus_pfp4"));
		mus_pfp4_branch->SetAddress(&mus_pfp4_);
	}
	mus_sta_p4_branch = 0;
	if (tree->GetAlias("mus_sta_p4") != 0) {
		mus_sta_p4_branch = tree->GetBranch(tree->GetAlias("mus_sta_p4"));
		mus_sta_p4_branch->SetAddress(&mus_sta_p4_);
	}
	mus_sta_vertex_p4_branch = 0;
	if (tree->GetAlias("mus_sta_vertex_p4") != 0) {
		mus_sta_vertex_p4_branch = tree->GetBranch(tree->GetAlias("mus_sta_vertex_p4"));
		mus_sta_vertex_p4_branch->SetAddress(&mus_sta_vertex_p4_);
	}
	mus_trk_p4_branch = 0;
	if (tree->GetAlias("mus_trk_p4") != 0) {
		mus_trk_p4_branch = tree->GetBranch(tree->GetAlias("mus_trk_p4"));
		mus_trk_p4_branch->SetAddress(&mus_trk_p4_);
	}
	mus_vertex_p4_branch = 0;
	if (tree->GetAlias("mus_vertex_p4") != 0) {
		mus_vertex_p4_branch = tree->GetBranch(tree->GetAlias("mus_vertex_p4"));
		mus_vertex_p4_branch->SetAddress(&mus_vertex_p4_);
	}
	pfcands_p4_branch = 0;
	if (tree->GetAlias("pfcands_p4") != 0) {
		pfcands_p4_branch = tree->GetBranch(tree->GetAlias("pfcands_p4"));
		pfcands_p4_branch->SetAddress(&pfcands_p4_);
	}
	pfjets_p4_branch = 0;
	if (tree->GetAlias("pfjets_p4") != 0) {
		pfjets_p4_branch = tree->GetBranch(tree->GetAlias("pfjets_p4"));
		pfjets_p4_branch->SetAddress(&pfjets_p4_);
	}
	taus_pf_lead_chargecand_p4_branch = 0;
	if (tree->GetAlias("taus_pf_lead_chargecand_p4") != 0) {
		taus_pf_lead_chargecand_p4_branch = tree->GetBranch(tree->GetAlias("taus_pf_lead_chargecand_p4"));
		taus_pf_lead_chargecand_p4_branch->SetAddress(&taus_pf_lead_chargecand_p4_);
	}
	taus_pf_lead_neutrcand_p4_branch = 0;
	if (tree->GetAlias("taus_pf_lead_neutrcand_p4") != 0) {
		taus_pf_lead_neutrcand_p4_branch = tree->GetBranch(tree->GetAlias("taus_pf_lead_neutrcand_p4"));
		taus_pf_lead_neutrcand_p4_branch->SetAddress(&taus_pf_lead_neutrcand_p4_);
	}
	taus_pf_p4_branch = 0;
	if (tree->GetAlias("taus_pf_p4") != 0) {
		taus_pf_p4_branch = tree->GetBranch(tree->GetAlias("taus_pf_p4"));
		taus_pf_p4_branch->SetAddress(&taus_pf_p4_);
	}
	photons_p4_branch = 0;
	if (tree->GetAlias("photons_p4") != 0) {
		photons_p4_branch = tree->GetBranch(tree->GetAlias("photons_p4"));
		photons_p4_branch->SetAddress(&photons_p4_);
	}
	vtxs_position_branch = 0;
	if (tree->GetAlias("vtxs_position") != 0) {
		vtxs_position_branch = tree->GetBranch(tree->GetAlias("vtxs_position"));
		vtxs_position_branch->SetAddress(&vtxs_position_);
	}
  tree->SetMakeClass(1);
	evt_CMS2tag_branch = 0;
	if (tree->GetAlias("evt_CMS2tag") != 0) {
		evt_CMS2tag_branch = tree->GetBranch(tree->GetAlias("evt_CMS2tag"));
		evt_CMS2tag_branch->SetAddress(&evt_CMS2tag_);
	}
	evt_dataset_branch = 0;
	if (tree->GetAlias("evt_dataset") != 0) {
		evt_dataset_branch = tree->GetBranch(tree->GetAlias("evt_dataset"));
		evt_dataset_branch->SetAddress(&evt_dataset_);
	}
	hlt_trigNames_branch = 0;
	if (tree->GetAlias("hlt_trigNames") != 0) {
		hlt_trigNames_branch = tree->GetBranch(tree->GetAlias("hlt_trigNames"));
		hlt_trigNames_branch->SetAddress(&hlt_trigNames_);
	}
	hcalnoise_HasBadRBXTS4TS5_branch = 0;
	if (tree->GetAlias("hcalnoise_HasBadRBXTS4TS5") != 0) {
		hcalnoise_HasBadRBXTS4TS5_branch = tree->GetBranch(tree->GetAlias("hcalnoise_HasBadRBXTS4TS5"));
		hcalnoise_HasBadRBXTS4TS5_branch->SetAddress(&hcalnoise_HasBadRBXTS4TS5_);
	}
	evt_cscTightHaloId_branch = 0;
	if (tree->GetAlias("evt_cscTightHaloId") != 0) {
		evt_cscTightHaloId_branch = tree->GetBranch(tree->GetAlias("evt_cscTightHaloId"));
		evt_cscTightHaloId_branch->SetAddress(&evt_cscTightHaloId_);
	}
	evt_hbheFilter_branch = 0;
	if (tree->GetAlias("evt_hbheFilter") != 0) {
		evt_hbheFilter_branch = tree->GetBranch(tree->GetAlias("evt_hbheFilter"));
		evt_hbheFilter_branch->SetAddress(&evt_hbheFilter_);
	}
	filt_cscBeamHalo_branch = 0;
	if (tree->GetAlias("filt_cscBeamHalo") != 0) {
		filt_cscBeamHalo_branch = tree->GetBranch(tree->GetAlias("filt_cscBeamHalo"));
		filt_cscBeamHalo_branch->SetAddress(&filt_cscBeamHalo_);
	}
	filt_ecalLaser_branch = 0;
	if (tree->GetAlias("filt_ecalLaser") != 0) {
		filt_ecalLaser_branch = tree->GetBranch(tree->GetAlias("filt_ecalLaser"));
		filt_ecalLaser_branch->SetAddress(&filt_ecalLaser_);
	}
	filt_ecalTP_branch = 0;
	if (tree->GetAlias("filt_ecalTP") != 0) {
		filt_ecalTP_branch = tree->GetBranch(tree->GetAlias("filt_ecalTP"));
		filt_ecalTP_branch->SetAddress(&filt_ecalTP_);
	}
	filt_eeBadSc_branch = 0;
	if (tree->GetAlias("filt_eeBadSc") != 0) {
		filt_eeBadSc_branch = tree->GetBranch(tree->GetAlias("filt_eeBadSc"));
		filt_eeBadSc_branch->SetAddress(&filt_eeBadSc_);
	}
	filt_goodVertices_branch = 0;
	if (tree->GetAlias("filt_goodVertices") != 0) {
		filt_goodVertices_branch = tree->GetBranch(tree->GetAlias("filt_goodVertices"));
		filt_goodVertices_branch->SetAddress(&filt_goodVertices_);
	}
	filt_hbheNoise_branch = 0;
	if (tree->GetAlias("filt_hbheNoise") != 0) {
		filt_hbheNoise_branch = tree->GetBranch(tree->GetAlias("filt_hbheNoise"));
		filt_hbheNoise_branch->SetAddress(&filt_hbheNoise_);
	}
	filt_hcalLaser_branch = 0;
	if (tree->GetAlias("filt_hcalLaser") != 0) {
		filt_hcalLaser_branch = tree->GetBranch(tree->GetAlias("filt_hcalLaser"));
		filt_hcalLaser_branch->SetAddress(&filt_hcalLaser_);
	}
	filt_metfilter_branch = 0;
	if (tree->GetAlias("filt_metfilter") != 0) {
		filt_metfilter_branch = tree->GetBranch(tree->GetAlias("filt_metfilter"));
		filt_metfilter_branch->SetAddress(&filt_metfilter_);
	}
	filt_trackingFailure_branch = 0;
	if (tree->GetAlias("filt_trackingFailure") != 0) {
		filt_trackingFailure_branch = tree->GetBranch(tree->GetAlias("filt_trackingFailure"));
		filt_trackingFailure_branch->SetAddress(&filt_trackingFailure_);
	}
	filt_trkPOGFilters_branch = 0;
	if (tree->GetAlias("filt_trkPOGFilters") != 0) {
		filt_trkPOGFilters_branch = tree->GetBranch(tree->GetAlias("filt_trkPOGFilters"));
		filt_trkPOGFilters_branch->SetAddress(&filt_trkPOGFilters_);
	}
	filt_trkPOG_logErrorTooManyClusters_branch = 0;
	if (tree->GetAlias("filt_trkPOG_logErrorTooManyClusters") != 0) {
		filt_trkPOG_logErrorTooManyClusters_branch = tree->GetBranch(tree->GetAlias("filt_trkPOG_logErrorTooManyClusters"));
		filt_trkPOG_logErrorTooManyClusters_branch->SetAddress(&filt_trkPOG_logErrorTooManyClusters_);
	}
	filt_trkPOG_manystripclus53X_branch = 0;
	if (tree->GetAlias("filt_trkPOG_manystripclus53X") != 0) {
		filt_trkPOG_manystripclus53X_branch = tree->GetBranch(tree->GetAlias("filt_trkPOG_manystripclus53X"));
		filt_trkPOG_manystripclus53X_branch->SetAddress(&filt_trkPOG_manystripclus53X_);
	}
	filt_trkPOG_toomanystripclus53X_branch = 0;
	if (tree->GetAlias("filt_trkPOG_toomanystripclus53X") != 0) {
		filt_trkPOG_toomanystripclus53X_branch = tree->GetBranch(tree->GetAlias("filt_trkPOG_toomanystripclus53X"));
		filt_trkPOG_toomanystripclus53X_branch->SetAddress(&filt_trkPOG_toomanystripclus53X_);
	}
	els_conv_vtx_flag_branch = 0;
	if (tree->GetAlias("els_conv_vtx_flag") != 0) {
		els_conv_vtx_flag_branch = tree->GetBranch(tree->GetAlias("els_conv_vtx_flag"));
		els_conv_vtx_flag_branch->SetAddress(&els_conv_vtx_flag_);
	}
	els_isGsfCtfScPixChargeConsistent_branch = 0;
	if (tree->GetAlias("els_isGsfCtfScPixChargeConsistent") != 0) {
		els_isGsfCtfScPixChargeConsistent_branch = tree->GetBranch(tree->GetAlias("els_isGsfCtfScPixChargeConsistent"));
		els_isGsfCtfScPixChargeConsistent_branch->SetAddress(&els_isGsfCtfScPixChargeConsistent_);
	}
	els_passingMvaPreselection_branch = 0;
	if (tree->GetAlias("els_passingMvaPreselection") != 0) {
		els_passingMvaPreselection_branch = tree->GetBranch(tree->GetAlias("els_passingMvaPreselection"));
		els_passingMvaPreselection_branch->SetAddress(&els_passingMvaPreselection_);
	}
	els_passingPflowPreselection_branch = 0;
	if (tree->GetAlias("els_passingPflowPreselection") != 0) {
		els_passingPflowPreselection_branch = tree->GetBranch(tree->GetAlias("els_passingPflowPreselection"));
		els_passingPflowPreselection_branch->SetAddress(&els_passingPflowPreselection_);
	}
	mus_isRPCMuon_branch = 0;
	if (tree->GetAlias("mus_isRPCMuon") != 0) {
		mus_isRPCMuon_branch = tree->GetBranch(tree->GetAlias("mus_isRPCMuon"));
		mus_isRPCMuon_branch->SetAddress(&mus_isRPCMuon_);
	}
	mus_tightMatch_branch = 0;
	if (tree->GetAlias("mus_tightMatch") != 0) {
		mus_tightMatch_branch = tree->GetBranch(tree->GetAlias("mus_tightMatch"));
		mus_tightMatch_branch->SetAddress(&mus_tightMatch_);
	}
	mus_updatedSta_branch = 0;
	if (tree->GetAlias("mus_updatedSta") != 0) {
		mus_updatedSta_branch = tree->GetBranch(tree->GetAlias("mus_updatedSta"));
		mus_updatedSta_branch->SetAddress(&mus_updatedSta_);
	}
	photons_haspixelSeed_branch = 0;
	if (tree->GetAlias("photons_haspixelSeed") != 0) {
		photons_haspixelSeed_branch = tree->GetBranch(tree->GetAlias("photons_haspixelSeed"));
		photons_haspixelSeed_branch->SetAddress(&photons_haspixelSeed_);
	}
	evt_bs_Xwidth_branch = 0;
	if (tree->GetAlias("evt_bs_Xwidth") != 0) {
		evt_bs_Xwidth_branch = tree->GetBranch(tree->GetAlias("evt_bs_Xwidth"));
		evt_bs_Xwidth_branch->SetAddress(&evt_bs_Xwidth_);
	}
	evt_bs_XwidthErr_branch = 0;
	if (tree->GetAlias("evt_bs_XwidthErr") != 0) {
		evt_bs_XwidthErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_XwidthErr"));
		evt_bs_XwidthErr_branch->SetAddress(&evt_bs_XwidthErr_);
	}
	evt_bs_Ywidth_branch = 0;
	if (tree->GetAlias("evt_bs_Ywidth") != 0) {
		evt_bs_Ywidth_branch = tree->GetBranch(tree->GetAlias("evt_bs_Ywidth"));
		evt_bs_Ywidth_branch->SetAddress(&evt_bs_Ywidth_);
	}
	evt_bs_YwidthErr_branch = 0;
	if (tree->GetAlias("evt_bs_YwidthErr") != 0) {
		evt_bs_YwidthErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_YwidthErr"));
		evt_bs_YwidthErr_branch->SetAddress(&evt_bs_YwidthErr_);
	}
	evt_bs_dxdz_branch = 0;
	if (tree->GetAlias("evt_bs_dxdz") != 0) {
		evt_bs_dxdz_branch = tree->GetBranch(tree->GetAlias("evt_bs_dxdz"));
		evt_bs_dxdz_branch->SetAddress(&evt_bs_dxdz_);
	}
	evt_bs_dxdzErr_branch = 0;
	if (tree->GetAlias("evt_bs_dxdzErr") != 0) {
		evt_bs_dxdzErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_dxdzErr"));
		evt_bs_dxdzErr_branch->SetAddress(&evt_bs_dxdzErr_);
	}
	evt_bs_dydz_branch = 0;
	if (tree->GetAlias("evt_bs_dydz") != 0) {
		evt_bs_dydz_branch = tree->GetBranch(tree->GetAlias("evt_bs_dydz"));
		evt_bs_dydz_branch->SetAddress(&evt_bs_dydz_);
	}
	evt_bs_dydzErr_branch = 0;
	if (tree->GetAlias("evt_bs_dydzErr") != 0) {
		evt_bs_dydzErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_dydzErr"));
		evt_bs_dydzErr_branch->SetAddress(&evt_bs_dydzErr_);
	}
	evt_bs_sigmaZ_branch = 0;
	if (tree->GetAlias("evt_bs_sigmaZ") != 0) {
		evt_bs_sigmaZ_branch = tree->GetBranch(tree->GetAlias("evt_bs_sigmaZ"));
		evt_bs_sigmaZ_branch->SetAddress(&evt_bs_sigmaZ_);
	}
	evt_bs_sigmaZErr_branch = 0;
	if (tree->GetAlias("evt_bs_sigmaZErr") != 0) {
		evt_bs_sigmaZErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_sigmaZErr"));
		evt_bs_sigmaZErr_branch->SetAddress(&evt_bs_sigmaZErr_);
	}
	evt_bs_xErr_branch = 0;
	if (tree->GetAlias("evt_bs_xErr") != 0) {
		evt_bs_xErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_xErr"));
		evt_bs_xErr_branch->SetAddress(&evt_bs_xErr_);
	}
	evt_bs_yErr_branch = 0;
	if (tree->GetAlias("evt_bs_yErr") != 0) {
		evt_bs_yErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_yErr"));
		evt_bs_yErr_branch->SetAddress(&evt_bs_yErr_);
	}
	evt_bs_zErr_branch = 0;
	if (tree->GetAlias("evt_bs_zErr") != 0) {
		evt_bs_zErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_zErr"));
		evt_bs_zErr_branch->SetAddress(&evt_bs_zErr_);
	}
	evt_bField_branch = 0;
	if (tree->GetAlias("evt_bField") != 0) {
		evt_bField_branch = tree->GetBranch(tree->GetAlias("evt_bField"));
		evt_bField_branch->SetAddress(&evt_bField_);
	}
	evt_fixgrid_all_rho_branch = 0;
	if (tree->GetAlias("evt_fixgrid_all_rho") != 0) {
		evt_fixgrid_all_rho_branch = tree->GetBranch(tree->GetAlias("evt_fixgrid_all_rho"));
		evt_fixgrid_all_rho_branch->SetAddress(&evt_fixgrid_all_rho_);
	}
	evt_fixgridfastjet_allcalo_rho_branch = 0;
	if (tree->GetAlias("evt_fixgridfastjet_allcalo_rho") != 0) {
		evt_fixgridfastjet_allcalo_rho_branch = tree->GetBranch(tree->GetAlias("evt_fixgridfastjet_allcalo_rho"));
		evt_fixgridfastjet_allcalo_rho_branch->SetAddress(&evt_fixgridfastjet_allcalo_rho_);
	}
	evt_fixgridfastjet_all_rho_branch = 0;
	if (tree->GetAlias("evt_fixgridfastjet_all_rho") != 0) {
		evt_fixgridfastjet_all_rho_branch = tree->GetBranch(tree->GetAlias("evt_fixgridfastjet_all_rho"));
		evt_fixgridfastjet_all_rho_branch->SetAddress(&evt_fixgridfastjet_all_rho_);
	}
	evt_fixgridfastjet_centralcalo_rho_branch = 0;
	if (tree->GetAlias("evt_fixgridfastjet_centralcalo_rho") != 0) {
		evt_fixgridfastjet_centralcalo_rho_branch = tree->GetBranch(tree->GetAlias("evt_fixgridfastjet_centralcalo_rho"));
		evt_fixgridfastjet_centralcalo_rho_branch->SetAddress(&evt_fixgridfastjet_centralcalo_rho_);
	}
	evt_fixgridfastjet_centralchargedpileup_rho_branch = 0;
	if (tree->GetAlias("evt_fixgridfastjet_centralchargedpileup_rho") != 0) {
		evt_fixgridfastjet_centralchargedpileup_rho_branch = tree->GetBranch(tree->GetAlias("evt_fixgridfastjet_centralchargedpileup_rho"));
		evt_fixgridfastjet_centralchargedpileup_rho_branch->SetAddress(&evt_fixgridfastjet_centralchargedpileup_rho_);
	}
	evt_fixgridfastjet_centralneutral_rho_branch = 0;
	if (tree->GetAlias("evt_fixgridfastjet_centralneutral_rho") != 0) {
		evt_fixgridfastjet_centralneutral_rho_branch = tree->GetBranch(tree->GetAlias("evt_fixgridfastjet_centralneutral_rho"));
		evt_fixgridfastjet_centralneutral_rho_branch->SetAddress(&evt_fixgridfastjet_centralneutral_rho_);
	}
	evt_kfactor_branch = 0;
	if (tree->GetAlias("evt_kfactor") != 0) {
		evt_kfactor_branch = tree->GetBranch(tree->GetAlias("evt_kfactor"));
		evt_kfactor_branch->SetAddress(&evt_kfactor_);
	}
	evt_scale1fb_branch = 0;
	if (tree->GetAlias("evt_scale1fb") != 0) {
		evt_scale1fb_branch = tree->GetBranch(tree->GetAlias("evt_scale1fb"));
		evt_scale1fb_branch->SetAddress(&evt_scale1fb_);
	}
	evt_xsec_excl_branch = 0;
	if (tree->GetAlias("evt_xsec_excl") != 0) {
		evt_xsec_excl_branch = tree->GetBranch(tree->GetAlias("evt_xsec_excl"));
		evt_xsec_excl_branch->SetAddress(&evt_xsec_excl_);
	}
	evt_xsec_incl_branch = 0;
	if (tree->GetAlias("evt_xsec_incl") != 0) {
		evt_xsec_incl_branch = tree->GetBranch(tree->GetAlias("evt_xsec_incl"));
		evt_xsec_incl_branch->SetAddress(&evt_xsec_incl_);
	}
	gen_met_branch = 0;
	if (tree->GetAlias("gen_met") != 0) {
		gen_met_branch = tree->GetBranch(tree->GetAlias("gen_met"));
		gen_met_branch->SetAddress(&gen_met_);
	}
	gen_metPhi_branch = 0;
	if (tree->GetAlias("gen_metPhi") != 0) {
		gen_metPhi_branch = tree->GetBranch(tree->GetAlias("gen_metPhi"));
		gen_metPhi_branch->SetAddress(&gen_metPhi_);
	}
	genps_alphaQCD_branch = 0;
	if (tree->GetAlias("genps_alphaQCD") != 0) {
		genps_alphaQCD_branch = tree->GetBranch(tree->GetAlias("genps_alphaQCD"));
		genps_alphaQCD_branch->SetAddress(&genps_alphaQCD_);
	}
	genps_pthat_branch = 0;
	if (tree->GetAlias("genps_pthat") != 0) {
		genps_pthat_branch = tree->GetBranch(tree->GetAlias("genps_pthat"));
		genps_pthat_branch->SetAddress(&genps_pthat_);
	}
	genps_qScale_branch = 0;
	if (tree->GetAlias("genps_qScale") != 0) {
		genps_qScale_branch = tree->GetBranch(tree->GetAlias("genps_qScale"));
		genps_qScale_branch->SetAddress(&genps_qScale_);
	}
	genps_weight_branch = 0;
	if (tree->GetAlias("genps_weight") != 0) {
		genps_weight_branch = tree->GetBranch(tree->GetAlias("genps_weight"));
		genps_weight_branch->SetAddress(&genps_weight_);
	}
	gen_sumEt_branch = 0;
	if (tree->GetAlias("gen_sumEt") != 0) {
		gen_sumEt_branch = tree->GetBranch(tree->GetAlias("gen_sumEt"));
		gen_sumEt_branch->SetAddress(&gen_sumEt_);
	}
	hcalnoise_GetRecHitEnergy_branch = 0;
	if (tree->GetAlias("hcalnoise_GetRecHitEnergy") != 0) {
		hcalnoise_GetRecHitEnergy_branch = tree->GetBranch(tree->GetAlias("hcalnoise_GetRecHitEnergy"));
		hcalnoise_GetRecHitEnergy_branch->SetAddress(&hcalnoise_GetRecHitEnergy_);
	}
	hcalnoise_GetRecHitEnergy15_branch = 0;
	if (tree->GetAlias("hcalnoise_GetRecHitEnergy15") != 0) {
		hcalnoise_GetRecHitEnergy15_branch = tree->GetBranch(tree->GetAlias("hcalnoise_GetRecHitEnergy15"));
		hcalnoise_GetRecHitEnergy15_branch->SetAddress(&hcalnoise_GetRecHitEnergy15_);
	}
	hcalnoise_GetTotalCalibCharge_branch = 0;
	if (tree->GetAlias("hcalnoise_GetTotalCalibCharge") != 0) {
		hcalnoise_GetTotalCalibCharge_branch = tree->GetBranch(tree->GetAlias("hcalnoise_GetTotalCalibCharge"));
		hcalnoise_GetTotalCalibCharge_branch->SetAddress(&hcalnoise_GetTotalCalibCharge_);
	}
	hcalnoise_TS4TS5NoiseSumE_branch = 0;
	if (tree->GetAlias("hcalnoise_TS4TS5NoiseSumE") != 0) {
		hcalnoise_TS4TS5NoiseSumE_branch = tree->GetBranch(tree->GetAlias("hcalnoise_TS4TS5NoiseSumE"));
		hcalnoise_TS4TS5NoiseSumE_branch->SetAddress(&hcalnoise_TS4TS5NoiseSumE_);
	}
	hcalnoise_TS4TS5NoiseSumEt_branch = 0;
	if (tree->GetAlias("hcalnoise_TS4TS5NoiseSumEt") != 0) {
		hcalnoise_TS4TS5NoiseSumEt_branch = tree->GetBranch(tree->GetAlias("hcalnoise_TS4TS5NoiseSumEt"));
		hcalnoise_TS4TS5NoiseSumEt_branch->SetAddress(&hcalnoise_TS4TS5NoiseSumEt_);
	}
	hcalnoise_eventChargeFraction_branch = 0;
	if (tree->GetAlias("hcalnoise_eventChargeFraction") != 0) {
		hcalnoise_eventChargeFraction_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventChargeFraction"));
		hcalnoise_eventChargeFraction_branch->SetAddress(&hcalnoise_eventChargeFraction_);
	}
	hcalnoise_eventEMEnergy_branch = 0;
	if (tree->GetAlias("hcalnoise_eventEMEnergy") != 0) {
		hcalnoise_eventEMEnergy_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventEMEnergy"));
		hcalnoise_eventEMEnergy_branch->SetAddress(&hcalnoise_eventEMEnergy_);
	}
	hcalnoise_eventEMFraction_branch = 0;
	if (tree->GetAlias("hcalnoise_eventEMFraction") != 0) {
		hcalnoise_eventEMFraction_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventEMFraction"));
		hcalnoise_eventEMFraction_branch->SetAddress(&hcalnoise_eventEMFraction_);
	}
	hcalnoise_eventHadEnergy_branch = 0;
	if (tree->GetAlias("hcalnoise_eventHadEnergy") != 0) {
		hcalnoise_eventHadEnergy_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventHadEnergy"));
		hcalnoise_eventHadEnergy_branch->SetAddress(&hcalnoise_eventHadEnergy_);
	}
	hcalnoise_eventTrackEnergy_branch = 0;
	if (tree->GetAlias("hcalnoise_eventTrackEnergy") != 0) {
		hcalnoise_eventTrackEnergy_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventTrackEnergy"));
		hcalnoise_eventTrackEnergy_branch->SetAddress(&hcalnoise_eventTrackEnergy_);
	}
	hcalnoise_flatNoiseSumE_branch = 0;
	if (tree->GetAlias("hcalnoise_flatNoiseSumE") != 0) {
		hcalnoise_flatNoiseSumE_branch = tree->GetBranch(tree->GetAlias("hcalnoise_flatNoiseSumE"));
		hcalnoise_flatNoiseSumE_branch->SetAddress(&hcalnoise_flatNoiseSumE_);
	}
	hcalnoise_flatNoiseSumEt_branch = 0;
	if (tree->GetAlias("hcalnoise_flatNoiseSumEt") != 0) {
		hcalnoise_flatNoiseSumEt_branch = tree->GetBranch(tree->GetAlias("hcalnoise_flatNoiseSumEt"));
		hcalnoise_flatNoiseSumEt_branch->SetAddress(&hcalnoise_flatNoiseSumEt_);
	}
	hcalnoise_isolatedNoiseSumE_branch = 0;
	if (tree->GetAlias("hcalnoise_isolatedNoiseSumE") != 0) {
		hcalnoise_isolatedNoiseSumE_branch = tree->GetBranch(tree->GetAlias("hcalnoise_isolatedNoiseSumE"));
		hcalnoise_isolatedNoiseSumE_branch->SetAddress(&hcalnoise_isolatedNoiseSumE_);
	}
	hcalnoise_isolatedNoiseSumEt_branch = 0;
	if (tree->GetAlias("hcalnoise_isolatedNoiseSumEt") != 0) {
		hcalnoise_isolatedNoiseSumEt_branch = tree->GetBranch(tree->GetAlias("hcalnoise_isolatedNoiseSumEt"));
		hcalnoise_isolatedNoiseSumEt_branch->SetAddress(&hcalnoise_isolatedNoiseSumEt_);
	}
	hcalnoise_max10GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_max10GeVHitTime") != 0) {
		hcalnoise_max10GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_max10GeVHitTime"));
		hcalnoise_max10GeVHitTime_branch->SetAddress(&hcalnoise_max10GeVHitTime_);
	}
	hcalnoise_max25GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_max25GeVHitTime") != 0) {
		hcalnoise_max25GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_max25GeVHitTime"));
		hcalnoise_max25GeVHitTime_branch->SetAddress(&hcalnoise_max25GeVHitTime_);
	}
	hcalnoise_maxE10TS_branch = 0;
	if (tree->GetAlias("hcalnoise_maxE10TS") != 0) {
		hcalnoise_maxE10TS_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxE10TS"));
		hcalnoise_maxE10TS_branch->SetAddress(&hcalnoise_maxE10TS_);
	}
	hcalnoise_maxE2Over10TS_branch = 0;
	if (tree->GetAlias("hcalnoise_maxE2Over10TS") != 0) {
		hcalnoise_maxE2Over10TS_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxE2Over10TS"));
		hcalnoise_maxE2Over10TS_branch->SetAddress(&hcalnoise_maxE2Over10TS_);
	}
	hcalnoise_maxE2TS_branch = 0;
	if (tree->GetAlias("hcalnoise_maxE2TS") != 0) {
		hcalnoise_maxE2TS_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxE2TS"));
		hcalnoise_maxE2TS_branch->SetAddress(&hcalnoise_maxE2TS_);
	}
	hcalnoise_min10GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_min10GeVHitTime") != 0) {
		hcalnoise_min10GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_min10GeVHitTime"));
		hcalnoise_min10GeVHitTime_branch->SetAddress(&hcalnoise_min10GeVHitTime_);
	}
	hcalnoise_min25GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_min25GeVHitTime") != 0) {
		hcalnoise_min25GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_min25GeVHitTime"));
		hcalnoise_min25GeVHitTime_branch->SetAddress(&hcalnoise_min25GeVHitTime_);
	}
	hcalnoise_minE10TS_branch = 0;
	if (tree->GetAlias("hcalnoise_minE10TS") != 0) {
		hcalnoise_minE10TS_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minE10TS"));
		hcalnoise_minE10TS_branch->SetAddress(&hcalnoise_minE10TS_);
	}
	hcalnoise_minE2Over10TS_branch = 0;
	if (tree->GetAlias("hcalnoise_minE2Over10TS") != 0) {
		hcalnoise_minE2Over10TS_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minE2Over10TS"));
		hcalnoise_minE2Over10TS_branch->SetAddress(&hcalnoise_minE2Over10TS_);
	}
	hcalnoise_minE2TS_branch = 0;
	if (tree->GetAlias("hcalnoise_minE2TS") != 0) {
		hcalnoise_minE2TS_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minE2TS"));
		hcalnoise_minE2TS_branch->SetAddress(&hcalnoise_minE2TS_);
	}
	hcalnoise_minHPDEMF_branch = 0;
	if (tree->GetAlias("hcalnoise_minHPDEMF") != 0) {
		hcalnoise_minHPDEMF_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minHPDEMF"));
		hcalnoise_minHPDEMF_branch->SetAddress(&hcalnoise_minHPDEMF_);
	}
	hcalnoise_minRBXEMF_branch = 0;
	if (tree->GetAlias("hcalnoise_minRBXEMF") != 0) {
		hcalnoise_minRBXEMF_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minRBXEMF"));
		hcalnoise_minRBXEMF_branch->SetAddress(&hcalnoise_minRBXEMF_);
	}
	hcalnoise_rms10GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_rms10GeVHitTime") != 0) {
		hcalnoise_rms10GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_rms10GeVHitTime"));
		hcalnoise_rms10GeVHitTime_branch->SetAddress(&hcalnoise_rms10GeVHitTime_);
	}
	hcalnoise_rms25GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_rms25GeVHitTime") != 0) {
		hcalnoise_rms25GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_rms25GeVHitTime"));
		hcalnoise_rms25GeVHitTime_branch->SetAddress(&hcalnoise_rms25GeVHitTime_);
	}
	hcalnoise_spikeNoiseSumE_branch = 0;
	if (tree->GetAlias("hcalnoise_spikeNoiseSumE") != 0) {
		hcalnoise_spikeNoiseSumE_branch = tree->GetBranch(tree->GetAlias("hcalnoise_spikeNoiseSumE"));
		hcalnoise_spikeNoiseSumE_branch->SetAddress(&hcalnoise_spikeNoiseSumE_);
	}
	hcalnoise_spikeNoiseSumEt_branch = 0;
	if (tree->GetAlias("hcalnoise_spikeNoiseSumEt") != 0) {
		hcalnoise_spikeNoiseSumEt_branch = tree->GetBranch(tree->GetAlias("hcalnoise_spikeNoiseSumEt"));
		hcalnoise_spikeNoiseSumEt_branch->SetAddress(&hcalnoise_spikeNoiseSumEt_);
	}
	hcalnoise_triangleNoiseSumE_branch = 0;
	if (tree->GetAlias("hcalnoise_triangleNoiseSumE") != 0) {
		hcalnoise_triangleNoiseSumE_branch = tree->GetBranch(tree->GetAlias("hcalnoise_triangleNoiseSumE"));
		hcalnoise_triangleNoiseSumE_branch->SetAddress(&hcalnoise_triangleNoiseSumE_);
	}
	hcalnoise_triangleNoiseSumEt_branch = 0;
	if (tree->GetAlias("hcalnoise_triangleNoiseSumEt") != 0) {
		hcalnoise_triangleNoiseSumEt_branch = tree->GetBranch(tree->GetAlias("hcalnoise_triangleNoiseSumEt"));
		hcalnoise_triangleNoiseSumEt_branch->SetAddress(&hcalnoise_triangleNoiseSumEt_);
	}
	pdfinfo_pdf1_branch = 0;
	if (tree->GetAlias("pdfinfo_pdf1") != 0) {
		pdfinfo_pdf1_branch = tree->GetBranch(tree->GetAlias("pdfinfo_pdf1"));
		pdfinfo_pdf1_branch->SetAddress(&pdfinfo_pdf1_);
	}
	pdfinfo_pdf2_branch = 0;
	if (tree->GetAlias("pdfinfo_pdf2") != 0) {
		pdfinfo_pdf2_branch = tree->GetBranch(tree->GetAlias("pdfinfo_pdf2"));
		pdfinfo_pdf2_branch->SetAddress(&pdfinfo_pdf2_);
	}
	pdfinfo_scale_branch = 0;
	if (tree->GetAlias("pdfinfo_scale") != 0) {
		pdfinfo_scale_branch = tree->GetBranch(tree->GetAlias("pdfinfo_scale"));
		pdfinfo_scale_branch->SetAddress(&pdfinfo_scale_);
	}
	pdfinfo_x1_branch = 0;
	if (tree->GetAlias("pdfinfo_x1") != 0) {
		pdfinfo_x1_branch = tree->GetBranch(tree->GetAlias("pdfinfo_x1"));
		pdfinfo_x1_branch->SetAddress(&pdfinfo_x1_);
	}
	pdfinfo_x2_branch = 0;
	if (tree->GetAlias("pdfinfo_x2") != 0) {
		pdfinfo_x2_branch = tree->GetBranch(tree->GetAlias("pdfinfo_x2"));
		pdfinfo_x2_branch->SetAddress(&pdfinfo_x2_);
	}
	evt_fixgrid_rho_all_branch = 0;
	if (tree->GetAlias("evt_fixgrid_rho_all") != 0) {
		evt_fixgrid_rho_all_branch = tree->GetBranch(tree->GetAlias("evt_fixgrid_rho_all"));
		evt_fixgrid_rho_all_branch->SetAddress(&evt_fixgrid_rho_all_);
	}
	evt_fixgrid_rho_ctr_branch = 0;
	if (tree->GetAlias("evt_fixgrid_rho_ctr") != 0) {
		evt_fixgrid_rho_ctr_branch = tree->GetBranch(tree->GetAlias("evt_fixgrid_rho_ctr"));
		evt_fixgrid_rho_ctr_branch->SetAddress(&evt_fixgrid_rho_ctr_);
	}
	evt_fixgrid_rho_fwd_branch = 0;
	if (tree->GetAlias("evt_fixgrid_rho_fwd") != 0) {
		evt_fixgrid_rho_fwd_branch = tree->GetBranch(tree->GetAlias("evt_fixgrid_rho_fwd"));
		evt_fixgrid_rho_fwd_branch->SetAddress(&evt_fixgrid_rho_fwd_);
	}
	evt_pfmet_branch = 0;
	if (tree->GetAlias("evt_pfmet") != 0) {
		evt_pfmet_branch = tree->GetBranch(tree->GetAlias("evt_pfmet"));
		evt_pfmet_branch->SetAddress(&evt_pfmet_);
	}
	evt_pfmetPhi_branch = 0;
	if (tree->GetAlias("evt_pfmetPhi") != 0) {
		evt_pfmetPhi_branch = tree->GetBranch(tree->GetAlias("evt_pfmetPhi"));
		evt_pfmetPhi_branch->SetAddress(&evt_pfmetPhi_);
	}
	evt_pfmetSig_branch = 0;
	if (tree->GetAlias("evt_pfmetSig") != 0) {
		evt_pfmetSig_branch = tree->GetBranch(tree->GetAlias("evt_pfmetSig"));
		evt_pfmetSig_branch->SetAddress(&evt_pfmetSig_);
	}
	evt_pfmetSignificance_branch = 0;
	if (tree->GetAlias("evt_pfmetSignificance") != 0) {
		evt_pfmetSignificance_branch = tree->GetBranch(tree->GetAlias("evt_pfmetSignificance"));
		evt_pfmetSignificance_branch->SetAddress(&evt_pfmetSignificance_);
	}
	evt_pfsumet_branch = 0;
	if (tree->GetAlias("evt_pfsumet") != 0) {
		evt_pfsumet_branch = tree->GetBranch(tree->GetAlias("evt_pfsumet"));
		evt_pfsumet_branch->SetAddress(&evt_pfsumet_);
	}
	genps_lepdaughter_p4_branch = 0;
	if (tree->GetAlias("genps_lepdaughter_p4") != 0) {
		genps_lepdaughter_p4_branch = tree->GetBranch(tree->GetAlias("genps_lepdaughter_p4"));
		genps_lepdaughter_p4_branch->SetAddress(&genps_lepdaughter_p4_);
	}
	hlt_trigObjs_p4_branch = 0;
	if (tree->GetAlias("hlt_trigObjs_p4") != 0) {
		hlt_trigObjs_p4_branch = tree->GetBranch(tree->GetAlias("hlt_trigObjs_p4"));
		hlt_trigObjs_p4_branch->SetAddress(&hlt_trigObjs_p4_);
	}
	taus_pf_isocands_p4_branch = 0;
	if (tree->GetAlias("taus_pf_isocands_p4") != 0) {
		taus_pf_isocands_p4_branch = tree->GetBranch(tree->GetAlias("taus_pf_isocands_p4"));
		taus_pf_isocands_p4_branch->SetAddress(&taus_pf_isocands_p4_);
	}
	taus_pf_signalcands_p4_branch = 0;
	if (tree->GetAlias("taus_pf_signalcands_p4") != 0) {
		taus_pf_signalcands_p4_branch = tree->GetBranch(tree->GetAlias("taus_pf_signalcands_p4"));
		taus_pf_signalcands_p4_branch->SetAddress(&taus_pf_signalcands_p4_);
	}
	evt_bs_covMatrix_branch = 0;
	if (tree->GetAlias("evt_bs_covMatrix") != 0) {
		evt_bs_covMatrix_branch = tree->GetBranch(tree->GetAlias("evt_bs_covMatrix"));
		evt_bs_covMatrix_branch->SetAddress(&evt_bs_covMatrix_);
	}
	els_mc3dr_branch = 0;
	if (tree->GetAlias("els_mc3dr") != 0) {
		els_mc3dr_branch = tree->GetBranch(tree->GetAlias("els_mc3dr"));
		els_mc3dr_branch->SetAddress(&els_mc3dr_);
	}
	els_mcdr_branch = 0;
	if (tree->GetAlias("els_mcdr") != 0) {
		els_mcdr_branch = tree->GetBranch(tree->GetAlias("els_mcdr"));
		els_mcdr_branch->SetAddress(&els_mcdr_);
	}
	mus_mc3dr_branch = 0;
	if (tree->GetAlias("mus_mc3dr") != 0) {
		mus_mc3dr_branch = tree->GetBranch(tree->GetAlias("mus_mc3dr"));
		mus_mc3dr_branch->SetAddress(&mus_mc3dr_);
	}
	mus_mcdr_branch = 0;
	if (tree->GetAlias("mus_mcdr") != 0) {
		mus_mcdr_branch = tree->GetBranch(tree->GetAlias("mus_mcdr"));
		mus_mcdr_branch->SetAddress(&mus_mcdr_);
	}
	pfjets_mc3dr_branch = 0;
	if (tree->GetAlias("pfjets_mc3dr") != 0) {
		pfjets_mc3dr_branch = tree->GetBranch(tree->GetAlias("pfjets_mc3dr"));
		pfjets_mc3dr_branch->SetAddress(&pfjets_mc3dr_);
	}
	pfjets_mcdr_branch = 0;
	if (tree->GetAlias("pfjets_mcdr") != 0) {
		pfjets_mcdr_branch = tree->GetBranch(tree->GetAlias("pfjets_mcdr"));
		pfjets_mcdr_branch->SetAddress(&pfjets_mcdr_);
	}
	pfjets_mc_emEnergy_branch = 0;
	if (tree->GetAlias("pfjets_mc_emEnergy") != 0) {
		pfjets_mc_emEnergy_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_emEnergy"));
		pfjets_mc_emEnergy_branch->SetAddress(&pfjets_mc_emEnergy_);
	}
	pfjets_mc_gpdr_branch = 0;
	if (tree->GetAlias("pfjets_mc_gpdr") != 0) {
		pfjets_mc_gpdr_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_gpdr"));
		pfjets_mc_gpdr_branch->SetAddress(&pfjets_mc_gpdr_);
	}
	pfjets_mc_hadEnergy_branch = 0;
	if (tree->GetAlias("pfjets_mc_hadEnergy") != 0) {
		pfjets_mc_hadEnergy_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_hadEnergy"));
		pfjets_mc_hadEnergy_branch->SetAddress(&pfjets_mc_hadEnergy_);
	}
	pfjets_mc_invEnergy_branch = 0;
	if (tree->GetAlias("pfjets_mc_invEnergy") != 0) {
		pfjets_mc_invEnergy_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_invEnergy"));
		pfjets_mc_invEnergy_branch->SetAddress(&pfjets_mc_invEnergy_);
	}
	pfjets_mc_otherEnergy_branch = 0;
	if (tree->GetAlias("pfjets_mc_otherEnergy") != 0) {
		pfjets_mc_otherEnergy_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_otherEnergy"));
		pfjets_mc_otherEnergy_branch->SetAddress(&pfjets_mc_otherEnergy_);
	}
	photons_mc3dr_branch = 0;
	if (tree->GetAlias("photons_mc3dr") != 0) {
		photons_mc3dr_branch = tree->GetBranch(tree->GetAlias("photons_mc3dr"));
		photons_mc3dr_branch->SetAddress(&photons_mc3dr_);
	}
	photons_mcdr_branch = 0;
	if (tree->GetAlias("photons_mcdr") != 0) {
		photons_mcdr_branch = tree->GetBranch(tree->GetAlias("photons_mcdr"));
		photons_mcdr_branch->SetAddress(&photons_mcdr_);
	}
	els_chi2_branch = 0;
	if (tree->GetAlias("els_chi2") != 0) {
		els_chi2_branch = tree->GetBranch(tree->GetAlias("els_chi2"));
		els_chi2_branch->SetAddress(&els_chi2_);
	}
	els_ckf_chi2_branch = 0;
	if (tree->GetAlias("els_ckf_chi2") != 0) {
		els_ckf_chi2_branch = tree->GetBranch(tree->GetAlias("els_ckf_chi2"));
		els_ckf_chi2_branch->SetAddress(&els_ckf_chi2_);
	}
	els_ckf_ndof_branch = 0;
	if (tree->GetAlias("els_ckf_ndof") != 0) {
		els_ckf_ndof_branch = tree->GetBranch(tree->GetAlias("els_ckf_ndof"));
		els_ckf_ndof_branch->SetAddress(&els_ckf_ndof_);
	}
	els_conv_dcot_branch = 0;
	if (tree->GetAlias("els_conv_dcot") != 0) {
		els_conv_dcot_branch = tree->GetBranch(tree->GetAlias("els_conv_dcot"));
		els_conv_dcot_branch->SetAddress(&els_conv_dcot_);
	}
	els_conv_dist_branch = 0;
	if (tree->GetAlias("els_conv_dist") != 0) {
		els_conv_dist_branch = tree->GetBranch(tree->GetAlias("els_conv_dist"));
		els_conv_dist_branch->SetAddress(&els_conv_dist_);
	}
	els_conv_old_dcot_branch = 0;
	if (tree->GetAlias("els_conv_old_dcot") != 0) {
		els_conv_old_dcot_branch = tree->GetBranch(tree->GetAlias("els_conv_old_dcot"));
		els_conv_old_dcot_branch->SetAddress(&els_conv_old_dcot_);
	}
	els_conv_old_dist_branch = 0;
	if (tree->GetAlias("els_conv_old_dist") != 0) {
		els_conv_old_dist_branch = tree->GetBranch(tree->GetAlias("els_conv_old_dist"));
		els_conv_old_dist_branch->SetAddress(&els_conv_old_dist_);
	}
	els_conv_old_radius_branch = 0;
	if (tree->GetAlias("els_conv_old_radius") != 0) {
		els_conv_old_radius_branch = tree->GetBranch(tree->GetAlias("els_conv_old_radius"));
		els_conv_old_radius_branch->SetAddress(&els_conv_old_radius_);
	}
	els_conv_radius_branch = 0;
	if (tree->GetAlias("els_conv_radius") != 0) {
		els_conv_radius_branch = tree->GetBranch(tree->GetAlias("els_conv_radius"));
		els_conv_radius_branch->SetAddress(&els_conv_radius_);
	}
	els_d0_branch = 0;
	if (tree->GetAlias("els_d0") != 0) {
		els_d0_branch = tree->GetBranch(tree->GetAlias("els_d0"));
		els_d0_branch->SetAddress(&els_d0_);
	}
	els_d0Err_branch = 0;
	if (tree->GetAlias("els_d0Err") != 0) {
		els_d0Err_branch = tree->GetBranch(tree->GetAlias("els_d0Err"));
		els_d0Err_branch->SetAddress(&els_d0Err_);
	}
	els_d0corr_branch = 0;
	if (tree->GetAlias("els_d0corr") != 0) {
		els_d0corr_branch = tree->GetBranch(tree->GetAlias("els_d0corr"));
		els_d0corr_branch->SetAddress(&els_d0corr_);
	}
	els_d0corrPhi_branch = 0;
	if (tree->GetAlias("els_d0corrPhi") != 0) {
		els_d0corrPhi_branch = tree->GetBranch(tree->GetAlias("els_d0corrPhi"));
		els_d0corrPhi_branch->SetAddress(&els_d0corrPhi_);
	}
	els_d0phiCov_branch = 0;
	if (tree->GetAlias("els_d0phiCov") != 0) {
		els_d0phiCov_branch = tree->GetBranch(tree->GetAlias("els_d0phiCov"));
		els_d0phiCov_branch->SetAddress(&els_d0phiCov_);
	}
	els_dEtaIn_branch = 0;
	if (tree->GetAlias("els_dEtaIn") != 0) {
		els_dEtaIn_branch = tree->GetBranch(tree->GetAlias("els_dEtaIn"));
		els_dEtaIn_branch->SetAddress(&els_dEtaIn_);
	}
	els_dEtaOut_branch = 0;
	if (tree->GetAlias("els_dEtaOut") != 0) {
		els_dEtaOut_branch = tree->GetBranch(tree->GetAlias("els_dEtaOut"));
		els_dEtaOut_branch->SetAddress(&els_dEtaOut_);
	}
	els_dPhiIn_branch = 0;
	if (tree->GetAlias("els_dPhiIn") != 0) {
		els_dPhiIn_branch = tree->GetBranch(tree->GetAlias("els_dPhiIn"));
		els_dPhiIn_branch->SetAddress(&els_dPhiIn_);
	}
	els_dPhiInPhiOut_branch = 0;
	if (tree->GetAlias("els_dPhiInPhiOut") != 0) {
		els_dPhiInPhiOut_branch = tree->GetBranch(tree->GetAlias("els_dPhiInPhiOut"));
		els_dPhiInPhiOut_branch->SetAddress(&els_dPhiInPhiOut_);
	}
	els_dPhiOut_branch = 0;
	if (tree->GetAlias("els_dPhiOut") != 0) {
		els_dPhiOut_branch = tree->GetBranch(tree->GetAlias("els_dPhiOut"));
		els_dPhiOut_branch->SetAddress(&els_dPhiOut_);
	}
	els_deltaEtaEleClusterTrackAtCalo_branch = 0;
	if (tree->GetAlias("els_deltaEtaEleClusterTrackAtCalo") != 0) {
		els_deltaEtaEleClusterTrackAtCalo_branch = tree->GetBranch(tree->GetAlias("els_deltaEtaEleClusterTrackAtCalo"));
		els_deltaEtaEleClusterTrackAtCalo_branch->SetAddress(&els_deltaEtaEleClusterTrackAtCalo_);
	}
	els_deltaPhiEleClusterTrackAtCalo_branch = 0;
	if (tree->GetAlias("els_deltaPhiEleClusterTrackAtCalo") != 0) {
		els_deltaPhiEleClusterTrackAtCalo_branch = tree->GetBranch(tree->GetAlias("els_deltaPhiEleClusterTrackAtCalo"));
		els_deltaPhiEleClusterTrackAtCalo_branch->SetAddress(&els_deltaPhiEleClusterTrackAtCalo_);
	}
	els_dxyPV_branch = 0;
	if (tree->GetAlias("els_dxyPV") != 0) {
		els_dxyPV_branch = tree->GetBranch(tree->GetAlias("els_dxyPV"));
		els_dxyPV_branch->SetAddress(&els_dxyPV_);
	}
	els_dzPV_branch = 0;
	if (tree->GetAlias("els_dzPV") != 0) {
		els_dzPV_branch = tree->GetBranch(tree->GetAlias("els_dzPV"));
		els_dzPV_branch->SetAddress(&els_dzPV_);
	}
	els_e1x5_branch = 0;
	if (tree->GetAlias("els_e1x5") != 0) {
		els_e1x5_branch = tree->GetBranch(tree->GetAlias("els_e1x5"));
		els_e1x5_branch->SetAddress(&els_e1x5_);
	}
	els_e1x5_full5x5_branch = 0;
	if (tree->GetAlias("els_e1x5_full5x5") != 0) {
		els_e1x5_full5x5_branch = tree->GetBranch(tree->GetAlias("els_e1x5_full5x5"));
		els_e1x5_full5x5_branch->SetAddress(&els_e1x5_full5x5_);
	}
	els_e2x5Max_branch = 0;
	if (tree->GetAlias("els_e2x5Max") != 0) {
		els_e2x5Max_branch = tree->GetBranch(tree->GetAlias("els_e2x5Max"));
		els_e2x5Max_branch->SetAddress(&els_e2x5Max_);
	}
	els_e2x5Max_full5x5_branch = 0;
	if (tree->GetAlias("els_e2x5Max_full5x5") != 0) {
		els_e2x5Max_full5x5_branch = tree->GetBranch(tree->GetAlias("els_e2x5Max_full5x5"));
		els_e2x5Max_full5x5_branch->SetAddress(&els_e2x5Max_full5x5_);
	}
	els_e3x3_branch = 0;
	if (tree->GetAlias("els_e3x3") != 0) {
		els_e3x3_branch = tree->GetBranch(tree->GetAlias("els_e3x3"));
		els_e3x3_branch->SetAddress(&els_e3x3_);
	}
	els_e5x5_branch = 0;
	if (tree->GetAlias("els_e5x5") != 0) {
		els_e5x5_branch = tree->GetBranch(tree->GetAlias("els_e5x5"));
		els_e5x5_branch->SetAddress(&els_e5x5_);
	}
	els_e5x5_full5x5_branch = 0;
	if (tree->GetAlias("els_e5x5_full5x5") != 0) {
		els_e5x5_full5x5_branch = tree->GetBranch(tree->GetAlias("els_e5x5_full5x5"));
		els_e5x5_full5x5_branch->SetAddress(&els_e5x5_full5x5_);
	}
	els_eMax_branch = 0;
	if (tree->GetAlias("els_eMax") != 0) {
		els_eMax_branch = tree->GetBranch(tree->GetAlias("els_eMax"));
		els_eMax_branch->SetAddress(&els_eMax_);
	}
	els_eOverPIn_branch = 0;
	if (tree->GetAlias("els_eOverPIn") != 0) {
		els_eOverPIn_branch = tree->GetBranch(tree->GetAlias("els_eOverPIn"));
		els_eOverPIn_branch->SetAddress(&els_eOverPIn_);
	}
	els_eOverPOut_branch = 0;
	if (tree->GetAlias("els_eOverPOut") != 0) {
		els_eOverPOut_branch = tree->GetBranch(tree->GetAlias("els_eOverPOut"));
		els_eOverPOut_branch->SetAddress(&els_eOverPOut_);
	}
	els_eSC_branch = 0;
	if (tree->GetAlias("els_eSC") != 0) {
		els_eSC_branch = tree->GetBranch(tree->GetAlias("els_eSC"));
		els_eSC_branch->SetAddress(&els_eSC_);
	}
	els_eSCPresh_branch = 0;
	if (tree->GetAlias("els_eSCPresh") != 0) {
		els_eSCPresh_branch = tree->GetBranch(tree->GetAlias("els_eSCPresh"));
		els_eSCPresh_branch->SetAddress(&els_eSCPresh_);
	}
	els_eSCRaw_branch = 0;
	if (tree->GetAlias("els_eSCRaw") != 0) {
		els_eSCRaw_branch = tree->GetBranch(tree->GetAlias("els_eSCRaw"));
		els_eSCRaw_branch->SetAddress(&els_eSCRaw_);
	}
	els_eSeed_branch = 0;
	if (tree->GetAlias("els_eSeed") != 0) {
		els_eSeed_branch = tree->GetBranch(tree->GetAlias("els_eSeed"));
		els_eSeed_branch->SetAddress(&els_eSeed_);
	}
	els_eSeedOverPIn_branch = 0;
	if (tree->GetAlias("els_eSeedOverPIn") != 0) {
		els_eSeedOverPIn_branch = tree->GetBranch(tree->GetAlias("els_eSeedOverPIn"));
		els_eSeedOverPIn_branch->SetAddress(&els_eSeedOverPIn_);
	}
	els_eSeedOverPOut_branch = 0;
	if (tree->GetAlias("els_eSeedOverPOut") != 0) {
		els_eSeedOverPOut_branch = tree->GetBranch(tree->GetAlias("els_eSeedOverPOut"));
		els_eSeedOverPOut_branch->SetAddress(&els_eSeedOverPOut_);
	}
	els_ecalEnergy_branch = 0;
	if (tree->GetAlias("els_ecalEnergy") != 0) {
		els_ecalEnergy_branch = tree->GetBranch(tree->GetAlias("els_ecalEnergy"));
		els_ecalEnergy_branch->SetAddress(&els_ecalEnergy_);
	}
	els_ecalEnergyError_branch = 0;
	if (tree->GetAlias("els_ecalEnergyError") != 0) {
		els_ecalEnergyError_branch = tree->GetBranch(tree->GetAlias("els_ecalEnergyError"));
		els_ecalEnergyError_branch->SetAddress(&els_ecalEnergyError_);
	}
	els_ecalIso_branch = 0;
	if (tree->GetAlias("els_ecalIso") != 0) {
		els_ecalIso_branch = tree->GetBranch(tree->GetAlias("els_ecalIso"));
		els_ecalIso_branch->SetAddress(&els_ecalIso_);
	}
	els_ecalIso04_branch = 0;
	if (tree->GetAlias("els_ecalIso04") != 0) {
		els_ecalIso04_branch = tree->GetBranch(tree->GetAlias("els_ecalIso04"));
		els_ecalIso04_branch->SetAddress(&els_ecalIso04_);
	}
	els_etaErr_branch = 0;
	if (tree->GetAlias("els_etaErr") != 0) {
		els_etaErr_branch = tree->GetBranch(tree->GetAlias("els_etaErr"));
		els_etaErr_branch->SetAddress(&els_etaErr_);
	}
	els_etaSC_branch = 0;
	if (tree->GetAlias("els_etaSC") != 0) {
		els_etaSC_branch = tree->GetBranch(tree->GetAlias("els_etaSC"));
		els_etaSC_branch->SetAddress(&els_etaSC_);
	}
	els_etaSCwidth_branch = 0;
	if (tree->GetAlias("els_etaSCwidth") != 0) {
		els_etaSCwidth_branch = tree->GetBranch(tree->GetAlias("els_etaSCwidth"));
		els_etaSCwidth_branch->SetAddress(&els_etaSCwidth_);
	}
	els_fbrem_branch = 0;
	if (tree->GetAlias("els_fbrem") != 0) {
		els_fbrem_branch = tree->GetBranch(tree->GetAlias("els_fbrem"));
		els_fbrem_branch->SetAddress(&els_fbrem_);
	}
	els_hOverE_branch = 0;
	if (tree->GetAlias("els_hOverE") != 0) {
		els_hOverE_branch = tree->GetBranch(tree->GetAlias("els_hOverE"));
		els_hOverE_branch->SetAddress(&els_hOverE_);
	}
	els_hOverEBC_branch = 0;
	if (tree->GetAlias("els_hOverEBC") != 0) {
		els_hOverEBC_branch = tree->GetBranch(tree->GetAlias("els_hOverEBC"));
		els_hOverEBC_branch->SetAddress(&els_hOverEBC_);
	}
	els_hcalDepth1OverEcal_branch = 0;
	if (tree->GetAlias("els_hcalDepth1OverEcal") != 0) {
		els_hcalDepth1OverEcal_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth1OverEcal"));
		els_hcalDepth1OverEcal_branch->SetAddress(&els_hcalDepth1OverEcal_);
	}
	els_hcalDepth1TowerSumEt_branch = 0;
	if (tree->GetAlias("els_hcalDepth1TowerSumEt") != 0) {
		els_hcalDepth1TowerSumEt_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth1TowerSumEt"));
		els_hcalDepth1TowerSumEt_branch->SetAddress(&els_hcalDepth1TowerSumEt_);
	}
	els_hcalDepth1TowerSumEt04_branch = 0;
	if (tree->GetAlias("els_hcalDepth1TowerSumEt04") != 0) {
		els_hcalDepth1TowerSumEt04_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth1TowerSumEt04"));
		els_hcalDepth1TowerSumEt04_branch->SetAddress(&els_hcalDepth1TowerSumEt04_);
	}
	els_hcalDepth2OverEcal_branch = 0;
	if (tree->GetAlias("els_hcalDepth2OverEcal") != 0) {
		els_hcalDepth2OverEcal_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth2OverEcal"));
		els_hcalDepth2OverEcal_branch->SetAddress(&els_hcalDepth2OverEcal_);
	}
	els_hcalDepth2TowerSumEt_branch = 0;
	if (tree->GetAlias("els_hcalDepth2TowerSumEt") != 0) {
		els_hcalDepth2TowerSumEt_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth2TowerSumEt"));
		els_hcalDepth2TowerSumEt_branch->SetAddress(&els_hcalDepth2TowerSumEt_);
	}
	els_hcalDepth2TowerSumEt04_branch = 0;
	if (tree->GetAlias("els_hcalDepth2TowerSumEt04") != 0) {
		els_hcalDepth2TowerSumEt04_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth2TowerSumEt04"));
		els_hcalDepth2TowerSumEt04_branch->SetAddress(&els_hcalDepth2TowerSumEt04_);
	}
	els_hcalIso_branch = 0;
	if (tree->GetAlias("els_hcalIso") != 0) {
		els_hcalIso_branch = tree->GetBranch(tree->GetAlias("els_hcalIso"));
		els_hcalIso_branch->SetAddress(&els_hcalIso_);
	}
	els_hcalIso04_branch = 0;
	if (tree->GetAlias("els_hcalIso04") != 0) {
		els_hcalIso04_branch = tree->GetBranch(tree->GetAlias("els_hcalIso04"));
		els_hcalIso04_branch->SetAddress(&els_hcalIso04_);
	}
	els_ip3d_branch = 0;
	if (tree->GetAlias("els_ip3d") != 0) {
		els_ip3d_branch = tree->GetBranch(tree->GetAlias("els_ip3d"));
		els_ip3d_branch->SetAddress(&els_ip3d_);
	}
	els_ip3derr_branch = 0;
	if (tree->GetAlias("els_ip3derr") != 0) {
		els_ip3derr_branch = tree->GetBranch(tree->GetAlias("els_ip3derr"));
		els_ip3derr_branch->SetAddress(&els_ip3derr_);
	}
	els_iso03_pf_branch = 0;
	if (tree->GetAlias("els_iso03_pf") != 0) {
		els_iso03_pf_branch = tree->GetBranch(tree->GetAlias("els_iso03_pf"));
		els_iso03_pf_branch->SetAddress(&els_iso03_pf_);
	}
	els_iso03_pf2012_ch_branch = 0;
	if (tree->GetAlias("els_iso03_pf2012_ch") != 0) {
		els_iso03_pf2012_ch_branch = tree->GetBranch(tree->GetAlias("els_iso03_pf2012_ch"));
		els_iso03_pf2012_ch_branch->SetAddress(&els_iso03_pf2012_ch_);
	}
	els_iso03_pf2012_em_branch = 0;
	if (tree->GetAlias("els_iso03_pf2012_em") != 0) {
		els_iso03_pf2012_em_branch = tree->GetBranch(tree->GetAlias("els_iso03_pf2012_em"));
		els_iso03_pf2012_em_branch->SetAddress(&els_iso03_pf2012_em_);
	}
	els_iso03_pf2012ext_ch_branch = 0;
	if (tree->GetAlias("els_iso03_pf2012ext_ch") != 0) {
		els_iso03_pf2012ext_ch_branch = tree->GetBranch(tree->GetAlias("els_iso03_pf2012ext_ch"));
		els_iso03_pf2012ext_ch_branch->SetAddress(&els_iso03_pf2012ext_ch_);
	}
	els_iso03_pf2012ext_em_branch = 0;
	if (tree->GetAlias("els_iso03_pf2012ext_em") != 0) {
		els_iso03_pf2012ext_em_branch = tree->GetBranch(tree->GetAlias("els_iso03_pf2012ext_em"));
		els_iso03_pf2012ext_em_branch->SetAddress(&els_iso03_pf2012ext_em_);
	}
	els_iso03_pf2012ext_nh_branch = 0;
	if (tree->GetAlias("els_iso03_pf2012ext_nh") != 0) {
		els_iso03_pf2012ext_nh_branch = tree->GetBranch(tree->GetAlias("els_iso03_pf2012ext_nh"));
		els_iso03_pf2012ext_nh_branch->SetAddress(&els_iso03_pf2012ext_nh_);
	}
	els_iso03_pf2012_nh_branch = 0;
	if (tree->GetAlias("els_iso03_pf2012_nh") != 0) {
		els_iso03_pf2012_nh_branch = tree->GetBranch(tree->GetAlias("els_iso03_pf2012_nh"));
		els_iso03_pf2012_nh_branch->SetAddress(&els_iso03_pf2012_nh_);
	}
	els_iso03_pf_ch_branch = 0;
	if (tree->GetAlias("els_iso03_pf_ch") != 0) {
		els_iso03_pf_ch_branch = tree->GetBranch(tree->GetAlias("els_iso03_pf_ch"));
		els_iso03_pf_ch_branch->SetAddress(&els_iso03_pf_ch_);
	}
	els_iso03_pf_gamma05_branch = 0;
	if (tree->GetAlias("els_iso03_pf_gamma05") != 0) {
		els_iso03_pf_gamma05_branch = tree->GetBranch(tree->GetAlias("els_iso03_pf_gamma05"));
		els_iso03_pf_gamma05_branch->SetAddress(&els_iso03_pf_gamma05_);
	}
	els_iso03_pf_nhad05_branch = 0;
	if (tree->GetAlias("els_iso03_pf_nhad05") != 0) {
		els_iso03_pf_nhad05_branch = tree->GetBranch(tree->GetAlias("els_iso03_pf_nhad05"));
		els_iso03_pf_nhad05_branch->SetAddress(&els_iso03_pf_nhad05_);
	}
	els_iso04_pf_branch = 0;
	if (tree->GetAlias("els_iso04_pf") != 0) {
		els_iso04_pf_branch = tree->GetBranch(tree->GetAlias("els_iso04_pf"));
		els_iso04_pf_branch->SetAddress(&els_iso04_pf_);
	}
	els_iso04_pf2012_ch_branch = 0;
	if (tree->GetAlias("els_iso04_pf2012_ch") != 0) {
		els_iso04_pf2012_ch_branch = tree->GetBranch(tree->GetAlias("els_iso04_pf2012_ch"));
		els_iso04_pf2012_ch_branch->SetAddress(&els_iso04_pf2012_ch_);
	}
	els_iso04_pf2012_em_branch = 0;
	if (tree->GetAlias("els_iso04_pf2012_em") != 0) {
		els_iso04_pf2012_em_branch = tree->GetBranch(tree->GetAlias("els_iso04_pf2012_em"));
		els_iso04_pf2012_em_branch->SetAddress(&els_iso04_pf2012_em_);
	}
	els_iso04_pf2012ext_ch_branch = 0;
	if (tree->GetAlias("els_iso04_pf2012ext_ch") != 0) {
		els_iso04_pf2012ext_ch_branch = tree->GetBranch(tree->GetAlias("els_iso04_pf2012ext_ch"));
		els_iso04_pf2012ext_ch_branch->SetAddress(&els_iso04_pf2012ext_ch_);
	}
	els_iso04_pf2012ext_em_branch = 0;
	if (tree->GetAlias("els_iso04_pf2012ext_em") != 0) {
		els_iso04_pf2012ext_em_branch = tree->GetBranch(tree->GetAlias("els_iso04_pf2012ext_em"));
		els_iso04_pf2012ext_em_branch->SetAddress(&els_iso04_pf2012ext_em_);
	}
	els_iso04_pf2012ext_nh_branch = 0;
	if (tree->GetAlias("els_iso04_pf2012ext_nh") != 0) {
		els_iso04_pf2012ext_nh_branch = tree->GetBranch(tree->GetAlias("els_iso04_pf2012ext_nh"));
		els_iso04_pf2012ext_nh_branch->SetAddress(&els_iso04_pf2012ext_nh_);
	}
	els_iso04_pf2012_nh_branch = 0;
	if (tree->GetAlias("els_iso04_pf2012_nh") != 0) {
		els_iso04_pf2012_nh_branch = tree->GetBranch(tree->GetAlias("els_iso04_pf2012_nh"));
		els_iso04_pf2012_nh_branch->SetAddress(&els_iso04_pf2012_nh_);
	}
	els_iso04_pf_ch_branch = 0;
	if (tree->GetAlias("els_iso04_pf_ch") != 0) {
		els_iso04_pf_ch_branch = tree->GetBranch(tree->GetAlias("els_iso04_pf_ch"));
		els_iso04_pf_ch_branch->SetAddress(&els_iso04_pf_ch_);
	}
	els_iso04_pf_gamma05_branch = 0;
	if (tree->GetAlias("els_iso04_pf_gamma05") != 0) {
		els_iso04_pf_gamma05_branch = tree->GetBranch(tree->GetAlias("els_iso04_pf_gamma05"));
		els_iso04_pf_gamma05_branch->SetAddress(&els_iso04_pf_gamma05_);
	}
	els_iso04_pf_nhad05_branch = 0;
	if (tree->GetAlias("els_iso04_pf_nhad05") != 0) {
		els_iso04_pf_nhad05_branch = tree->GetBranch(tree->GetAlias("els_iso04_pf_nhad05"));
		els_iso04_pf_nhad05_branch->SetAddress(&els_iso04_pf_nhad05_);
	}
	els_layer1_charge_branch = 0;
	if (tree->GetAlias("els_layer1_charge") != 0) {
		els_layer1_charge_branch = tree->GetBranch(tree->GetAlias("els_layer1_charge"));
		els_layer1_charge_branch->SetAddress(&els_layer1_charge_);
	}
	els_lh_branch = 0;
	if (tree->GetAlias("els_lh") != 0) {
		els_lh_branch = tree->GetBranch(tree->GetAlias("els_lh"));
		els_lh_branch->SetAddress(&els_lh_);
	}
	els_mc_patMatch_dr_branch = 0;
	if (tree->GetAlias("els_mc_patMatch_dr") != 0) {
		els_mc_patMatch_dr_branch = tree->GetBranch(tree->GetAlias("els_mc_patMatch_dr"));
		els_mc_patMatch_dr_branch->SetAddress(&els_mc_patMatch_dr_);
	}
	els_mva_branch = 0;
	if (tree->GetAlias("els_mva") != 0) {
		els_mva_branch = tree->GetBranch(tree->GetAlias("els_mva"));
		els_mva_branch->SetAddress(&els_mva_);
	}
	els_ndof_branch = 0;
	if (tree->GetAlias("els_ndof") != 0) {
		els_ndof_branch = tree->GetBranch(tree->GetAlias("els_ndof"));
		els_ndof_branch->SetAddress(&els_ndof_);
	}
	els_pfChargedHadronIso_branch = 0;
	if (tree->GetAlias("els_pfChargedHadronIso") != 0) {
		els_pfChargedHadronIso_branch = tree->GetBranch(tree->GetAlias("els_pfChargedHadronIso"));
		els_pfChargedHadronIso_branch->SetAddress(&els_pfChargedHadronIso_);
	}
	els_pfNeutralHadronIso_branch = 0;
	if (tree->GetAlias("els_pfNeutralHadronIso") != 0) {
		els_pfNeutralHadronIso_branch = tree->GetBranch(tree->GetAlias("els_pfNeutralHadronIso"));
		els_pfNeutralHadronIso_branch->SetAddress(&els_pfNeutralHadronIso_);
	}
	els_pfPUIso_branch = 0;
	if (tree->GetAlias("els_pfPUIso") != 0) {
		els_pfPUIso_branch = tree->GetBranch(tree->GetAlias("els_pfPUIso"));
		els_pfPUIso_branch->SetAddress(&els_pfPUIso_);
	}
	els_pfPhotonIso_branch = 0;
	if (tree->GetAlias("els_pfPhotonIso") != 0) {
		els_pfPhotonIso_branch = tree->GetBranch(tree->GetAlias("els_pfPhotonIso"));
		els_pfPhotonIso_branch->SetAddress(&els_pfPhotonIso_);
	}
	els_phiErr_branch = 0;
	if (tree->GetAlias("els_phiErr") != 0) {
		els_phiErr_branch = tree->GetBranch(tree->GetAlias("els_phiErr"));
		els_phiErr_branch->SetAddress(&els_phiErr_);
	}
	els_phiSC_branch = 0;
	if (tree->GetAlias("els_phiSC") != 0) {
		els_phiSC_branch = tree->GetBranch(tree->GetAlias("els_phiSC"));
		els_phiSC_branch->SetAddress(&els_phiSC_);
	}
	els_phiSCwidth_branch = 0;
	if (tree->GetAlias("els_phiSCwidth") != 0) {
		els_phiSCwidth_branch = tree->GetBranch(tree->GetAlias("els_phiSCwidth"));
		els_phiSCwidth_branch->SetAddress(&els_phiSCwidth_);
	}
	els_ptErr_branch = 0;
	if (tree->GetAlias("els_ptErr") != 0) {
		els_ptErr_branch = tree->GetBranch(tree->GetAlias("els_ptErr"));
		els_ptErr_branch->SetAddress(&els_ptErr_);
	}
	els_ptErrGsf_branch = 0;
	if (tree->GetAlias("els_ptErrGsf") != 0) {
		els_ptErrGsf_branch = tree->GetBranch(tree->GetAlias("els_ptErrGsf"));
		els_ptErrGsf_branch->SetAddress(&els_ptErrGsf_);
	}
	els_r9_branch = 0;
	if (tree->GetAlias("els_r9") != 0) {
		els_r9_branch = tree->GetBranch(tree->GetAlias("els_r9"));
		els_r9_branch->SetAddress(&els_r9_);
	}
	els_r9_full5x5_branch = 0;
	if (tree->GetAlias("els_r9_full5x5") != 0) {
		els_r9_full5x5_branch = tree->GetBranch(tree->GetAlias("els_r9_full5x5"));
		els_r9_full5x5_branch->SetAddress(&els_r9_full5x5_);
	}
	els_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("els_sigmaEtaEta") != 0) {
		els_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("els_sigmaEtaEta"));
		els_sigmaEtaEta_branch->SetAddress(&els_sigmaEtaEta_);
	}
	els_sigmaEtaEta_full5x5_branch = 0;
	if (tree->GetAlias("els_sigmaEtaEta_full5x5") != 0) {
		els_sigmaEtaEta_full5x5_branch = tree->GetBranch(tree->GetAlias("els_sigmaEtaEta_full5x5"));
		els_sigmaEtaEta_full5x5_branch->SetAddress(&els_sigmaEtaEta_full5x5_);
	}
	els_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("els_sigmaIEtaIEta") != 0) {
		els_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("els_sigmaIEtaIEta"));
		els_sigmaIEtaIEta_branch->SetAddress(&els_sigmaIEtaIEta_);
	}
	els_sigmaIEtaIEtaSC_branch = 0;
	if (tree->GetAlias("els_sigmaIEtaIEtaSC") != 0) {
		els_sigmaIEtaIEtaSC_branch = tree->GetBranch(tree->GetAlias("els_sigmaIEtaIEtaSC"));
		els_sigmaIEtaIEtaSC_branch->SetAddress(&els_sigmaIEtaIEtaSC_);
	}
	els_sigmaIEtaIEta_full5x5_branch = 0;
	if (tree->GetAlias("els_sigmaIEtaIEta_full5x5") != 0) {
		els_sigmaIEtaIEta_full5x5_branch = tree->GetBranch(tree->GetAlias("els_sigmaIEtaIEta_full5x5"));
		els_sigmaIEtaIEta_full5x5_branch->SetAddress(&els_sigmaIEtaIEta_full5x5_);
	}
	els_sigmaIEtaIPhi_branch = 0;
	if (tree->GetAlias("els_sigmaIEtaIPhi") != 0) {
		els_sigmaIEtaIPhi_branch = tree->GetBranch(tree->GetAlias("els_sigmaIEtaIPhi"));
		els_sigmaIEtaIPhi_branch->SetAddress(&els_sigmaIEtaIPhi_);
	}
	els_sigmaIPhiIPhi_branch = 0;
	if (tree->GetAlias("els_sigmaIPhiIPhi") != 0) {
		els_sigmaIPhiIPhi_branch = tree->GetBranch(tree->GetAlias("els_sigmaIPhiIPhi"));
		els_sigmaIPhiIPhi_branch->SetAddress(&els_sigmaIPhiIPhi_);
	}
	els_sigmaIPhiIPhiSC_branch = 0;
	if (tree->GetAlias("els_sigmaIPhiIPhiSC") != 0) {
		els_sigmaIPhiIPhiSC_branch = tree->GetBranch(tree->GetAlias("els_sigmaIPhiIPhiSC"));
		els_sigmaIPhiIPhiSC_branch->SetAddress(&els_sigmaIPhiIPhiSC_);
	}
	els_sigmaIPhiIPhi_full5x5_branch = 0;
	if (tree->GetAlias("els_sigmaIPhiIPhi_full5x5") != 0) {
		els_sigmaIPhiIPhi_full5x5_branch = tree->GetBranch(tree->GetAlias("els_sigmaIPhiIPhi_full5x5"));
		els_sigmaIPhiIPhi_full5x5_branch->SetAddress(&els_sigmaIPhiIPhi_full5x5_);
	}
	els_sigmaIphiIphi_branch = 0;
	if (tree->GetAlias("els_sigmaIphiIphi") != 0) {
		els_sigmaIphiIphi_branch = tree->GetBranch(tree->GetAlias("els_sigmaIphiIphi"));
		els_sigmaIphiIphi_branch->SetAddress(&els_sigmaIphiIphi_);
	}
	els_sigmaPhiPhi_branch = 0;
	if (tree->GetAlias("els_sigmaPhiPhi") != 0) {
		els_sigmaPhiPhi_branch = tree->GetBranch(tree->GetAlias("els_sigmaPhiPhi"));
		els_sigmaPhiPhi_branch->SetAddress(&els_sigmaPhiPhi_);
	}
	els_tkIso_branch = 0;
	if (tree->GetAlias("els_tkIso") != 0) {
		els_tkIso_branch = tree->GetBranch(tree->GetAlias("els_tkIso"));
		els_tkIso_branch->SetAddress(&els_tkIso_);
	}
	els_tkIso04_branch = 0;
	if (tree->GetAlias("els_tkIso04") != 0) {
		els_tkIso04_branch = tree->GetBranch(tree->GetAlias("els_tkIso04"));
		els_tkIso04_branch->SetAddress(&els_tkIso04_);
	}
	els_trackMomentumError_branch = 0;
	if (tree->GetAlias("els_trackMomentumError") != 0) {
		els_trackMomentumError_branch = tree->GetBranch(tree->GetAlias("els_trackMomentumError"));
		els_trackMomentumError_branch->SetAddress(&els_trackMomentumError_);
	}
	els_trkdr_branch = 0;
	if (tree->GetAlias("els_trkdr") != 0) {
		els_trkdr_branch = tree->GetBranch(tree->GetAlias("els_trkdr"));
		els_trkdr_branch->SetAddress(&els_trkdr_);
	}
	els_trkshFrac_branch = 0;
	if (tree->GetAlias("els_trkshFrac") != 0) {
		els_trkshFrac_branch = tree->GetBranch(tree->GetAlias("els_trkshFrac"));
		els_trkshFrac_branch->SetAddress(&els_trkshFrac_);
	}
	els_z0_branch = 0;
	if (tree->GetAlias("els_z0") != 0) {
		els_z0_branch = tree->GetBranch(tree->GetAlias("els_z0"));
		els_z0_branch->SetAddress(&els_z0_);
	}
	els_z0Err_branch = 0;
	if (tree->GetAlias("els_z0Err") != 0) {
		els_z0Err_branch = tree->GetBranch(tree->GetAlias("els_z0Err"));
		els_z0Err_branch->SetAddress(&els_z0Err_);
	}
	els_z0corr_branch = 0;
	if (tree->GetAlias("els_z0corr") != 0) {
		els_z0corr_branch = tree->GetBranch(tree->GetAlias("els_z0corr"));
		els_z0corr_branch->SetAddress(&els_z0corr_);
	}
	mus_backToBackCompat_branch = 0;
	if (tree->GetAlias("mus_backToBackCompat") != 0) {
		mus_backToBackCompat_branch = tree->GetBranch(tree->GetAlias("mus_backToBackCompat"));
		mus_backToBackCompat_branch->SetAddress(&mus_backToBackCompat_);
	}
	mus_best_dxyPV_branch = 0;
	if (tree->GetAlias("mus_best_dxyPV") != 0) {
		mus_best_dxyPV_branch = tree->GetBranch(tree->GetAlias("mus_best_dxyPV"));
		mus_best_dxyPV_branch->SetAddress(&mus_best_dxyPV_);
	}
	mus_best_dzPV_branch = 0;
	if (tree->GetAlias("mus_best_dzPV") != 0) {
		mus_best_dzPV_branch = tree->GetBranch(tree->GetAlias("mus_best_dzPV"));
		mus_best_dzPV_branch->SetAddress(&mus_best_dzPV_);
	}
	mus_caloCompatibility_branch = 0;
	if (tree->GetAlias("mus_caloCompatibility") != 0) {
		mus_caloCompatibility_branch = tree->GetBranch(tree->GetAlias("mus_caloCompatibility"));
		mus_caloCompatibility_branch->SetAddress(&mus_caloCompatibility_);
	}
	mus_chi2_branch = 0;
	if (tree->GetAlias("mus_chi2") != 0) {
		mus_chi2_branch = tree->GetBranch(tree->GetAlias("mus_chi2"));
		mus_chi2_branch->SetAddress(&mus_chi2_);
	}
	mus_chi2LocalMomentum_branch = 0;
	if (tree->GetAlias("mus_chi2LocalMomentum") != 0) {
		mus_chi2LocalMomentum_branch = tree->GetBranch(tree->GetAlias("mus_chi2LocalMomentum"));
		mus_chi2LocalMomentum_branch->SetAddress(&mus_chi2LocalMomentum_);
	}
	mus_chi2LocalPosition_branch = 0;
	if (tree->GetAlias("mus_chi2LocalPosition") != 0) {
		mus_chi2LocalPosition_branch = tree->GetBranch(tree->GetAlias("mus_chi2LocalPosition"));
		mus_chi2LocalPosition_branch->SetAddress(&mus_chi2LocalPosition_);
	}
	mus_cosmicCompat_branch = 0;
	if (tree->GetAlias("mus_cosmicCompat") != 0) {
		mus_cosmicCompat_branch = tree->GetBranch(tree->GetAlias("mus_cosmicCompat"));
		mus_cosmicCompat_branch->SetAddress(&mus_cosmicCompat_);
	}
	mus_d0_branch = 0;
	if (tree->GetAlias("mus_d0") != 0) {
		mus_d0_branch = tree->GetBranch(tree->GetAlias("mus_d0"));
		mus_d0_branch->SetAddress(&mus_d0_);
	}
	mus_d0Err_branch = 0;
	if (tree->GetAlias("mus_d0Err") != 0) {
		mus_d0Err_branch = tree->GetBranch(tree->GetAlias("mus_d0Err"));
		mus_d0Err_branch->SetAddress(&mus_d0Err_);
	}
	mus_d0corr_branch = 0;
	if (tree->GetAlias("mus_d0corr") != 0) {
		mus_d0corr_branch = tree->GetBranch(tree->GetAlias("mus_d0corr"));
		mus_d0corr_branch->SetAddress(&mus_d0corr_);
	}
	mus_d0corrPhi_branch = 0;
	if (tree->GetAlias("mus_d0corrPhi") != 0) {
		mus_d0corrPhi_branch = tree->GetBranch(tree->GetAlias("mus_d0corrPhi"));
		mus_d0corrPhi_branch->SetAddress(&mus_d0corrPhi_);
	}
	mus_d0phiCov_branch = 0;
	if (tree->GetAlias("mus_d0phiCov") != 0) {
		mus_d0phiCov_branch = tree->GetBranch(tree->GetAlias("mus_d0phiCov"));
		mus_d0phiCov_branch->SetAddress(&mus_d0phiCov_);
	}
	mus_dxyPV_branch = 0;
	if (tree->GetAlias("mus_dxyPV") != 0) {
		mus_dxyPV_branch = tree->GetBranch(tree->GetAlias("mus_dxyPV"));
		mus_dxyPV_branch->SetAddress(&mus_dxyPV_);
	}
	mus_dzPV_branch = 0;
	if (tree->GetAlias("mus_dzPV") != 0) {
		mus_dzPV_branch = tree->GetBranch(tree->GetAlias("mus_dzPV"));
		mus_dzPV_branch->SetAddress(&mus_dzPV_);
	}
	mus_ecal_time_branch = 0;
	if (tree->GetAlias("mus_ecal_time") != 0) {
		mus_ecal_time_branch = tree->GetBranch(tree->GetAlias("mus_ecal_time"));
		mus_ecal_time_branch->SetAddress(&mus_ecal_time_);
	}
	mus_ecal_timeError_branch = 0;
	if (tree->GetAlias("mus_ecal_timeError") != 0) {
		mus_ecal_timeError_branch = tree->GetBranch(tree->GetAlias("mus_ecal_timeError"));
		mus_ecal_timeError_branch->SetAddress(&mus_ecal_timeError_);
	}
	mus_e_em_branch = 0;
	if (tree->GetAlias("mus_e_em") != 0) {
		mus_e_em_branch = tree->GetBranch(tree->GetAlias("mus_e_em"));
		mus_e_em_branch->SetAddress(&mus_e_em_);
	}
	mus_e_emS9_branch = 0;
	if (tree->GetAlias("mus_e_emS9") != 0) {
		mus_e_emS9_branch = tree->GetBranch(tree->GetAlias("mus_e_emS9"));
		mus_e_emS9_branch->SetAddress(&mus_e_emS9_);
	}
	mus_e_had_branch = 0;
	if (tree->GetAlias("mus_e_had") != 0) {
		mus_e_had_branch = tree->GetBranch(tree->GetAlias("mus_e_had"));
		mus_e_had_branch->SetAddress(&mus_e_had_);
	}
	mus_e_hadS9_branch = 0;
	if (tree->GetAlias("mus_e_hadS9") != 0) {
		mus_e_hadS9_branch = tree->GetBranch(tree->GetAlias("mus_e_hadS9"));
		mus_e_hadS9_branch->SetAddress(&mus_e_hadS9_);
	}
	mus_e_ho_branch = 0;
	if (tree->GetAlias("mus_e_ho") != 0) {
		mus_e_ho_branch = tree->GetBranch(tree->GetAlias("mus_e_ho"));
		mus_e_ho_branch->SetAddress(&mus_e_ho_);
	}
	mus_e_hoS9_branch = 0;
	if (tree->GetAlias("mus_e_hoS9") != 0) {
		mus_e_hoS9_branch = tree->GetBranch(tree->GetAlias("mus_e_hoS9"));
		mus_e_hoS9_branch->SetAddress(&mus_e_hoS9_);
	}
	mus_emMax_branch = 0;
	if (tree->GetAlias("mus_emMax") != 0) {
		mus_emMax_branch = tree->GetBranch(tree->GetAlias("mus_emMax"));
		mus_emMax_branch->SetAddress(&mus_emMax_);
	}
	mus_emS25_branch = 0;
	if (tree->GetAlias("mus_emS25") != 0) {
		mus_emS25_branch = tree->GetBranch(tree->GetAlias("mus_emS25"));
		mus_emS25_branch->SetAddress(&mus_emS25_);
	}
	mus_etaErr_branch = 0;
	if (tree->GetAlias("mus_etaErr") != 0) {
		mus_etaErr_branch = tree->GetBranch(tree->GetAlias("mus_etaErr"));
		mus_etaErr_branch->SetAddress(&mus_etaErr_);
	}
	mus_gfit_chi2_branch = 0;
	if (tree->GetAlias("mus_gfit_chi2") != 0) {
		mus_gfit_chi2_branch = tree->GetBranch(tree->GetAlias("mus_gfit_chi2"));
		mus_gfit_chi2_branch->SetAddress(&mus_gfit_chi2_);
	}
	mus_gfit_d0_branch = 0;
	if (tree->GetAlias("mus_gfit_d0") != 0) {
		mus_gfit_d0_branch = tree->GetBranch(tree->GetAlias("mus_gfit_d0"));
		mus_gfit_d0_branch->SetAddress(&mus_gfit_d0_);
	}
	mus_gfit_d0Err_branch = 0;
	if (tree->GetAlias("mus_gfit_d0Err") != 0) {
		mus_gfit_d0Err_branch = tree->GetBranch(tree->GetAlias("mus_gfit_d0Err"));
		mus_gfit_d0Err_branch->SetAddress(&mus_gfit_d0Err_);
	}
	mus_gfit_d0corr_branch = 0;
	if (tree->GetAlias("mus_gfit_d0corr") != 0) {
		mus_gfit_d0corr_branch = tree->GetBranch(tree->GetAlias("mus_gfit_d0corr"));
		mus_gfit_d0corr_branch->SetAddress(&mus_gfit_d0corr_);
	}
	mus_gfit_d0corrPhi_branch = 0;
	if (tree->GetAlias("mus_gfit_d0corrPhi") != 0) {
		mus_gfit_d0corrPhi_branch = tree->GetBranch(tree->GetAlias("mus_gfit_d0corrPhi"));
		mus_gfit_d0corrPhi_branch->SetAddress(&mus_gfit_d0corrPhi_);
	}
	mus_gfit_d0phiCov_branch = 0;
	if (tree->GetAlias("mus_gfit_d0phiCov") != 0) {
		mus_gfit_d0phiCov_branch = tree->GetBranch(tree->GetAlias("mus_gfit_d0phiCov"));
		mus_gfit_d0phiCov_branch->SetAddress(&mus_gfit_d0phiCov_);
	}
	mus_gfit_ndof_branch = 0;
	if (tree->GetAlias("mus_gfit_ndof") != 0) {
		mus_gfit_ndof_branch = tree->GetBranch(tree->GetAlias("mus_gfit_ndof"));
		mus_gfit_ndof_branch->SetAddress(&mus_gfit_ndof_);
	}
	mus_gfit_qoverp_branch = 0;
	if (tree->GetAlias("mus_gfit_qoverp") != 0) {
		mus_gfit_qoverp_branch = tree->GetBranch(tree->GetAlias("mus_gfit_qoverp"));
		mus_gfit_qoverp_branch->SetAddress(&mus_gfit_qoverp_);
	}
	mus_gfit_qoverpError_branch = 0;
	if (tree->GetAlias("mus_gfit_qoverpError") != 0) {
		mus_gfit_qoverpError_branch = tree->GetBranch(tree->GetAlias("mus_gfit_qoverpError"));
		mus_gfit_qoverpError_branch->SetAddress(&mus_gfit_qoverpError_);
	}
	mus_gfit_z0_branch = 0;
	if (tree->GetAlias("mus_gfit_z0") != 0) {
		mus_gfit_z0_branch = tree->GetBranch(tree->GetAlias("mus_gfit_z0"));
		mus_gfit_z0_branch->SetAddress(&mus_gfit_z0_);
	}
	mus_gfit_z0Err_branch = 0;
	if (tree->GetAlias("mus_gfit_z0Err") != 0) {
		mus_gfit_z0Err_branch = tree->GetBranch(tree->GetAlias("mus_gfit_z0Err"));
		mus_gfit_z0Err_branch->SetAddress(&mus_gfit_z0Err_);
	}
	mus_gfit_z0corr_branch = 0;
	if (tree->GetAlias("mus_gfit_z0corr") != 0) {
		mus_gfit_z0corr_branch = tree->GetBranch(tree->GetAlias("mus_gfit_z0corr"));
		mus_gfit_z0corr_branch->SetAddress(&mus_gfit_z0corr_);
	}
	mus_glbKink_branch = 0;
	if (tree->GetAlias("mus_glbKink") != 0) {
		mus_glbKink_branch = tree->GetBranch(tree->GetAlias("mus_glbKink"));
		mus_glbKink_branch->SetAddress(&mus_glbKink_);
	}
	mus_glbTrackProbability_branch = 0;
	if (tree->GetAlias("mus_glbTrackProbability") != 0) {
		mus_glbTrackProbability_branch = tree->GetBranch(tree->GetAlias("mus_glbTrackProbability"));
		mus_glbTrackProbability_branch->SetAddress(&mus_glbTrackProbability_);
	}
	mus_globalDeltaEtaPhi_branch = 0;
	if (tree->GetAlias("mus_globalDeltaEtaPhi") != 0) {
		mus_globalDeltaEtaPhi_branch = tree->GetBranch(tree->GetAlias("mus_globalDeltaEtaPhi"));
		mus_globalDeltaEtaPhi_branch->SetAddress(&mus_globalDeltaEtaPhi_);
	}
	mus_hadMax_branch = 0;
	if (tree->GetAlias("mus_hadMax") != 0) {
		mus_hadMax_branch = tree->GetBranch(tree->GetAlias("mus_hadMax"));
		mus_hadMax_branch->SetAddress(&mus_hadMax_);
	}
	mus_hcal_time_branch = 0;
	if (tree->GetAlias("mus_hcal_time") != 0) {
		mus_hcal_time_branch = tree->GetBranch(tree->GetAlias("mus_hcal_time"));
		mus_hcal_time_branch->SetAddress(&mus_hcal_time_);
	}
	mus_hcal_timeError_branch = 0;
	if (tree->GetAlias("mus_hcal_timeError") != 0) {
		mus_hcal_timeError_branch = tree->GetBranch(tree->GetAlias("mus_hcal_timeError"));
		mus_hcal_timeError_branch->SetAddress(&mus_hcal_timeError_);
	}
	mus_ip3d_branch = 0;
	if (tree->GetAlias("mus_ip3d") != 0) {
		mus_ip3d_branch = tree->GetBranch(tree->GetAlias("mus_ip3d"));
		mus_ip3d_branch->SetAddress(&mus_ip3d_);
	}
	mus_ip3derr_branch = 0;
	if (tree->GetAlias("mus_ip3derr") != 0) {
		mus_ip3derr_branch = tree->GetBranch(tree->GetAlias("mus_ip3derr"));
		mus_ip3derr_branch->SetAddress(&mus_ip3derr_);
	}
	mus_iso03_emEt_branch = 0;
	if (tree->GetAlias("mus_iso03_emEt") != 0) {
		mus_iso03_emEt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_emEt"));
		mus_iso03_emEt_branch->SetAddress(&mus_iso03_emEt_);
	}
	mus_iso03_hadEt_branch = 0;
	if (tree->GetAlias("mus_iso03_hadEt") != 0) {
		mus_iso03_hadEt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_hadEt"));
		mus_iso03_hadEt_branch->SetAddress(&mus_iso03_hadEt_);
	}
	mus_iso03_hoEt_branch = 0;
	if (tree->GetAlias("mus_iso03_hoEt") != 0) {
		mus_iso03_hoEt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_hoEt"));
		mus_iso03_hoEt_branch->SetAddress(&mus_iso03_hoEt_);
	}
	mus_iso03_sumPt_branch = 0;
	if (tree->GetAlias("mus_iso03_sumPt") != 0) {
		mus_iso03_sumPt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_sumPt"));
		mus_iso03_sumPt_branch->SetAddress(&mus_iso03_sumPt_);
	}
	mus_iso05_emEt_branch = 0;
	if (tree->GetAlias("mus_iso05_emEt") != 0) {
		mus_iso05_emEt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_emEt"));
		mus_iso05_emEt_branch->SetAddress(&mus_iso05_emEt_);
	}
	mus_iso05_hadEt_branch = 0;
	if (tree->GetAlias("mus_iso05_hadEt") != 0) {
		mus_iso05_hadEt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_hadEt"));
		mus_iso05_hadEt_branch->SetAddress(&mus_iso05_hadEt_);
	}
	mus_iso05_hoEt_branch = 0;
	if (tree->GetAlias("mus_iso05_hoEt") != 0) {
		mus_iso05_hoEt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_hoEt"));
		mus_iso05_hoEt_branch->SetAddress(&mus_iso05_hoEt_);
	}
	mus_iso05_sumPt_branch = 0;
	if (tree->GetAlias("mus_iso05_sumPt") != 0) {
		mus_iso05_sumPt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_sumPt"));
		mus_iso05_sumPt_branch->SetAddress(&mus_iso05_sumPt_);
	}
	mus_isoMeanDRR03_pf_ChargedHadronPt_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR03_pf_ChargedHadronPt") != 0) {
		mus_isoMeanDRR03_pf_ChargedHadronPt_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR03_pf_ChargedHadronPt"));
		mus_isoMeanDRR03_pf_ChargedHadronPt_branch->SetAddress(&mus_isoMeanDRR03_pf_ChargedHadronPt_);
	}
	mus_isoMeanDRR03_pf_ChargedParticlePt_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR03_pf_ChargedParticlePt") != 0) {
		mus_isoMeanDRR03_pf_ChargedParticlePt_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR03_pf_ChargedParticlePt"));
		mus_isoMeanDRR03_pf_ChargedParticlePt_branch->SetAddress(&mus_isoMeanDRR03_pf_ChargedParticlePt_);
	}
	mus_isoMeanDRR03_pf_NeutralHadronEt_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR03_pf_NeutralHadronEt") != 0) {
		mus_isoMeanDRR03_pf_NeutralHadronEt_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR03_pf_NeutralHadronEt"));
		mus_isoMeanDRR03_pf_NeutralHadronEt_branch->SetAddress(&mus_isoMeanDRR03_pf_NeutralHadronEt_);
	}
	mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold") != 0) {
		mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold"));
		mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_branch->SetAddress(&mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_);
	}
	mus_isoMeanDRR03_pf_PUPt_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR03_pf_PUPt") != 0) {
		mus_isoMeanDRR03_pf_PUPt_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR03_pf_PUPt"));
		mus_isoMeanDRR03_pf_PUPt_branch->SetAddress(&mus_isoMeanDRR03_pf_PUPt_);
	}
	mus_isoMeanDRR03_pf_PhotonEt_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR03_pf_PhotonEt") != 0) {
		mus_isoMeanDRR03_pf_PhotonEt_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR03_pf_PhotonEt"));
		mus_isoMeanDRR03_pf_PhotonEt_branch->SetAddress(&mus_isoMeanDRR03_pf_PhotonEt_);
	}
	mus_isoMeanDRR03_pf_PhotonEtHighThreshold_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR03_pf_PhotonEtHighThreshold") != 0) {
		mus_isoMeanDRR03_pf_PhotonEtHighThreshold_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR03_pf_PhotonEtHighThreshold"));
		mus_isoMeanDRR03_pf_PhotonEtHighThreshold_branch->SetAddress(&mus_isoMeanDRR03_pf_PhotonEtHighThreshold_);
	}
	mus_isoMeanDRR04_pf_ChargedHadronPt_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR04_pf_ChargedHadronPt") != 0) {
		mus_isoMeanDRR04_pf_ChargedHadronPt_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR04_pf_ChargedHadronPt"));
		mus_isoMeanDRR04_pf_ChargedHadronPt_branch->SetAddress(&mus_isoMeanDRR04_pf_ChargedHadronPt_);
	}
	mus_isoMeanDRR04_pf_ChargedParticlePt_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR04_pf_ChargedParticlePt") != 0) {
		mus_isoMeanDRR04_pf_ChargedParticlePt_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR04_pf_ChargedParticlePt"));
		mus_isoMeanDRR04_pf_ChargedParticlePt_branch->SetAddress(&mus_isoMeanDRR04_pf_ChargedParticlePt_);
	}
	mus_isoMeanDRR04_pf_NeutralHadronEt_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR04_pf_NeutralHadronEt") != 0) {
		mus_isoMeanDRR04_pf_NeutralHadronEt_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR04_pf_NeutralHadronEt"));
		mus_isoMeanDRR04_pf_NeutralHadronEt_branch->SetAddress(&mus_isoMeanDRR04_pf_NeutralHadronEt_);
	}
	mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold") != 0) {
		mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold"));
		mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_branch->SetAddress(&mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_);
	}
	mus_isoMeanDRR04_pf_PUPt_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR04_pf_PUPt") != 0) {
		mus_isoMeanDRR04_pf_PUPt_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR04_pf_PUPt"));
		mus_isoMeanDRR04_pf_PUPt_branch->SetAddress(&mus_isoMeanDRR04_pf_PUPt_);
	}
	mus_isoMeanDRR04_pf_PhotonEt_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR04_pf_PhotonEt") != 0) {
		mus_isoMeanDRR04_pf_PhotonEt_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR04_pf_PhotonEt"));
		mus_isoMeanDRR04_pf_PhotonEt_branch->SetAddress(&mus_isoMeanDRR04_pf_PhotonEt_);
	}
	mus_isoMeanDRR04_pf_PhotonEtHighThreshold_branch = 0;
	if (tree->GetAlias("mus_isoMeanDRR04_pf_PhotonEtHighThreshold") != 0) {
		mus_isoMeanDRR04_pf_PhotonEtHighThreshold_branch = tree->GetBranch(tree->GetAlias("mus_isoMeanDRR04_pf_PhotonEtHighThreshold"));
		mus_isoMeanDRR04_pf_PhotonEtHighThreshold_branch->SetAddress(&mus_isoMeanDRR04_pf_PhotonEtHighThreshold_);
	}
	mus_isoR03_pf_ChargedHadronPt_branch = 0;
	if (tree->GetAlias("mus_isoR03_pf_ChargedHadronPt") != 0) {
		mus_isoR03_pf_ChargedHadronPt_branch = tree->GetBranch(tree->GetAlias("mus_isoR03_pf_ChargedHadronPt"));
		mus_isoR03_pf_ChargedHadronPt_branch->SetAddress(&mus_isoR03_pf_ChargedHadronPt_);
	}
	mus_isoR03_pf_ChargedParticlePt_branch = 0;
	if (tree->GetAlias("mus_isoR03_pf_ChargedParticlePt") != 0) {
		mus_isoR03_pf_ChargedParticlePt_branch = tree->GetBranch(tree->GetAlias("mus_isoR03_pf_ChargedParticlePt"));
		mus_isoR03_pf_ChargedParticlePt_branch->SetAddress(&mus_isoR03_pf_ChargedParticlePt_);
	}
	mus_isoR03_pf_NeutralHadronEt_branch = 0;
	if (tree->GetAlias("mus_isoR03_pf_NeutralHadronEt") != 0) {
		mus_isoR03_pf_NeutralHadronEt_branch = tree->GetBranch(tree->GetAlias("mus_isoR03_pf_NeutralHadronEt"));
		mus_isoR03_pf_NeutralHadronEt_branch->SetAddress(&mus_isoR03_pf_NeutralHadronEt_);
	}
	mus_isoR03_pf_NeutralHadronEtHighThreshold_branch = 0;
	if (tree->GetAlias("mus_isoR03_pf_NeutralHadronEtHighThreshold") != 0) {
		mus_isoR03_pf_NeutralHadronEtHighThreshold_branch = tree->GetBranch(tree->GetAlias("mus_isoR03_pf_NeutralHadronEtHighThreshold"));
		mus_isoR03_pf_NeutralHadronEtHighThreshold_branch->SetAddress(&mus_isoR03_pf_NeutralHadronEtHighThreshold_);
	}
	mus_isoR03_pf_PUPt_branch = 0;
	if (tree->GetAlias("mus_isoR03_pf_PUPt") != 0) {
		mus_isoR03_pf_PUPt_branch = tree->GetBranch(tree->GetAlias("mus_isoR03_pf_PUPt"));
		mus_isoR03_pf_PUPt_branch->SetAddress(&mus_isoR03_pf_PUPt_);
	}
	mus_isoR03_pf_PhotonEt_branch = 0;
	if (tree->GetAlias("mus_isoR03_pf_PhotonEt") != 0) {
		mus_isoR03_pf_PhotonEt_branch = tree->GetBranch(tree->GetAlias("mus_isoR03_pf_PhotonEt"));
		mus_isoR03_pf_PhotonEt_branch->SetAddress(&mus_isoR03_pf_PhotonEt_);
	}
	mus_isoR03_pf_PhotonEtHighThreshold_branch = 0;
	if (tree->GetAlias("mus_isoR03_pf_PhotonEtHighThreshold") != 0) {
		mus_isoR03_pf_PhotonEtHighThreshold_branch = tree->GetBranch(tree->GetAlias("mus_isoR03_pf_PhotonEtHighThreshold"));
		mus_isoR03_pf_PhotonEtHighThreshold_branch->SetAddress(&mus_isoR03_pf_PhotonEtHighThreshold_);
	}
	mus_isoR04_pf_ChargedHadronPt_branch = 0;
	if (tree->GetAlias("mus_isoR04_pf_ChargedHadronPt") != 0) {
		mus_isoR04_pf_ChargedHadronPt_branch = tree->GetBranch(tree->GetAlias("mus_isoR04_pf_ChargedHadronPt"));
		mus_isoR04_pf_ChargedHadronPt_branch->SetAddress(&mus_isoR04_pf_ChargedHadronPt_);
	}
	mus_isoR04_pf_ChargedParticlePt_branch = 0;
	if (tree->GetAlias("mus_isoR04_pf_ChargedParticlePt") != 0) {
		mus_isoR04_pf_ChargedParticlePt_branch = tree->GetBranch(tree->GetAlias("mus_isoR04_pf_ChargedParticlePt"));
		mus_isoR04_pf_ChargedParticlePt_branch->SetAddress(&mus_isoR04_pf_ChargedParticlePt_);
	}
	mus_isoR04_pf_NeutralHadronEt_branch = 0;
	if (tree->GetAlias("mus_isoR04_pf_NeutralHadronEt") != 0) {
		mus_isoR04_pf_NeutralHadronEt_branch = tree->GetBranch(tree->GetAlias("mus_isoR04_pf_NeutralHadronEt"));
		mus_isoR04_pf_NeutralHadronEt_branch->SetAddress(&mus_isoR04_pf_NeutralHadronEt_);
	}
	mus_isoR04_pf_NeutralHadronEtHighThreshold_branch = 0;
	if (tree->GetAlias("mus_isoR04_pf_NeutralHadronEtHighThreshold") != 0) {
		mus_isoR04_pf_NeutralHadronEtHighThreshold_branch = tree->GetBranch(tree->GetAlias("mus_isoR04_pf_NeutralHadronEtHighThreshold"));
		mus_isoR04_pf_NeutralHadronEtHighThreshold_branch->SetAddress(&mus_isoR04_pf_NeutralHadronEtHighThreshold_);
	}
	mus_isoR04_pf_PUPt_branch = 0;
	if (tree->GetAlias("mus_isoR04_pf_PUPt") != 0) {
		mus_isoR04_pf_PUPt_branch = tree->GetBranch(tree->GetAlias("mus_isoR04_pf_PUPt"));
		mus_isoR04_pf_PUPt_branch->SetAddress(&mus_isoR04_pf_PUPt_);
	}
	mus_isoR04_pf_PhotonEt_branch = 0;
	if (tree->GetAlias("mus_isoR04_pf_PhotonEt") != 0) {
		mus_isoR04_pf_PhotonEt_branch = tree->GetBranch(tree->GetAlias("mus_isoR04_pf_PhotonEt"));
		mus_isoR04_pf_PhotonEt_branch->SetAddress(&mus_isoR04_pf_PhotonEt_);
	}
	mus_isoR04_pf_PhotonEtHighThreshold_branch = 0;
	if (tree->GetAlias("mus_isoR04_pf_PhotonEtHighThreshold") != 0) {
		mus_isoR04_pf_PhotonEtHighThreshold_branch = tree->GetBranch(tree->GetAlias("mus_isoR04_pf_PhotonEtHighThreshold"));
		mus_isoR04_pf_PhotonEtHighThreshold_branch->SetAddress(&mus_isoR04_pf_PhotonEtHighThreshold_);
	}
	mus_isoSumDRR03_pf_ChargedHadronPt_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR03_pf_ChargedHadronPt") != 0) {
		mus_isoSumDRR03_pf_ChargedHadronPt_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR03_pf_ChargedHadronPt"));
		mus_isoSumDRR03_pf_ChargedHadronPt_branch->SetAddress(&mus_isoSumDRR03_pf_ChargedHadronPt_);
	}
	mus_isoSumDRR03_pf_ChargedParticlePt_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR03_pf_ChargedParticlePt") != 0) {
		mus_isoSumDRR03_pf_ChargedParticlePt_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR03_pf_ChargedParticlePt"));
		mus_isoSumDRR03_pf_ChargedParticlePt_branch->SetAddress(&mus_isoSumDRR03_pf_ChargedParticlePt_);
	}
	mus_isoSumDRR03_pf_NeutralHadronEt_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR03_pf_NeutralHadronEt") != 0) {
		mus_isoSumDRR03_pf_NeutralHadronEt_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR03_pf_NeutralHadronEt"));
		mus_isoSumDRR03_pf_NeutralHadronEt_branch->SetAddress(&mus_isoSumDRR03_pf_NeutralHadronEt_);
	}
	mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold") != 0) {
		mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold"));
		mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_branch->SetAddress(&mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_);
	}
	mus_isoSumDRR03_pf_PUPt_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR03_pf_PUPt") != 0) {
		mus_isoSumDRR03_pf_PUPt_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR03_pf_PUPt"));
		mus_isoSumDRR03_pf_PUPt_branch->SetAddress(&mus_isoSumDRR03_pf_PUPt_);
	}
	mus_isoSumDRR03_pf_PhotonEt_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR03_pf_PhotonEt") != 0) {
		mus_isoSumDRR03_pf_PhotonEt_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR03_pf_PhotonEt"));
		mus_isoSumDRR03_pf_PhotonEt_branch->SetAddress(&mus_isoSumDRR03_pf_PhotonEt_);
	}
	mus_isoSumDRR03_pf_PhotonEtHighThreshold_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR03_pf_PhotonEtHighThreshold") != 0) {
		mus_isoSumDRR03_pf_PhotonEtHighThreshold_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR03_pf_PhotonEtHighThreshold"));
		mus_isoSumDRR03_pf_PhotonEtHighThreshold_branch->SetAddress(&mus_isoSumDRR03_pf_PhotonEtHighThreshold_);
	}
	mus_isoSumDRR04_pf_ChargedHadronPt_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR04_pf_ChargedHadronPt") != 0) {
		mus_isoSumDRR04_pf_ChargedHadronPt_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR04_pf_ChargedHadronPt"));
		mus_isoSumDRR04_pf_ChargedHadronPt_branch->SetAddress(&mus_isoSumDRR04_pf_ChargedHadronPt_);
	}
	mus_isoSumDRR04_pf_ChargedParticlePt_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR04_pf_ChargedParticlePt") != 0) {
		mus_isoSumDRR04_pf_ChargedParticlePt_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR04_pf_ChargedParticlePt"));
		mus_isoSumDRR04_pf_ChargedParticlePt_branch->SetAddress(&mus_isoSumDRR04_pf_ChargedParticlePt_);
	}
	mus_isoSumDRR04_pf_NeutralHadronEt_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR04_pf_NeutralHadronEt") != 0) {
		mus_isoSumDRR04_pf_NeutralHadronEt_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR04_pf_NeutralHadronEt"));
		mus_isoSumDRR04_pf_NeutralHadronEt_branch->SetAddress(&mus_isoSumDRR04_pf_NeutralHadronEt_);
	}
	mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold") != 0) {
		mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold"));
		mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_branch->SetAddress(&mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_);
	}
	mus_isoSumDRR04_pf_PUPt_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR04_pf_PUPt") != 0) {
		mus_isoSumDRR04_pf_PUPt_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR04_pf_PUPt"));
		mus_isoSumDRR04_pf_PUPt_branch->SetAddress(&mus_isoSumDRR04_pf_PUPt_);
	}
	mus_isoSumDRR04_pf_PhotonEt_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR04_pf_PhotonEt") != 0) {
		mus_isoSumDRR04_pf_PhotonEt_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR04_pf_PhotonEt"));
		mus_isoSumDRR04_pf_PhotonEt_branch->SetAddress(&mus_isoSumDRR04_pf_PhotonEt_);
	}
	mus_isoSumDRR04_pf_PhotonEtHighThreshold_branch = 0;
	if (tree->GetAlias("mus_isoSumDRR04_pf_PhotonEtHighThreshold") != 0) {
		mus_isoSumDRR04_pf_PhotonEtHighThreshold_branch = tree->GetBranch(tree->GetAlias("mus_isoSumDRR04_pf_PhotonEtHighThreshold"));
		mus_isoSumDRR04_pf_PhotonEtHighThreshold_branch->SetAddress(&mus_isoSumDRR04_pf_PhotonEtHighThreshold_);
	}
	mus_iso_ecalvetoDep_branch = 0;
	if (tree->GetAlias("mus_iso_ecalvetoDep") != 0) {
		mus_iso_ecalvetoDep_branch = tree->GetBranch(tree->GetAlias("mus_iso_ecalvetoDep"));
		mus_iso_ecalvetoDep_branch->SetAddress(&mus_iso_ecalvetoDep_);
	}
	mus_iso_hcalvetoDep_branch = 0;
	if (tree->GetAlias("mus_iso_hcalvetoDep") != 0) {
		mus_iso_hcalvetoDep_branch = tree->GetBranch(tree->GetAlias("mus_iso_hcalvetoDep"));
		mus_iso_hcalvetoDep_branch->SetAddress(&mus_iso_hcalvetoDep_);
	}
	mus_iso_hovetoDep_branch = 0;
	if (tree->GetAlias("mus_iso_hovetoDep") != 0) {
		mus_iso_hovetoDep_branch = tree->GetBranch(tree->GetAlias("mus_iso_hovetoDep"));
		mus_iso_hovetoDep_branch->SetAddress(&mus_iso_hovetoDep_);
	}
	mus_iso_trckvetoDep_branch = 0;
	if (tree->GetAlias("mus_iso_trckvetoDep") != 0) {
		mus_iso_trckvetoDep_branch = tree->GetBranch(tree->GetAlias("mus_iso_trckvetoDep"));
		mus_iso_trckvetoDep_branch->SetAddress(&mus_iso_trckvetoDep_);
	}
	mus_localDistance_branch = 0;
	if (tree->GetAlias("mus_localDistance") != 0) {
		mus_localDistance_branch = tree->GetBranch(tree->GetAlias("mus_localDistance"));
		mus_localDistance_branch->SetAddress(&mus_localDistance_);
	}
	mus_mc_patMatch_dr_branch = 0;
	if (tree->GetAlias("mus_mc_patMatch_dr") != 0) {
		mus_mc_patMatch_dr_branch = tree->GetBranch(tree->GetAlias("mus_mc_patMatch_dr"));
		mus_mc_patMatch_dr_branch->SetAddress(&mus_mc_patMatch_dr_);
	}
	mus_ndof_branch = 0;
	if (tree->GetAlias("mus_ndof") != 0) {
		mus_ndof_branch = tree->GetBranch(tree->GetAlias("mus_ndof"));
		mus_ndof_branch->SetAddress(&mus_ndof_);
	}
	mus_overlapCompat_branch = 0;
	if (tree->GetAlias("mus_overlapCompat") != 0) {
		mus_overlapCompat_branch = tree->GetBranch(tree->GetAlias("mus_overlapCompat"));
		mus_overlapCompat_branch->SetAddress(&mus_overlapCompat_);
	}
	mus_phiErr_branch = 0;
	if (tree->GetAlias("mus_phiErr") != 0) {
		mus_phiErr_branch = tree->GetBranch(tree->GetAlias("mus_phiErr"));
		mus_phiErr_branch->SetAddress(&mus_phiErr_);
	}
	mus_ptErr_branch = 0;
	if (tree->GetAlias("mus_ptErr") != 0) {
		mus_ptErr_branch = tree->GetBranch(tree->GetAlias("mus_ptErr"));
		mus_ptErr_branch->SetAddress(&mus_ptErr_);
	}
	mus_qoverp_branch = 0;
	if (tree->GetAlias("mus_qoverp") != 0) {
		mus_qoverp_branch = tree->GetBranch(tree->GetAlias("mus_qoverp"));
		mus_qoverp_branch->SetAddress(&mus_qoverp_);
	}
	mus_qoverpError_branch = 0;
	if (tree->GetAlias("mus_qoverpError") != 0) {
		mus_qoverpError_branch = tree->GetBranch(tree->GetAlias("mus_qoverpError"));
		mus_qoverpError_branch->SetAddress(&mus_qoverpError_);
	}
	mus_segmCompatibility_branch = 0;
	if (tree->GetAlias("mus_segmCompatibility") != 0) {
		mus_segmCompatibility_branch = tree->GetBranch(tree->GetAlias("mus_segmCompatibility"));
		mus_segmCompatibility_branch->SetAddress(&mus_segmCompatibility_);
	}
	mus_staRelChi2_branch = 0;
	if (tree->GetAlias("mus_staRelChi2") != 0) {
		mus_staRelChi2_branch = tree->GetBranch(tree->GetAlias("mus_staRelChi2"));
		mus_staRelChi2_branch->SetAddress(&mus_staRelChi2_);
	}
	mus_sta_chi2_branch = 0;
	if (tree->GetAlias("mus_sta_chi2") != 0) {
		mus_sta_chi2_branch = tree->GetBranch(tree->GetAlias("mus_sta_chi2"));
		mus_sta_chi2_branch->SetAddress(&mus_sta_chi2_);
	}
	mus_sta_d0_branch = 0;
	if (tree->GetAlias("mus_sta_d0") != 0) {
		mus_sta_d0_branch = tree->GetBranch(tree->GetAlias("mus_sta_d0"));
		mus_sta_d0_branch->SetAddress(&mus_sta_d0_);
	}
	mus_sta_d0Err_branch = 0;
	if (tree->GetAlias("mus_sta_d0Err") != 0) {
		mus_sta_d0Err_branch = tree->GetBranch(tree->GetAlias("mus_sta_d0Err"));
		mus_sta_d0Err_branch->SetAddress(&mus_sta_d0Err_);
	}
	mus_sta_d0corr_branch = 0;
	if (tree->GetAlias("mus_sta_d0corr") != 0) {
		mus_sta_d0corr_branch = tree->GetBranch(tree->GetAlias("mus_sta_d0corr"));
		mus_sta_d0corr_branch->SetAddress(&mus_sta_d0corr_);
	}
	mus_sta_d0corrPhi_branch = 0;
	if (tree->GetAlias("mus_sta_d0corrPhi") != 0) {
		mus_sta_d0corrPhi_branch = tree->GetBranch(tree->GetAlias("mus_sta_d0corrPhi"));
		mus_sta_d0corrPhi_branch->SetAddress(&mus_sta_d0corrPhi_);
	}
	mus_sta_d0phiCov_branch = 0;
	if (tree->GetAlias("mus_sta_d0phiCov") != 0) {
		mus_sta_d0phiCov_branch = tree->GetBranch(tree->GetAlias("mus_sta_d0phiCov"));
		mus_sta_d0phiCov_branch->SetAddress(&mus_sta_d0phiCov_);
	}
	mus_sta_ndof_branch = 0;
	if (tree->GetAlias("mus_sta_ndof") != 0) {
		mus_sta_ndof_branch = tree->GetBranch(tree->GetAlias("mus_sta_ndof"));
		mus_sta_ndof_branch->SetAddress(&mus_sta_ndof_);
	}
	mus_sta_qoverp_branch = 0;
	if (tree->GetAlias("mus_sta_qoverp") != 0) {
		mus_sta_qoverp_branch = tree->GetBranch(tree->GetAlias("mus_sta_qoverp"));
		mus_sta_qoverp_branch->SetAddress(&mus_sta_qoverp_);
	}
	mus_sta_qoverpError_branch = 0;
	if (tree->GetAlias("mus_sta_qoverpError") != 0) {
		mus_sta_qoverpError_branch = tree->GetBranch(tree->GetAlias("mus_sta_qoverpError"));
		mus_sta_qoverpError_branch->SetAddress(&mus_sta_qoverpError_);
	}
	mus_sta_z0_branch = 0;
	if (tree->GetAlias("mus_sta_z0") != 0) {
		mus_sta_z0_branch = tree->GetBranch(tree->GetAlias("mus_sta_z0"));
		mus_sta_z0_branch->SetAddress(&mus_sta_z0_);
	}
	mus_sta_z0Err_branch = 0;
	if (tree->GetAlias("mus_sta_z0Err") != 0) {
		mus_sta_z0Err_branch = tree->GetBranch(tree->GetAlias("mus_sta_z0Err"));
		mus_sta_z0Err_branch->SetAddress(&mus_sta_z0Err_);
	}
	mus_sta_z0corr_branch = 0;
	if (tree->GetAlias("mus_sta_z0corr") != 0) {
		mus_sta_z0corr_branch = tree->GetBranch(tree->GetAlias("mus_sta_z0corr"));
		mus_sta_z0corr_branch->SetAddress(&mus_sta_z0corr_);
	}
	mus_timeAtIpInOut_branch = 0;
	if (tree->GetAlias("mus_timeAtIpInOut") != 0) {
		mus_timeAtIpInOut_branch = tree->GetBranch(tree->GetAlias("mus_timeAtIpInOut"));
		mus_timeAtIpInOut_branch->SetAddress(&mus_timeAtIpInOut_);
	}
	mus_timeAtIpInOutErr_branch = 0;
	if (tree->GetAlias("mus_timeAtIpInOutErr") != 0) {
		mus_timeAtIpInOutErr_branch = tree->GetBranch(tree->GetAlias("mus_timeAtIpInOutErr"));
		mus_timeAtIpInOutErr_branch->SetAddress(&mus_timeAtIpInOutErr_);
	}
	mus_timeAtIpOutIn_branch = 0;
	if (tree->GetAlias("mus_timeAtIpOutIn") != 0) {
		mus_timeAtIpOutIn_branch = tree->GetBranch(tree->GetAlias("mus_timeAtIpOutIn"));
		mus_timeAtIpOutIn_branch->SetAddress(&mus_timeAtIpOutIn_);
	}
	mus_timeAtIpOutInErr_branch = 0;
	if (tree->GetAlias("mus_timeAtIpOutInErr") != 0) {
		mus_timeAtIpOutInErr_branch = tree->GetBranch(tree->GetAlias("mus_timeAtIpOutInErr"));
		mus_timeAtIpOutInErr_branch->SetAddress(&mus_timeAtIpOutInErr_);
	}
	mus_timeCompat_branch = 0;
	if (tree->GetAlias("mus_timeCompat") != 0) {
		mus_timeCompat_branch = tree->GetBranch(tree->GetAlias("mus_timeCompat"));
		mus_timeCompat_branch->SetAddress(&mus_timeCompat_);
	}
	mus_tower_branch = 0;
	if (tree->GetAlias("mus_tower") != 0) {
		mus_tower_branch = tree->GetBranch(tree->GetAlias("mus_tower"));
		mus_tower_branch->SetAddress(&mus_tower_);
	}
	mus_towerS9_branch = 0;
	if (tree->GetAlias("mus_towerS9") != 0) {
		mus_towerS9_branch = tree->GetBranch(tree->GetAlias("mus_towerS9"));
		mus_towerS9_branch->SetAddress(&mus_towerS9_);
	}
	mus_trkKink_branch = 0;
	if (tree->GetAlias("mus_trkKink") != 0) {
		mus_trkKink_branch = tree->GetBranch(tree->GetAlias("mus_trkKink"));
		mus_trkKink_branch->SetAddress(&mus_trkKink_);
	}
	mus_trkRelChi2_branch = 0;
	if (tree->GetAlias("mus_trkRelChi2") != 0) {
		mus_trkRelChi2_branch = tree->GetBranch(tree->GetAlias("mus_trkRelChi2"));
		mus_trkRelChi2_branch->SetAddress(&mus_trkRelChi2_);
	}
	mus_vertexCompat_branch = 0;
	if (tree->GetAlias("mus_vertexCompat") != 0) {
		mus_vertexCompat_branch = tree->GetBranch(tree->GetAlias("mus_vertexCompat"));
		mus_vertexCompat_branch->SetAddress(&mus_vertexCompat_);
	}
	mus_vertexphi_branch = 0;
	if (tree->GetAlias("mus_vertexphi") != 0) {
		mus_vertexphi_branch = tree->GetBranch(tree->GetAlias("mus_vertexphi"));
		mus_vertexphi_branch->SetAddress(&mus_vertexphi_);
	}
	mus_z0_branch = 0;
	if (tree->GetAlias("mus_z0") != 0) {
		mus_z0_branch = tree->GetBranch(tree->GetAlias("mus_z0"));
		mus_z0_branch->SetAddress(&mus_z0_);
	}
	mus_z0Err_branch = 0;
	if (tree->GetAlias("mus_z0Err") != 0) {
		mus_z0Err_branch = tree->GetBranch(tree->GetAlias("mus_z0Err"));
		mus_z0Err_branch->SetAddress(&mus_z0Err_);
	}
	mus_z0corr_branch = 0;
	if (tree->GetAlias("mus_z0corr") != 0) {
		mus_z0corr_branch = tree->GetBranch(tree->GetAlias("mus_z0corr"));
		mus_z0corr_branch->SetAddress(&mus_z0corr_);
	}
	pfcands_dz_branch = 0;
	if (tree->GetAlias("pfcands_dz") != 0) {
		pfcands_dz_branch = tree->GetBranch(tree->GetAlias("pfcands_dz"));
		pfcands_dz_branch->SetAddress(&pfcands_dz_);
	}
	pfjets_area_branch = 0;
	if (tree->GetAlias("pfjets_area") != 0) {
		pfjets_area_branch = tree->GetBranch(tree->GetAlias("pfjets_area"));
		pfjets_area_branch->SetAddress(&pfjets_area_);
	}
	pfjets_chargedEmE_branch = 0;
	if (tree->GetAlias("pfjets_chargedEmE") != 0) {
		pfjets_chargedEmE_branch = tree->GetBranch(tree->GetAlias("pfjets_chargedEmE"));
		pfjets_chargedEmE_branch->SetAddress(&pfjets_chargedEmE_);
	}
	pfjets_chargedHadronE_branch = 0;
	if (tree->GetAlias("pfjets_chargedHadronE") != 0) {
		pfjets_chargedHadronE_branch = tree->GetBranch(tree->GetAlias("pfjets_chargedHadronE"));
		pfjets_chargedHadronE_branch->SetAddress(&pfjets_chargedHadronE_);
	}
	pfjets_combinedSecondaryVertexBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_combinedSecondaryVertexBJetTag") != 0) {
		pfjets_combinedSecondaryVertexBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_combinedSecondaryVertexBJetTag"));
		pfjets_combinedSecondaryVertexBJetTag_branch->SetAddress(&pfjets_combinedSecondaryVertexBJetTag_);
	}
	pfjets_cor_branch = 0;
	if (tree->GetAlias("pfjets_cor") != 0) {
		pfjets_cor_branch = tree->GetBranch(tree->GetAlias("pfjets_cor"));
		pfjets_cor_branch->SetAddress(&pfjets_cor_);
	}
	pfjets_corL1Fast_branch = 0;
	if (tree->GetAlias("pfjets_corL1Fast") != 0) {
		pfjets_corL1Fast_branch = tree->GetBranch(tree->GetAlias("pfjets_corL1Fast"));
		pfjets_corL1Fast_branch->SetAddress(&pfjets_corL1Fast_);
	}
	pfjets_corL1FastL2L3_branch = 0;
	if (tree->GetAlias("pfjets_corL1FastL2L3") != 0) {
		pfjets_corL1FastL2L3_branch = tree->GetBranch(tree->GetAlias("pfjets_corL1FastL2L3"));
		pfjets_corL1FastL2L3_branch->SetAddress(&pfjets_corL1FastL2L3_);
	}
	pfjets_corL1FastL2L3residual_branch = 0;
	if (tree->GetAlias("pfjets_corL1FastL2L3residual") != 0) {
		pfjets_corL1FastL2L3residual_branch = tree->GetBranch(tree->GetAlias("pfjets_corL1FastL2L3residual"));
		pfjets_corL1FastL2L3residual_branch->SetAddress(&pfjets_corL1FastL2L3residual_);
	}
	pfjets_corL1L2L3_branch = 0;
	if (tree->GetAlias("pfjets_corL1L2L3") != 0) {
		pfjets_corL1L2L3_branch = tree->GetBranch(tree->GetAlias("pfjets_corL1L2L3"));
		pfjets_corL1L2L3_branch->SetAddress(&pfjets_corL1L2L3_);
	}
	pfjets_electronE_branch = 0;
	if (tree->GetAlias("pfjets_electronE") != 0) {
		pfjets_electronE_branch = tree->GetBranch(tree->GetAlias("pfjets_electronE"));
		pfjets_electronE_branch->SetAddress(&pfjets_electronE_);
	}
	pfjets_hfEmE_branch = 0;
	if (tree->GetAlias("pfjets_hfEmE") != 0) {
		pfjets_hfEmE_branch = tree->GetBranch(tree->GetAlias("pfjets_hfEmE"));
		pfjets_hfEmE_branch->SetAddress(&pfjets_hfEmE_);
	}
	pfjets_hfHadronE_branch = 0;
	if (tree->GetAlias("pfjets_hfHadronE") != 0) {
		pfjets_hfHadronE_branch = tree->GetBranch(tree->GetAlias("pfjets_hfHadronE"));
		pfjets_hfHadronE_branch->SetAddress(&pfjets_hfHadronE_);
	}
	pfjets_jetBProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_jetBProbabilityBJetTag") != 0) {
		pfjets_jetBProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_jetBProbabilityBJetTag"));
		pfjets_jetBProbabilityBJetTag_branch->SetAddress(&pfjets_jetBProbabilityBJetTag_);
	}
	pfjets_jetProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_jetProbabilityBJetTag") != 0) {
		pfjets_jetProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_jetProbabilityBJetTag"));
		pfjets_jetProbabilityBJetTag_branch->SetAddress(&pfjets_jetProbabilityBJetTag_);
	}
	pfjets_muonE_branch = 0;
	if (tree->GetAlias("pfjets_muonE") != 0) {
		pfjets_muonE_branch = tree->GetBranch(tree->GetAlias("pfjets_muonE"));
		pfjets_muonE_branch->SetAddress(&pfjets_muonE_);
	}
	pfjets_neutralEmE_branch = 0;
	if (tree->GetAlias("pfjets_neutralEmE") != 0) {
		pfjets_neutralEmE_branch = tree->GetBranch(tree->GetAlias("pfjets_neutralEmE"));
		pfjets_neutralEmE_branch->SetAddress(&pfjets_neutralEmE_);
	}
	pfjets_neutralHadronE_branch = 0;
	if (tree->GetAlias("pfjets_neutralHadronE") != 0) {
		pfjets_neutralHadronE_branch = tree->GetBranch(tree->GetAlias("pfjets_neutralHadronE"));
		pfjets_neutralHadronE_branch->SetAddress(&pfjets_neutralHadronE_);
	}
	pfjets_photonE_branch = 0;
	if (tree->GetAlias("pfjets_photonE") != 0) {
		pfjets_photonE_branch = tree->GetBranch(tree->GetAlias("pfjets_photonE"));
		pfjets_photonE_branch->SetAddress(&pfjets_photonE_);
	}
	pfjets_pileupJetId_branch = 0;
	if (tree->GetAlias("pfjets_pileupJetId") != 0) {
		pfjets_pileupJetId_branch = tree->GetBranch(tree->GetAlias("pfjets_pileupJetId"));
		pfjets_pileupJetId_branch->SetAddress(&pfjets_pileupJetId_);
	}
	pfjets_simpleSecondaryVertexHighEffBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_simpleSecondaryVertexHighEffBJetTag") != 0) {
		pfjets_simpleSecondaryVertexHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_simpleSecondaryVertexHighEffBJetTag"));
		pfjets_simpleSecondaryVertexHighEffBJetTag_branch->SetAddress(&pfjets_simpleSecondaryVertexHighEffBJetTag_);
	}
	pfjets_simpleSecondaryVertexHighPurBJetTags_branch = 0;
	if (tree->GetAlias("pfjets_simpleSecondaryVertexHighPurBJetTags") != 0) {
		pfjets_simpleSecondaryVertexHighPurBJetTags_branch = tree->GetBranch(tree->GetAlias("pfjets_simpleSecondaryVertexHighPurBJetTags"));
		pfjets_simpleSecondaryVertexHighPurBJetTags_branch->SetAddress(&pfjets_simpleSecondaryVertexHighPurBJetTags_);
	}
	pfjets_trackCountingHighEffBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_trackCountingHighEffBJetTag") != 0) {
		pfjets_trackCountingHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_trackCountingHighEffBJetTag"));
		pfjets_trackCountingHighEffBJetTag_branch->SetAddress(&pfjets_trackCountingHighEffBJetTag_);
	}
	pfjets_trackCountingHighPurBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_trackCountingHighPurBJetTag") != 0) {
		pfjets_trackCountingHighPurBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_trackCountingHighPurBJetTag"));
		pfjets_trackCountingHighPurBJetTag_branch->SetAddress(&pfjets_trackCountingHighPurBJetTag_);
	}
	taus_pf_againstElectronDeadECAL_branch = 0;
	if (tree->GetAlias("taus_pf_againstElectronDeadECAL") != 0) {
		taus_pf_againstElectronDeadECAL_branch = tree->GetBranch(tree->GetAlias("taus_pf_againstElectronDeadECAL"));
		taus_pf_againstElectronDeadECAL_branch->SetAddress(&taus_pf_againstElectronDeadECAL_);
	}
	taus_pf_againstElectronLoose_branch = 0;
	if (tree->GetAlias("taus_pf_againstElectronLoose") != 0) {
		taus_pf_againstElectronLoose_branch = tree->GetBranch(tree->GetAlias("taus_pf_againstElectronLoose"));
		taus_pf_againstElectronLoose_branch->SetAddress(&taus_pf_againstElectronLoose_);
	}
	taus_pf_againstElectronMedium_branch = 0;
	if (tree->GetAlias("taus_pf_againstElectronMedium") != 0) {
		taus_pf_againstElectronMedium_branch = tree->GetBranch(tree->GetAlias("taus_pf_againstElectronMedium"));
		taus_pf_againstElectronMedium_branch->SetAddress(&taus_pf_againstElectronMedium_);
	}
	taus_pf_againstElectronTight_branch = 0;
	if (tree->GetAlias("taus_pf_againstElectronTight") != 0) {
		taus_pf_againstElectronTight_branch = tree->GetBranch(tree->GetAlias("taus_pf_againstElectronTight"));
		taus_pf_againstElectronTight_branch->SetAddress(&taus_pf_againstElectronTight_);
	}
	taus_pf_againstMuonLoose_branch = 0;
	if (tree->GetAlias("taus_pf_againstMuonLoose") != 0) {
		taus_pf_againstMuonLoose_branch = tree->GetBranch(tree->GetAlias("taus_pf_againstMuonLoose"));
		taus_pf_againstMuonLoose_branch->SetAddress(&taus_pf_againstMuonLoose_);
	}
	taus_pf_againstMuonLoose2_branch = 0;
	if (tree->GetAlias("taus_pf_againstMuonLoose2") != 0) {
		taus_pf_againstMuonLoose2_branch = tree->GetBranch(tree->GetAlias("taus_pf_againstMuonLoose2"));
		taus_pf_againstMuonLoose2_branch->SetAddress(&taus_pf_againstMuonLoose2_);
	}
	taus_pf_againstMuonLoose3_branch = 0;
	if (tree->GetAlias("taus_pf_againstMuonLoose3") != 0) {
		taus_pf_againstMuonLoose3_branch = tree->GetBranch(tree->GetAlias("taus_pf_againstMuonLoose3"));
		taus_pf_againstMuonLoose3_branch->SetAddress(&taus_pf_againstMuonLoose3_);
	}
	taus_pf_againstMuonMedium_branch = 0;
	if (tree->GetAlias("taus_pf_againstMuonMedium") != 0) {
		taus_pf_againstMuonMedium_branch = tree->GetBranch(tree->GetAlias("taus_pf_againstMuonMedium"));
		taus_pf_againstMuonMedium_branch->SetAddress(&taus_pf_againstMuonMedium_);
	}
	taus_pf_againstMuonMedium2_branch = 0;
	if (tree->GetAlias("taus_pf_againstMuonMedium2") != 0) {
		taus_pf_againstMuonMedium2_branch = tree->GetBranch(tree->GetAlias("taus_pf_againstMuonMedium2"));
		taus_pf_againstMuonMedium2_branch->SetAddress(&taus_pf_againstMuonMedium2_);
	}
	taus_pf_againstMuonTight_branch = 0;
	if (tree->GetAlias("taus_pf_againstMuonTight") != 0) {
		taus_pf_againstMuonTight_branch = tree->GetBranch(tree->GetAlias("taus_pf_againstMuonTight"));
		taus_pf_againstMuonTight_branch->SetAddress(&taus_pf_againstMuonTight_);
	}
	taus_pf_againstMuonTight2_branch = 0;
	if (tree->GetAlias("taus_pf_againstMuonTight2") != 0) {
		taus_pf_againstMuonTight2_branch = tree->GetBranch(tree->GetAlias("taus_pf_againstMuonTight2"));
		taus_pf_againstMuonTight2_branch->SetAddress(&taus_pf_againstMuonTight2_);
	}
	taus_pf_againstMuonTight3_branch = 0;
	if (tree->GetAlias("taus_pf_againstMuonTight3") != 0) {
		taus_pf_againstMuonTight3_branch = tree->GetBranch(tree->GetAlias("taus_pf_againstMuonTight3"));
		taus_pf_againstMuonTight3_branch->SetAddress(&taus_pf_againstMuonTight3_);
	}
	taus_pf_byCombinedIsolationDeltaBetaCorrRaw_branch = 0;
	if (tree->GetAlias("taus_pf_byCombinedIsolationDeltaBetaCorrRaw") != 0) {
		taus_pf_byCombinedIsolationDeltaBetaCorrRaw_branch = tree->GetBranch(tree->GetAlias("taus_pf_byCombinedIsolationDeltaBetaCorrRaw"));
		taus_pf_byCombinedIsolationDeltaBetaCorrRaw_branch->SetAddress(&taus_pf_byCombinedIsolationDeltaBetaCorrRaw_);
	}
	taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_branch = 0;
	if (tree->GetAlias("taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits") != 0) {
		taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_branch = tree->GetBranch(tree->GetAlias("taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits"));
		taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_branch->SetAddress(&taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_);
	}
	taus_pf_byDecayModeFinding_branch = 0;
	if (tree->GetAlias("taus_pf_byDecayModeFinding") != 0) {
		taus_pf_byDecayModeFinding_branch = tree->GetBranch(tree->GetAlias("taus_pf_byDecayModeFinding"));
		taus_pf_byDecayModeFinding_branch->SetAddress(&taus_pf_byDecayModeFinding_);
	}
	taus_pf_byLooseCombinedIsolationDeltaBetaCorr_branch = 0;
	if (tree->GetAlias("taus_pf_byLooseCombinedIsolationDeltaBetaCorr") != 0) {
		taus_pf_byLooseCombinedIsolationDeltaBetaCorr_branch = tree->GetBranch(tree->GetAlias("taus_pf_byLooseCombinedIsolationDeltaBetaCorr"));
		taus_pf_byLooseCombinedIsolationDeltaBetaCorr_branch->SetAddress(&taus_pf_byLooseCombinedIsolationDeltaBetaCorr_);
	}
	taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_branch = 0;
	if (tree->GetAlias("taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits") != 0) {
		taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_branch = tree->GetBranch(tree->GetAlias("taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits"));
		taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_branch->SetAddress(&taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_);
	}
	taus_pf_byMediumCombinedIsolationDeltaBetaCorr_branch = 0;
	if (tree->GetAlias("taus_pf_byMediumCombinedIsolationDeltaBetaCorr") != 0) {
		taus_pf_byMediumCombinedIsolationDeltaBetaCorr_branch = tree->GetBranch(tree->GetAlias("taus_pf_byMediumCombinedIsolationDeltaBetaCorr"));
		taus_pf_byMediumCombinedIsolationDeltaBetaCorr_branch->SetAddress(&taus_pf_byMediumCombinedIsolationDeltaBetaCorr_);
	}
	taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_branch = 0;
	if (tree->GetAlias("taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits") != 0) {
		taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_branch = tree->GetBranch(tree->GetAlias("taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits"));
		taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_branch->SetAddress(&taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_);
	}
	taus_pf_byTightCombinedIsolationDeltaBetaCorr_branch = 0;
	if (tree->GetAlias("taus_pf_byTightCombinedIsolationDeltaBetaCorr") != 0) {
		taus_pf_byTightCombinedIsolationDeltaBetaCorr_branch = tree->GetBranch(tree->GetAlias("taus_pf_byTightCombinedIsolationDeltaBetaCorr"));
		taus_pf_byTightCombinedIsolationDeltaBetaCorr_branch->SetAddress(&taus_pf_byTightCombinedIsolationDeltaBetaCorr_);
	}
	taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_branch = 0;
	if (tree->GetAlias("taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits") != 0) {
		taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_branch = tree->GetBranch(tree->GetAlias("taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits"));
		taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_branch->SetAddress(&taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_);
	}
	taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_branch = 0;
	if (tree->GetAlias("taus_pf_byVLooseCombinedIsolationDeltaBetaCorr") != 0) {
		taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_branch = tree->GetBranch(tree->GetAlias("taus_pf_byVLooseCombinedIsolationDeltaBetaCorr"));
		taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_branch->SetAddress(&taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_);
	}
	photons_e1x5_branch = 0;
	if (tree->GetAlias("photons_e1x5") != 0) {
		photons_e1x5_branch = tree->GetBranch(tree->GetAlias("photons_e1x5"));
		photons_e1x5_branch->SetAddress(&photons_e1x5_);
	}
	photons_e2x5Max_branch = 0;
	if (tree->GetAlias("photons_e2x5Max") != 0) {
		photons_e2x5Max_branch = tree->GetBranch(tree->GetAlias("photons_e2x5Max"));
		photons_e2x5Max_branch->SetAddress(&photons_e2x5Max_);
	}
	photons_e3x3_branch = 0;
	if (tree->GetAlias("photons_e3x3") != 0) {
		photons_e3x3_branch = tree->GetBranch(tree->GetAlias("photons_e3x3"));
		photons_e3x3_branch->SetAddress(&photons_e3x3_);
	}
	photons_e5x5_branch = 0;
	if (tree->GetAlias("photons_e5x5") != 0) {
		photons_e5x5_branch = tree->GetBranch(tree->GetAlias("photons_e5x5"));
		photons_e5x5_branch->SetAddress(&photons_e5x5_);
	}
	photons_ecalIso03_branch = 0;
	if (tree->GetAlias("photons_ecalIso03") != 0) {
		photons_ecalIso03_branch = tree->GetBranch(tree->GetAlias("photons_ecalIso03"));
		photons_ecalIso03_branch->SetAddress(&photons_ecalIso03_);
	}
	photons_ecalIso04_branch = 0;
	if (tree->GetAlias("photons_ecalIso04") != 0) {
		photons_ecalIso04_branch = tree->GetBranch(tree->GetAlias("photons_ecalIso04"));
		photons_ecalIso04_branch->SetAddress(&photons_ecalIso04_);
	}
	photons_full3x3_e3x3_branch = 0;
	if (tree->GetAlias("photons_full3x3_e3x3") != 0) {
		photons_full3x3_e3x3_branch = tree->GetBranch(tree->GetAlias("photons_full3x3_e3x3"));
		photons_full3x3_e3x3_branch->SetAddress(&photons_full3x3_e3x3_);
	}
	photons_full5x5_e1x5_branch = 0;
	if (tree->GetAlias("photons_full5x5_e1x5") != 0) {
		photons_full5x5_e1x5_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_e1x5"));
		photons_full5x5_e1x5_branch->SetAddress(&photons_full5x5_e1x5_);
	}
	photons_full5x5_e2x5Max_branch = 0;
	if (tree->GetAlias("photons_full5x5_e2x5Max") != 0) {
		photons_full5x5_e2x5Max_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_e2x5Max"));
		photons_full5x5_e2x5Max_branch->SetAddress(&photons_full5x5_e2x5Max_);
	}
	photons_full5x5_e5x5_branch = 0;
	if (tree->GetAlias("photons_full5x5_e5x5") != 0) {
		photons_full5x5_e5x5_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_e5x5"));
		photons_full5x5_e5x5_branch->SetAddress(&photons_full5x5_e5x5_);
	}
	photons_full5x5_hOverE_branch = 0;
	if (tree->GetAlias("photons_full5x5_hOverE") != 0) {
		photons_full5x5_hOverE_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_hOverE"));
		photons_full5x5_hOverE_branch->SetAddress(&photons_full5x5_hOverE_);
	}
	photons_full5x5_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("photons_full5x5_sigmaEtaEta") != 0) {
		photons_full5x5_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_sigmaEtaEta"));
		photons_full5x5_sigmaEtaEta_branch->SetAddress(&photons_full5x5_sigmaEtaEta_);
	}
	photons_full5x5_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("photons_full5x5_sigmaIEtaIEta") != 0) {
		photons_full5x5_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_sigmaIEtaIEta"));
		photons_full5x5_sigmaIEtaIEta_branch->SetAddress(&photons_full5x5_sigmaIEtaIEta_);
	}
	photons_hOverE_branch = 0;
	if (tree->GetAlias("photons_hOverE") != 0) {
		photons_hOverE_branch = tree->GetBranch(tree->GetAlias("photons_hOverE"));
		photons_hOverE_branch->SetAddress(&photons_hOverE_);
	}
	photons_hcalDepth1TowerSumEtBcConeDR03_branch = 0;
	if (tree->GetAlias("photons_hcalDepth1TowerSumEtBcConeDR03") != 0) {
		photons_hcalDepth1TowerSumEtBcConeDR03_branch = tree->GetBranch(tree->GetAlias("photons_hcalDepth1TowerSumEtBcConeDR03"));
		photons_hcalDepth1TowerSumEtBcConeDR03_branch->SetAddress(&photons_hcalDepth1TowerSumEtBcConeDR03_);
	}
	photons_hcalDepth1TowerSumEtBcConeDR04_branch = 0;
	if (tree->GetAlias("photons_hcalDepth1TowerSumEtBcConeDR04") != 0) {
		photons_hcalDepth1TowerSumEtBcConeDR04_branch = tree->GetBranch(tree->GetAlias("photons_hcalDepth1TowerSumEtBcConeDR04"));
		photons_hcalDepth1TowerSumEtBcConeDR04_branch->SetAddress(&photons_hcalDepth1TowerSumEtBcConeDR04_);
	}
	photons_hcalDepth2TowerSumEtBcConeDR03_branch = 0;
	if (tree->GetAlias("photons_hcalDepth2TowerSumEtBcConeDR03") != 0) {
		photons_hcalDepth2TowerSumEtBcConeDR03_branch = tree->GetBranch(tree->GetAlias("photons_hcalDepth2TowerSumEtBcConeDR03"));
		photons_hcalDepth2TowerSumEtBcConeDR03_branch->SetAddress(&photons_hcalDepth2TowerSumEtBcConeDR03_);
	}
	photons_hcalDepth2TowerSumEtBcConeDR04_branch = 0;
	if (tree->GetAlias("photons_hcalDepth2TowerSumEtBcConeDR04") != 0) {
		photons_hcalDepth2TowerSumEtBcConeDR04_branch = tree->GetBranch(tree->GetAlias("photons_hcalDepth2TowerSumEtBcConeDR04"));
		photons_hcalDepth2TowerSumEtBcConeDR04_branch->SetAddress(&photons_hcalDepth2TowerSumEtBcConeDR04_);
	}
	photons_hcalIso03_branch = 0;
	if (tree->GetAlias("photons_hcalIso03") != 0) {
		photons_hcalIso03_branch = tree->GetBranch(tree->GetAlias("photons_hcalIso03"));
		photons_hcalIso03_branch->SetAddress(&photons_hcalIso03_);
	}
	photons_hcalIso04_branch = 0;
	if (tree->GetAlias("photons_hcalIso04") != 0) {
		photons_hcalIso04_branch = tree->GetBranch(tree->GetAlias("photons_hcalIso04"));
		photons_hcalIso04_branch->SetAddress(&photons_hcalIso04_);
	}
	photons_hcalTowerSumEtBcConeDR03_branch = 0;
	if (tree->GetAlias("photons_hcalTowerSumEtBcConeDR03") != 0) {
		photons_hcalTowerSumEtBcConeDR03_branch = tree->GetBranch(tree->GetAlias("photons_hcalTowerSumEtBcConeDR03"));
		photons_hcalTowerSumEtBcConeDR03_branch->SetAddress(&photons_hcalTowerSumEtBcConeDR03_);
	}
	photons_hcalTowerSumEtBcConeDR04_branch = 0;
	if (tree->GetAlias("photons_hcalTowerSumEtBcConeDR04") != 0) {
		photons_hcalTowerSumEtBcConeDR04_branch = tree->GetBranch(tree->GetAlias("photons_hcalTowerSumEtBcConeDR04"));
		photons_hcalTowerSumEtBcConeDR04_branch->SetAddress(&photons_hcalTowerSumEtBcConeDR04_);
	}
	photons_ntkIsoHollow03_branch = 0;
	if (tree->GetAlias("photons_ntkIsoHollow03") != 0) {
		photons_ntkIsoHollow03_branch = tree->GetBranch(tree->GetAlias("photons_ntkIsoHollow03"));
		photons_ntkIsoHollow03_branch->SetAddress(&photons_ntkIsoHollow03_);
	}
	photons_ntkIsoHollow04_branch = 0;
	if (tree->GetAlias("photons_ntkIsoHollow04") != 0) {
		photons_ntkIsoHollow04_branch = tree->GetBranch(tree->GetAlias("photons_ntkIsoHollow04"));
		photons_ntkIsoHollow04_branch->SetAddress(&photons_ntkIsoHollow04_);
	}
	photons_ntkIsoSolid03_branch = 0;
	if (tree->GetAlias("photons_ntkIsoSolid03") != 0) {
		photons_ntkIsoSolid03_branch = tree->GetBranch(tree->GetAlias("photons_ntkIsoSolid03"));
		photons_ntkIsoSolid03_branch->SetAddress(&photons_ntkIsoSolid03_);
	}
	photons_ntkIsoSolid04_branch = 0;
	if (tree->GetAlias("photons_ntkIsoSolid04") != 0) {
		photons_ntkIsoSolid04_branch = tree->GetBranch(tree->GetAlias("photons_ntkIsoSolid04"));
		photons_ntkIsoSolid04_branch->SetAddress(&photons_ntkIsoSolid04_);
	}
	photons_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("photons_sigmaEtaEta") != 0) {
		photons_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("photons_sigmaEtaEta"));
		photons_sigmaEtaEta_branch->SetAddress(&photons_sigmaEtaEta_);
	}
	photons_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("photons_sigmaIEtaIEta") != 0) {
		photons_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("photons_sigmaIEtaIEta"));
		photons_sigmaIEtaIEta_branch->SetAddress(&photons_sigmaIEtaIEta_);
	}
	photons_tkIsoHollow03_branch = 0;
	if (tree->GetAlias("photons_tkIsoHollow03") != 0) {
		photons_tkIsoHollow03_branch = tree->GetBranch(tree->GetAlias("photons_tkIsoHollow03"));
		photons_tkIsoHollow03_branch->SetAddress(&photons_tkIsoHollow03_);
	}
	photons_tkIsoHollow04_branch = 0;
	if (tree->GetAlias("photons_tkIsoHollow04") != 0) {
		photons_tkIsoHollow04_branch = tree->GetBranch(tree->GetAlias("photons_tkIsoHollow04"));
		photons_tkIsoHollow04_branch->SetAddress(&photons_tkIsoHollow04_);
	}
	photons_tkIsoSolid03_branch = 0;
	if (tree->GetAlias("photons_tkIsoSolid03") != 0) {
		photons_tkIsoSolid03_branch = tree->GetBranch(tree->GetAlias("photons_tkIsoSolid03"));
		photons_tkIsoSolid03_branch->SetAddress(&photons_tkIsoSolid03_);
	}
	photons_tkIsoSolid04_branch = 0;
	if (tree->GetAlias("photons_tkIsoSolid04") != 0) {
		photons_tkIsoSolid04_branch = tree->GetBranch(tree->GetAlias("photons_tkIsoSolid04"));
		photons_tkIsoSolid04_branch->SetAddress(&photons_tkIsoSolid04_);
	}
	puInfo_trueNumInteractions_branch = 0;
	if (tree->GetAlias("puInfo_trueNumInteractions") != 0) {
		puInfo_trueNumInteractions_branch = tree->GetBranch(tree->GetAlias("puInfo_trueNumInteractions"));
		puInfo_trueNumInteractions_branch->SetAddress(&puInfo_trueNumInteractions_);
	}
	convs_chi2_branch = 0;
	if (tree->GetAlias("convs_chi2") != 0) {
		convs_chi2_branch = tree->GetBranch(tree->GetAlias("convs_chi2"));
		convs_chi2_branch->SetAddress(&convs_chi2_);
	}
	convs_dl_branch = 0;
	if (tree->GetAlias("convs_dl") != 0) {
		convs_dl_branch = tree->GetBranch(tree->GetAlias("convs_dl"));
		convs_dl_branch->SetAddress(&convs_dl_);
	}
	convs_ndof_branch = 0;
	if (tree->GetAlias("convs_ndof") != 0) {
		convs_ndof_branch = tree->GetBranch(tree->GetAlias("convs_ndof"));
		convs_ndof_branch->SetAddress(&convs_ndof_);
	}
	vtxs_chi2_branch = 0;
	if (tree->GetAlias("vtxs_chi2") != 0) {
		vtxs_chi2_branch = tree->GetBranch(tree->GetAlias("vtxs_chi2"));
		vtxs_chi2_branch->SetAddress(&vtxs_chi2_);
	}
	vtxs_ndof_branch = 0;
	if (tree->GetAlias("vtxs_ndof") != 0) {
		vtxs_ndof_branch = tree->GetBranch(tree->GetAlias("vtxs_ndof"));
		vtxs_ndof_branch->SetAddress(&vtxs_ndof_);
	}
	vtxs_sumpt_branch = 0;
	if (tree->GetAlias("vtxs_sumpt") != 0) {
		vtxs_sumpt_branch = tree->GetBranch(tree->GetAlias("vtxs_sumpt"));
		vtxs_sumpt_branch->SetAddress(&vtxs_sumpt_);
	}
	vtxs_xError_branch = 0;
	if (tree->GetAlias("vtxs_xError") != 0) {
		vtxs_xError_branch = tree->GetBranch(tree->GetAlias("vtxs_xError"));
		vtxs_xError_branch->SetAddress(&vtxs_xError_);
	}
	vtxs_yError_branch = 0;
	if (tree->GetAlias("vtxs_yError") != 0) {
		vtxs_yError_branch = tree->GetBranch(tree->GetAlias("vtxs_yError"));
		vtxs_yError_branch->SetAddress(&vtxs_yError_);
	}
	vtxs_zError_branch = 0;
	if (tree->GetAlias("vtxs_zError") != 0) {
		vtxs_zError_branch = tree->GetBranch(tree->GetAlias("vtxs_zError"));
		vtxs_zError_branch->SetAddress(&vtxs_zError_);
	}
	els_convs_dcot_branch = 0;
	if (tree->GetAlias("els_convs_dcot") != 0) {
		els_convs_dcot_branch = tree->GetBranch(tree->GetAlias("els_convs_dcot"));
		els_convs_dcot_branch->SetAddress(&els_convs_dcot_);
	}
	els_convs_dist_branch = 0;
	if (tree->GetAlias("els_convs_dist") != 0) {
		els_convs_dist_branch = tree->GetBranch(tree->GetAlias("els_convs_dist"));
		els_convs_dist_branch->SetAddress(&els_convs_dist_);
	}
	els_convs_radius_branch = 0;
	if (tree->GetAlias("els_convs_radius") != 0) {
		els_convs_radius_branch = tree->GetBranch(tree->GetAlias("els_convs_radius"));
		els_convs_radius_branch->SetAddress(&els_convs_radius_);
	}
	mus_stationShowerDeltaR_branch = 0;
	if (tree->GetAlias("mus_stationShowerDeltaR") != 0) {
		mus_stationShowerDeltaR_branch = tree->GetBranch(tree->GetAlias("mus_stationShowerDeltaR"));
		mus_stationShowerDeltaR_branch->SetAddress(&mus_stationShowerDeltaR_);
	}
	mus_stationShowerSizeT_branch = 0;
	if (tree->GetAlias("mus_stationShowerSizeT") != 0) {
		mus_stationShowerSizeT_branch = tree->GetBranch(tree->GetAlias("mus_stationShowerSizeT"));
		mus_stationShowerSizeT_branch->SetAddress(&mus_stationShowerSizeT_);
	}
	puInfo_instLumi_branch = 0;
	if (tree->GetAlias("puInfo_instLumi") != 0) {
		puInfo_instLumi_branch = tree->GetBranch(tree->GetAlias("puInfo_instLumi"));
		puInfo_instLumi_branch->SetAddress(&puInfo_instLumi_);
	}
	vtxs_covMatrix_branch = 0;
	if (tree->GetAlias("vtxs_covMatrix") != 0) {
		vtxs_covMatrix_branch = tree->GetBranch(tree->GetAlias("vtxs_covMatrix"));
		vtxs_covMatrix_branch->SetAddress(&vtxs_covMatrix_);
	}
	evt_bsType_branch = 0;
	if (tree->GetAlias("evt_bsType") != 0) {
		evt_bsType_branch = tree->GetBranch(tree->GetAlias("evt_bsType"));
		evt_bsType_branch->SetAddress(&evt_bsType_);
	}
	evt_bunchCrossing_branch = 0;
	if (tree->GetAlias("evt_bunchCrossing") != 0) {
		evt_bunchCrossing_branch = tree->GetBranch(tree->GetAlias("evt_bunchCrossing"));
		evt_bunchCrossing_branch->SetAddress(&evt_bunchCrossing_);
	}
	evt_experimentType_branch = 0;
	if (tree->GetAlias("evt_experimentType") != 0) {
		evt_experimentType_branch = tree->GetBranch(tree->GetAlias("evt_experimentType"));
		evt_experimentType_branch->SetAddress(&evt_experimentType_);
	}
	evt_isRealData_branch = 0;
	if (tree->GetAlias("evt_isRealData") != 0) {
		evt_isRealData_branch = tree->GetBranch(tree->GetAlias("evt_isRealData"));
		evt_isRealData_branch->SetAddress(&evt_isRealData_);
	}
	evt_orbitNumber_branch = 0;
	if (tree->GetAlias("evt_orbitNumber") != 0) {
		evt_orbitNumber_branch = tree->GetBranch(tree->GetAlias("evt_orbitNumber"));
		evt_orbitNumber_branch->SetAddress(&evt_orbitNumber_);
	}
	evt_storeNumber_branch = 0;
	if (tree->GetAlias("evt_storeNumber") != 0) {
		evt_storeNumber_branch = tree->GetBranch(tree->GetAlias("evt_storeNumber"));
		evt_storeNumber_branch->SetAddress(&evt_storeNumber_);
	}
	hcalnoise_GetRecHitCount_branch = 0;
	if (tree->GetAlias("hcalnoise_GetRecHitCount") != 0) {
		hcalnoise_GetRecHitCount_branch = tree->GetBranch(tree->GetAlias("hcalnoise_GetRecHitCount"));
		hcalnoise_GetRecHitCount_branch->SetAddress(&hcalnoise_GetRecHitCount_);
	}
	hcalnoise_GetRecHitCount15_branch = 0;
	if (tree->GetAlias("hcalnoise_GetRecHitCount15") != 0) {
		hcalnoise_GetRecHitCount15_branch = tree->GetBranch(tree->GetAlias("hcalnoise_GetRecHitCount15"));
		hcalnoise_GetRecHitCount15_branch->SetAddress(&hcalnoise_GetRecHitCount15_);
	}
	hcalnoise_maxHPDHits_branch = 0;
	if (tree->GetAlias("hcalnoise_maxHPDHits") != 0) {
		hcalnoise_maxHPDHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxHPDHits"));
		hcalnoise_maxHPDHits_branch->SetAddress(&hcalnoise_maxHPDHits_);
	}
	hcalnoise_maxHPDNoOtherHits_branch = 0;
	if (tree->GetAlias("hcalnoise_maxHPDNoOtherHits") != 0) {
		hcalnoise_maxHPDNoOtherHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxHPDNoOtherHits"));
		hcalnoise_maxHPDNoOtherHits_branch->SetAddress(&hcalnoise_maxHPDNoOtherHits_);
	}
	hcalnoise_maxRBXHits_branch = 0;
	if (tree->GetAlias("hcalnoise_maxRBXHits") != 0) {
		hcalnoise_maxRBXHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxRBXHits"));
		hcalnoise_maxRBXHits_branch->SetAddress(&hcalnoise_maxRBXHits_);
	}
	hcalnoise_maxZeros_branch = 0;
	if (tree->GetAlias("hcalnoise_maxZeros") != 0) {
		hcalnoise_maxZeros_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxZeros"));
		hcalnoise_maxZeros_branch->SetAddress(&hcalnoise_maxZeros_);
	}
	hcalnoise_noiseFilterStatus_branch = 0;
	if (tree->GetAlias("hcalnoise_noiseFilterStatus") != 0) {
		hcalnoise_noiseFilterStatus_branch = tree->GetBranch(tree->GetAlias("hcalnoise_noiseFilterStatus"));
		hcalnoise_noiseFilterStatus_branch->SetAddress(&hcalnoise_noiseFilterStatus_);
	}
	hcalnoise_noiseType_branch = 0;
	if (tree->GetAlias("hcalnoise_noiseType") != 0) {
		hcalnoise_noiseType_branch = tree->GetBranch(tree->GetAlias("hcalnoise_noiseType"));
		hcalnoise_noiseType_branch->SetAddress(&hcalnoise_noiseType_);
	}
	hcalnoise_num10GeVHits_branch = 0;
	if (tree->GetAlias("hcalnoise_num10GeVHits") != 0) {
		hcalnoise_num10GeVHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_num10GeVHits"));
		hcalnoise_num10GeVHits_branch->SetAddress(&hcalnoise_num10GeVHits_);
	}
	hcalnoise_num25GeVHits_branch = 0;
	if (tree->GetAlias("hcalnoise_num25GeVHits") != 0) {
		hcalnoise_num25GeVHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_num25GeVHits"));
		hcalnoise_num25GeVHits_branch->SetAddress(&hcalnoise_num25GeVHits_);
	}
	hcalnoise_numFlatNoiseChannels_branch = 0;
	if (tree->GetAlias("hcalnoise_numFlatNoiseChannels") != 0) {
		hcalnoise_numFlatNoiseChannels_branch = tree->GetBranch(tree->GetAlias("hcalnoise_numFlatNoiseChannels"));
		hcalnoise_numFlatNoiseChannels_branch->SetAddress(&hcalnoise_numFlatNoiseChannels_);
	}
	hcalnoise_numIsolatedNoiseChannels_branch = 0;
	if (tree->GetAlias("hcalnoise_numIsolatedNoiseChannels") != 0) {
		hcalnoise_numIsolatedNoiseChannels_branch = tree->GetBranch(tree->GetAlias("hcalnoise_numIsolatedNoiseChannels"));
		hcalnoise_numIsolatedNoiseChannels_branch->SetAddress(&hcalnoise_numIsolatedNoiseChannels_);
	}
	hcalnoise_numProblematicRBXs_branch = 0;
	if (tree->GetAlias("hcalnoise_numProblematicRBXs") != 0) {
		hcalnoise_numProblematicRBXs_branch = tree->GetBranch(tree->GetAlias("hcalnoise_numProblematicRBXs"));
		hcalnoise_numProblematicRBXs_branch->SetAddress(&hcalnoise_numProblematicRBXs_);
	}
	hcalnoise_numSpikeNoiseChannels_branch = 0;
	if (tree->GetAlias("hcalnoise_numSpikeNoiseChannels") != 0) {
		hcalnoise_numSpikeNoiseChannels_branch = tree->GetBranch(tree->GetAlias("hcalnoise_numSpikeNoiseChannels"));
		hcalnoise_numSpikeNoiseChannels_branch->SetAddress(&hcalnoise_numSpikeNoiseChannels_);
	}
	hcalnoise_numTS4TS5NoiseChannels_branch = 0;
	if (tree->GetAlias("hcalnoise_numTS4TS5NoiseChannels") != 0) {
		hcalnoise_numTS4TS5NoiseChannels_branch = tree->GetBranch(tree->GetAlias("hcalnoise_numTS4TS5NoiseChannels"));
		hcalnoise_numTS4TS5NoiseChannels_branch->SetAddress(&hcalnoise_numTS4TS5NoiseChannels_);
	}
	hcalnoise_numTriangleNoiseChannels_branch = 0;
	if (tree->GetAlias("hcalnoise_numTriangleNoiseChannels") != 0) {
		hcalnoise_numTriangleNoiseChannels_branch = tree->GetBranch(tree->GetAlias("hcalnoise_numTriangleNoiseChannels"));
		hcalnoise_numTriangleNoiseChannels_branch->SetAddress(&hcalnoise_numTriangleNoiseChannels_);
	}
	hcalnoise_passHighLevelNoiseFilter_branch = 0;
	if (tree->GetAlias("hcalnoise_passHighLevelNoiseFilter") != 0) {
		hcalnoise_passHighLevelNoiseFilter_branch = tree->GetBranch(tree->GetAlias("hcalnoise_passHighLevelNoiseFilter"));
		hcalnoise_passHighLevelNoiseFilter_branch->SetAddress(&hcalnoise_passHighLevelNoiseFilter_);
	}
	hcalnoise_passLooseNoiseFilter_branch = 0;
	if (tree->GetAlias("hcalnoise_passLooseNoiseFilter") != 0) {
		hcalnoise_passLooseNoiseFilter_branch = tree->GetBranch(tree->GetAlias("hcalnoise_passLooseNoiseFilter"));
		hcalnoise_passLooseNoiseFilter_branch->SetAddress(&hcalnoise_passLooseNoiseFilter_);
	}
	hcalnoise_passTightNoiseFilter_branch = 0;
	if (tree->GetAlias("hcalnoise_passTightNoiseFilter") != 0) {
		hcalnoise_passTightNoiseFilter_branch = tree->GetBranch(tree->GetAlias("hcalnoise_passTightNoiseFilter"));
		hcalnoise_passTightNoiseFilter_branch->SetAddress(&hcalnoise_passTightNoiseFilter_);
	}
	pdfinfo_id1_branch = 0;
	if (tree->GetAlias("pdfinfo_id1") != 0) {
		pdfinfo_id1_branch = tree->GetBranch(tree->GetAlias("pdfinfo_id1"));
		pdfinfo_id1_branch->SetAddress(&pdfinfo_id1_);
	}
	pdfinfo_id2_branch = 0;
	if (tree->GetAlias("pdfinfo_id2") != 0) {
		pdfinfo_id2_branch = tree->GetBranch(tree->GetAlias("pdfinfo_id2"));
		pdfinfo_id2_branch->SetAddress(&pdfinfo_id2_);
	}
	els_mc3_id_branch = 0;
	if (tree->GetAlias("els_mc3_id") != 0) {
		els_mc3_id_branch = tree->GetBranch(tree->GetAlias("els_mc3_id"));
		els_mc3_id_branch->SetAddress(&els_mc3_id_);
	}
	els_mc3idx_branch = 0;
	if (tree->GetAlias("els_mc3idx") != 0) {
		els_mc3idx_branch = tree->GetBranch(tree->GetAlias("els_mc3idx"));
		els_mc3idx_branch->SetAddress(&els_mc3idx_);
	}
	els_mc3_motherid_branch = 0;
	if (tree->GetAlias("els_mc3_motherid") != 0) {
		els_mc3_motherid_branch = tree->GetBranch(tree->GetAlias("els_mc3_motherid"));
		els_mc3_motherid_branch->SetAddress(&els_mc3_motherid_);
	}
	els_mc3_motheridx_branch = 0;
	if (tree->GetAlias("els_mc3_motheridx") != 0) {
		els_mc3_motheridx_branch = tree->GetBranch(tree->GetAlias("els_mc3_motheridx"));
		els_mc3_motheridx_branch->SetAddress(&els_mc3_motheridx_);
	}
	els_mc_id_branch = 0;
	if (tree->GetAlias("els_mc_id") != 0) {
		els_mc_id_branch = tree->GetBranch(tree->GetAlias("els_mc_id"));
		els_mc_id_branch->SetAddress(&els_mc_id_);
	}
	els_mcidx_branch = 0;
	if (tree->GetAlias("els_mcidx") != 0) {
		els_mcidx_branch = tree->GetBranch(tree->GetAlias("els_mcidx"));
		els_mcidx_branch->SetAddress(&els_mcidx_);
	}
	els_mc_motherid_branch = 0;
	if (tree->GetAlias("els_mc_motherid") != 0) {
		els_mc_motherid_branch = tree->GetBranch(tree->GetAlias("els_mc_motherid"));
		els_mc_motherid_branch->SetAddress(&els_mc_motherid_);
	}
	mus_mc3_id_branch = 0;
	if (tree->GetAlias("mus_mc3_id") != 0) {
		mus_mc3_id_branch = tree->GetBranch(tree->GetAlias("mus_mc3_id"));
		mus_mc3_id_branch->SetAddress(&mus_mc3_id_);
	}
	mus_mc3idx_branch = 0;
	if (tree->GetAlias("mus_mc3idx") != 0) {
		mus_mc3idx_branch = tree->GetBranch(tree->GetAlias("mus_mc3idx"));
		mus_mc3idx_branch->SetAddress(&mus_mc3idx_);
	}
	mus_mc3_motherid_branch = 0;
	if (tree->GetAlias("mus_mc3_motherid") != 0) {
		mus_mc3_motherid_branch = tree->GetBranch(tree->GetAlias("mus_mc3_motherid"));
		mus_mc3_motherid_branch->SetAddress(&mus_mc3_motherid_);
	}
	mus_mc3_motheridx_branch = 0;
	if (tree->GetAlias("mus_mc3_motheridx") != 0) {
		mus_mc3_motheridx_branch = tree->GetBranch(tree->GetAlias("mus_mc3_motheridx"));
		mus_mc3_motheridx_branch->SetAddress(&mus_mc3_motheridx_);
	}
	mus_mc_id_branch = 0;
	if (tree->GetAlias("mus_mc_id") != 0) {
		mus_mc_id_branch = tree->GetBranch(tree->GetAlias("mus_mc_id"));
		mus_mc_id_branch->SetAddress(&mus_mc_id_);
	}
	mus_mcidx_branch = 0;
	if (tree->GetAlias("mus_mcidx") != 0) {
		mus_mcidx_branch = tree->GetBranch(tree->GetAlias("mus_mcidx"));
		mus_mcidx_branch->SetAddress(&mus_mcidx_);
	}
	mus_mc_motherid_branch = 0;
	if (tree->GetAlias("mus_mc_motherid") != 0) {
		mus_mc_motherid_branch = tree->GetBranch(tree->GetAlias("mus_mc_motherid"));
		mus_mc_motherid_branch->SetAddress(&mus_mc_motherid_);
	}
	pfjets_mc3_id_branch = 0;
	if (tree->GetAlias("pfjets_mc3_id") != 0) {
		pfjets_mc3_id_branch = tree->GetBranch(tree->GetAlias("pfjets_mc3_id"));
		pfjets_mc3_id_branch->SetAddress(&pfjets_mc3_id_);
	}
	pfjets_mc3idx_branch = 0;
	if (tree->GetAlias("pfjets_mc3idx") != 0) {
		pfjets_mc3idx_branch = tree->GetBranch(tree->GetAlias("pfjets_mc3idx"));
		pfjets_mc3idx_branch->SetAddress(&pfjets_mc3idx_);
	}
	pfjets_mc_gpidx_branch = 0;
	if (tree->GetAlias("pfjets_mc_gpidx") != 0) {
		pfjets_mc_gpidx_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_gpidx"));
		pfjets_mc_gpidx_branch->SetAddress(&pfjets_mc_gpidx_);
	}
	pfjets_mc_id_branch = 0;
	if (tree->GetAlias("pfjets_mc_id") != 0) {
		pfjets_mc_id_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_id"));
		pfjets_mc_id_branch->SetAddress(&pfjets_mc_id_);
	}
	pfjets_mcidx_branch = 0;
	if (tree->GetAlias("pfjets_mcidx") != 0) {
		pfjets_mcidx_branch = tree->GetBranch(tree->GetAlias("pfjets_mcidx"));
		pfjets_mcidx_branch->SetAddress(&pfjets_mcidx_);
	}
	pfjets_mc_motherid_branch = 0;
	if (tree->GetAlias("pfjets_mc_motherid") != 0) {
		pfjets_mc_motherid_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_motherid"));
		pfjets_mc_motherid_branch->SetAddress(&pfjets_mc_motherid_);
	}
	photons_mc3_id_branch = 0;
	if (tree->GetAlias("photons_mc3_id") != 0) {
		photons_mc3_id_branch = tree->GetBranch(tree->GetAlias("photons_mc3_id"));
		photons_mc3_id_branch->SetAddress(&photons_mc3_id_);
	}
	photons_mc3idx_branch = 0;
	if (tree->GetAlias("photons_mc3idx") != 0) {
		photons_mc3idx_branch = tree->GetBranch(tree->GetAlias("photons_mc3idx"));
		photons_mc3idx_branch->SetAddress(&photons_mc3idx_);
	}
	photons_mc3_motherid_branch = 0;
	if (tree->GetAlias("photons_mc3_motherid") != 0) {
		photons_mc3_motherid_branch = tree->GetBranch(tree->GetAlias("photons_mc3_motherid"));
		photons_mc3_motherid_branch->SetAddress(&photons_mc3_motherid_);
	}
	photons_mc3_motheridx_branch = 0;
	if (tree->GetAlias("photons_mc3_motheridx") != 0) {
		photons_mc3_motheridx_branch = tree->GetBranch(tree->GetAlias("photons_mc3_motheridx"));
		photons_mc3_motheridx_branch->SetAddress(&photons_mc3_motheridx_);
	}
	photons_mc_id_branch = 0;
	if (tree->GetAlias("photons_mc_id") != 0) {
		photons_mc_id_branch = tree->GetBranch(tree->GetAlias("photons_mc_id"));
		photons_mc_id_branch->SetAddress(&photons_mc_id_);
	}
	photons_mcidx_branch = 0;
	if (tree->GetAlias("photons_mcidx") != 0) {
		photons_mcidx_branch = tree->GetBranch(tree->GetAlias("photons_mcidx"));
		photons_mcidx_branch->SetAddress(&photons_mcidx_);
	}
	photons_mc_motherid_branch = 0;
	if (tree->GetAlias("photons_mc_motherid") != 0) {
		photons_mc_motherid_branch = tree->GetBranch(tree->GetAlias("photons_mc_motherid"));
		photons_mc_motherid_branch->SetAddress(&photons_mc_motherid_);
	}
	els_category_branch = 0;
	if (tree->GetAlias("els_category") != 0) {
		els_category_branch = tree->GetBranch(tree->GetAlias("els_category"));
		els_category_branch->SetAddress(&els_category_);
	}
	els_charge_branch = 0;
	if (tree->GetAlias("els_charge") != 0) {
		els_charge_branch = tree->GetBranch(tree->GetAlias("els_charge"));
		els_charge_branch->SetAddress(&els_charge_);
	}
	els_ckf_charge_branch = 0;
	if (tree->GetAlias("els_ckf_charge") != 0) {
		els_ckf_charge_branch = tree->GetBranch(tree->GetAlias("els_ckf_charge"));
		els_ckf_charge_branch->SetAddress(&els_ckf_charge_);
	}
	els_ckf_laywithmeas_branch = 0;
	if (tree->GetAlias("els_ckf_laywithmeas") != 0) {
		els_ckf_laywithmeas_branch = tree->GetBranch(tree->GetAlias("els_ckf_laywithmeas"));
		els_ckf_laywithmeas_branch->SetAddress(&els_ckf_laywithmeas_);
	}
	els_class_branch = 0;
	if (tree->GetAlias("els_class") != 0) {
		els_class_branch = tree->GetBranch(tree->GetAlias("els_class"));
		els_class_branch->SetAddress(&els_class_);
	}
	els_conv_delMissHits_branch = 0;
	if (tree->GetAlias("els_conv_delMissHits") != 0) {
		els_conv_delMissHits_branch = tree->GetBranch(tree->GetAlias("els_conv_delMissHits"));
		els_conv_delMissHits_branch->SetAddress(&els_conv_delMissHits_);
	}
	els_conv_flag_branch = 0;
	if (tree->GetAlias("els_conv_flag") != 0) {
		els_conv_flag_branch = tree->GetBranch(tree->GetAlias("els_conv_flag"));
		els_conv_flag_branch->SetAddress(&els_conv_flag_);
	}
	els_conv_gsftkidx_branch = 0;
	if (tree->GetAlias("els_conv_gsftkidx") != 0) {
		els_conv_gsftkidx_branch = tree->GetBranch(tree->GetAlias("els_conv_gsftkidx"));
		els_conv_gsftkidx_branch->SetAddress(&els_conv_gsftkidx_);
	}
	els_conv_old_delMissHits_branch = 0;
	if (tree->GetAlias("els_conv_old_delMissHits") != 0) {
		els_conv_old_delMissHits_branch = tree->GetBranch(tree->GetAlias("els_conv_old_delMissHits"));
		els_conv_old_delMissHits_branch->SetAddress(&els_conv_old_delMissHits_);
	}
	els_conv_old_flag_branch = 0;
	if (tree->GetAlias("els_conv_old_flag") != 0) {
		els_conv_old_flag_branch = tree->GetBranch(tree->GetAlias("els_conv_old_flag"));
		els_conv_old_flag_branch->SetAddress(&els_conv_old_flag_);
	}
	els_conv_old_gsftkidx_branch = 0;
	if (tree->GetAlias("els_conv_old_gsftkidx") != 0) {
		els_conv_old_gsftkidx_branch = tree->GetBranch(tree->GetAlias("els_conv_old_gsftkidx"));
		els_conv_old_gsftkidx_branch->SetAddress(&els_conv_old_gsftkidx_);
	}
	els_conv_old_tkidx_branch = 0;
	if (tree->GetAlias("els_conv_old_tkidx") != 0) {
		els_conv_old_tkidx_branch = tree->GetBranch(tree->GetAlias("els_conv_old_tkidx"));
		els_conv_old_tkidx_branch->SetAddress(&els_conv_old_tkidx_);
	}
	els_conv_tkidx_branch = 0;
	if (tree->GetAlias("els_conv_tkidx") != 0) {
		els_conv_tkidx_branch = tree->GetBranch(tree->GetAlias("els_conv_tkidx"));
		els_conv_tkidx_branch->SetAddress(&els_conv_tkidx_);
	}
	els_exp_innerlayers_branch = 0;
	if (tree->GetAlias("els_exp_innerlayers") != 0) {
		els_exp_innerlayers_branch = tree->GetBranch(tree->GetAlias("els_exp_innerlayers"));
		els_exp_innerlayers_branch->SetAddress(&els_exp_innerlayers_);
	}
	els_exp_outerlayers_branch = 0;
	if (tree->GetAlias("els_exp_outerlayers") != 0) {
		els_exp_outerlayers_branch = tree->GetBranch(tree->GetAlias("els_exp_outerlayers"));
		els_exp_outerlayers_branch->SetAddress(&els_exp_outerlayers_);
	}
	els_fiduciality_branch = 0;
	if (tree->GetAlias("els_fiduciality") != 0) {
		els_fiduciality_branch = tree->GetBranch(tree->GetAlias("els_fiduciality"));
		els_fiduciality_branch->SetAddress(&els_fiduciality_);
	}
	els_gsftrkidx_branch = 0;
	if (tree->GetAlias("els_gsftrkidx") != 0) {
		els_gsftrkidx_branch = tree->GetBranch(tree->GetAlias("els_gsftrkidx"));
		els_gsftrkidx_branch->SetAddress(&els_gsftrkidx_);
	}
	els_layer1_det_branch = 0;
	if (tree->GetAlias("els_layer1_det") != 0) {
		els_layer1_det_branch = tree->GetBranch(tree->GetAlias("els_layer1_det"));
		els_layer1_det_branch->SetAddress(&els_layer1_det_);
	}
	els_layer1_layer_branch = 0;
	if (tree->GetAlias("els_layer1_layer") != 0) {
		els_layer1_layer_branch = tree->GetBranch(tree->GetAlias("els_layer1_layer"));
		els_layer1_layer_branch->SetAddress(&els_layer1_layer_);
	}
	els_layer1_sizerphi_branch = 0;
	if (tree->GetAlias("els_layer1_sizerphi") != 0) {
		els_layer1_sizerphi_branch = tree->GetBranch(tree->GetAlias("els_layer1_sizerphi"));
		els_layer1_sizerphi_branch->SetAddress(&els_layer1_sizerphi_);
	}
	els_layer1_sizerz_branch = 0;
	if (tree->GetAlias("els_layer1_sizerz") != 0) {
		els_layer1_sizerz_branch = tree->GetBranch(tree->GetAlias("els_layer1_sizerz"));
		els_layer1_sizerz_branch->SetAddress(&els_layer1_sizerz_);
	}
	els_lostHits_branch = 0;
	if (tree->GetAlias("els_lostHits") != 0) {
		els_lostHits_branch = tree->GetBranch(tree->GetAlias("els_lostHits"));
		els_lostHits_branch->SetAddress(&els_lostHits_);
	}
	els_lost_pixelhits_branch = 0;
	if (tree->GetAlias("els_lost_pixelhits") != 0) {
		els_lost_pixelhits_branch = tree->GetBranch(tree->GetAlias("els_lost_pixelhits"));
		els_lost_pixelhits_branch->SetAddress(&els_lost_pixelhits_);
	}
	els_mc_patMatch_id_branch = 0;
	if (tree->GetAlias("els_mc_patMatch_id") != 0) {
		els_mc_patMatch_id_branch = tree->GetBranch(tree->GetAlias("els_mc_patMatch_id"));
		els_mc_patMatch_id_branch->SetAddress(&els_mc_patMatch_id_);
	}
	els_nSeed_branch = 0;
	if (tree->GetAlias("els_nSeed") != 0) {
		els_nSeed_branch = tree->GetBranch(tree->GetAlias("els_nSeed"));
		els_nSeed_branch->SetAddress(&els_nSeed_);
	}
	els_nlayers_branch = 0;
	if (tree->GetAlias("els_nlayers") != 0) {
		els_nlayers_branch = tree->GetBranch(tree->GetAlias("els_nlayers"));
		els_nlayers_branch->SetAddress(&els_nlayers_);
	}
	els_nlayers3D_branch = 0;
	if (tree->GetAlias("els_nlayers3D") != 0) {
		els_nlayers3D_branch = tree->GetBranch(tree->GetAlias("els_nlayers3D"));
		els_nlayers3D_branch->SetAddress(&els_nlayers3D_);
	}
	els_nlayersLost_branch = 0;
	if (tree->GetAlias("els_nlayersLost") != 0) {
		els_nlayersLost_branch = tree->GetBranch(tree->GetAlias("els_nlayersLost"));
		els_nlayersLost_branch->SetAddress(&els_nlayersLost_);
	}
	els_sccharge_branch = 0;
	if (tree->GetAlias("els_sccharge") != 0) {
		els_sccharge_branch = tree->GetBranch(tree->GetAlias("els_sccharge"));
		els_sccharge_branch->SetAddress(&els_sccharge_);
	}
	els_scindex_branch = 0;
	if (tree->GetAlias("els_scindex") != 0) {
		els_scindex_branch = tree->GetBranch(tree->GetAlias("els_scindex"));
		els_scindex_branch->SetAddress(&els_scindex_);
	}
	els_trk_charge_branch = 0;
	if (tree->GetAlias("els_trk_charge") != 0) {
		els_trk_charge_branch = tree->GetBranch(tree->GetAlias("els_trk_charge"));
		els_trk_charge_branch->SetAddress(&els_trk_charge_);
	}
	els_trkidx_branch = 0;
	if (tree->GetAlias("els_trkidx") != 0) {
		els_trkidx_branch = tree->GetBranch(tree->GetAlias("els_trkidx"));
		els_trkidx_branch->SetAddress(&els_trkidx_);
	}
	els_type_branch = 0;
	if (tree->GetAlias("els_type") != 0) {
		els_type_branch = tree->GetBranch(tree->GetAlias("els_type"));
		els_type_branch->SetAddress(&els_type_);
	}
	els_validHits_branch = 0;
	if (tree->GetAlias("els_validHits") != 0) {
		els_validHits_branch = tree->GetBranch(tree->GetAlias("els_validHits"));
		els_validHits_branch->SetAddress(&els_validHits_);
	}
	els_valid_pixelhits_branch = 0;
	if (tree->GetAlias("els_valid_pixelhits") != 0) {
		els_valid_pixelhits_branch = tree->GetBranch(tree->GetAlias("els_valid_pixelhits"));
		els_valid_pixelhits_branch->SetAddress(&els_valid_pixelhits_);
	}
	genps_id_branch = 0;
	if (tree->GetAlias("genps_id") != 0) {
		genps_id_branch = tree->GetBranch(tree->GetAlias("genps_id"));
		genps_id_branch->SetAddress(&genps_id_);
	}
	genps_id_mother_branch = 0;
	if (tree->GetAlias("genps_id_mother") != 0) {
		genps_id_mother_branch = tree->GetBranch(tree->GetAlias("genps_id_mother"));
		genps_id_mother_branch->SetAddress(&genps_id_mother_);
	}
	genps_status_branch = 0;
	if (tree->GetAlias("genps_status") != 0) {
		genps_status_branch = tree->GetBranch(tree->GetAlias("genps_status"));
		genps_status_branch->SetAddress(&genps_status_);
	}
	mus_algo_branch = 0;
	if (tree->GetAlias("mus_algo") != 0) {
		mus_algo_branch = tree->GetBranch(tree->GetAlias("mus_algo"));
		mus_algo_branch->SetAddress(&mus_algo_);
	}
	mus_charge_branch = 0;
	if (tree->GetAlias("mus_charge") != 0) {
		mus_charge_branch = tree->GetBranch(tree->GetAlias("mus_charge"));
		mus_charge_branch->SetAddress(&mus_charge_);
	}
	mus_ecal_rawId_branch = 0;
	if (tree->GetAlias("mus_ecal_rawId") != 0) {
		mus_ecal_rawId_branch = tree->GetBranch(tree->GetAlias("mus_ecal_rawId"));
		mus_ecal_rawId_branch->SetAddress(&mus_ecal_rawId_);
	}
	mus_exp_innerlayers_branch = 0;
	if (tree->GetAlias("mus_exp_innerlayers") != 0) {
		mus_exp_innerlayers_branch = tree->GetBranch(tree->GetAlias("mus_exp_innerlayers"));
		mus_exp_innerlayers_branch->SetAddress(&mus_exp_innerlayers_);
	}
	mus_exp_outerlayers_branch = 0;
	if (tree->GetAlias("mus_exp_outerlayers") != 0) {
		mus_exp_outerlayers_branch = tree->GetBranch(tree->GetAlias("mus_exp_outerlayers"));
		mus_exp_outerlayers_branch->SetAddress(&mus_exp_outerlayers_);
	}
	mus_gfit_algo_branch = 0;
	if (tree->GetAlias("mus_gfit_algo") != 0) {
		mus_gfit_algo_branch = tree->GetBranch(tree->GetAlias("mus_gfit_algo"));
		mus_gfit_algo_branch->SetAddress(&mus_gfit_algo_);
	}
	mus_gfit_etaErr_branch = 0;
	if (tree->GetAlias("mus_gfit_etaErr") != 0) {
		mus_gfit_etaErr_branch = tree->GetBranch(tree->GetAlias("mus_gfit_etaErr"));
		mus_gfit_etaErr_branch->SetAddress(&mus_gfit_etaErr_);
	}
	mus_gfit_exp_innerlayer_branch = 0;
	if (tree->GetAlias("mus_gfit_exp_innerlayer") != 0) {
		mus_gfit_exp_innerlayer_branch = tree->GetBranch(tree->GetAlias("mus_gfit_exp_innerlayer"));
		mus_gfit_exp_innerlayer_branch->SetAddress(&mus_gfit_exp_innerlayer_);
	}
	mus_gfit_exp_outerlayer_branch = 0;
	if (tree->GetAlias("mus_gfit_exp_outerlayer") != 0) {
		mus_gfit_exp_outerlayer_branch = tree->GetBranch(tree->GetAlias("mus_gfit_exp_outerlayer"));
		mus_gfit_exp_outerlayer_branch->SetAddress(&mus_gfit_exp_outerlayer_);
	}
	mus_gfit_lostHits_branch = 0;
	if (tree->GetAlias("mus_gfit_lostHits") != 0) {
		mus_gfit_lostHits_branch = tree->GetBranch(tree->GetAlias("mus_gfit_lostHits"));
		mus_gfit_lostHits_branch->SetAddress(&mus_gfit_lostHits_);
	}
	mus_gfit_lostPixelHits_branch = 0;
	if (tree->GetAlias("mus_gfit_lostPixelHits") != 0) {
		mus_gfit_lostPixelHits_branch = tree->GetBranch(tree->GetAlias("mus_gfit_lostPixelHits"));
		mus_gfit_lostPixelHits_branch->SetAddress(&mus_gfit_lostPixelHits_);
	}
	mus_gfit_nlayers_branch = 0;
	if (tree->GetAlias("mus_gfit_nlayers") != 0) {
		mus_gfit_nlayers_branch = tree->GetBranch(tree->GetAlias("mus_gfit_nlayers"));
		mus_gfit_nlayers_branch->SetAddress(&mus_gfit_nlayers_);
	}
	mus_gfit_nlayers3D_branch = 0;
	if (tree->GetAlias("mus_gfit_nlayers3D") != 0) {
		mus_gfit_nlayers3D_branch = tree->GetBranch(tree->GetAlias("mus_gfit_nlayers3D"));
		mus_gfit_nlayers3D_branch->SetAddress(&mus_gfit_nlayers3D_);
	}
	mus_gfit_nlayersLost_branch = 0;
	if (tree->GetAlias("mus_gfit_nlayersLost") != 0) {
		mus_gfit_nlayersLost_branch = tree->GetBranch(tree->GetAlias("mus_gfit_nlayersLost"));
		mus_gfit_nlayersLost_branch->SetAddress(&mus_gfit_nlayersLost_);
	}
	mus_gfit_phiErr_branch = 0;
	if (tree->GetAlias("mus_gfit_phiErr") != 0) {
		mus_gfit_phiErr_branch = tree->GetBranch(tree->GetAlias("mus_gfit_phiErr"));
		mus_gfit_phiErr_branch->SetAddress(&mus_gfit_phiErr_);
	}
	mus_gfit_ptErr_branch = 0;
	if (tree->GetAlias("mus_gfit_ptErr") != 0) {
		mus_gfit_ptErr_branch = tree->GetBranch(tree->GetAlias("mus_gfit_ptErr"));
		mus_gfit_ptErr_branch->SetAddress(&mus_gfit_ptErr_);
	}
	mus_gfit_qualityMask_branch = 0;
	if (tree->GetAlias("mus_gfit_qualityMask") != 0) {
		mus_gfit_qualityMask_branch = tree->GetBranch(tree->GetAlias("mus_gfit_qualityMask"));
		mus_gfit_qualityMask_branch->SetAddress(&mus_gfit_qualityMask_);
	}
	mus_gfit_trk_charge_branch = 0;
	if (tree->GetAlias("mus_gfit_trk_charge") != 0) {
		mus_gfit_trk_charge_branch = tree->GetBranch(tree->GetAlias("mus_gfit_trk_charge"));
		mus_gfit_trk_charge_branch->SetAddress(&mus_gfit_trk_charge_);
	}
	mus_gfit_validHits_branch = 0;
	if (tree->GetAlias("mus_gfit_validHits") != 0) {
		mus_gfit_validHits_branch = tree->GetBranch(tree->GetAlias("mus_gfit_validHits"));
		mus_gfit_validHits_branch->SetAddress(&mus_gfit_validHits_);
	}
	mus_gfit_validPixelHits_branch = 0;
	if (tree->GetAlias("mus_gfit_validPixelHits") != 0) {
		mus_gfit_validPixelHits_branch = tree->GetBranch(tree->GetAlias("mus_gfit_validPixelHits"));
		mus_gfit_validPixelHits_branch->SetAddress(&mus_gfit_validPixelHits_);
	}
	mus_gfit_validSTAHits_branch = 0;
	if (tree->GetAlias("mus_gfit_validSTAHits") != 0) {
		mus_gfit_validSTAHits_branch = tree->GetBranch(tree->GetAlias("mus_gfit_validSTAHits"));
		mus_gfit_validSTAHits_branch->SetAddress(&mus_gfit_validSTAHits_);
	}
	mus_gfit_validSiHits_branch = 0;
	if (tree->GetAlias("mus_gfit_validSiHits") != 0) {
		mus_gfit_validSiHits_branch = tree->GetBranch(tree->GetAlias("mus_gfit_validSiHits"));
		mus_gfit_validSiHits_branch->SetAddress(&mus_gfit_validSiHits_);
	}
	mus_goodmask_branch = 0;
	if (tree->GetAlias("mus_goodmask") != 0) {
		mus_goodmask_branch = tree->GetBranch(tree->GetAlias("mus_goodmask"));
		mus_goodmask_branch->SetAddress(&mus_goodmask_);
	}
	mus_hcal_rawId_branch = 0;
	if (tree->GetAlias("mus_hcal_rawId") != 0) {
		mus_hcal_rawId_branch = tree->GetBranch(tree->GetAlias("mus_hcal_rawId"));
		mus_hcal_rawId_branch->SetAddress(&mus_hcal_rawId_);
	}
	mus_iso03_njets_branch = 0;
	if (tree->GetAlias("mus_iso03_njets") != 0) {
		mus_iso03_njets_branch = tree->GetBranch(tree->GetAlias("mus_iso03_njets"));
		mus_iso03_njets_branch->SetAddress(&mus_iso03_njets_);
	}
	mus_iso03_ntrk_branch = 0;
	if (tree->GetAlias("mus_iso03_ntrk") != 0) {
		mus_iso03_ntrk_branch = tree->GetBranch(tree->GetAlias("mus_iso03_ntrk"));
		mus_iso03_ntrk_branch->SetAddress(&mus_iso03_ntrk_);
	}
	mus_iso05_ntrk_branch = 0;
	if (tree->GetAlias("mus_iso05_ntrk") != 0) {
		mus_iso05_ntrk_branch = tree->GetBranch(tree->GetAlias("mus_iso05_ntrk"));
		mus_iso05_ntrk_branch->SetAddress(&mus_iso05_ntrk_);
	}
	mus_lostHits_branch = 0;
	if (tree->GetAlias("mus_lostHits") != 0) {
		mus_lostHits_branch = tree->GetBranch(tree->GetAlias("mus_lostHits"));
		mus_lostHits_branch->SetAddress(&mus_lostHits_);
	}
	mus_lostPixelHits_branch = 0;
	if (tree->GetAlias("mus_lostPixelHits") != 0) {
		mus_lostPixelHits_branch = tree->GetBranch(tree->GetAlias("mus_lostPixelHits"));
		mus_lostPixelHits_branch->SetAddress(&mus_lostPixelHits_);
	}
	mus_mc_patMatch_id_branch = 0;
	if (tree->GetAlias("mus_mc_patMatch_id") != 0) {
		mus_mc_patMatch_id_branch = tree->GetBranch(tree->GetAlias("mus_mc_patMatch_id"));
		mus_mc_patMatch_id_branch->SetAddress(&mus_mc_patMatch_id_);
	}
	mus_muonBestTrackType_branch = 0;
	if (tree->GetAlias("mus_muonBestTrackType") != 0) {
		mus_muonBestTrackType_branch = tree->GetBranch(tree->GetAlias("mus_muonBestTrackType"));
		mus_muonBestTrackType_branch->SetAddress(&mus_muonBestTrackType_);
	}
	mus_nOverlaps_branch = 0;
	if (tree->GetAlias("mus_nOverlaps") != 0) {
		mus_nOverlaps_branch = tree->GetBranch(tree->GetAlias("mus_nOverlaps"));
		mus_nOverlaps_branch->SetAddress(&mus_nOverlaps_);
	}
	mus_nlayers_branch = 0;
	if (tree->GetAlias("mus_nlayers") != 0) {
		mus_nlayers_branch = tree->GetBranch(tree->GetAlias("mus_nlayers"));
		mus_nlayers_branch->SetAddress(&mus_nlayers_);
	}
	mus_nlayers3D_branch = 0;
	if (tree->GetAlias("mus_nlayers3D") != 0) {
		mus_nlayers3D_branch = tree->GetBranch(tree->GetAlias("mus_nlayers3D"));
		mus_nlayers3D_branch->SetAddress(&mus_nlayers3D_);
	}
	mus_nlayersLost_branch = 0;
	if (tree->GetAlias("mus_nlayersLost") != 0) {
		mus_nlayersLost_branch = tree->GetBranch(tree->GetAlias("mus_nlayersLost"));
		mus_nlayersLost_branch->SetAddress(&mus_nlayersLost_);
	}
	mus_nmatches_branch = 0;
	if (tree->GetAlias("mus_nmatches") != 0) {
		mus_nmatches_branch = tree->GetBranch(tree->GetAlias("mus_nmatches"));
		mus_nmatches_branch->SetAddress(&mus_nmatches_);
	}
	mus_numberOfMatchedStations_branch = 0;
	if (tree->GetAlias("mus_numberOfMatchedStations") != 0) {
		mus_numberOfMatchedStations_branch = tree->GetBranch(tree->GetAlias("mus_numberOfMatchedStations"));
		mus_numberOfMatchedStations_branch->SetAddress(&mus_numberOfMatchedStations_);
	}
	mus_overlap0_branch = 0;
	if (tree->GetAlias("mus_overlap0") != 0) {
		mus_overlap0_branch = tree->GetBranch(tree->GetAlias("mus_overlap0"));
		mus_overlap0_branch->SetAddress(&mus_overlap0_);
	}
	mus_overlap1_branch = 0;
	if (tree->GetAlias("mus_overlap1") != 0) {
		mus_overlap1_branch = tree->GetBranch(tree->GetAlias("mus_overlap1"));
		mus_overlap1_branch->SetAddress(&mus_overlap1_);
	}
	mus_pfcharge_branch = 0;
	if (tree->GetAlias("mus_pfcharge") != 0) {
		mus_pfcharge_branch = tree->GetBranch(tree->GetAlias("mus_pfcharge"));
		mus_pfcharge_branch->SetAddress(&mus_pfcharge_);
	}
	mus_pfparticleId_branch = 0;
	if (tree->GetAlias("mus_pfparticleId") != 0) {
		mus_pfparticleId_branch = tree->GetBranch(tree->GetAlias("mus_pfparticleId"));
		mus_pfparticleId_branch->SetAddress(&mus_pfparticleId_);
	}
	mus_pid_PFMuon_branch = 0;
	if (tree->GetAlias("mus_pid_PFMuon") != 0) {
		mus_pid_PFMuon_branch = tree->GetBranch(tree->GetAlias("mus_pid_PFMuon"));
		mus_pid_PFMuon_branch->SetAddress(&mus_pid_PFMuon_);
	}
	mus_pid_TM2DCompatibilityLoose_branch = 0;
	if (tree->GetAlias("mus_pid_TM2DCompatibilityLoose") != 0) {
		mus_pid_TM2DCompatibilityLoose_branch = tree->GetBranch(tree->GetAlias("mus_pid_TM2DCompatibilityLoose"));
		mus_pid_TM2DCompatibilityLoose_branch->SetAddress(&mus_pid_TM2DCompatibilityLoose_);
	}
	mus_pid_TM2DCompatibilityTight_branch = 0;
	if (tree->GetAlias("mus_pid_TM2DCompatibilityTight") != 0) {
		mus_pid_TM2DCompatibilityTight_branch = tree->GetBranch(tree->GetAlias("mus_pid_TM2DCompatibilityTight"));
		mus_pid_TM2DCompatibilityTight_branch->SetAddress(&mus_pid_TM2DCompatibilityTight_);
	}
	mus_pid_TMLastStationLoose_branch = 0;
	if (tree->GetAlias("mus_pid_TMLastStationLoose") != 0) {
		mus_pid_TMLastStationLoose_branch = tree->GetBranch(tree->GetAlias("mus_pid_TMLastStationLoose"));
		mus_pid_TMLastStationLoose_branch->SetAddress(&mus_pid_TMLastStationLoose_);
	}
	mus_pid_TMLastStationTight_branch = 0;
	if (tree->GetAlias("mus_pid_TMLastStationTight") != 0) {
		mus_pid_TMLastStationTight_branch = tree->GetBranch(tree->GetAlias("mus_pid_TMLastStationTight"));
		mus_pid_TMLastStationTight_branch->SetAddress(&mus_pid_TMLastStationTight_);
	}
	mus_qualityMask_branch = 0;
	if (tree->GetAlias("mus_qualityMask") != 0) {
		mus_qualityMask_branch = tree->GetBranch(tree->GetAlias("mus_qualityMask"));
		mus_qualityMask_branch->SetAddress(&mus_qualityMask_);
	}
	mus_sta_algo_branch = 0;
	if (tree->GetAlias("mus_sta_algo") != 0) {
		mus_sta_algo_branch = tree->GetBranch(tree->GetAlias("mus_sta_algo"));
		mus_sta_algo_branch->SetAddress(&mus_sta_algo_);
	}
	mus_sta_etaErr_branch = 0;
	if (tree->GetAlias("mus_sta_etaErr") != 0) {
		mus_sta_etaErr_branch = tree->GetBranch(tree->GetAlias("mus_sta_etaErr"));
		mus_sta_etaErr_branch->SetAddress(&mus_sta_etaErr_);
	}
	mus_sta_exp_innerlayer_branch = 0;
	if (tree->GetAlias("mus_sta_exp_innerlayer") != 0) {
		mus_sta_exp_innerlayer_branch = tree->GetBranch(tree->GetAlias("mus_sta_exp_innerlayer"));
		mus_sta_exp_innerlayer_branch->SetAddress(&mus_sta_exp_innerlayer_);
	}
	mus_sta_exp_outerlayer_branch = 0;
	if (tree->GetAlias("mus_sta_exp_outerlayer") != 0) {
		mus_sta_exp_outerlayer_branch = tree->GetBranch(tree->GetAlias("mus_sta_exp_outerlayer"));
		mus_sta_exp_outerlayer_branch->SetAddress(&mus_sta_exp_outerlayer_);
	}
	mus_sta_lostHits_branch = 0;
	if (tree->GetAlias("mus_sta_lostHits") != 0) {
		mus_sta_lostHits_branch = tree->GetBranch(tree->GetAlias("mus_sta_lostHits"));
		mus_sta_lostHits_branch->SetAddress(&mus_sta_lostHits_);
	}
	mus_sta_lostPixelHits_branch = 0;
	if (tree->GetAlias("mus_sta_lostPixelHits") != 0) {
		mus_sta_lostPixelHits_branch = tree->GetBranch(tree->GetAlias("mus_sta_lostPixelHits"));
		mus_sta_lostPixelHits_branch->SetAddress(&mus_sta_lostPixelHits_);
	}
	mus_sta_nlayers_branch = 0;
	if (tree->GetAlias("mus_sta_nlayers") != 0) {
		mus_sta_nlayers_branch = tree->GetBranch(tree->GetAlias("mus_sta_nlayers"));
		mus_sta_nlayers_branch->SetAddress(&mus_sta_nlayers_);
	}
	mus_sta_nlayers3D_branch = 0;
	if (tree->GetAlias("mus_sta_nlayers3D") != 0) {
		mus_sta_nlayers3D_branch = tree->GetBranch(tree->GetAlias("mus_sta_nlayers3D"));
		mus_sta_nlayers3D_branch->SetAddress(&mus_sta_nlayers3D_);
	}
	mus_sta_nlayersLost_branch = 0;
	if (tree->GetAlias("mus_sta_nlayersLost") != 0) {
		mus_sta_nlayersLost_branch = tree->GetBranch(tree->GetAlias("mus_sta_nlayersLost"));
		mus_sta_nlayersLost_branch->SetAddress(&mus_sta_nlayersLost_);
	}
	mus_sta_phiErr_branch = 0;
	if (tree->GetAlias("mus_sta_phiErr") != 0) {
		mus_sta_phiErr_branch = tree->GetBranch(tree->GetAlias("mus_sta_phiErr"));
		mus_sta_phiErr_branch->SetAddress(&mus_sta_phiErr_);
	}
	mus_sta_ptErr_branch = 0;
	if (tree->GetAlias("mus_sta_ptErr") != 0) {
		mus_sta_ptErr_branch = tree->GetBranch(tree->GetAlias("mus_sta_ptErr"));
		mus_sta_ptErr_branch->SetAddress(&mus_sta_ptErr_);
	}
	mus_sta_qualityMask_branch = 0;
	if (tree->GetAlias("mus_sta_qualityMask") != 0) {
		mus_sta_qualityMask_branch = tree->GetBranch(tree->GetAlias("mus_sta_qualityMask"));
		mus_sta_qualityMask_branch->SetAddress(&mus_sta_qualityMask_);
	}
	mus_sta_trk_charge_branch = 0;
	if (tree->GetAlias("mus_sta_trk_charge") != 0) {
		mus_sta_trk_charge_branch = tree->GetBranch(tree->GetAlias("mus_sta_trk_charge"));
		mus_sta_trk_charge_branch->SetAddress(&mus_sta_trk_charge_);
	}
	mus_sta_validHits_branch = 0;
	if (tree->GetAlias("mus_sta_validHits") != 0) {
		mus_sta_validHits_branch = tree->GetBranch(tree->GetAlias("mus_sta_validHits"));
		mus_sta_validHits_branch->SetAddress(&mus_sta_validHits_);
	}
	mus_sta_validPixelHits_branch = 0;
	if (tree->GetAlias("mus_sta_validPixelHits") != 0) {
		mus_sta_validPixelHits_branch = tree->GetBranch(tree->GetAlias("mus_sta_validPixelHits"));
		mus_sta_validPixelHits_branch->SetAddress(&mus_sta_validPixelHits_);
	}
	mus_timeDirection_branch = 0;
	if (tree->GetAlias("mus_timeDirection") != 0) {
		mus_timeDirection_branch = tree->GetBranch(tree->GetAlias("mus_timeDirection"));
		mus_timeDirection_branch->SetAddress(&mus_timeDirection_);
	}
	mus_timeNumStationsUsed_branch = 0;
	if (tree->GetAlias("mus_timeNumStationsUsed") != 0) {
		mus_timeNumStationsUsed_branch = tree->GetBranch(tree->GetAlias("mus_timeNumStationsUsed"));
		mus_timeNumStationsUsed_branch->SetAddress(&mus_timeNumStationsUsed_);
	}
	mus_trk_charge_branch = 0;
	if (tree->GetAlias("mus_trk_charge") != 0) {
		mus_trk_charge_branch = tree->GetBranch(tree->GetAlias("mus_trk_charge"));
		mus_trk_charge_branch->SetAddress(&mus_trk_charge_);
	}
	mus_trkidx_branch = 0;
	if (tree->GetAlias("mus_trkidx") != 0) {
		mus_trkidx_branch = tree->GetBranch(tree->GetAlias("mus_trkidx"));
		mus_trkidx_branch->SetAddress(&mus_trkidx_);
	}
	mus_type_branch = 0;
	if (tree->GetAlias("mus_type") != 0) {
		mus_type_branch = tree->GetBranch(tree->GetAlias("mus_type"));
		mus_type_branch->SetAddress(&mus_type_);
	}
	mus_validHits_branch = 0;
	if (tree->GetAlias("mus_validHits") != 0) {
		mus_validHits_branch = tree->GetBranch(tree->GetAlias("mus_validHits"));
		mus_validHits_branch->SetAddress(&mus_validHits_);
	}
	mus_validPixelHits_branch = 0;
	if (tree->GetAlias("mus_validPixelHits") != 0) {
		mus_validPixelHits_branch = tree->GetBranch(tree->GetAlias("mus_validPixelHits"));
		mus_validPixelHits_branch->SetAddress(&mus_validPixelHits_);
	}
	pfcands_charge_branch = 0;
	if (tree->GetAlias("pfcands_charge") != 0) {
		pfcands_charge_branch = tree->GetBranch(tree->GetAlias("pfcands_charge"));
		pfcands_charge_branch->SetAddress(&pfcands_charge_);
	}
	pfcands_particleId_branch = 0;
	if (tree->GetAlias("pfcands_particleId") != 0) {
		pfcands_particleId_branch = tree->GetBranch(tree->GetAlias("pfcands_particleId"));
		pfcands_particleId_branch->SetAddress(&pfcands_particleId_);
	}
	pfjets_chargedHadronMultiplicity_branch = 0;
	if (tree->GetAlias("pfjets_chargedHadronMultiplicity") != 0) {
		pfjets_chargedHadronMultiplicity_branch = tree->GetBranch(tree->GetAlias("pfjets_chargedHadronMultiplicity"));
		pfjets_chargedHadronMultiplicity_branch->SetAddress(&pfjets_chargedHadronMultiplicity_);
	}
	pfjets_chargedMultiplicity_branch = 0;
	if (tree->GetAlias("pfjets_chargedMultiplicity") != 0) {
		pfjets_chargedMultiplicity_branch = tree->GetBranch(tree->GetAlias("pfjets_chargedMultiplicity"));
		pfjets_chargedMultiplicity_branch->SetAddress(&pfjets_chargedMultiplicity_);
	}
	pfjets_electronMultiplicity_branch = 0;
	if (tree->GetAlias("pfjets_electronMultiplicity") != 0) {
		pfjets_electronMultiplicity_branch = tree->GetBranch(tree->GetAlias("pfjets_electronMultiplicity"));
		pfjets_electronMultiplicity_branch->SetAddress(&pfjets_electronMultiplicity_);
	}
	pfjets_hfEmMultiplicity_branch = 0;
	if (tree->GetAlias("pfjets_hfEmMultiplicity") != 0) {
		pfjets_hfEmMultiplicity_branch = tree->GetBranch(tree->GetAlias("pfjets_hfEmMultiplicity"));
		pfjets_hfEmMultiplicity_branch->SetAddress(&pfjets_hfEmMultiplicity_);
	}
	pfjets_hfHadronMultiplicity_branch = 0;
	if (tree->GetAlias("pfjets_hfHadronMultiplicity") != 0) {
		pfjets_hfHadronMultiplicity_branch = tree->GetBranch(tree->GetAlias("pfjets_hfHadronMultiplicity"));
		pfjets_hfHadronMultiplicity_branch->SetAddress(&pfjets_hfHadronMultiplicity_);
	}
	pfjets_muonMultiplicity_branch = 0;
	if (tree->GetAlias("pfjets_muonMultiplicity") != 0) {
		pfjets_muonMultiplicity_branch = tree->GetBranch(tree->GetAlias("pfjets_muonMultiplicity"));
		pfjets_muonMultiplicity_branch->SetAddress(&pfjets_muonMultiplicity_);
	}
	pfjets_neutralHadronMultiplicity_branch = 0;
	if (tree->GetAlias("pfjets_neutralHadronMultiplicity") != 0) {
		pfjets_neutralHadronMultiplicity_branch = tree->GetBranch(tree->GetAlias("pfjets_neutralHadronMultiplicity"));
		pfjets_neutralHadronMultiplicity_branch->SetAddress(&pfjets_neutralHadronMultiplicity_);
	}
	pfjets_neutralMultiplicity_branch = 0;
	if (tree->GetAlias("pfjets_neutralMultiplicity") != 0) {
		pfjets_neutralMultiplicity_branch = tree->GetBranch(tree->GetAlias("pfjets_neutralMultiplicity"));
		pfjets_neutralMultiplicity_branch->SetAddress(&pfjets_neutralMultiplicity_);
	}
	pfjets_partonFlavour_branch = 0;
	if (tree->GetAlias("pfjets_partonFlavour") != 0) {
		pfjets_partonFlavour_branch = tree->GetBranch(tree->GetAlias("pfjets_partonFlavour"));
		pfjets_partonFlavour_branch->SetAddress(&pfjets_partonFlavour_);
	}
	pfjets_photonMultiplicity_branch = 0;
	if (tree->GetAlias("pfjets_photonMultiplicity") != 0) {
		pfjets_photonMultiplicity_branch = tree->GetBranch(tree->GetAlias("pfjets_photonMultiplicity"));
		pfjets_photonMultiplicity_branch->SetAddress(&pfjets_photonMultiplicity_);
	}
	taus_pf_charge_branch = 0;
	if (tree->GetAlias("taus_pf_charge") != 0) {
		taus_pf_charge_branch = tree->GetBranch(tree->GetAlias("taus_pf_charge"));
		taus_pf_charge_branch->SetAddress(&taus_pf_charge_);
	}
	photons_fiduciality_branch = 0;
	if (tree->GetAlias("photons_fiduciality") != 0) {
		photons_fiduciality_branch = tree->GetBranch(tree->GetAlias("photons_fiduciality"));
		photons_fiduciality_branch->SetAddress(&photons_fiduciality_);
	}
	puInfo_bunchCrossing_branch = 0;
	if (tree->GetAlias("puInfo_bunchCrossing") != 0) {
		puInfo_bunchCrossing_branch = tree->GetBranch(tree->GetAlias("puInfo_bunchCrossing"));
		puInfo_bunchCrossing_branch->SetAddress(&puInfo_bunchCrossing_);
	}
	puInfo_nPUvertices_branch = 0;
	if (tree->GetAlias("puInfo_nPUvertices") != 0) {
		puInfo_nPUvertices_branch = tree->GetBranch(tree->GetAlias("puInfo_nPUvertices"));
		puInfo_nPUvertices_branch->SetAddress(&puInfo_nPUvertices_);
	}
	convs_algo_branch = 0;
	if (tree->GetAlias("convs_algo") != 0) {
		convs_algo_branch = tree->GetBranch(tree->GetAlias("convs_algo"));
		convs_algo_branch->SetAddress(&convs_algo_);
	}
	convs_isConverted_branch = 0;
	if (tree->GetAlias("convs_isConverted") != 0) {
		convs_isConverted_branch = tree->GetBranch(tree->GetAlias("convs_isConverted"));
		convs_isConverted_branch->SetAddress(&convs_isConverted_);
	}
	convs_quality_branch = 0;
	if (tree->GetAlias("convs_quality") != 0) {
		convs_quality_branch = tree->GetBranch(tree->GetAlias("convs_quality"));
		convs_quality_branch->SetAddress(&convs_quality_);
	}
	vtxs_isFake_branch = 0;
	if (tree->GetAlias("vtxs_isFake") != 0) {
		vtxs_isFake_branch = tree->GetBranch(tree->GetAlias("vtxs_isFake"));
		vtxs_isFake_branch->SetAddress(&vtxs_isFake_);
	}
	vtxs_isValid_branch = 0;
	if (tree->GetAlias("vtxs_isValid") != 0) {
		vtxs_isValid_branch = tree->GetBranch(tree->GetAlias("vtxs_isValid"));
		vtxs_isValid_branch->SetAddress(&vtxs_isValid_);
	}
	vtxs_tracksSize_branch = 0;
	if (tree->GetAlias("vtxs_tracksSize") != 0) {
		vtxs_tracksSize_branch = tree->GetBranch(tree->GetAlias("vtxs_tracksSize"));
		vtxs_tracksSize_branch->SetAddress(&vtxs_tracksSize_);
	}
	els_convs_delMissHits_branch = 0;
	if (tree->GetAlias("els_convs_delMissHits") != 0) {
		els_convs_delMissHits_branch = tree->GetBranch(tree->GetAlias("els_convs_delMissHits"));
		els_convs_delMissHits_branch->SetAddress(&els_convs_delMissHits_);
	}
	els_convs_flag_branch = 0;
	if (tree->GetAlias("els_convs_flag") != 0) {
		els_convs_flag_branch = tree->GetBranch(tree->GetAlias("els_convs_flag"));
		els_convs_flag_branch->SetAddress(&els_convs_flag_);
	}
	els_convs_gsftkidx_branch = 0;
	if (tree->GetAlias("els_convs_gsftkidx") != 0) {
		els_convs_gsftkidx_branch = tree->GetBranch(tree->GetAlias("els_convs_gsftkidx"));
		els_convs_gsftkidx_branch->SetAddress(&els_convs_gsftkidx_);
	}
	els_convs_tkidx_branch = 0;
	if (tree->GetAlias("els_convs_tkidx") != 0) {
		els_convs_tkidx_branch = tree->GetBranch(tree->GetAlias("els_convs_tkidx"));
		els_convs_tkidx_branch->SetAddress(&els_convs_tkidx_);
	}
	els_PFCand_idx_branch = 0;
	if (tree->GetAlias("els_PFCand_idx") != 0) {
		els_PFCand_idx_branch = tree->GetBranch(tree->GetAlias("els_PFCand_idx"));
		els_PFCand_idx_branch->SetAddress(&els_PFCand_idx_);
	}
	genps_lepdaughter_id_branch = 0;
	if (tree->GetAlias("genps_lepdaughter_id") != 0) {
		genps_lepdaughter_id_branch = tree->GetBranch(tree->GetAlias("genps_lepdaughter_id"));
		genps_lepdaughter_id_branch->SetAddress(&genps_lepdaughter_id_);
	}
	genps_lepdaughter_idx_branch = 0;
	if (tree->GetAlias("genps_lepdaughter_idx") != 0) {
		genps_lepdaughter_idx_branch = tree->GetBranch(tree->GetAlias("genps_lepdaughter_idx"));
		genps_lepdaughter_idx_branch->SetAddress(&genps_lepdaughter_idx_);
	}
	hlt_trigObjs_id_branch = 0;
	if (tree->GetAlias("hlt_trigObjs_id") != 0) {
		hlt_trigObjs_id_branch = tree->GetBranch(tree->GetAlias("hlt_trigObjs_id"));
		hlt_trigObjs_id_branch->SetAddress(&hlt_trigObjs_id_);
	}
	mus_nStationCorrelatedHits_branch = 0;
	if (tree->GetAlias("mus_nStationCorrelatedHits") != 0) {
		mus_nStationCorrelatedHits_branch = tree->GetBranch(tree->GetAlias("mus_nStationCorrelatedHits"));
		mus_nStationCorrelatedHits_branch->SetAddress(&mus_nStationCorrelatedHits_);
	}
	mus_nStationHits_branch = 0;
	if (tree->GetAlias("mus_nStationHits") != 0) {
		mus_nStationHits_branch = tree->GetBranch(tree->GetAlias("mus_nStationHits"));
		mus_nStationHits_branch->SetAddress(&mus_nStationHits_);
	}
	pfjets_pfcandIndicies_branch = 0;
	if (tree->GetAlias("pfjets_pfcandIndicies") != 0) {
		pfjets_pfcandIndicies_branch = tree->GetBranch(tree->GetAlias("pfjets_pfcandIndicies"));
		pfjets_pfcandIndicies_branch->SetAddress(&pfjets_pfcandIndicies_);
	}
	convs_nHitsBeforeVtx_branch = 0;
	if (tree->GetAlias("convs_nHitsBeforeVtx") != 0) {
		convs_nHitsBeforeVtx_branch = tree->GetBranch(tree->GetAlias("convs_nHitsBeforeVtx"));
		convs_nHitsBeforeVtx_branch->SetAddress(&convs_nHitsBeforeVtx_);
	}
	convs_tkalgo_branch = 0;
	if (tree->GetAlias("convs_tkalgo") != 0) {
		convs_tkalgo_branch = tree->GetBranch(tree->GetAlias("convs_tkalgo"));
		convs_tkalgo_branch->SetAddress(&convs_tkalgo_);
	}
	convs_tkidx_branch = 0;
	if (tree->GetAlias("convs_tkidx") != 0) {
		convs_tkidx_branch = tree->GetBranch(tree->GetAlias("convs_tkidx"));
		convs_tkidx_branch->SetAddress(&convs_tkidx_);
	}
	els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch = 0;
	if (tree->GetAlias("els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version") != 0) {
		els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version"));
		els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch->SetAddress(&els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_);
	}
	els_HLT_Ele17_Ele8_LeadingLeg_version_branch = 0;
	if (tree->GetAlias("els_HLT_Ele17_Ele8_LeadingLeg_version") != 0) {
		els_HLT_Ele17_Ele8_LeadingLeg_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele17_Ele8_LeadingLeg_version"));
		els_HLT_Ele17_Ele8_LeadingLeg_version_branch->SetAddress(&els_HLT_Ele17_Ele8_LeadingLeg_version_);
	}
	els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch = 0;
	if (tree->GetAlias("els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version") != 0) {
		els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version"));
		els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch->SetAddress(&els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_);
	}
	els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch = 0;
	if (tree->GetAlias("els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version") != 0) {
		els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version"));
		els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch->SetAddress(&els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_);
	}
	els_HLT_Ele17_Ele8_TrailingLeg_version_branch = 0;
	if (tree->GetAlias("els_HLT_Ele17_Ele8_TrailingLeg_version") != 0) {
		els_HLT_Ele17_Ele8_TrailingLeg_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele17_Ele8_TrailingLeg_version"));
		els_HLT_Ele17_Ele8_TrailingLeg_version_branch->SetAddress(&els_HLT_Ele17_Ele8_TrailingLeg_version_);
	}
	els_HLT_Ele17_Ele8_version_branch = 0;
	if (tree->GetAlias("els_HLT_Ele17_Ele8_version") != 0) {
		els_HLT_Ele17_Ele8_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele17_Ele8_version"));
		els_HLT_Ele17_Ele8_version_branch->SetAddress(&els_HLT_Ele17_Ele8_version_);
	}
	els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch = 0;
	if (tree->GetAlias("els_HLT_Ele20_SC4_Mass50_LeadingLeg_version") != 0) {
		els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele20_SC4_Mass50_LeadingLeg_version"));
		els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch->SetAddress(&els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_);
	}
	els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch = 0;
	if (tree->GetAlias("els_HLT_Ele20_SC4_Mass50_TrailingLeg_version") != 0) {
		els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele20_SC4_Mass50_TrailingLeg_version"));
		els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch->SetAddress(&els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_);
	}
	els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch = 0;
	if (tree->GetAlias("els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version") != 0) {
		els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version"));
		els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch->SetAddress(&els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_);
	}
	els_HLT_Ele27_WP80_version_branch = 0;
	if (tree->GetAlias("els_HLT_Ele27_WP80_version") != 0) {
		els_HLT_Ele27_WP80_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele27_WP80_version"));
		els_HLT_Ele27_WP80_version_branch->SetAddress(&els_HLT_Ele27_WP80_version_);
	}
	els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch = 0;
	if (tree->GetAlias("els_HLT_Ele32_SC17_Mass50_LeadingLeg_version") != 0) {
		els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele32_SC17_Mass50_LeadingLeg_version"));
		els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch->SetAddress(&els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_);
	}
	els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch = 0;
	if (tree->GetAlias("els_HLT_Ele32_SC17_Mass50_TrailingLeg_version") != 0) {
		els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele32_SC17_Mass50_TrailingLeg_version"));
		els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch->SetAddress(&els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_);
	}
	els_HLT_Mu17_Ele8_TrailingLeg_version_branch = 0;
	if (tree->GetAlias("els_HLT_Mu17_Ele8_TrailingLeg_version") != 0) {
		els_HLT_Mu17_Ele8_TrailingLeg_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Mu17_Ele8_TrailingLeg_version"));
		els_HLT_Mu17_Ele8_TrailingLeg_version_branch->SetAddress(&els_HLT_Mu17_Ele8_TrailingLeg_version_);
	}
	els_HLT_Mu17_Ele8_version_branch = 0;
	if (tree->GetAlias("els_HLT_Mu17_Ele8_version") != 0) {
		els_HLT_Mu17_Ele8_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Mu17_Ele8_version"));
		els_HLT_Mu17_Ele8_version_branch->SetAddress(&els_HLT_Mu17_Ele8_version_);
	}
	els_HLT_Mu8_Ele17_version_branch = 0;
	if (tree->GetAlias("els_HLT_Mu8_Ele17_version") != 0) {
		els_HLT_Mu8_Ele17_version_branch = tree->GetBranch(tree->GetAlias("els_HLT_Mu8_Ele17_version"));
		els_HLT_Mu8_Ele17_version_branch->SetAddress(&els_HLT_Mu8_Ele17_version_);
	}
	evt_nels_branch = 0;
	if (tree->GetAlias("evt_nels") != 0) {
		evt_nels_branch = tree->GetBranch(tree->GetAlias("evt_nels"));
		evt_nels_branch->SetAddress(&evt_nels_);
	}
	evt_detectorStatus_branch = 0;
	if (tree->GetAlias("evt_detectorStatus") != 0) {
		evt_detectorStatus_branch = tree->GetBranch(tree->GetAlias("evt_detectorStatus"));
		evt_detectorStatus_branch->SetAddress(&evt_detectorStatus_);
	}
	evt_event_branch = 0;
	if (tree->GetAlias("evt_event") != 0) {
		evt_event_branch = tree->GetBranch(tree->GetAlias("evt_event"));
		evt_event_branch->SetAddress(&evt_event_);
	}
	evt_lumiBlock_branch = 0;
	if (tree->GetAlias("evt_lumiBlock") != 0) {
		evt_lumiBlock_branch = tree->GetBranch(tree->GetAlias("evt_lumiBlock"));
		evt_lumiBlock_branch->SetAddress(&evt_lumiBlock_);
	}
	evt_run_branch = 0;
	if (tree->GetAlias("evt_run") != 0) {
		evt_run_branch = tree->GetBranch(tree->GetAlias("evt_run"));
		evt_run_branch->SetAddress(&evt_run_);
	}
	evt_ngenjetsNoMuNoNu_branch = 0;
	if (tree->GetAlias("evt_ngenjetsNoMuNoNu") != 0) {
		evt_ngenjetsNoMuNoNu_branch = tree->GetBranch(tree->GetAlias("evt_ngenjetsNoMuNoNu"));
		evt_ngenjetsNoMuNoNu_branch->SetAddress(&evt_ngenjetsNoMuNoNu_);
	}
	genps_signalProcessID_branch = 0;
	if (tree->GetAlias("genps_signalProcessID") != 0) {
		genps_signalProcessID_branch = tree->GetBranch(tree->GetAlias("genps_signalProcessID"));
		genps_signalProcessID_branch->SetAddress(&genps_signalProcessID_);
	}
	mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch = 0;
	if (tree->GetAlias("mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version") != 0) {
		mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version"));
		mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch->SetAddress(&mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_);
	}
	mus_HLT_IsoMu24_eta2p1_version_branch = 0;
	if (tree->GetAlias("mus_HLT_IsoMu24_eta2p1_version") != 0) {
		mus_HLT_IsoMu24_eta2p1_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_IsoMu24_eta2p1_version"));
		mus_HLT_IsoMu24_eta2p1_version_branch->SetAddress(&mus_HLT_IsoMu24_eta2p1_version_);
	}
	mus_HLT_Mu17_Ele8_LeadingLeg_version_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_Ele8_LeadingLeg_version") != 0) {
		mus_HLT_Mu17_Ele8_LeadingLeg_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_Ele8_LeadingLeg_version"));
		mus_HLT_Mu17_Ele8_LeadingLeg_version_branch->SetAddress(&mus_HLT_Mu17_Ele8_LeadingLeg_version_);
	}
	mus_HLT_Mu17_Ele8_version_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_Ele8_version") != 0) {
		mus_HLT_Mu17_Ele8_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_Ele8_version"));
		mus_HLT_Mu17_Ele8_version_branch->SetAddress(&mus_HLT_Mu17_Ele8_version_);
	}
	mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version") != 0) {
		mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version"));
		mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch->SetAddress(&mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_);
	}
	mus_HLT_Mu17_Mu8_LeadingLeg_version_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_Mu8_LeadingLeg_version") != 0) {
		mus_HLT_Mu17_Mu8_LeadingLeg_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_Mu8_LeadingLeg_version"));
		mus_HLT_Mu17_Mu8_LeadingLeg_version_branch->SetAddress(&mus_HLT_Mu17_Mu8_LeadingLeg_version_);
	}
	mus_HLT_Mu17_Mu8_TrailingLeg_version_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_Mu8_TrailingLeg_version") != 0) {
		mus_HLT_Mu17_Mu8_TrailingLeg_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_Mu8_TrailingLeg_version"));
		mus_HLT_Mu17_Mu8_TrailingLeg_version_branch->SetAddress(&mus_HLT_Mu17_Mu8_TrailingLeg_version_);
	}
	mus_HLT_Mu17_Mu8_version_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_Mu8_version") != 0) {
		mus_HLT_Mu17_Mu8_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_Mu8_version"));
		mus_HLT_Mu17_Mu8_version_branch->SetAddress(&mus_HLT_Mu17_Mu8_version_);
	}
	mus_HLT_Mu17_TkMu8_LeadingLeg_version_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_TkMu8_LeadingLeg_version") != 0) {
		mus_HLT_Mu17_TkMu8_LeadingLeg_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_TkMu8_LeadingLeg_version"));
		mus_HLT_Mu17_TkMu8_LeadingLeg_version_branch->SetAddress(&mus_HLT_Mu17_TkMu8_LeadingLeg_version_);
	}
	mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version") != 0) {
		mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version"));
		mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch->SetAddress(&mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_);
	}
	mus_HLT_Mu17_TkMu8_TrailingLeg_version_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_TkMu8_TrailingLeg_version") != 0) {
		mus_HLT_Mu17_TkMu8_TrailingLeg_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_TkMu8_TrailingLeg_version"));
		mus_HLT_Mu17_TkMu8_TrailingLeg_version_branch->SetAddress(&mus_HLT_Mu17_TkMu8_TrailingLeg_version_);
	}
	mus_HLT_Mu17_TkMu8_version_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_TkMu8_version") != 0) {
		mus_HLT_Mu17_TkMu8_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_TkMu8_version"));
		mus_HLT_Mu17_TkMu8_version_branch->SetAddress(&mus_HLT_Mu17_TkMu8_version_);
	}
	mus_HLT_Mu8_Ele17_TrailingLeg_version_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu8_Ele17_TrailingLeg_version") != 0) {
		mus_HLT_Mu8_Ele17_TrailingLeg_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu8_Ele17_TrailingLeg_version"));
		mus_HLT_Mu8_Ele17_TrailingLeg_version_branch->SetAddress(&mus_HLT_Mu8_Ele17_TrailingLeg_version_);
	}
	mus_HLT_Mu8_Ele17_version_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu8_Ele17_version") != 0) {
		mus_HLT_Mu8_Ele17_version_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu8_Ele17_version"));
		mus_HLT_Mu8_Ele17_version_branch->SetAddress(&mus_HLT_Mu8_Ele17_version_);
	}
	evt_nphotons_branch = 0;
	if (tree->GetAlias("evt_nphotons") != 0) {
		evt_nphotons_branch = tree->GetBranch(tree->GetAlias("evt_nphotons"));
		evt_nphotons_branch->SetAddress(&evt_nphotons_);
	}
	evt_nvtxs_branch = 0;
	if (tree->GetAlias("evt_nvtxs") != 0) {
		evt_nvtxs_branch = tree->GetBranch(tree->GetAlias("evt_nvtxs"));
		evt_nvtxs_branch->SetAddress(&evt_nvtxs_);
	}
	els_HLT_Ele17_Ele8_branch = 0;
	if (tree->GetAlias("els_HLT_Ele17_Ele8") != 0) {
		els_HLT_Ele17_Ele8_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele17_Ele8"));
		els_HLT_Ele17_Ele8_branch->SetAddress(&els_HLT_Ele17_Ele8_);
	}
	els_HLT_Ele17_Ele8_L1sL1DoubleEG137_branch = 0;
	if (tree->GetAlias("els_HLT_Ele17_Ele8_L1sL1DoubleEG137") != 0) {
		els_HLT_Ele17_Ele8_L1sL1DoubleEG137_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele17_Ele8_L1sL1DoubleEG137"));
		els_HLT_Ele17_Ele8_L1sL1DoubleEG137_branch->SetAddress(&els_HLT_Ele17_Ele8_L1sL1DoubleEG137_);
	}
	els_HLT_Ele17_Ele8_LeadingLeg_branch = 0;
	if (tree->GetAlias("els_HLT_Ele17_Ele8_LeadingLeg") != 0) {
		els_HLT_Ele17_Ele8_LeadingLeg_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele17_Ele8_LeadingLeg"));
		els_HLT_Ele17_Ele8_LeadingLeg_branch->SetAddress(&els_HLT_Ele17_Ele8_LeadingLeg_);
	}
	els_HLT_Ele17_Ele8_Mass50_LeadingLeg_branch = 0;
	if (tree->GetAlias("els_HLT_Ele17_Ele8_Mass50_LeadingLeg") != 0) {
		els_HLT_Ele17_Ele8_Mass50_LeadingLeg_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele17_Ele8_Mass50_LeadingLeg"));
		els_HLT_Ele17_Ele8_Mass50_LeadingLeg_branch->SetAddress(&els_HLT_Ele17_Ele8_Mass50_LeadingLeg_);
	}
	els_HLT_Ele17_Ele8_Mass50_TrailingLeg_branch = 0;
	if (tree->GetAlias("els_HLT_Ele17_Ele8_Mass50_TrailingLeg") != 0) {
		els_HLT_Ele17_Ele8_Mass50_TrailingLeg_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele17_Ele8_Mass50_TrailingLeg"));
		els_HLT_Ele17_Ele8_Mass50_TrailingLeg_branch->SetAddress(&els_HLT_Ele17_Ele8_Mass50_TrailingLeg_);
	}
	els_HLT_Ele17_Ele8_TrailingLeg_branch = 0;
	if (tree->GetAlias("els_HLT_Ele17_Ele8_TrailingLeg") != 0) {
		els_HLT_Ele17_Ele8_TrailingLeg_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele17_Ele8_TrailingLeg"));
		els_HLT_Ele17_Ele8_TrailingLeg_branch->SetAddress(&els_HLT_Ele17_Ele8_TrailingLeg_);
	}
	els_HLT_Ele20_SC4_Mass50_LeadingLeg_branch = 0;
	if (tree->GetAlias("els_HLT_Ele20_SC4_Mass50_LeadingLeg") != 0) {
		els_HLT_Ele20_SC4_Mass50_LeadingLeg_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele20_SC4_Mass50_LeadingLeg"));
		els_HLT_Ele20_SC4_Mass50_LeadingLeg_branch->SetAddress(&els_HLT_Ele20_SC4_Mass50_LeadingLeg_);
	}
	els_HLT_Ele20_SC4_Mass50_TrailingLeg_branch = 0;
	if (tree->GetAlias("els_HLT_Ele20_SC4_Mass50_TrailingLeg") != 0) {
		els_HLT_Ele20_SC4_Mass50_TrailingLeg_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele20_SC4_Mass50_TrailingLeg"));
		els_HLT_Ele20_SC4_Mass50_TrailingLeg_branch->SetAddress(&els_HLT_Ele20_SC4_Mass50_TrailingLeg_);
	}
	els_HLT_Ele27_WP80_branch = 0;
	if (tree->GetAlias("els_HLT_Ele27_WP80") != 0) {
		els_HLT_Ele27_WP80_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele27_WP80"));
		els_HLT_Ele27_WP80_branch->SetAddress(&els_HLT_Ele27_WP80_);
	}
	els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_branch = 0;
	if (tree->GetAlias("els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22") != 0) {
		els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22"));
		els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_branch->SetAddress(&els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_);
	}
	els_HLT_Ele32_SC17_Mass50_LeadingLeg_branch = 0;
	if (tree->GetAlias("els_HLT_Ele32_SC17_Mass50_LeadingLeg") != 0) {
		els_HLT_Ele32_SC17_Mass50_LeadingLeg_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele32_SC17_Mass50_LeadingLeg"));
		els_HLT_Ele32_SC17_Mass50_LeadingLeg_branch->SetAddress(&els_HLT_Ele32_SC17_Mass50_LeadingLeg_);
	}
	els_HLT_Ele32_SC17_Mass50_TrailingLeg_branch = 0;
	if (tree->GetAlias("els_HLT_Ele32_SC17_Mass50_TrailingLeg") != 0) {
		els_HLT_Ele32_SC17_Mass50_TrailingLeg_branch = tree->GetBranch(tree->GetAlias("els_HLT_Ele32_SC17_Mass50_TrailingLeg"));
		els_HLT_Ele32_SC17_Mass50_TrailingLeg_branch->SetAddress(&els_HLT_Ele32_SC17_Mass50_TrailingLeg_);
	}
	els_HLT_Mu17_Ele8_branch = 0;
	if (tree->GetAlias("els_HLT_Mu17_Ele8") != 0) {
		els_HLT_Mu17_Ele8_branch = tree->GetBranch(tree->GetAlias("els_HLT_Mu17_Ele8"));
		els_HLT_Mu17_Ele8_branch->SetAddress(&els_HLT_Mu17_Ele8_);
	}
	els_HLT_Mu17_Ele8_TrailingLeg_branch = 0;
	if (tree->GetAlias("els_HLT_Mu17_Ele8_TrailingLeg") != 0) {
		els_HLT_Mu17_Ele8_TrailingLeg_branch = tree->GetBranch(tree->GetAlias("els_HLT_Mu17_Ele8_TrailingLeg"));
		els_HLT_Mu17_Ele8_TrailingLeg_branch->SetAddress(&els_HLT_Mu17_Ele8_TrailingLeg_);
	}
	els_HLT_Mu8_Ele17_branch = 0;
	if (tree->GetAlias("els_HLT_Mu8_Ele17") != 0) {
		els_HLT_Mu8_Ele17_branch = tree->GetBranch(tree->GetAlias("els_HLT_Mu8_Ele17"));
		els_HLT_Mu8_Ele17_branch->SetAddress(&els_HLT_Mu8_Ele17_);
	}
	els_id2012ext_loose_branch = 0;
	if (tree->GetAlias("els_id2012ext_loose") != 0) {
		els_id2012ext_loose_branch = tree->GetBranch(tree->GetAlias("els_id2012ext_loose"));
		els_id2012ext_loose_branch->SetAddress(&els_id2012ext_loose_);
	}
	els_id2012ext_medium_branch = 0;
	if (tree->GetAlias("els_id2012ext_medium") != 0) {
		els_id2012ext_medium_branch = tree->GetBranch(tree->GetAlias("els_id2012ext_medium"));
		els_id2012ext_medium_branch->SetAddress(&els_id2012ext_medium_);
	}
	els_id2012ext_tight_branch = 0;
	if (tree->GetAlias("els_id2012ext_tight") != 0) {
		els_id2012ext_tight_branch = tree->GetBranch(tree->GetAlias("els_id2012ext_tight"));
		els_id2012ext_tight_branch->SetAddress(&els_id2012ext_tight_);
	}
	els_id2012ext_veto_branch = 0;
	if (tree->GetAlias("els_id2012ext_veto") != 0) {
		els_id2012ext_veto_branch = tree->GetBranch(tree->GetAlias("els_id2012ext_veto"));
		els_id2012ext_veto_branch->SetAddress(&els_id2012ext_veto_);
	}
	els_id2012_loose_branch = 0;
	if (tree->GetAlias("els_id2012_loose") != 0) {
		els_id2012_loose_branch = tree->GetBranch(tree->GetAlias("els_id2012_loose"));
		els_id2012_loose_branch->SetAddress(&els_id2012_loose_);
	}
	els_id2012_medium_branch = 0;
	if (tree->GetAlias("els_id2012_medium") != 0) {
		els_id2012_medium_branch = tree->GetBranch(tree->GetAlias("els_id2012_medium"));
		els_id2012_medium_branch->SetAddress(&els_id2012_medium_);
	}
	els_id2012_tight_branch = 0;
	if (tree->GetAlias("els_id2012_tight") != 0) {
		els_id2012_tight_branch = tree->GetBranch(tree->GetAlias("els_id2012_tight"));
		els_id2012_tight_branch->SetAddress(&els_id2012_tight_);
	}
	els_id2012_veto_branch = 0;
	if (tree->GetAlias("els_id2012_veto") != 0) {
		els_id2012_veto_branch = tree->GetBranch(tree->GetAlias("els_id2012_veto"));
		els_id2012_veto_branch->SetAddress(&els_id2012_veto_);
	}
	hlt_prescales_branch = 0;
	if (tree->GetAlias("hlt_prescales") != 0) {
		hlt_prescales_branch = tree->GetBranch(tree->GetAlias("hlt_prescales"));
		hlt_prescales_branch->SetAddress(&hlt_prescales_);
	}
	mus_HLT_IsoMu24_eta2p1_branch = 0;
	if (tree->GetAlias("mus_HLT_IsoMu24_eta2p1") != 0) {
		mus_HLT_IsoMu24_eta2p1_branch = tree->GetBranch(tree->GetAlias("mus_HLT_IsoMu24_eta2p1"));
		mus_HLT_IsoMu24_eta2p1_branch->SetAddress(&mus_HLT_IsoMu24_eta2p1_);
	}
	mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_branch = 0;
	if (tree->GetAlias("mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1") != 0) {
		mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_branch = tree->GetBranch(tree->GetAlias("mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1"));
		mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_branch->SetAddress(&mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_);
	}
	mus_HLT_Mu17_Ele8_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_Ele8") != 0) {
		mus_HLT_Mu17_Ele8_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_Ele8"));
		mus_HLT_Mu17_Ele8_branch->SetAddress(&mus_HLT_Mu17_Ele8_);
	}
	mus_HLT_Mu17_Ele8_LeadingLeg_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_Ele8_LeadingLeg") != 0) {
		mus_HLT_Mu17_Ele8_LeadingLeg_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_Ele8_LeadingLeg"));
		mus_HLT_Mu17_Ele8_LeadingLeg_branch->SetAddress(&mus_HLT_Mu17_Ele8_LeadingLeg_);
	}
	mus_HLT_Mu17_Mu8_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_Mu8") != 0) {
		mus_HLT_Mu17_Mu8_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_Mu8"));
		mus_HLT_Mu17_Mu8_branch->SetAddress(&mus_HLT_Mu17_Mu8_);
	}
	mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen") != 0) {
		mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen"));
		mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_branch->SetAddress(&mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_);
	}
	mus_HLT_Mu17_Mu8_LeadingLeg_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_Mu8_LeadingLeg") != 0) {
		mus_HLT_Mu17_Mu8_LeadingLeg_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_Mu8_LeadingLeg"));
		mus_HLT_Mu17_Mu8_LeadingLeg_branch->SetAddress(&mus_HLT_Mu17_Mu8_LeadingLeg_);
	}
	mus_HLT_Mu17_Mu8_TrailingLeg_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_Mu8_TrailingLeg") != 0) {
		mus_HLT_Mu17_Mu8_TrailingLeg_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_Mu8_TrailingLeg"));
		mus_HLT_Mu17_Mu8_TrailingLeg_branch->SetAddress(&mus_HLT_Mu17_Mu8_TrailingLeg_);
	}
	mus_HLT_Mu17_TkMu8_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_TkMu8") != 0) {
		mus_HLT_Mu17_TkMu8_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_TkMu8"));
		mus_HLT_Mu17_TkMu8_branch->SetAddress(&mus_HLT_Mu17_TkMu8_);
	}
	mus_HLT_Mu17_TkMu8_LeadingLeg_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_TkMu8_LeadingLeg") != 0) {
		mus_HLT_Mu17_TkMu8_LeadingLeg_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_TkMu8_LeadingLeg"));
		mus_HLT_Mu17_TkMu8_LeadingLeg_branch->SetAddress(&mus_HLT_Mu17_TkMu8_LeadingLeg_);
	}
	mus_HLT_Mu17_TkMu8_TrailingLeg_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_TkMu8_TrailingLeg") != 0) {
		mus_HLT_Mu17_TkMu8_TrailingLeg_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_TkMu8_TrailingLeg"));
		mus_HLT_Mu17_TkMu8_TrailingLeg_branch->SetAddress(&mus_HLT_Mu17_TkMu8_TrailingLeg_);
	}
	mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered") != 0) {
		mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered"));
		mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_branch->SetAddress(&mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_);
	}
	mus_HLT_Mu8_Ele17_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu8_Ele17") != 0) {
		mus_HLT_Mu8_Ele17_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu8_Ele17"));
		mus_HLT_Mu8_Ele17_branch->SetAddress(&mus_HLT_Mu8_Ele17_);
	}
	mus_HLT_Mu8_Ele17_TrailingLeg_branch = 0;
	if (tree->GetAlias("mus_HLT_Mu8_Ele17_TrailingLeg") != 0) {
		mus_HLT_Mu8_Ele17_TrailingLeg_branch = tree->GetBranch(tree->GetAlias("mus_HLT_Mu8_Ele17_TrailingLeg"));
		mus_HLT_Mu8_Ele17_TrailingLeg_branch->SetAddress(&mus_HLT_Mu8_Ele17_TrailingLeg_);
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		hlt_bits_isLoaded = false;
		evt_CMS2tag_isLoaded = false;
		evt_dataset_isLoaded = false;
		hlt_trigNames_isLoaded = false;
		hcalnoise_HasBadRBXTS4TS5_isLoaded = false;
		evt_cscTightHaloId_isLoaded = false;
		evt_hbheFilter_isLoaded = false;
		filt_cscBeamHalo_isLoaded = false;
		filt_ecalLaser_isLoaded = false;
		filt_ecalTP_isLoaded = false;
		filt_eeBadSc_isLoaded = false;
		filt_goodVertices_isLoaded = false;
		filt_hbheNoise_isLoaded = false;
		filt_hcalLaser_isLoaded = false;
		filt_metfilter_isLoaded = false;
		filt_trackingFailure_isLoaded = false;
		filt_trkPOGFilters_isLoaded = false;
		filt_trkPOG_logErrorTooManyClusters_isLoaded = false;
		filt_trkPOG_manystripclus53X_isLoaded = false;
		filt_trkPOG_toomanystripclus53X_isLoaded = false;
		els_conv_vtx_flag_isLoaded = false;
		els_isGsfCtfScPixChargeConsistent_isLoaded = false;
		els_passingMvaPreselection_isLoaded = false;
		els_passingPflowPreselection_isLoaded = false;
		mus_isRPCMuon_isLoaded = false;
		mus_tightMatch_isLoaded = false;
		mus_updatedSta_isLoaded = false;
		photons_haspixelSeed_isLoaded = false;
		evt_bs_Xwidth_isLoaded = false;
		evt_bs_XwidthErr_isLoaded = false;
		evt_bs_Ywidth_isLoaded = false;
		evt_bs_YwidthErr_isLoaded = false;
		evt_bs_dxdz_isLoaded = false;
		evt_bs_dxdzErr_isLoaded = false;
		evt_bs_dydz_isLoaded = false;
		evt_bs_dydzErr_isLoaded = false;
		evt_bs_sigmaZ_isLoaded = false;
		evt_bs_sigmaZErr_isLoaded = false;
		evt_bs_xErr_isLoaded = false;
		evt_bs_yErr_isLoaded = false;
		evt_bs_zErr_isLoaded = false;
		evt_bField_isLoaded = false;
		evt_fixgrid_all_rho_isLoaded = false;
		evt_fixgridfastjet_allcalo_rho_isLoaded = false;
		evt_fixgridfastjet_all_rho_isLoaded = false;
		evt_fixgridfastjet_centralcalo_rho_isLoaded = false;
		evt_fixgridfastjet_centralchargedpileup_rho_isLoaded = false;
		evt_fixgridfastjet_centralneutral_rho_isLoaded = false;
		evt_kfactor_isLoaded = false;
		evt_scale1fb_isLoaded = false;
		evt_xsec_excl_isLoaded = false;
		evt_xsec_incl_isLoaded = false;
		gen_met_isLoaded = false;
		gen_metPhi_isLoaded = false;
		genps_alphaQCD_isLoaded = false;
		genps_pthat_isLoaded = false;
		genps_qScale_isLoaded = false;
		genps_weight_isLoaded = false;
		gen_sumEt_isLoaded = false;
		hcalnoise_GetRecHitEnergy_isLoaded = false;
		hcalnoise_GetRecHitEnergy15_isLoaded = false;
		hcalnoise_GetTotalCalibCharge_isLoaded = false;
		hcalnoise_TS4TS5NoiseSumE_isLoaded = false;
		hcalnoise_TS4TS5NoiseSumEt_isLoaded = false;
		hcalnoise_eventChargeFraction_isLoaded = false;
		hcalnoise_eventEMEnergy_isLoaded = false;
		hcalnoise_eventEMFraction_isLoaded = false;
		hcalnoise_eventHadEnergy_isLoaded = false;
		hcalnoise_eventTrackEnergy_isLoaded = false;
		hcalnoise_flatNoiseSumE_isLoaded = false;
		hcalnoise_flatNoiseSumEt_isLoaded = false;
		hcalnoise_isolatedNoiseSumE_isLoaded = false;
		hcalnoise_isolatedNoiseSumEt_isLoaded = false;
		hcalnoise_max10GeVHitTime_isLoaded = false;
		hcalnoise_max25GeVHitTime_isLoaded = false;
		hcalnoise_maxE10TS_isLoaded = false;
		hcalnoise_maxE2Over10TS_isLoaded = false;
		hcalnoise_maxE2TS_isLoaded = false;
		hcalnoise_min10GeVHitTime_isLoaded = false;
		hcalnoise_min25GeVHitTime_isLoaded = false;
		hcalnoise_minE10TS_isLoaded = false;
		hcalnoise_minE2Over10TS_isLoaded = false;
		hcalnoise_minE2TS_isLoaded = false;
		hcalnoise_minHPDEMF_isLoaded = false;
		hcalnoise_minRBXEMF_isLoaded = false;
		hcalnoise_rms10GeVHitTime_isLoaded = false;
		hcalnoise_rms25GeVHitTime_isLoaded = false;
		hcalnoise_spikeNoiseSumE_isLoaded = false;
		hcalnoise_spikeNoiseSumEt_isLoaded = false;
		hcalnoise_triangleNoiseSumE_isLoaded = false;
		hcalnoise_triangleNoiseSumEt_isLoaded = false;
		pdfinfo_pdf1_isLoaded = false;
		pdfinfo_pdf2_isLoaded = false;
		pdfinfo_scale_isLoaded = false;
		pdfinfo_x1_isLoaded = false;
		pdfinfo_x2_isLoaded = false;
		evt_fixgrid_rho_all_isLoaded = false;
		evt_fixgrid_rho_ctr_isLoaded = false;
		evt_fixgrid_rho_fwd_isLoaded = false;
		evt_pfmet_isLoaded = false;
		evt_pfmetPhi_isLoaded = false;
		evt_pfmetSig_isLoaded = false;
		evt_pfmetSignificance_isLoaded = false;
		evt_pfsumet_isLoaded = false;
		evt_bsp4_isLoaded = false;
		els_mc_motherp4_isLoaded = false;
		els_mc_p4_isLoaded = false;
		mus_mc_motherp4_isLoaded = false;
		mus_mc_p4_isLoaded = false;
		pfjets_mc_gp_p4_isLoaded = false;
		pfjets_mc_motherp4_isLoaded = false;
		pfjets_mc_p4_isLoaded = false;
		photons_mc_motherp4_isLoaded = false;
		photons_mc_p4_isLoaded = false;
		els_conv_pos_p4_isLoaded = false;
		els_inner_position_isLoaded = false;
		els_mc_patMatch_p4_isLoaded = false;
		els_outer_position_isLoaded = false;
		els_p4_isLoaded = false;
		els_p4In_isLoaded = false;
		els_p4Out_isLoaded = false;
		els_trk_p4_isLoaded = false;
		els_trk_vertex_p4_isLoaded = false;
		els_vertex_p4_isLoaded = false;
		genjets_p4NoMuNoNu_isLoaded = false;
		genps_p4_isLoaded = false;
		genps_prod_vtx_isLoaded = false;
		mus_gfit_p4_isLoaded = false;
		mus_gfit_vertex_p4_isLoaded = false;
		mus_mc_patMatch_p4_isLoaded = false;
		mus_p4_isLoaded = false;
		mus_pfp4_isLoaded = false;
		mus_sta_p4_isLoaded = false;
		mus_sta_vertex_p4_isLoaded = false;
		mus_trk_p4_isLoaded = false;
		mus_vertex_p4_isLoaded = false;
		pfcands_p4_isLoaded = false;
		pfjets_p4_isLoaded = false;
		taus_pf_lead_chargecand_p4_isLoaded = false;
		taus_pf_lead_neutrcand_p4_isLoaded = false;
		taus_pf_p4_isLoaded = false;
		photons_p4_isLoaded = false;
		vtxs_position_isLoaded = false;
		genps_lepdaughter_p4_isLoaded = false;
		hlt_trigObjs_p4_isLoaded = false;
		taus_pf_isocands_p4_isLoaded = false;
		taus_pf_signalcands_p4_isLoaded = false;
		evt_bs_covMatrix_isLoaded = false;
		els_mc3dr_isLoaded = false;
		els_mcdr_isLoaded = false;
		mus_mc3dr_isLoaded = false;
		mus_mcdr_isLoaded = false;
		pfjets_mc3dr_isLoaded = false;
		pfjets_mcdr_isLoaded = false;
		pfjets_mc_emEnergy_isLoaded = false;
		pfjets_mc_gpdr_isLoaded = false;
		pfjets_mc_hadEnergy_isLoaded = false;
		pfjets_mc_invEnergy_isLoaded = false;
		pfjets_mc_otherEnergy_isLoaded = false;
		photons_mc3dr_isLoaded = false;
		photons_mcdr_isLoaded = false;
		els_chi2_isLoaded = false;
		els_ckf_chi2_isLoaded = false;
		els_ckf_ndof_isLoaded = false;
		els_conv_dcot_isLoaded = false;
		els_conv_dist_isLoaded = false;
		els_conv_old_dcot_isLoaded = false;
		els_conv_old_dist_isLoaded = false;
		els_conv_old_radius_isLoaded = false;
		els_conv_radius_isLoaded = false;
		els_d0_isLoaded = false;
		els_d0Err_isLoaded = false;
		els_d0corr_isLoaded = false;
		els_d0corrPhi_isLoaded = false;
		els_d0phiCov_isLoaded = false;
		els_dEtaIn_isLoaded = false;
		els_dEtaOut_isLoaded = false;
		els_dPhiIn_isLoaded = false;
		els_dPhiInPhiOut_isLoaded = false;
		els_dPhiOut_isLoaded = false;
		els_deltaEtaEleClusterTrackAtCalo_isLoaded = false;
		els_deltaPhiEleClusterTrackAtCalo_isLoaded = false;
		els_dxyPV_isLoaded = false;
		els_dzPV_isLoaded = false;
		els_e1x5_isLoaded = false;
		els_e1x5_full5x5_isLoaded = false;
		els_e2x5Max_isLoaded = false;
		els_e2x5Max_full5x5_isLoaded = false;
		els_e3x3_isLoaded = false;
		els_e5x5_isLoaded = false;
		els_e5x5_full5x5_isLoaded = false;
		els_eMax_isLoaded = false;
		els_eOverPIn_isLoaded = false;
		els_eOverPOut_isLoaded = false;
		els_eSC_isLoaded = false;
		els_eSCPresh_isLoaded = false;
		els_eSCRaw_isLoaded = false;
		els_eSeed_isLoaded = false;
		els_eSeedOverPIn_isLoaded = false;
		els_eSeedOverPOut_isLoaded = false;
		els_ecalEnergy_isLoaded = false;
		els_ecalEnergyError_isLoaded = false;
		els_ecalIso_isLoaded = false;
		els_ecalIso04_isLoaded = false;
		els_etaErr_isLoaded = false;
		els_etaSC_isLoaded = false;
		els_etaSCwidth_isLoaded = false;
		els_fbrem_isLoaded = false;
		els_hOverE_isLoaded = false;
		els_hOverEBC_isLoaded = false;
		els_hcalDepth1OverEcal_isLoaded = false;
		els_hcalDepth1TowerSumEt_isLoaded = false;
		els_hcalDepth1TowerSumEt04_isLoaded = false;
		els_hcalDepth2OverEcal_isLoaded = false;
		els_hcalDepth2TowerSumEt_isLoaded = false;
		els_hcalDepth2TowerSumEt04_isLoaded = false;
		els_hcalIso_isLoaded = false;
		els_hcalIso04_isLoaded = false;
		els_ip3d_isLoaded = false;
		els_ip3derr_isLoaded = false;
		els_iso03_pf_isLoaded = false;
		els_iso03_pf2012_ch_isLoaded = false;
		els_iso03_pf2012_em_isLoaded = false;
		els_iso03_pf2012ext_ch_isLoaded = false;
		els_iso03_pf2012ext_em_isLoaded = false;
		els_iso03_pf2012ext_nh_isLoaded = false;
		els_iso03_pf2012_nh_isLoaded = false;
		els_iso03_pf_ch_isLoaded = false;
		els_iso03_pf_gamma05_isLoaded = false;
		els_iso03_pf_nhad05_isLoaded = false;
		els_iso04_pf_isLoaded = false;
		els_iso04_pf2012_ch_isLoaded = false;
		els_iso04_pf2012_em_isLoaded = false;
		els_iso04_pf2012ext_ch_isLoaded = false;
		els_iso04_pf2012ext_em_isLoaded = false;
		els_iso04_pf2012ext_nh_isLoaded = false;
		els_iso04_pf2012_nh_isLoaded = false;
		els_iso04_pf_ch_isLoaded = false;
		els_iso04_pf_gamma05_isLoaded = false;
		els_iso04_pf_nhad05_isLoaded = false;
		els_layer1_charge_isLoaded = false;
		els_lh_isLoaded = false;
		els_mc_patMatch_dr_isLoaded = false;
		els_mva_isLoaded = false;
		els_ndof_isLoaded = false;
		els_pfChargedHadronIso_isLoaded = false;
		els_pfNeutralHadronIso_isLoaded = false;
		els_pfPUIso_isLoaded = false;
		els_pfPhotonIso_isLoaded = false;
		els_phiErr_isLoaded = false;
		els_phiSC_isLoaded = false;
		els_phiSCwidth_isLoaded = false;
		els_ptErr_isLoaded = false;
		els_ptErrGsf_isLoaded = false;
		els_r9_isLoaded = false;
		els_r9_full5x5_isLoaded = false;
		els_sigmaEtaEta_isLoaded = false;
		els_sigmaEtaEta_full5x5_isLoaded = false;
		els_sigmaIEtaIEta_isLoaded = false;
		els_sigmaIEtaIEtaSC_isLoaded = false;
		els_sigmaIEtaIEta_full5x5_isLoaded = false;
		els_sigmaIEtaIPhi_isLoaded = false;
		els_sigmaIPhiIPhi_isLoaded = false;
		els_sigmaIPhiIPhiSC_isLoaded = false;
		els_sigmaIPhiIPhi_full5x5_isLoaded = false;
		els_sigmaIphiIphi_isLoaded = false;
		els_sigmaPhiPhi_isLoaded = false;
		els_tkIso_isLoaded = false;
		els_tkIso04_isLoaded = false;
		els_trackMomentumError_isLoaded = false;
		els_trkdr_isLoaded = false;
		els_trkshFrac_isLoaded = false;
		els_z0_isLoaded = false;
		els_z0Err_isLoaded = false;
		els_z0corr_isLoaded = false;
		mus_backToBackCompat_isLoaded = false;
		mus_best_dxyPV_isLoaded = false;
		mus_best_dzPV_isLoaded = false;
		mus_caloCompatibility_isLoaded = false;
		mus_chi2_isLoaded = false;
		mus_chi2LocalMomentum_isLoaded = false;
		mus_chi2LocalPosition_isLoaded = false;
		mus_cosmicCompat_isLoaded = false;
		mus_d0_isLoaded = false;
		mus_d0Err_isLoaded = false;
		mus_d0corr_isLoaded = false;
		mus_d0corrPhi_isLoaded = false;
		mus_d0phiCov_isLoaded = false;
		mus_dxyPV_isLoaded = false;
		mus_dzPV_isLoaded = false;
		mus_ecal_time_isLoaded = false;
		mus_ecal_timeError_isLoaded = false;
		mus_e_em_isLoaded = false;
		mus_e_emS9_isLoaded = false;
		mus_e_had_isLoaded = false;
		mus_e_hadS9_isLoaded = false;
		mus_e_ho_isLoaded = false;
		mus_e_hoS9_isLoaded = false;
		mus_emMax_isLoaded = false;
		mus_emS25_isLoaded = false;
		mus_etaErr_isLoaded = false;
		mus_gfit_chi2_isLoaded = false;
		mus_gfit_d0_isLoaded = false;
		mus_gfit_d0Err_isLoaded = false;
		mus_gfit_d0corr_isLoaded = false;
		mus_gfit_d0corrPhi_isLoaded = false;
		mus_gfit_d0phiCov_isLoaded = false;
		mus_gfit_ndof_isLoaded = false;
		mus_gfit_qoverp_isLoaded = false;
		mus_gfit_qoverpError_isLoaded = false;
		mus_gfit_z0_isLoaded = false;
		mus_gfit_z0Err_isLoaded = false;
		mus_gfit_z0corr_isLoaded = false;
		mus_glbKink_isLoaded = false;
		mus_glbTrackProbability_isLoaded = false;
		mus_globalDeltaEtaPhi_isLoaded = false;
		mus_hadMax_isLoaded = false;
		mus_hcal_time_isLoaded = false;
		mus_hcal_timeError_isLoaded = false;
		mus_ip3d_isLoaded = false;
		mus_ip3derr_isLoaded = false;
		mus_iso03_emEt_isLoaded = false;
		mus_iso03_hadEt_isLoaded = false;
		mus_iso03_hoEt_isLoaded = false;
		mus_iso03_sumPt_isLoaded = false;
		mus_iso05_emEt_isLoaded = false;
		mus_iso05_hadEt_isLoaded = false;
		mus_iso05_hoEt_isLoaded = false;
		mus_iso05_sumPt_isLoaded = false;
		mus_isoMeanDRR03_pf_ChargedHadronPt_isLoaded = false;
		mus_isoMeanDRR03_pf_ChargedParticlePt_isLoaded = false;
		mus_isoMeanDRR03_pf_NeutralHadronEt_isLoaded = false;
		mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_isLoaded = false;
		mus_isoMeanDRR03_pf_PUPt_isLoaded = false;
		mus_isoMeanDRR03_pf_PhotonEt_isLoaded = false;
		mus_isoMeanDRR03_pf_PhotonEtHighThreshold_isLoaded = false;
		mus_isoMeanDRR04_pf_ChargedHadronPt_isLoaded = false;
		mus_isoMeanDRR04_pf_ChargedParticlePt_isLoaded = false;
		mus_isoMeanDRR04_pf_NeutralHadronEt_isLoaded = false;
		mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_isLoaded = false;
		mus_isoMeanDRR04_pf_PUPt_isLoaded = false;
		mus_isoMeanDRR04_pf_PhotonEt_isLoaded = false;
		mus_isoMeanDRR04_pf_PhotonEtHighThreshold_isLoaded = false;
		mus_isoR03_pf_ChargedHadronPt_isLoaded = false;
		mus_isoR03_pf_ChargedParticlePt_isLoaded = false;
		mus_isoR03_pf_NeutralHadronEt_isLoaded = false;
		mus_isoR03_pf_NeutralHadronEtHighThreshold_isLoaded = false;
		mus_isoR03_pf_PUPt_isLoaded = false;
		mus_isoR03_pf_PhotonEt_isLoaded = false;
		mus_isoR03_pf_PhotonEtHighThreshold_isLoaded = false;
		mus_isoR04_pf_ChargedHadronPt_isLoaded = false;
		mus_isoR04_pf_ChargedParticlePt_isLoaded = false;
		mus_isoR04_pf_NeutralHadronEt_isLoaded = false;
		mus_isoR04_pf_NeutralHadronEtHighThreshold_isLoaded = false;
		mus_isoR04_pf_PUPt_isLoaded = false;
		mus_isoR04_pf_PhotonEt_isLoaded = false;
		mus_isoR04_pf_PhotonEtHighThreshold_isLoaded = false;
		mus_isoSumDRR03_pf_ChargedHadronPt_isLoaded = false;
		mus_isoSumDRR03_pf_ChargedParticlePt_isLoaded = false;
		mus_isoSumDRR03_pf_NeutralHadronEt_isLoaded = false;
		mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_isLoaded = false;
		mus_isoSumDRR03_pf_PUPt_isLoaded = false;
		mus_isoSumDRR03_pf_PhotonEt_isLoaded = false;
		mus_isoSumDRR03_pf_PhotonEtHighThreshold_isLoaded = false;
		mus_isoSumDRR04_pf_ChargedHadronPt_isLoaded = false;
		mus_isoSumDRR04_pf_ChargedParticlePt_isLoaded = false;
		mus_isoSumDRR04_pf_NeutralHadronEt_isLoaded = false;
		mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_isLoaded = false;
		mus_isoSumDRR04_pf_PUPt_isLoaded = false;
		mus_isoSumDRR04_pf_PhotonEt_isLoaded = false;
		mus_isoSumDRR04_pf_PhotonEtHighThreshold_isLoaded = false;
		mus_iso_ecalvetoDep_isLoaded = false;
		mus_iso_hcalvetoDep_isLoaded = false;
		mus_iso_hovetoDep_isLoaded = false;
		mus_iso_trckvetoDep_isLoaded = false;
		mus_localDistance_isLoaded = false;
		mus_mc_patMatch_dr_isLoaded = false;
		mus_ndof_isLoaded = false;
		mus_overlapCompat_isLoaded = false;
		mus_phiErr_isLoaded = false;
		mus_ptErr_isLoaded = false;
		mus_qoverp_isLoaded = false;
		mus_qoverpError_isLoaded = false;
		mus_segmCompatibility_isLoaded = false;
		mus_staRelChi2_isLoaded = false;
		mus_sta_chi2_isLoaded = false;
		mus_sta_d0_isLoaded = false;
		mus_sta_d0Err_isLoaded = false;
		mus_sta_d0corr_isLoaded = false;
		mus_sta_d0corrPhi_isLoaded = false;
		mus_sta_d0phiCov_isLoaded = false;
		mus_sta_ndof_isLoaded = false;
		mus_sta_qoverp_isLoaded = false;
		mus_sta_qoverpError_isLoaded = false;
		mus_sta_z0_isLoaded = false;
		mus_sta_z0Err_isLoaded = false;
		mus_sta_z0corr_isLoaded = false;
		mus_timeAtIpInOut_isLoaded = false;
		mus_timeAtIpInOutErr_isLoaded = false;
		mus_timeAtIpOutIn_isLoaded = false;
		mus_timeAtIpOutInErr_isLoaded = false;
		mus_timeCompat_isLoaded = false;
		mus_tower_isLoaded = false;
		mus_towerS9_isLoaded = false;
		mus_trkKink_isLoaded = false;
		mus_trkRelChi2_isLoaded = false;
		mus_vertexCompat_isLoaded = false;
		mus_vertexphi_isLoaded = false;
		mus_z0_isLoaded = false;
		mus_z0Err_isLoaded = false;
		mus_z0corr_isLoaded = false;
		pfcands_dz_isLoaded = false;
		pfjets_area_isLoaded = false;
		pfjets_chargedEmE_isLoaded = false;
		pfjets_chargedHadronE_isLoaded = false;
		pfjets_combinedSecondaryVertexBJetTag_isLoaded = false;
		pfjets_cor_isLoaded = false;
		pfjets_corL1Fast_isLoaded = false;
		pfjets_corL1FastL2L3_isLoaded = false;
		pfjets_corL1FastL2L3residual_isLoaded = false;
		pfjets_corL1L2L3_isLoaded = false;
		pfjets_electronE_isLoaded = false;
		pfjets_hfEmE_isLoaded = false;
		pfjets_hfHadronE_isLoaded = false;
		pfjets_jetBProbabilityBJetTag_isLoaded = false;
		pfjets_jetProbabilityBJetTag_isLoaded = false;
		pfjets_muonE_isLoaded = false;
		pfjets_neutralEmE_isLoaded = false;
		pfjets_neutralHadronE_isLoaded = false;
		pfjets_photonE_isLoaded = false;
		pfjets_pileupJetId_isLoaded = false;
		pfjets_simpleSecondaryVertexHighEffBJetTag_isLoaded = false;
		pfjets_simpleSecondaryVertexHighPurBJetTags_isLoaded = false;
		pfjets_trackCountingHighEffBJetTag_isLoaded = false;
		pfjets_trackCountingHighPurBJetTag_isLoaded = false;
		taus_pf_againstElectronDeadECAL_isLoaded = false;
		taus_pf_againstElectronLoose_isLoaded = false;
		taus_pf_againstElectronMedium_isLoaded = false;
		taus_pf_againstElectronTight_isLoaded = false;
		taus_pf_againstMuonLoose_isLoaded = false;
		taus_pf_againstMuonLoose2_isLoaded = false;
		taus_pf_againstMuonLoose3_isLoaded = false;
		taus_pf_againstMuonMedium_isLoaded = false;
		taus_pf_againstMuonMedium2_isLoaded = false;
		taus_pf_againstMuonTight_isLoaded = false;
		taus_pf_againstMuonTight2_isLoaded = false;
		taus_pf_againstMuonTight3_isLoaded = false;
		taus_pf_byCombinedIsolationDeltaBetaCorrRaw_isLoaded = false;
		taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_isLoaded = false;
		taus_pf_byDecayModeFinding_isLoaded = false;
		taus_pf_byLooseCombinedIsolationDeltaBetaCorr_isLoaded = false;
		taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_isLoaded = false;
		taus_pf_byMediumCombinedIsolationDeltaBetaCorr_isLoaded = false;
		taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_isLoaded = false;
		taus_pf_byTightCombinedIsolationDeltaBetaCorr_isLoaded = false;
		taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_isLoaded = false;
		taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_isLoaded = false;
		photons_e1x5_isLoaded = false;
		photons_e2x5Max_isLoaded = false;
		photons_e3x3_isLoaded = false;
		photons_e5x5_isLoaded = false;
		photons_ecalIso03_isLoaded = false;
		photons_ecalIso04_isLoaded = false;
		photons_full3x3_e3x3_isLoaded = false;
		photons_full5x5_e1x5_isLoaded = false;
		photons_full5x5_e2x5Max_isLoaded = false;
		photons_full5x5_e5x5_isLoaded = false;
		photons_full5x5_hOverE_isLoaded = false;
		photons_full5x5_sigmaEtaEta_isLoaded = false;
		photons_full5x5_sigmaIEtaIEta_isLoaded = false;
		photons_hOverE_isLoaded = false;
		photons_hcalDepth1TowerSumEtBcConeDR03_isLoaded = false;
		photons_hcalDepth1TowerSumEtBcConeDR04_isLoaded = false;
		photons_hcalDepth2TowerSumEtBcConeDR03_isLoaded = false;
		photons_hcalDepth2TowerSumEtBcConeDR04_isLoaded = false;
		photons_hcalIso03_isLoaded = false;
		photons_hcalIso04_isLoaded = false;
		photons_hcalTowerSumEtBcConeDR03_isLoaded = false;
		photons_hcalTowerSumEtBcConeDR04_isLoaded = false;
		photons_ntkIsoHollow03_isLoaded = false;
		photons_ntkIsoHollow04_isLoaded = false;
		photons_ntkIsoSolid03_isLoaded = false;
		photons_ntkIsoSolid04_isLoaded = false;
		photons_sigmaEtaEta_isLoaded = false;
		photons_sigmaIEtaIEta_isLoaded = false;
		photons_tkIsoHollow03_isLoaded = false;
		photons_tkIsoHollow04_isLoaded = false;
		photons_tkIsoSolid03_isLoaded = false;
		photons_tkIsoSolid04_isLoaded = false;
		puInfo_trueNumInteractions_isLoaded = false;
		convs_chi2_isLoaded = false;
		convs_dl_isLoaded = false;
		convs_ndof_isLoaded = false;
		vtxs_chi2_isLoaded = false;
		vtxs_ndof_isLoaded = false;
		vtxs_sumpt_isLoaded = false;
		vtxs_xError_isLoaded = false;
		vtxs_yError_isLoaded = false;
		vtxs_zError_isLoaded = false;
		els_convs_dcot_isLoaded = false;
		els_convs_dist_isLoaded = false;
		els_convs_radius_isLoaded = false;
		mus_stationShowerDeltaR_isLoaded = false;
		mus_stationShowerSizeT_isLoaded = false;
		puInfo_instLumi_isLoaded = false;
		vtxs_covMatrix_isLoaded = false;
		evt_bsType_isLoaded = false;
		evt_bunchCrossing_isLoaded = false;
		evt_experimentType_isLoaded = false;
		evt_isRealData_isLoaded = false;
		evt_orbitNumber_isLoaded = false;
		evt_storeNumber_isLoaded = false;
		hcalnoise_GetRecHitCount_isLoaded = false;
		hcalnoise_GetRecHitCount15_isLoaded = false;
		hcalnoise_maxHPDHits_isLoaded = false;
		hcalnoise_maxHPDNoOtherHits_isLoaded = false;
		hcalnoise_maxRBXHits_isLoaded = false;
		hcalnoise_maxZeros_isLoaded = false;
		hcalnoise_noiseFilterStatus_isLoaded = false;
		hcalnoise_noiseType_isLoaded = false;
		hcalnoise_num10GeVHits_isLoaded = false;
		hcalnoise_num25GeVHits_isLoaded = false;
		hcalnoise_numFlatNoiseChannels_isLoaded = false;
		hcalnoise_numIsolatedNoiseChannels_isLoaded = false;
		hcalnoise_numProblematicRBXs_isLoaded = false;
		hcalnoise_numSpikeNoiseChannels_isLoaded = false;
		hcalnoise_numTS4TS5NoiseChannels_isLoaded = false;
		hcalnoise_numTriangleNoiseChannels_isLoaded = false;
		hcalnoise_passHighLevelNoiseFilter_isLoaded = false;
		hcalnoise_passLooseNoiseFilter_isLoaded = false;
		hcalnoise_passTightNoiseFilter_isLoaded = false;
		pdfinfo_id1_isLoaded = false;
		pdfinfo_id2_isLoaded = false;
		els_mc3_id_isLoaded = false;
		els_mc3idx_isLoaded = false;
		els_mc3_motherid_isLoaded = false;
		els_mc3_motheridx_isLoaded = false;
		els_mc_id_isLoaded = false;
		els_mcidx_isLoaded = false;
		els_mc_motherid_isLoaded = false;
		mus_mc3_id_isLoaded = false;
		mus_mc3idx_isLoaded = false;
		mus_mc3_motherid_isLoaded = false;
		mus_mc3_motheridx_isLoaded = false;
		mus_mc_id_isLoaded = false;
		mus_mcidx_isLoaded = false;
		mus_mc_motherid_isLoaded = false;
		pfjets_mc3_id_isLoaded = false;
		pfjets_mc3idx_isLoaded = false;
		pfjets_mc_gpidx_isLoaded = false;
		pfjets_mc_id_isLoaded = false;
		pfjets_mcidx_isLoaded = false;
		pfjets_mc_motherid_isLoaded = false;
		photons_mc3_id_isLoaded = false;
		photons_mc3idx_isLoaded = false;
		photons_mc3_motherid_isLoaded = false;
		photons_mc3_motheridx_isLoaded = false;
		photons_mc_id_isLoaded = false;
		photons_mcidx_isLoaded = false;
		photons_mc_motherid_isLoaded = false;
		els_category_isLoaded = false;
		els_charge_isLoaded = false;
		els_ckf_charge_isLoaded = false;
		els_ckf_laywithmeas_isLoaded = false;
		els_class_isLoaded = false;
		els_conv_delMissHits_isLoaded = false;
		els_conv_flag_isLoaded = false;
		els_conv_gsftkidx_isLoaded = false;
		els_conv_old_delMissHits_isLoaded = false;
		els_conv_old_flag_isLoaded = false;
		els_conv_old_gsftkidx_isLoaded = false;
		els_conv_old_tkidx_isLoaded = false;
		els_conv_tkidx_isLoaded = false;
		els_exp_innerlayers_isLoaded = false;
		els_exp_outerlayers_isLoaded = false;
		els_fiduciality_isLoaded = false;
		els_gsftrkidx_isLoaded = false;
		els_layer1_det_isLoaded = false;
		els_layer1_layer_isLoaded = false;
		els_layer1_sizerphi_isLoaded = false;
		els_layer1_sizerz_isLoaded = false;
		els_lostHits_isLoaded = false;
		els_lost_pixelhits_isLoaded = false;
		els_mc_patMatch_id_isLoaded = false;
		els_nSeed_isLoaded = false;
		els_nlayers_isLoaded = false;
		els_nlayers3D_isLoaded = false;
		els_nlayersLost_isLoaded = false;
		els_sccharge_isLoaded = false;
		els_scindex_isLoaded = false;
		els_trk_charge_isLoaded = false;
		els_trkidx_isLoaded = false;
		els_type_isLoaded = false;
		els_validHits_isLoaded = false;
		els_valid_pixelhits_isLoaded = false;
		genps_id_isLoaded = false;
		genps_id_mother_isLoaded = false;
		genps_status_isLoaded = false;
		mus_algo_isLoaded = false;
		mus_charge_isLoaded = false;
		mus_ecal_rawId_isLoaded = false;
		mus_exp_innerlayers_isLoaded = false;
		mus_exp_outerlayers_isLoaded = false;
		mus_gfit_algo_isLoaded = false;
		mus_gfit_etaErr_isLoaded = false;
		mus_gfit_exp_innerlayer_isLoaded = false;
		mus_gfit_exp_outerlayer_isLoaded = false;
		mus_gfit_lostHits_isLoaded = false;
		mus_gfit_lostPixelHits_isLoaded = false;
		mus_gfit_nlayers_isLoaded = false;
		mus_gfit_nlayers3D_isLoaded = false;
		mus_gfit_nlayersLost_isLoaded = false;
		mus_gfit_phiErr_isLoaded = false;
		mus_gfit_ptErr_isLoaded = false;
		mus_gfit_qualityMask_isLoaded = false;
		mus_gfit_trk_charge_isLoaded = false;
		mus_gfit_validHits_isLoaded = false;
		mus_gfit_validPixelHits_isLoaded = false;
		mus_gfit_validSTAHits_isLoaded = false;
		mus_gfit_validSiHits_isLoaded = false;
		mus_goodmask_isLoaded = false;
		mus_hcal_rawId_isLoaded = false;
		mus_iso03_njets_isLoaded = false;
		mus_iso03_ntrk_isLoaded = false;
		mus_iso05_ntrk_isLoaded = false;
		mus_lostHits_isLoaded = false;
		mus_lostPixelHits_isLoaded = false;
		mus_mc_patMatch_id_isLoaded = false;
		mus_muonBestTrackType_isLoaded = false;
		mus_nOverlaps_isLoaded = false;
		mus_nlayers_isLoaded = false;
		mus_nlayers3D_isLoaded = false;
		mus_nlayersLost_isLoaded = false;
		mus_nmatches_isLoaded = false;
		mus_numberOfMatchedStations_isLoaded = false;
		mus_overlap0_isLoaded = false;
		mus_overlap1_isLoaded = false;
		mus_pfcharge_isLoaded = false;
		mus_pfparticleId_isLoaded = false;
		mus_pid_PFMuon_isLoaded = false;
		mus_pid_TM2DCompatibilityLoose_isLoaded = false;
		mus_pid_TM2DCompatibilityTight_isLoaded = false;
		mus_pid_TMLastStationLoose_isLoaded = false;
		mus_pid_TMLastStationTight_isLoaded = false;
		mus_qualityMask_isLoaded = false;
		mus_sta_algo_isLoaded = false;
		mus_sta_etaErr_isLoaded = false;
		mus_sta_exp_innerlayer_isLoaded = false;
		mus_sta_exp_outerlayer_isLoaded = false;
		mus_sta_lostHits_isLoaded = false;
		mus_sta_lostPixelHits_isLoaded = false;
		mus_sta_nlayers_isLoaded = false;
		mus_sta_nlayers3D_isLoaded = false;
		mus_sta_nlayersLost_isLoaded = false;
		mus_sta_phiErr_isLoaded = false;
		mus_sta_ptErr_isLoaded = false;
		mus_sta_qualityMask_isLoaded = false;
		mus_sta_trk_charge_isLoaded = false;
		mus_sta_validHits_isLoaded = false;
		mus_sta_validPixelHits_isLoaded = false;
		mus_timeDirection_isLoaded = false;
		mus_timeNumStationsUsed_isLoaded = false;
		mus_trk_charge_isLoaded = false;
		mus_trkidx_isLoaded = false;
		mus_type_isLoaded = false;
		mus_validHits_isLoaded = false;
		mus_validPixelHits_isLoaded = false;
		pfcands_charge_isLoaded = false;
		pfcands_particleId_isLoaded = false;
		pfjets_chargedHadronMultiplicity_isLoaded = false;
		pfjets_chargedMultiplicity_isLoaded = false;
		pfjets_electronMultiplicity_isLoaded = false;
		pfjets_hfEmMultiplicity_isLoaded = false;
		pfjets_hfHadronMultiplicity_isLoaded = false;
		pfjets_muonMultiplicity_isLoaded = false;
		pfjets_neutralHadronMultiplicity_isLoaded = false;
		pfjets_neutralMultiplicity_isLoaded = false;
		pfjets_partonFlavour_isLoaded = false;
		pfjets_photonMultiplicity_isLoaded = false;
		taus_pf_charge_isLoaded = false;
		photons_fiduciality_isLoaded = false;
		puInfo_bunchCrossing_isLoaded = false;
		puInfo_nPUvertices_isLoaded = false;
		convs_algo_isLoaded = false;
		convs_isConverted_isLoaded = false;
		convs_quality_isLoaded = false;
		vtxs_isFake_isLoaded = false;
		vtxs_isValid_isLoaded = false;
		vtxs_tracksSize_isLoaded = false;
		els_convs_delMissHits_isLoaded = false;
		els_convs_flag_isLoaded = false;
		els_convs_gsftkidx_isLoaded = false;
		els_convs_tkidx_isLoaded = false;
		els_PFCand_idx_isLoaded = false;
		genps_lepdaughter_id_isLoaded = false;
		genps_lepdaughter_idx_isLoaded = false;
		hlt_trigObjs_id_isLoaded = false;
		mus_nStationCorrelatedHits_isLoaded = false;
		mus_nStationHits_isLoaded = false;
		pfjets_pfcandIndicies_isLoaded = false;
		convs_nHitsBeforeVtx_isLoaded = false;
		convs_tkalgo_isLoaded = false;
		convs_tkidx_isLoaded = false;
		els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_isLoaded = false;
		els_HLT_Ele17_Ele8_LeadingLeg_version_isLoaded = false;
		els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_isLoaded = false;
		els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_isLoaded = false;
		els_HLT_Ele17_Ele8_TrailingLeg_version_isLoaded = false;
		els_HLT_Ele17_Ele8_version_isLoaded = false;
		els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_isLoaded = false;
		els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_isLoaded = false;
		els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_isLoaded = false;
		els_HLT_Ele27_WP80_version_isLoaded = false;
		els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_isLoaded = false;
		els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_isLoaded = false;
		els_HLT_Mu17_Ele8_TrailingLeg_version_isLoaded = false;
		els_HLT_Mu17_Ele8_version_isLoaded = false;
		els_HLT_Mu8_Ele17_version_isLoaded = false;
		evt_nels_isLoaded = false;
		evt_detectorStatus_isLoaded = false;
		evt_event_isLoaded = false;
		evt_lumiBlock_isLoaded = false;
		evt_run_isLoaded = false;
		evt_ngenjetsNoMuNoNu_isLoaded = false;
		genps_signalProcessID_isLoaded = false;
		mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_isLoaded = false;
		mus_HLT_IsoMu24_eta2p1_version_isLoaded = false;
		mus_HLT_Mu17_Ele8_LeadingLeg_version_isLoaded = false;
		mus_HLT_Mu17_Ele8_version_isLoaded = false;
		mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_isLoaded = false;
		mus_HLT_Mu17_Mu8_LeadingLeg_version_isLoaded = false;
		mus_HLT_Mu17_Mu8_TrailingLeg_version_isLoaded = false;
		mus_HLT_Mu17_Mu8_version_isLoaded = false;
		mus_HLT_Mu17_TkMu8_LeadingLeg_version_isLoaded = false;
		mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_isLoaded = false;
		mus_HLT_Mu17_TkMu8_TrailingLeg_version_isLoaded = false;
		mus_HLT_Mu17_TkMu8_version_isLoaded = false;
		mus_HLT_Mu8_Ele17_TrailingLeg_version_isLoaded = false;
		mus_HLT_Mu8_Ele17_version_isLoaded = false;
		evt_nphotons_isLoaded = false;
		evt_nvtxs_isLoaded = false;
		els_HLT_Ele17_Ele8_isLoaded = false;
		els_HLT_Ele17_Ele8_L1sL1DoubleEG137_isLoaded = false;
		els_HLT_Ele17_Ele8_LeadingLeg_isLoaded = false;
		els_HLT_Ele17_Ele8_Mass50_LeadingLeg_isLoaded = false;
		els_HLT_Ele17_Ele8_Mass50_TrailingLeg_isLoaded = false;
		els_HLT_Ele17_Ele8_TrailingLeg_isLoaded = false;
		els_HLT_Ele20_SC4_Mass50_LeadingLeg_isLoaded = false;
		els_HLT_Ele20_SC4_Mass50_TrailingLeg_isLoaded = false;
		els_HLT_Ele27_WP80_isLoaded = false;
		els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_isLoaded = false;
		els_HLT_Ele32_SC17_Mass50_LeadingLeg_isLoaded = false;
		els_HLT_Ele32_SC17_Mass50_TrailingLeg_isLoaded = false;
		els_HLT_Mu17_Ele8_isLoaded = false;
		els_HLT_Mu17_Ele8_TrailingLeg_isLoaded = false;
		els_HLT_Mu8_Ele17_isLoaded = false;
		els_id2012ext_loose_isLoaded = false;
		els_id2012ext_medium_isLoaded = false;
		els_id2012ext_tight_isLoaded = false;
		els_id2012ext_veto_isLoaded = false;
		els_id2012_loose_isLoaded = false;
		els_id2012_medium_isLoaded = false;
		els_id2012_tight_isLoaded = false;
		els_id2012_veto_isLoaded = false;
		hlt_prescales_isLoaded = false;
		mus_HLT_IsoMu24_eta2p1_isLoaded = false;
		mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_isLoaded = false;
		mus_HLT_Mu17_Ele8_isLoaded = false;
		mus_HLT_Mu17_Ele8_LeadingLeg_isLoaded = false;
		mus_HLT_Mu17_Mu8_isLoaded = false;
		mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_isLoaded = false;
		mus_HLT_Mu17_Mu8_LeadingLeg_isLoaded = false;
		mus_HLT_Mu17_Mu8_TrailingLeg_isLoaded = false;
		mus_HLT_Mu17_TkMu8_isLoaded = false;
		mus_HLT_Mu17_TkMu8_LeadingLeg_isLoaded = false;
		mus_HLT_Mu17_TkMu8_TrailingLeg_isLoaded = false;
		mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_isLoaded = false;
		mus_HLT_Mu8_Ele17_isLoaded = false;
		mus_HLT_Mu8_Ele17_TrailingLeg_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (hlt_bits_branch != 0) hlt_bits();
	if (evt_CMS2tag_branch != 0) evt_CMS2tag();
	if (evt_dataset_branch != 0) evt_dataset();
	if (hlt_trigNames_branch != 0) hlt_trigNames();
	if (hcalnoise_HasBadRBXTS4TS5_branch != 0) hcalnoise_HasBadRBXTS4TS5();
	if (evt_cscTightHaloId_branch != 0) evt_cscTightHaloId();
	if (evt_hbheFilter_branch != 0) evt_hbheFilter();
	if (filt_cscBeamHalo_branch != 0) filt_cscBeamHalo();
	if (filt_ecalLaser_branch != 0) filt_ecalLaser();
	if (filt_ecalTP_branch != 0) filt_ecalTP();
	if (filt_eeBadSc_branch != 0) filt_eeBadSc();
	if (filt_goodVertices_branch != 0) filt_goodVertices();
	if (filt_hbheNoise_branch != 0) filt_hbheNoise();
	if (filt_hcalLaser_branch != 0) filt_hcalLaser();
	if (filt_metfilter_branch != 0) filt_metfilter();
	if (filt_trackingFailure_branch != 0) filt_trackingFailure();
	if (filt_trkPOGFilters_branch != 0) filt_trkPOGFilters();
	if (filt_trkPOG_logErrorTooManyClusters_branch != 0) filt_trkPOG_logErrorTooManyClusters();
	if (filt_trkPOG_manystripclus53X_branch != 0) filt_trkPOG_manystripclus53X();
	if (filt_trkPOG_toomanystripclus53X_branch != 0) filt_trkPOG_toomanystripclus53X();
	if (els_conv_vtx_flag_branch != 0) els_conv_vtx_flag();
	if (els_isGsfCtfScPixChargeConsistent_branch != 0) els_isGsfCtfScPixChargeConsistent();
	if (els_passingMvaPreselection_branch != 0) els_passingMvaPreselection();
	if (els_passingPflowPreselection_branch != 0) els_passingPflowPreselection();
	if (mus_isRPCMuon_branch != 0) mus_isRPCMuon();
	if (mus_tightMatch_branch != 0) mus_tightMatch();
	if (mus_updatedSta_branch != 0) mus_updatedSta();
	if (photons_haspixelSeed_branch != 0) photons_haspixelSeed();
	if (evt_bs_Xwidth_branch != 0) evt_bs_Xwidth();
	if (evt_bs_XwidthErr_branch != 0) evt_bs_XwidthErr();
	if (evt_bs_Ywidth_branch != 0) evt_bs_Ywidth();
	if (evt_bs_YwidthErr_branch != 0) evt_bs_YwidthErr();
	if (evt_bs_dxdz_branch != 0) evt_bs_dxdz();
	if (evt_bs_dxdzErr_branch != 0) evt_bs_dxdzErr();
	if (evt_bs_dydz_branch != 0) evt_bs_dydz();
	if (evt_bs_dydzErr_branch != 0) evt_bs_dydzErr();
	if (evt_bs_sigmaZ_branch != 0) evt_bs_sigmaZ();
	if (evt_bs_sigmaZErr_branch != 0) evt_bs_sigmaZErr();
	if (evt_bs_xErr_branch != 0) evt_bs_xErr();
	if (evt_bs_yErr_branch != 0) evt_bs_yErr();
	if (evt_bs_zErr_branch != 0) evt_bs_zErr();
	if (evt_bField_branch != 0) evt_bField();
	if (evt_fixgrid_all_rho_branch != 0) evt_fixgrid_all_rho();
	if (evt_fixgridfastjet_allcalo_rho_branch != 0) evt_fixgridfastjet_allcalo_rho();
	if (evt_fixgridfastjet_all_rho_branch != 0) evt_fixgridfastjet_all_rho();
	if (evt_fixgridfastjet_centralcalo_rho_branch != 0) evt_fixgridfastjet_centralcalo_rho();
	if (evt_fixgridfastjet_centralchargedpileup_rho_branch != 0) evt_fixgridfastjet_centralchargedpileup_rho();
	if (evt_fixgridfastjet_centralneutral_rho_branch != 0) evt_fixgridfastjet_centralneutral_rho();
	if (evt_kfactor_branch != 0) evt_kfactor();
	if (evt_scale1fb_branch != 0) evt_scale1fb();
	if (evt_xsec_excl_branch != 0) evt_xsec_excl();
	if (evt_xsec_incl_branch != 0) evt_xsec_incl();
	if (gen_met_branch != 0) gen_met();
	if (gen_metPhi_branch != 0) gen_metPhi();
	if (genps_alphaQCD_branch != 0) genps_alphaQCD();
	if (genps_pthat_branch != 0) genps_pthat();
	if (genps_qScale_branch != 0) genps_qScale();
	if (genps_weight_branch != 0) genps_weight();
	if (gen_sumEt_branch != 0) gen_sumEt();
	if (hcalnoise_GetRecHitEnergy_branch != 0) hcalnoise_GetRecHitEnergy();
	if (hcalnoise_GetRecHitEnergy15_branch != 0) hcalnoise_GetRecHitEnergy15();
	if (hcalnoise_GetTotalCalibCharge_branch != 0) hcalnoise_GetTotalCalibCharge();
	if (hcalnoise_TS4TS5NoiseSumE_branch != 0) hcalnoise_TS4TS5NoiseSumE();
	if (hcalnoise_TS4TS5NoiseSumEt_branch != 0) hcalnoise_TS4TS5NoiseSumEt();
	if (hcalnoise_eventChargeFraction_branch != 0) hcalnoise_eventChargeFraction();
	if (hcalnoise_eventEMEnergy_branch != 0) hcalnoise_eventEMEnergy();
	if (hcalnoise_eventEMFraction_branch != 0) hcalnoise_eventEMFraction();
	if (hcalnoise_eventHadEnergy_branch != 0) hcalnoise_eventHadEnergy();
	if (hcalnoise_eventTrackEnergy_branch != 0) hcalnoise_eventTrackEnergy();
	if (hcalnoise_flatNoiseSumE_branch != 0) hcalnoise_flatNoiseSumE();
	if (hcalnoise_flatNoiseSumEt_branch != 0) hcalnoise_flatNoiseSumEt();
	if (hcalnoise_isolatedNoiseSumE_branch != 0) hcalnoise_isolatedNoiseSumE();
	if (hcalnoise_isolatedNoiseSumEt_branch != 0) hcalnoise_isolatedNoiseSumEt();
	if (hcalnoise_max10GeVHitTime_branch != 0) hcalnoise_max10GeVHitTime();
	if (hcalnoise_max25GeVHitTime_branch != 0) hcalnoise_max25GeVHitTime();
	if (hcalnoise_maxE10TS_branch != 0) hcalnoise_maxE10TS();
	if (hcalnoise_maxE2Over10TS_branch != 0) hcalnoise_maxE2Over10TS();
	if (hcalnoise_maxE2TS_branch != 0) hcalnoise_maxE2TS();
	if (hcalnoise_min10GeVHitTime_branch != 0) hcalnoise_min10GeVHitTime();
	if (hcalnoise_min25GeVHitTime_branch != 0) hcalnoise_min25GeVHitTime();
	if (hcalnoise_minE10TS_branch != 0) hcalnoise_minE10TS();
	if (hcalnoise_minE2Over10TS_branch != 0) hcalnoise_minE2Over10TS();
	if (hcalnoise_minE2TS_branch != 0) hcalnoise_minE2TS();
	if (hcalnoise_minHPDEMF_branch != 0) hcalnoise_minHPDEMF();
	if (hcalnoise_minRBXEMF_branch != 0) hcalnoise_minRBXEMF();
	if (hcalnoise_rms10GeVHitTime_branch != 0) hcalnoise_rms10GeVHitTime();
	if (hcalnoise_rms25GeVHitTime_branch != 0) hcalnoise_rms25GeVHitTime();
	if (hcalnoise_spikeNoiseSumE_branch != 0) hcalnoise_spikeNoiseSumE();
	if (hcalnoise_spikeNoiseSumEt_branch != 0) hcalnoise_spikeNoiseSumEt();
	if (hcalnoise_triangleNoiseSumE_branch != 0) hcalnoise_triangleNoiseSumE();
	if (hcalnoise_triangleNoiseSumEt_branch != 0) hcalnoise_triangleNoiseSumEt();
	if (pdfinfo_pdf1_branch != 0) pdfinfo_pdf1();
	if (pdfinfo_pdf2_branch != 0) pdfinfo_pdf2();
	if (pdfinfo_scale_branch != 0) pdfinfo_scale();
	if (pdfinfo_x1_branch != 0) pdfinfo_x1();
	if (pdfinfo_x2_branch != 0) pdfinfo_x2();
	if (evt_fixgrid_rho_all_branch != 0) evt_fixgrid_rho_all();
	if (evt_fixgrid_rho_ctr_branch != 0) evt_fixgrid_rho_ctr();
	if (evt_fixgrid_rho_fwd_branch != 0) evt_fixgrid_rho_fwd();
	if (evt_pfmet_branch != 0) evt_pfmet();
	if (evt_pfmetPhi_branch != 0) evt_pfmetPhi();
	if (evt_pfmetSig_branch != 0) evt_pfmetSig();
	if (evt_pfmetSignificance_branch != 0) evt_pfmetSignificance();
	if (evt_pfsumet_branch != 0) evt_pfsumet();
	if (evt_bsp4_branch != 0) evt_bsp4();
	if (els_mc_motherp4_branch != 0) els_mc_motherp4();
	if (els_mc_p4_branch != 0) els_mc_p4();
	if (mus_mc_motherp4_branch != 0) mus_mc_motherp4();
	if (mus_mc_p4_branch != 0) mus_mc_p4();
	if (pfjets_mc_gp_p4_branch != 0) pfjets_mc_gp_p4();
	if (pfjets_mc_motherp4_branch != 0) pfjets_mc_motherp4();
	if (pfjets_mc_p4_branch != 0) pfjets_mc_p4();
	if (photons_mc_motherp4_branch != 0) photons_mc_motherp4();
	if (photons_mc_p4_branch != 0) photons_mc_p4();
	if (els_conv_pos_p4_branch != 0) els_conv_pos_p4();
	if (els_inner_position_branch != 0) els_inner_position();
	if (els_mc_patMatch_p4_branch != 0) els_mc_patMatch_p4();
	if (els_outer_position_branch != 0) els_outer_position();
	if (els_p4_branch != 0) els_p4();
	if (els_p4In_branch != 0) els_p4In();
	if (els_p4Out_branch != 0) els_p4Out();
	if (els_trk_p4_branch != 0) els_trk_p4();
	if (els_trk_vertex_p4_branch != 0) els_trk_vertex_p4();
	if (els_vertex_p4_branch != 0) els_vertex_p4();
	if (genjets_p4NoMuNoNu_branch != 0) genjets_p4NoMuNoNu();
	if (genps_p4_branch != 0) genps_p4();
	if (genps_prod_vtx_branch != 0) genps_prod_vtx();
	if (mus_gfit_p4_branch != 0) mus_gfit_p4();
	if (mus_gfit_vertex_p4_branch != 0) mus_gfit_vertex_p4();
	if (mus_mc_patMatch_p4_branch != 0) mus_mc_patMatch_p4();
	if (mus_p4_branch != 0) mus_p4();
	if (mus_pfp4_branch != 0) mus_pfp4();
	if (mus_sta_p4_branch != 0) mus_sta_p4();
	if (mus_sta_vertex_p4_branch != 0) mus_sta_vertex_p4();
	if (mus_trk_p4_branch != 0) mus_trk_p4();
	if (mus_vertex_p4_branch != 0) mus_vertex_p4();
	if (pfcands_p4_branch != 0) pfcands_p4();
	if (pfjets_p4_branch != 0) pfjets_p4();
	if (taus_pf_lead_chargecand_p4_branch != 0) taus_pf_lead_chargecand_p4();
	if (taus_pf_lead_neutrcand_p4_branch != 0) taus_pf_lead_neutrcand_p4();
	if (taus_pf_p4_branch != 0) taus_pf_p4();
	if (photons_p4_branch != 0) photons_p4();
	if (vtxs_position_branch != 0) vtxs_position();
	if (genps_lepdaughter_p4_branch != 0) genps_lepdaughter_p4();
	if (hlt_trigObjs_p4_branch != 0) hlt_trigObjs_p4();
	if (taus_pf_isocands_p4_branch != 0) taus_pf_isocands_p4();
	if (taus_pf_signalcands_p4_branch != 0) taus_pf_signalcands_p4();
	if (evt_bs_covMatrix_branch != 0) evt_bs_covMatrix();
	if (els_mc3dr_branch != 0) els_mc3dr();
	if (els_mcdr_branch != 0) els_mcdr();
	if (mus_mc3dr_branch != 0) mus_mc3dr();
	if (mus_mcdr_branch != 0) mus_mcdr();
	if (pfjets_mc3dr_branch != 0) pfjets_mc3dr();
	if (pfjets_mcdr_branch != 0) pfjets_mcdr();
	if (pfjets_mc_emEnergy_branch != 0) pfjets_mc_emEnergy();
	if (pfjets_mc_gpdr_branch != 0) pfjets_mc_gpdr();
	if (pfjets_mc_hadEnergy_branch != 0) pfjets_mc_hadEnergy();
	if (pfjets_mc_invEnergy_branch != 0) pfjets_mc_invEnergy();
	if (pfjets_mc_otherEnergy_branch != 0) pfjets_mc_otherEnergy();
	if (photons_mc3dr_branch != 0) photons_mc3dr();
	if (photons_mcdr_branch != 0) photons_mcdr();
	if (els_chi2_branch != 0) els_chi2();
	if (els_ckf_chi2_branch != 0) els_ckf_chi2();
	if (els_ckf_ndof_branch != 0) els_ckf_ndof();
	if (els_conv_dcot_branch != 0) els_conv_dcot();
	if (els_conv_dist_branch != 0) els_conv_dist();
	if (els_conv_old_dcot_branch != 0) els_conv_old_dcot();
	if (els_conv_old_dist_branch != 0) els_conv_old_dist();
	if (els_conv_old_radius_branch != 0) els_conv_old_radius();
	if (els_conv_radius_branch != 0) els_conv_radius();
	if (els_d0_branch != 0) els_d0();
	if (els_d0Err_branch != 0) els_d0Err();
	if (els_d0corr_branch != 0) els_d0corr();
	if (els_d0corrPhi_branch != 0) els_d0corrPhi();
	if (els_d0phiCov_branch != 0) els_d0phiCov();
	if (els_dEtaIn_branch != 0) els_dEtaIn();
	if (els_dEtaOut_branch != 0) els_dEtaOut();
	if (els_dPhiIn_branch != 0) els_dPhiIn();
	if (els_dPhiInPhiOut_branch != 0) els_dPhiInPhiOut();
	if (els_dPhiOut_branch != 0) els_dPhiOut();
	if (els_deltaEtaEleClusterTrackAtCalo_branch != 0) els_deltaEtaEleClusterTrackAtCalo();
	if (els_deltaPhiEleClusterTrackAtCalo_branch != 0) els_deltaPhiEleClusterTrackAtCalo();
	if (els_dxyPV_branch != 0) els_dxyPV();
	if (els_dzPV_branch != 0) els_dzPV();
	if (els_e1x5_branch != 0) els_e1x5();
	if (els_e1x5_full5x5_branch != 0) els_e1x5_full5x5();
	if (els_e2x5Max_branch != 0) els_e2x5Max();
	if (els_e2x5Max_full5x5_branch != 0) els_e2x5Max_full5x5();
	if (els_e3x3_branch != 0) els_e3x3();
	if (els_e5x5_branch != 0) els_e5x5();
	if (els_e5x5_full5x5_branch != 0) els_e5x5_full5x5();
	if (els_eMax_branch != 0) els_eMax();
	if (els_eOverPIn_branch != 0) els_eOverPIn();
	if (els_eOverPOut_branch != 0) els_eOverPOut();
	if (els_eSC_branch != 0) els_eSC();
	if (els_eSCPresh_branch != 0) els_eSCPresh();
	if (els_eSCRaw_branch != 0) els_eSCRaw();
	if (els_eSeed_branch != 0) els_eSeed();
	if (els_eSeedOverPIn_branch != 0) els_eSeedOverPIn();
	if (els_eSeedOverPOut_branch != 0) els_eSeedOverPOut();
	if (els_ecalEnergy_branch != 0) els_ecalEnergy();
	if (els_ecalEnergyError_branch != 0) els_ecalEnergyError();
	if (els_ecalIso_branch != 0) els_ecalIso();
	if (els_ecalIso04_branch != 0) els_ecalIso04();
	if (els_etaErr_branch != 0) els_etaErr();
	if (els_etaSC_branch != 0) els_etaSC();
	if (els_etaSCwidth_branch != 0) els_etaSCwidth();
	if (els_fbrem_branch != 0) els_fbrem();
	if (els_hOverE_branch != 0) els_hOverE();
	if (els_hOverEBC_branch != 0) els_hOverEBC();
	if (els_hcalDepth1OverEcal_branch != 0) els_hcalDepth1OverEcal();
	if (els_hcalDepth1TowerSumEt_branch != 0) els_hcalDepth1TowerSumEt();
	if (els_hcalDepth1TowerSumEt04_branch != 0) els_hcalDepth1TowerSumEt04();
	if (els_hcalDepth2OverEcal_branch != 0) els_hcalDepth2OverEcal();
	if (els_hcalDepth2TowerSumEt_branch != 0) els_hcalDepth2TowerSumEt();
	if (els_hcalDepth2TowerSumEt04_branch != 0) els_hcalDepth2TowerSumEt04();
	if (els_hcalIso_branch != 0) els_hcalIso();
	if (els_hcalIso04_branch != 0) els_hcalIso04();
	if (els_ip3d_branch != 0) els_ip3d();
	if (els_ip3derr_branch != 0) els_ip3derr();
	if (els_iso03_pf_branch != 0) els_iso03_pf();
	if (els_iso03_pf2012_ch_branch != 0) els_iso03_pf2012_ch();
	if (els_iso03_pf2012_em_branch != 0) els_iso03_pf2012_em();
	if (els_iso03_pf2012ext_ch_branch != 0) els_iso03_pf2012ext_ch();
	if (els_iso03_pf2012ext_em_branch != 0) els_iso03_pf2012ext_em();
	if (els_iso03_pf2012ext_nh_branch != 0) els_iso03_pf2012ext_nh();
	if (els_iso03_pf2012_nh_branch != 0) els_iso03_pf2012_nh();
	if (els_iso03_pf_ch_branch != 0) els_iso03_pf_ch();
	if (els_iso03_pf_gamma05_branch != 0) els_iso03_pf_gamma05();
	if (els_iso03_pf_nhad05_branch != 0) els_iso03_pf_nhad05();
	if (els_iso04_pf_branch != 0) els_iso04_pf();
	if (els_iso04_pf2012_ch_branch != 0) els_iso04_pf2012_ch();
	if (els_iso04_pf2012_em_branch != 0) els_iso04_pf2012_em();
	if (els_iso04_pf2012ext_ch_branch != 0) els_iso04_pf2012ext_ch();
	if (els_iso04_pf2012ext_em_branch != 0) els_iso04_pf2012ext_em();
	if (els_iso04_pf2012ext_nh_branch != 0) els_iso04_pf2012ext_nh();
	if (els_iso04_pf2012_nh_branch != 0) els_iso04_pf2012_nh();
	if (els_iso04_pf_ch_branch != 0) els_iso04_pf_ch();
	if (els_iso04_pf_gamma05_branch != 0) els_iso04_pf_gamma05();
	if (els_iso04_pf_nhad05_branch != 0) els_iso04_pf_nhad05();
	if (els_layer1_charge_branch != 0) els_layer1_charge();
	if (els_lh_branch != 0) els_lh();
	if (els_mc_patMatch_dr_branch != 0) els_mc_patMatch_dr();
	if (els_mva_branch != 0) els_mva();
	if (els_ndof_branch != 0) els_ndof();
	if (els_pfChargedHadronIso_branch != 0) els_pfChargedHadronIso();
	if (els_pfNeutralHadronIso_branch != 0) els_pfNeutralHadronIso();
	if (els_pfPUIso_branch != 0) els_pfPUIso();
	if (els_pfPhotonIso_branch != 0) els_pfPhotonIso();
	if (els_phiErr_branch != 0) els_phiErr();
	if (els_phiSC_branch != 0) els_phiSC();
	if (els_phiSCwidth_branch != 0) els_phiSCwidth();
	if (els_ptErr_branch != 0) els_ptErr();
	if (els_ptErrGsf_branch != 0) els_ptErrGsf();
	if (els_r9_branch != 0) els_r9();
	if (els_r9_full5x5_branch != 0) els_r9_full5x5();
	if (els_sigmaEtaEta_branch != 0) els_sigmaEtaEta();
	if (els_sigmaEtaEta_full5x5_branch != 0) els_sigmaEtaEta_full5x5();
	if (els_sigmaIEtaIEta_branch != 0) els_sigmaIEtaIEta();
	if (els_sigmaIEtaIEtaSC_branch != 0) els_sigmaIEtaIEtaSC();
	if (els_sigmaIEtaIEta_full5x5_branch != 0) els_sigmaIEtaIEta_full5x5();
	if (els_sigmaIEtaIPhi_branch != 0) els_sigmaIEtaIPhi();
	if (els_sigmaIPhiIPhi_branch != 0) els_sigmaIPhiIPhi();
	if (els_sigmaIPhiIPhiSC_branch != 0) els_sigmaIPhiIPhiSC();
	if (els_sigmaIPhiIPhi_full5x5_branch != 0) els_sigmaIPhiIPhi_full5x5();
	if (els_sigmaIphiIphi_branch != 0) els_sigmaIphiIphi();
	if (els_sigmaPhiPhi_branch != 0) els_sigmaPhiPhi();
	if (els_tkIso_branch != 0) els_tkIso();
	if (els_tkIso04_branch != 0) els_tkIso04();
	if (els_trackMomentumError_branch != 0) els_trackMomentumError();
	if (els_trkdr_branch != 0) els_trkdr();
	if (els_trkshFrac_branch != 0) els_trkshFrac();
	if (els_z0_branch != 0) els_z0();
	if (els_z0Err_branch != 0) els_z0Err();
	if (els_z0corr_branch != 0) els_z0corr();
	if (mus_backToBackCompat_branch != 0) mus_backToBackCompat();
	if (mus_best_dxyPV_branch != 0) mus_best_dxyPV();
	if (mus_best_dzPV_branch != 0) mus_best_dzPV();
	if (mus_caloCompatibility_branch != 0) mus_caloCompatibility();
	if (mus_chi2_branch != 0) mus_chi2();
	if (mus_chi2LocalMomentum_branch != 0) mus_chi2LocalMomentum();
	if (mus_chi2LocalPosition_branch != 0) mus_chi2LocalPosition();
	if (mus_cosmicCompat_branch != 0) mus_cosmicCompat();
	if (mus_d0_branch != 0) mus_d0();
	if (mus_d0Err_branch != 0) mus_d0Err();
	if (mus_d0corr_branch != 0) mus_d0corr();
	if (mus_d0corrPhi_branch != 0) mus_d0corrPhi();
	if (mus_d0phiCov_branch != 0) mus_d0phiCov();
	if (mus_dxyPV_branch != 0) mus_dxyPV();
	if (mus_dzPV_branch != 0) mus_dzPV();
	if (mus_ecal_time_branch != 0) mus_ecal_time();
	if (mus_ecal_timeError_branch != 0) mus_ecal_timeError();
	if (mus_e_em_branch != 0) mus_e_em();
	if (mus_e_emS9_branch != 0) mus_e_emS9();
	if (mus_e_had_branch != 0) mus_e_had();
	if (mus_e_hadS9_branch != 0) mus_e_hadS9();
	if (mus_e_ho_branch != 0) mus_e_ho();
	if (mus_e_hoS9_branch != 0) mus_e_hoS9();
	if (mus_emMax_branch != 0) mus_emMax();
	if (mus_emS25_branch != 0) mus_emS25();
	if (mus_etaErr_branch != 0) mus_etaErr();
	if (mus_gfit_chi2_branch != 0) mus_gfit_chi2();
	if (mus_gfit_d0_branch != 0) mus_gfit_d0();
	if (mus_gfit_d0Err_branch != 0) mus_gfit_d0Err();
	if (mus_gfit_d0corr_branch != 0) mus_gfit_d0corr();
	if (mus_gfit_d0corrPhi_branch != 0) mus_gfit_d0corrPhi();
	if (mus_gfit_d0phiCov_branch != 0) mus_gfit_d0phiCov();
	if (mus_gfit_ndof_branch != 0) mus_gfit_ndof();
	if (mus_gfit_qoverp_branch != 0) mus_gfit_qoverp();
	if (mus_gfit_qoverpError_branch != 0) mus_gfit_qoverpError();
	if (mus_gfit_z0_branch != 0) mus_gfit_z0();
	if (mus_gfit_z0Err_branch != 0) mus_gfit_z0Err();
	if (mus_gfit_z0corr_branch != 0) mus_gfit_z0corr();
	if (mus_glbKink_branch != 0) mus_glbKink();
	if (mus_glbTrackProbability_branch != 0) mus_glbTrackProbability();
	if (mus_globalDeltaEtaPhi_branch != 0) mus_globalDeltaEtaPhi();
	if (mus_hadMax_branch != 0) mus_hadMax();
	if (mus_hcal_time_branch != 0) mus_hcal_time();
	if (mus_hcal_timeError_branch != 0) mus_hcal_timeError();
	if (mus_ip3d_branch != 0) mus_ip3d();
	if (mus_ip3derr_branch != 0) mus_ip3derr();
	if (mus_iso03_emEt_branch != 0) mus_iso03_emEt();
	if (mus_iso03_hadEt_branch != 0) mus_iso03_hadEt();
	if (mus_iso03_hoEt_branch != 0) mus_iso03_hoEt();
	if (mus_iso03_sumPt_branch != 0) mus_iso03_sumPt();
	if (mus_iso05_emEt_branch != 0) mus_iso05_emEt();
	if (mus_iso05_hadEt_branch != 0) mus_iso05_hadEt();
	if (mus_iso05_hoEt_branch != 0) mus_iso05_hoEt();
	if (mus_iso05_sumPt_branch != 0) mus_iso05_sumPt();
	if (mus_isoMeanDRR03_pf_ChargedHadronPt_branch != 0) mus_isoMeanDRR03_pf_ChargedHadronPt();
	if (mus_isoMeanDRR03_pf_ChargedParticlePt_branch != 0) mus_isoMeanDRR03_pf_ChargedParticlePt();
	if (mus_isoMeanDRR03_pf_NeutralHadronEt_branch != 0) mus_isoMeanDRR03_pf_NeutralHadronEt();
	if (mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_branch != 0) mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold();
	if (mus_isoMeanDRR03_pf_PUPt_branch != 0) mus_isoMeanDRR03_pf_PUPt();
	if (mus_isoMeanDRR03_pf_PhotonEt_branch != 0) mus_isoMeanDRR03_pf_PhotonEt();
	if (mus_isoMeanDRR03_pf_PhotonEtHighThreshold_branch != 0) mus_isoMeanDRR03_pf_PhotonEtHighThreshold();
	if (mus_isoMeanDRR04_pf_ChargedHadronPt_branch != 0) mus_isoMeanDRR04_pf_ChargedHadronPt();
	if (mus_isoMeanDRR04_pf_ChargedParticlePt_branch != 0) mus_isoMeanDRR04_pf_ChargedParticlePt();
	if (mus_isoMeanDRR04_pf_NeutralHadronEt_branch != 0) mus_isoMeanDRR04_pf_NeutralHadronEt();
	if (mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_branch != 0) mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold();
	if (mus_isoMeanDRR04_pf_PUPt_branch != 0) mus_isoMeanDRR04_pf_PUPt();
	if (mus_isoMeanDRR04_pf_PhotonEt_branch != 0) mus_isoMeanDRR04_pf_PhotonEt();
	if (mus_isoMeanDRR04_pf_PhotonEtHighThreshold_branch != 0) mus_isoMeanDRR04_pf_PhotonEtHighThreshold();
	if (mus_isoR03_pf_ChargedHadronPt_branch != 0) mus_isoR03_pf_ChargedHadronPt();
	if (mus_isoR03_pf_ChargedParticlePt_branch != 0) mus_isoR03_pf_ChargedParticlePt();
	if (mus_isoR03_pf_NeutralHadronEt_branch != 0) mus_isoR03_pf_NeutralHadronEt();
	if (mus_isoR03_pf_NeutralHadronEtHighThreshold_branch != 0) mus_isoR03_pf_NeutralHadronEtHighThreshold();
	if (mus_isoR03_pf_PUPt_branch != 0) mus_isoR03_pf_PUPt();
	if (mus_isoR03_pf_PhotonEt_branch != 0) mus_isoR03_pf_PhotonEt();
	if (mus_isoR03_pf_PhotonEtHighThreshold_branch != 0) mus_isoR03_pf_PhotonEtHighThreshold();
	if (mus_isoR04_pf_ChargedHadronPt_branch != 0) mus_isoR04_pf_ChargedHadronPt();
	if (mus_isoR04_pf_ChargedParticlePt_branch != 0) mus_isoR04_pf_ChargedParticlePt();
	if (mus_isoR04_pf_NeutralHadronEt_branch != 0) mus_isoR04_pf_NeutralHadronEt();
	if (mus_isoR04_pf_NeutralHadronEtHighThreshold_branch != 0) mus_isoR04_pf_NeutralHadronEtHighThreshold();
	if (mus_isoR04_pf_PUPt_branch != 0) mus_isoR04_pf_PUPt();
	if (mus_isoR04_pf_PhotonEt_branch != 0) mus_isoR04_pf_PhotonEt();
	if (mus_isoR04_pf_PhotonEtHighThreshold_branch != 0) mus_isoR04_pf_PhotonEtHighThreshold();
	if (mus_isoSumDRR03_pf_ChargedHadronPt_branch != 0) mus_isoSumDRR03_pf_ChargedHadronPt();
	if (mus_isoSumDRR03_pf_ChargedParticlePt_branch != 0) mus_isoSumDRR03_pf_ChargedParticlePt();
	if (mus_isoSumDRR03_pf_NeutralHadronEt_branch != 0) mus_isoSumDRR03_pf_NeutralHadronEt();
	if (mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_branch != 0) mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold();
	if (mus_isoSumDRR03_pf_PUPt_branch != 0) mus_isoSumDRR03_pf_PUPt();
	if (mus_isoSumDRR03_pf_PhotonEt_branch != 0) mus_isoSumDRR03_pf_PhotonEt();
	if (mus_isoSumDRR03_pf_PhotonEtHighThreshold_branch != 0) mus_isoSumDRR03_pf_PhotonEtHighThreshold();
	if (mus_isoSumDRR04_pf_ChargedHadronPt_branch != 0) mus_isoSumDRR04_pf_ChargedHadronPt();
	if (mus_isoSumDRR04_pf_ChargedParticlePt_branch != 0) mus_isoSumDRR04_pf_ChargedParticlePt();
	if (mus_isoSumDRR04_pf_NeutralHadronEt_branch != 0) mus_isoSumDRR04_pf_NeutralHadronEt();
	if (mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_branch != 0) mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold();
	if (mus_isoSumDRR04_pf_PUPt_branch != 0) mus_isoSumDRR04_pf_PUPt();
	if (mus_isoSumDRR04_pf_PhotonEt_branch != 0) mus_isoSumDRR04_pf_PhotonEt();
	if (mus_isoSumDRR04_pf_PhotonEtHighThreshold_branch != 0) mus_isoSumDRR04_pf_PhotonEtHighThreshold();
	if (mus_iso_ecalvetoDep_branch != 0) mus_iso_ecalvetoDep();
	if (mus_iso_hcalvetoDep_branch != 0) mus_iso_hcalvetoDep();
	if (mus_iso_hovetoDep_branch != 0) mus_iso_hovetoDep();
	if (mus_iso_trckvetoDep_branch != 0) mus_iso_trckvetoDep();
	if (mus_localDistance_branch != 0) mus_localDistance();
	if (mus_mc_patMatch_dr_branch != 0) mus_mc_patMatch_dr();
	if (mus_ndof_branch != 0) mus_ndof();
	if (mus_overlapCompat_branch != 0) mus_overlapCompat();
	if (mus_phiErr_branch != 0) mus_phiErr();
	if (mus_ptErr_branch != 0) mus_ptErr();
	if (mus_qoverp_branch != 0) mus_qoverp();
	if (mus_qoverpError_branch != 0) mus_qoverpError();
	if (mus_segmCompatibility_branch != 0) mus_segmCompatibility();
	if (mus_staRelChi2_branch != 0) mus_staRelChi2();
	if (mus_sta_chi2_branch != 0) mus_sta_chi2();
	if (mus_sta_d0_branch != 0) mus_sta_d0();
	if (mus_sta_d0Err_branch != 0) mus_sta_d0Err();
	if (mus_sta_d0corr_branch != 0) mus_sta_d0corr();
	if (mus_sta_d0corrPhi_branch != 0) mus_sta_d0corrPhi();
	if (mus_sta_d0phiCov_branch != 0) mus_sta_d0phiCov();
	if (mus_sta_ndof_branch != 0) mus_sta_ndof();
	if (mus_sta_qoverp_branch != 0) mus_sta_qoverp();
	if (mus_sta_qoverpError_branch != 0) mus_sta_qoverpError();
	if (mus_sta_z0_branch != 0) mus_sta_z0();
	if (mus_sta_z0Err_branch != 0) mus_sta_z0Err();
	if (mus_sta_z0corr_branch != 0) mus_sta_z0corr();
	if (mus_timeAtIpInOut_branch != 0) mus_timeAtIpInOut();
	if (mus_timeAtIpInOutErr_branch != 0) mus_timeAtIpInOutErr();
	if (mus_timeAtIpOutIn_branch != 0) mus_timeAtIpOutIn();
	if (mus_timeAtIpOutInErr_branch != 0) mus_timeAtIpOutInErr();
	if (mus_timeCompat_branch != 0) mus_timeCompat();
	if (mus_tower_branch != 0) mus_tower();
	if (mus_towerS9_branch != 0) mus_towerS9();
	if (mus_trkKink_branch != 0) mus_trkKink();
	if (mus_trkRelChi2_branch != 0) mus_trkRelChi2();
	if (mus_vertexCompat_branch != 0) mus_vertexCompat();
	if (mus_vertexphi_branch != 0) mus_vertexphi();
	if (mus_z0_branch != 0) mus_z0();
	if (mus_z0Err_branch != 0) mus_z0Err();
	if (mus_z0corr_branch != 0) mus_z0corr();
	if (pfcands_dz_branch != 0) pfcands_dz();
	if (pfjets_area_branch != 0) pfjets_area();
	if (pfjets_chargedEmE_branch != 0) pfjets_chargedEmE();
	if (pfjets_chargedHadronE_branch != 0) pfjets_chargedHadronE();
	if (pfjets_combinedSecondaryVertexBJetTag_branch != 0) pfjets_combinedSecondaryVertexBJetTag();
	if (pfjets_cor_branch != 0) pfjets_cor();
	if (pfjets_corL1Fast_branch != 0) pfjets_corL1Fast();
	if (pfjets_corL1FastL2L3_branch != 0) pfjets_corL1FastL2L3();
	if (pfjets_corL1FastL2L3residual_branch != 0) pfjets_corL1FastL2L3residual();
	if (pfjets_corL1L2L3_branch != 0) pfjets_corL1L2L3();
	if (pfjets_electronE_branch != 0) pfjets_electronE();
	if (pfjets_hfEmE_branch != 0) pfjets_hfEmE();
	if (pfjets_hfHadronE_branch != 0) pfjets_hfHadronE();
	if (pfjets_jetBProbabilityBJetTag_branch != 0) pfjets_jetBProbabilityBJetTag();
	if (pfjets_jetProbabilityBJetTag_branch != 0) pfjets_jetProbabilityBJetTag();
	if (pfjets_muonE_branch != 0) pfjets_muonE();
	if (pfjets_neutralEmE_branch != 0) pfjets_neutralEmE();
	if (pfjets_neutralHadronE_branch != 0) pfjets_neutralHadronE();
	if (pfjets_photonE_branch != 0) pfjets_photonE();
	if (pfjets_pileupJetId_branch != 0) pfjets_pileupJetId();
	if (pfjets_simpleSecondaryVertexHighEffBJetTag_branch != 0) pfjets_simpleSecondaryVertexHighEffBJetTag();
	if (pfjets_simpleSecondaryVertexHighPurBJetTags_branch != 0) pfjets_simpleSecondaryVertexHighPurBJetTags();
	if (pfjets_trackCountingHighEffBJetTag_branch != 0) pfjets_trackCountingHighEffBJetTag();
	if (pfjets_trackCountingHighPurBJetTag_branch != 0) pfjets_trackCountingHighPurBJetTag();
	if (taus_pf_againstElectronDeadECAL_branch != 0) taus_pf_againstElectronDeadECAL();
	if (taus_pf_againstElectronLoose_branch != 0) taus_pf_againstElectronLoose();
	if (taus_pf_againstElectronMedium_branch != 0) taus_pf_againstElectronMedium();
	if (taus_pf_againstElectronTight_branch != 0) taus_pf_againstElectronTight();
	if (taus_pf_againstMuonLoose_branch != 0) taus_pf_againstMuonLoose();
	if (taus_pf_againstMuonLoose2_branch != 0) taus_pf_againstMuonLoose2();
	if (taus_pf_againstMuonLoose3_branch != 0) taus_pf_againstMuonLoose3();
	if (taus_pf_againstMuonMedium_branch != 0) taus_pf_againstMuonMedium();
	if (taus_pf_againstMuonMedium2_branch != 0) taus_pf_againstMuonMedium2();
	if (taus_pf_againstMuonTight_branch != 0) taus_pf_againstMuonTight();
	if (taus_pf_againstMuonTight2_branch != 0) taus_pf_againstMuonTight2();
	if (taus_pf_againstMuonTight3_branch != 0) taus_pf_againstMuonTight3();
	if (taus_pf_byCombinedIsolationDeltaBetaCorrRaw_branch != 0) taus_pf_byCombinedIsolationDeltaBetaCorrRaw();
	if (taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_branch != 0) taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits();
	if (taus_pf_byDecayModeFinding_branch != 0) taus_pf_byDecayModeFinding();
	if (taus_pf_byLooseCombinedIsolationDeltaBetaCorr_branch != 0) taus_pf_byLooseCombinedIsolationDeltaBetaCorr();
	if (taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_branch != 0) taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits();
	if (taus_pf_byMediumCombinedIsolationDeltaBetaCorr_branch != 0) taus_pf_byMediumCombinedIsolationDeltaBetaCorr();
	if (taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_branch != 0) taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits();
	if (taus_pf_byTightCombinedIsolationDeltaBetaCorr_branch != 0) taus_pf_byTightCombinedIsolationDeltaBetaCorr();
	if (taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_branch != 0) taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits();
	if (taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_branch != 0) taus_pf_byVLooseCombinedIsolationDeltaBetaCorr();
	if (photons_e1x5_branch != 0) photons_e1x5();
	if (photons_e2x5Max_branch != 0) photons_e2x5Max();
	if (photons_e3x3_branch != 0) photons_e3x3();
	if (photons_e5x5_branch != 0) photons_e5x5();
	if (photons_ecalIso03_branch != 0) photons_ecalIso03();
	if (photons_ecalIso04_branch != 0) photons_ecalIso04();
	if (photons_full3x3_e3x3_branch != 0) photons_full3x3_e3x3();
	if (photons_full5x5_e1x5_branch != 0) photons_full5x5_e1x5();
	if (photons_full5x5_e2x5Max_branch != 0) photons_full5x5_e2x5Max();
	if (photons_full5x5_e5x5_branch != 0) photons_full5x5_e5x5();
	if (photons_full5x5_hOverE_branch != 0) photons_full5x5_hOverE();
	if (photons_full5x5_sigmaEtaEta_branch != 0) photons_full5x5_sigmaEtaEta();
	if (photons_full5x5_sigmaIEtaIEta_branch != 0) photons_full5x5_sigmaIEtaIEta();
	if (photons_hOverE_branch != 0) photons_hOverE();
	if (photons_hcalDepth1TowerSumEtBcConeDR03_branch != 0) photons_hcalDepth1TowerSumEtBcConeDR03();
	if (photons_hcalDepth1TowerSumEtBcConeDR04_branch != 0) photons_hcalDepth1TowerSumEtBcConeDR04();
	if (photons_hcalDepth2TowerSumEtBcConeDR03_branch != 0) photons_hcalDepth2TowerSumEtBcConeDR03();
	if (photons_hcalDepth2TowerSumEtBcConeDR04_branch != 0) photons_hcalDepth2TowerSumEtBcConeDR04();
	if (photons_hcalIso03_branch != 0) photons_hcalIso03();
	if (photons_hcalIso04_branch != 0) photons_hcalIso04();
	if (photons_hcalTowerSumEtBcConeDR03_branch != 0) photons_hcalTowerSumEtBcConeDR03();
	if (photons_hcalTowerSumEtBcConeDR04_branch != 0) photons_hcalTowerSumEtBcConeDR04();
	if (photons_ntkIsoHollow03_branch != 0) photons_ntkIsoHollow03();
	if (photons_ntkIsoHollow04_branch != 0) photons_ntkIsoHollow04();
	if (photons_ntkIsoSolid03_branch != 0) photons_ntkIsoSolid03();
	if (photons_ntkIsoSolid04_branch != 0) photons_ntkIsoSolid04();
	if (photons_sigmaEtaEta_branch != 0) photons_sigmaEtaEta();
	if (photons_sigmaIEtaIEta_branch != 0) photons_sigmaIEtaIEta();
	if (photons_tkIsoHollow03_branch != 0) photons_tkIsoHollow03();
	if (photons_tkIsoHollow04_branch != 0) photons_tkIsoHollow04();
	if (photons_tkIsoSolid03_branch != 0) photons_tkIsoSolid03();
	if (photons_tkIsoSolid04_branch != 0) photons_tkIsoSolid04();
	if (puInfo_trueNumInteractions_branch != 0) puInfo_trueNumInteractions();
	if (convs_chi2_branch != 0) convs_chi2();
	if (convs_dl_branch != 0) convs_dl();
	if (convs_ndof_branch != 0) convs_ndof();
	if (vtxs_chi2_branch != 0) vtxs_chi2();
	if (vtxs_ndof_branch != 0) vtxs_ndof();
	if (vtxs_sumpt_branch != 0) vtxs_sumpt();
	if (vtxs_xError_branch != 0) vtxs_xError();
	if (vtxs_yError_branch != 0) vtxs_yError();
	if (vtxs_zError_branch != 0) vtxs_zError();
	if (els_convs_dcot_branch != 0) els_convs_dcot();
	if (els_convs_dist_branch != 0) els_convs_dist();
	if (els_convs_radius_branch != 0) els_convs_radius();
	if (mus_stationShowerDeltaR_branch != 0) mus_stationShowerDeltaR();
	if (mus_stationShowerSizeT_branch != 0) mus_stationShowerSizeT();
	if (puInfo_instLumi_branch != 0) puInfo_instLumi();
	if (vtxs_covMatrix_branch != 0) vtxs_covMatrix();
	if (evt_bsType_branch != 0) evt_bsType();
	if (evt_bunchCrossing_branch != 0) evt_bunchCrossing();
	if (evt_experimentType_branch != 0) evt_experimentType();
	if (evt_isRealData_branch != 0) evt_isRealData();
	if (evt_orbitNumber_branch != 0) evt_orbitNumber();
	if (evt_storeNumber_branch != 0) evt_storeNumber();
	if (hcalnoise_GetRecHitCount_branch != 0) hcalnoise_GetRecHitCount();
	if (hcalnoise_GetRecHitCount15_branch != 0) hcalnoise_GetRecHitCount15();
	if (hcalnoise_maxHPDHits_branch != 0) hcalnoise_maxHPDHits();
	if (hcalnoise_maxHPDNoOtherHits_branch != 0) hcalnoise_maxHPDNoOtherHits();
	if (hcalnoise_maxRBXHits_branch != 0) hcalnoise_maxRBXHits();
	if (hcalnoise_maxZeros_branch != 0) hcalnoise_maxZeros();
	if (hcalnoise_noiseFilterStatus_branch != 0) hcalnoise_noiseFilterStatus();
	if (hcalnoise_noiseType_branch != 0) hcalnoise_noiseType();
	if (hcalnoise_num10GeVHits_branch != 0) hcalnoise_num10GeVHits();
	if (hcalnoise_num25GeVHits_branch != 0) hcalnoise_num25GeVHits();
	if (hcalnoise_numFlatNoiseChannels_branch != 0) hcalnoise_numFlatNoiseChannels();
	if (hcalnoise_numIsolatedNoiseChannels_branch != 0) hcalnoise_numIsolatedNoiseChannels();
	if (hcalnoise_numProblematicRBXs_branch != 0) hcalnoise_numProblematicRBXs();
	if (hcalnoise_numSpikeNoiseChannels_branch != 0) hcalnoise_numSpikeNoiseChannels();
	if (hcalnoise_numTS4TS5NoiseChannels_branch != 0) hcalnoise_numTS4TS5NoiseChannels();
	if (hcalnoise_numTriangleNoiseChannels_branch != 0) hcalnoise_numTriangleNoiseChannels();
	if (hcalnoise_passHighLevelNoiseFilter_branch != 0) hcalnoise_passHighLevelNoiseFilter();
	if (hcalnoise_passLooseNoiseFilter_branch != 0) hcalnoise_passLooseNoiseFilter();
	if (hcalnoise_passTightNoiseFilter_branch != 0) hcalnoise_passTightNoiseFilter();
	if (pdfinfo_id1_branch != 0) pdfinfo_id1();
	if (pdfinfo_id2_branch != 0) pdfinfo_id2();
	if (els_mc3_id_branch != 0) els_mc3_id();
	if (els_mc3idx_branch != 0) els_mc3idx();
	if (els_mc3_motherid_branch != 0) els_mc3_motherid();
	if (els_mc3_motheridx_branch != 0) els_mc3_motheridx();
	if (els_mc_id_branch != 0) els_mc_id();
	if (els_mcidx_branch != 0) els_mcidx();
	if (els_mc_motherid_branch != 0) els_mc_motherid();
	if (mus_mc3_id_branch != 0) mus_mc3_id();
	if (mus_mc3idx_branch != 0) mus_mc3idx();
	if (mus_mc3_motherid_branch != 0) mus_mc3_motherid();
	if (mus_mc3_motheridx_branch != 0) mus_mc3_motheridx();
	if (mus_mc_id_branch != 0) mus_mc_id();
	if (mus_mcidx_branch != 0) mus_mcidx();
	if (mus_mc_motherid_branch != 0) mus_mc_motherid();
	if (pfjets_mc3_id_branch != 0) pfjets_mc3_id();
	if (pfjets_mc3idx_branch != 0) pfjets_mc3idx();
	if (pfjets_mc_gpidx_branch != 0) pfjets_mc_gpidx();
	if (pfjets_mc_id_branch != 0) pfjets_mc_id();
	if (pfjets_mcidx_branch != 0) pfjets_mcidx();
	if (pfjets_mc_motherid_branch != 0) pfjets_mc_motherid();
	if (photons_mc3_id_branch != 0) photons_mc3_id();
	if (photons_mc3idx_branch != 0) photons_mc3idx();
	if (photons_mc3_motherid_branch != 0) photons_mc3_motherid();
	if (photons_mc3_motheridx_branch != 0) photons_mc3_motheridx();
	if (photons_mc_id_branch != 0) photons_mc_id();
	if (photons_mcidx_branch != 0) photons_mcidx();
	if (photons_mc_motherid_branch != 0) photons_mc_motherid();
	if (els_category_branch != 0) els_category();
	if (els_charge_branch != 0) els_charge();
	if (els_ckf_charge_branch != 0) els_ckf_charge();
	if (els_ckf_laywithmeas_branch != 0) els_ckf_laywithmeas();
	if (els_class_branch != 0) els_class();
	if (els_conv_delMissHits_branch != 0) els_conv_delMissHits();
	if (els_conv_flag_branch != 0) els_conv_flag();
	if (els_conv_gsftkidx_branch != 0) els_conv_gsftkidx();
	if (els_conv_old_delMissHits_branch != 0) els_conv_old_delMissHits();
	if (els_conv_old_flag_branch != 0) els_conv_old_flag();
	if (els_conv_old_gsftkidx_branch != 0) els_conv_old_gsftkidx();
	if (els_conv_old_tkidx_branch != 0) els_conv_old_tkidx();
	if (els_conv_tkidx_branch != 0) els_conv_tkidx();
	if (els_exp_innerlayers_branch != 0) els_exp_innerlayers();
	if (els_exp_outerlayers_branch != 0) els_exp_outerlayers();
	if (els_fiduciality_branch != 0) els_fiduciality();
	if (els_gsftrkidx_branch != 0) els_gsftrkidx();
	if (els_layer1_det_branch != 0) els_layer1_det();
	if (els_layer1_layer_branch != 0) els_layer1_layer();
	if (els_layer1_sizerphi_branch != 0) els_layer1_sizerphi();
	if (els_layer1_sizerz_branch != 0) els_layer1_sizerz();
	if (els_lostHits_branch != 0) els_lostHits();
	if (els_lost_pixelhits_branch != 0) els_lost_pixelhits();
	if (els_mc_patMatch_id_branch != 0) els_mc_patMatch_id();
	if (els_nSeed_branch != 0) els_nSeed();
	if (els_nlayers_branch != 0) els_nlayers();
	if (els_nlayers3D_branch != 0) els_nlayers3D();
	if (els_nlayersLost_branch != 0) els_nlayersLost();
	if (els_sccharge_branch != 0) els_sccharge();
	if (els_scindex_branch != 0) els_scindex();
	if (els_trk_charge_branch != 0) els_trk_charge();
	if (els_trkidx_branch != 0) els_trkidx();
	if (els_type_branch != 0) els_type();
	if (els_validHits_branch != 0) els_validHits();
	if (els_valid_pixelhits_branch != 0) els_valid_pixelhits();
	if (genps_id_branch != 0) genps_id();
	if (genps_id_mother_branch != 0) genps_id_mother();
	if (genps_status_branch != 0) genps_status();
	if (mus_algo_branch != 0) mus_algo();
	if (mus_charge_branch != 0) mus_charge();
	if (mus_ecal_rawId_branch != 0) mus_ecal_rawId();
	if (mus_exp_innerlayers_branch != 0) mus_exp_innerlayers();
	if (mus_exp_outerlayers_branch != 0) mus_exp_outerlayers();
	if (mus_gfit_algo_branch != 0) mus_gfit_algo();
	if (mus_gfit_etaErr_branch != 0) mus_gfit_etaErr();
	if (mus_gfit_exp_innerlayer_branch != 0) mus_gfit_exp_innerlayer();
	if (mus_gfit_exp_outerlayer_branch != 0) mus_gfit_exp_outerlayer();
	if (mus_gfit_lostHits_branch != 0) mus_gfit_lostHits();
	if (mus_gfit_lostPixelHits_branch != 0) mus_gfit_lostPixelHits();
	if (mus_gfit_nlayers_branch != 0) mus_gfit_nlayers();
	if (mus_gfit_nlayers3D_branch != 0) mus_gfit_nlayers3D();
	if (mus_gfit_nlayersLost_branch != 0) mus_gfit_nlayersLost();
	if (mus_gfit_phiErr_branch != 0) mus_gfit_phiErr();
	if (mus_gfit_ptErr_branch != 0) mus_gfit_ptErr();
	if (mus_gfit_qualityMask_branch != 0) mus_gfit_qualityMask();
	if (mus_gfit_trk_charge_branch != 0) mus_gfit_trk_charge();
	if (mus_gfit_validHits_branch != 0) mus_gfit_validHits();
	if (mus_gfit_validPixelHits_branch != 0) mus_gfit_validPixelHits();
	if (mus_gfit_validSTAHits_branch != 0) mus_gfit_validSTAHits();
	if (mus_gfit_validSiHits_branch != 0) mus_gfit_validSiHits();
	if (mus_goodmask_branch != 0) mus_goodmask();
	if (mus_hcal_rawId_branch != 0) mus_hcal_rawId();
	if (mus_iso03_njets_branch != 0) mus_iso03_njets();
	if (mus_iso03_ntrk_branch != 0) mus_iso03_ntrk();
	if (mus_iso05_ntrk_branch != 0) mus_iso05_ntrk();
	if (mus_lostHits_branch != 0) mus_lostHits();
	if (mus_lostPixelHits_branch != 0) mus_lostPixelHits();
	if (mus_mc_patMatch_id_branch != 0) mus_mc_patMatch_id();
	if (mus_muonBestTrackType_branch != 0) mus_muonBestTrackType();
	if (mus_nOverlaps_branch != 0) mus_nOverlaps();
	if (mus_nlayers_branch != 0) mus_nlayers();
	if (mus_nlayers3D_branch != 0) mus_nlayers3D();
	if (mus_nlayersLost_branch != 0) mus_nlayersLost();
	if (mus_nmatches_branch != 0) mus_nmatches();
	if (mus_numberOfMatchedStations_branch != 0) mus_numberOfMatchedStations();
	if (mus_overlap0_branch != 0) mus_overlap0();
	if (mus_overlap1_branch != 0) mus_overlap1();
	if (mus_pfcharge_branch != 0) mus_pfcharge();
	if (mus_pfparticleId_branch != 0) mus_pfparticleId();
	if (mus_pid_PFMuon_branch != 0) mus_pid_PFMuon();
	if (mus_pid_TM2DCompatibilityLoose_branch != 0) mus_pid_TM2DCompatibilityLoose();
	if (mus_pid_TM2DCompatibilityTight_branch != 0) mus_pid_TM2DCompatibilityTight();
	if (mus_pid_TMLastStationLoose_branch != 0) mus_pid_TMLastStationLoose();
	if (mus_pid_TMLastStationTight_branch != 0) mus_pid_TMLastStationTight();
	if (mus_qualityMask_branch != 0) mus_qualityMask();
	if (mus_sta_algo_branch != 0) mus_sta_algo();
	if (mus_sta_etaErr_branch != 0) mus_sta_etaErr();
	if (mus_sta_exp_innerlayer_branch != 0) mus_sta_exp_innerlayer();
	if (mus_sta_exp_outerlayer_branch != 0) mus_sta_exp_outerlayer();
	if (mus_sta_lostHits_branch != 0) mus_sta_lostHits();
	if (mus_sta_lostPixelHits_branch != 0) mus_sta_lostPixelHits();
	if (mus_sta_nlayers_branch != 0) mus_sta_nlayers();
	if (mus_sta_nlayers3D_branch != 0) mus_sta_nlayers3D();
	if (mus_sta_nlayersLost_branch != 0) mus_sta_nlayersLost();
	if (mus_sta_phiErr_branch != 0) mus_sta_phiErr();
	if (mus_sta_ptErr_branch != 0) mus_sta_ptErr();
	if (mus_sta_qualityMask_branch != 0) mus_sta_qualityMask();
	if (mus_sta_trk_charge_branch != 0) mus_sta_trk_charge();
	if (mus_sta_validHits_branch != 0) mus_sta_validHits();
	if (mus_sta_validPixelHits_branch != 0) mus_sta_validPixelHits();
	if (mus_timeDirection_branch != 0) mus_timeDirection();
	if (mus_timeNumStationsUsed_branch != 0) mus_timeNumStationsUsed();
	if (mus_trk_charge_branch != 0) mus_trk_charge();
	if (mus_trkidx_branch != 0) mus_trkidx();
	if (mus_type_branch != 0) mus_type();
	if (mus_validHits_branch != 0) mus_validHits();
	if (mus_validPixelHits_branch != 0) mus_validPixelHits();
	if (pfcands_charge_branch != 0) pfcands_charge();
	if (pfcands_particleId_branch != 0) pfcands_particleId();
	if (pfjets_chargedHadronMultiplicity_branch != 0) pfjets_chargedHadronMultiplicity();
	if (pfjets_chargedMultiplicity_branch != 0) pfjets_chargedMultiplicity();
	if (pfjets_electronMultiplicity_branch != 0) pfjets_electronMultiplicity();
	if (pfjets_hfEmMultiplicity_branch != 0) pfjets_hfEmMultiplicity();
	if (pfjets_hfHadronMultiplicity_branch != 0) pfjets_hfHadronMultiplicity();
	if (pfjets_muonMultiplicity_branch != 0) pfjets_muonMultiplicity();
	if (pfjets_neutralHadronMultiplicity_branch != 0) pfjets_neutralHadronMultiplicity();
	if (pfjets_neutralMultiplicity_branch != 0) pfjets_neutralMultiplicity();
	if (pfjets_partonFlavour_branch != 0) pfjets_partonFlavour();
	if (pfjets_photonMultiplicity_branch != 0) pfjets_photonMultiplicity();
	if (taus_pf_charge_branch != 0) taus_pf_charge();
	if (photons_fiduciality_branch != 0) photons_fiduciality();
	if (puInfo_bunchCrossing_branch != 0) puInfo_bunchCrossing();
	if (puInfo_nPUvertices_branch != 0) puInfo_nPUvertices();
	if (convs_algo_branch != 0) convs_algo();
	if (convs_isConverted_branch != 0) convs_isConverted();
	if (convs_quality_branch != 0) convs_quality();
	if (vtxs_isFake_branch != 0) vtxs_isFake();
	if (vtxs_isValid_branch != 0) vtxs_isValid();
	if (vtxs_tracksSize_branch != 0) vtxs_tracksSize();
	if (els_convs_delMissHits_branch != 0) els_convs_delMissHits();
	if (els_convs_flag_branch != 0) els_convs_flag();
	if (els_convs_gsftkidx_branch != 0) els_convs_gsftkidx();
	if (els_convs_tkidx_branch != 0) els_convs_tkidx();
	if (els_PFCand_idx_branch != 0) els_PFCand_idx();
	if (genps_lepdaughter_id_branch != 0) genps_lepdaughter_id();
	if (genps_lepdaughter_idx_branch != 0) genps_lepdaughter_idx();
	if (hlt_trigObjs_id_branch != 0) hlt_trigObjs_id();
	if (mus_nStationCorrelatedHits_branch != 0) mus_nStationCorrelatedHits();
	if (mus_nStationHits_branch != 0) mus_nStationHits();
	if (pfjets_pfcandIndicies_branch != 0) pfjets_pfcandIndicies();
	if (convs_nHitsBeforeVtx_branch != 0) convs_nHitsBeforeVtx();
	if (convs_tkalgo_branch != 0) convs_tkalgo();
	if (convs_tkidx_branch != 0) convs_tkidx();
	if (els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch != 0) els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version();
	if (els_HLT_Ele17_Ele8_LeadingLeg_version_branch != 0) els_HLT_Ele17_Ele8_LeadingLeg_version();
	if (els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch != 0) els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version();
	if (els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch != 0) els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version();
	if (els_HLT_Ele17_Ele8_TrailingLeg_version_branch != 0) els_HLT_Ele17_Ele8_TrailingLeg_version();
	if (els_HLT_Ele17_Ele8_version_branch != 0) els_HLT_Ele17_Ele8_version();
	if (els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch != 0) els_HLT_Ele20_SC4_Mass50_LeadingLeg_version();
	if (els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch != 0) els_HLT_Ele20_SC4_Mass50_TrailingLeg_version();
	if (els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch != 0) els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version();
	if (els_HLT_Ele27_WP80_version_branch != 0) els_HLT_Ele27_WP80_version();
	if (els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch != 0) els_HLT_Ele32_SC17_Mass50_LeadingLeg_version();
	if (els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch != 0) els_HLT_Ele32_SC17_Mass50_TrailingLeg_version();
	if (els_HLT_Mu17_Ele8_TrailingLeg_version_branch != 0) els_HLT_Mu17_Ele8_TrailingLeg_version();
	if (els_HLT_Mu17_Ele8_version_branch != 0) els_HLT_Mu17_Ele8_version();
	if (els_HLT_Mu8_Ele17_version_branch != 0) els_HLT_Mu8_Ele17_version();
	if (evt_nels_branch != 0) evt_nels();
	if (evt_detectorStatus_branch != 0) evt_detectorStatus();
	if (evt_event_branch != 0) evt_event();
	if (evt_lumiBlock_branch != 0) evt_lumiBlock();
	if (evt_run_branch != 0) evt_run();
	if (evt_ngenjetsNoMuNoNu_branch != 0) evt_ngenjetsNoMuNoNu();
	if (genps_signalProcessID_branch != 0) genps_signalProcessID();
	if (mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch != 0) mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version();
	if (mus_HLT_IsoMu24_eta2p1_version_branch != 0) mus_HLT_IsoMu24_eta2p1_version();
	if (mus_HLT_Mu17_Ele8_LeadingLeg_version_branch != 0) mus_HLT_Mu17_Ele8_LeadingLeg_version();
	if (mus_HLT_Mu17_Ele8_version_branch != 0) mus_HLT_Mu17_Ele8_version();
	if (mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch != 0) mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version();
	if (mus_HLT_Mu17_Mu8_LeadingLeg_version_branch != 0) mus_HLT_Mu17_Mu8_LeadingLeg_version();
	if (mus_HLT_Mu17_Mu8_TrailingLeg_version_branch != 0) mus_HLT_Mu17_Mu8_TrailingLeg_version();
	if (mus_HLT_Mu17_Mu8_version_branch != 0) mus_HLT_Mu17_Mu8_version();
	if (mus_HLT_Mu17_TkMu8_LeadingLeg_version_branch != 0) mus_HLT_Mu17_TkMu8_LeadingLeg_version();
	if (mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch != 0) mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version();
	if (mus_HLT_Mu17_TkMu8_TrailingLeg_version_branch != 0) mus_HLT_Mu17_TkMu8_TrailingLeg_version();
	if (mus_HLT_Mu17_TkMu8_version_branch != 0) mus_HLT_Mu17_TkMu8_version();
	if (mus_HLT_Mu8_Ele17_TrailingLeg_version_branch != 0) mus_HLT_Mu8_Ele17_TrailingLeg_version();
	if (mus_HLT_Mu8_Ele17_version_branch != 0) mus_HLT_Mu8_Ele17_version();
	if (evt_nphotons_branch != 0) evt_nphotons();
	if (evt_nvtxs_branch != 0) evt_nvtxs();
	if (els_HLT_Ele17_Ele8_branch != 0) els_HLT_Ele17_Ele8();
	if (els_HLT_Ele17_Ele8_L1sL1DoubleEG137_branch != 0) els_HLT_Ele17_Ele8_L1sL1DoubleEG137();
	if (els_HLT_Ele17_Ele8_LeadingLeg_branch != 0) els_HLT_Ele17_Ele8_LeadingLeg();
	if (els_HLT_Ele17_Ele8_Mass50_LeadingLeg_branch != 0) els_HLT_Ele17_Ele8_Mass50_LeadingLeg();
	if (els_HLT_Ele17_Ele8_Mass50_TrailingLeg_branch != 0) els_HLT_Ele17_Ele8_Mass50_TrailingLeg();
	if (els_HLT_Ele17_Ele8_TrailingLeg_branch != 0) els_HLT_Ele17_Ele8_TrailingLeg();
	if (els_HLT_Ele20_SC4_Mass50_LeadingLeg_branch != 0) els_HLT_Ele20_SC4_Mass50_LeadingLeg();
	if (els_HLT_Ele20_SC4_Mass50_TrailingLeg_branch != 0) els_HLT_Ele20_SC4_Mass50_TrailingLeg();
	if (els_HLT_Ele27_WP80_branch != 0) els_HLT_Ele27_WP80();
	if (els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_branch != 0) els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22();
	if (els_HLT_Ele32_SC17_Mass50_LeadingLeg_branch != 0) els_HLT_Ele32_SC17_Mass50_LeadingLeg();
	if (els_HLT_Ele32_SC17_Mass50_TrailingLeg_branch != 0) els_HLT_Ele32_SC17_Mass50_TrailingLeg();
	if (els_HLT_Mu17_Ele8_branch != 0) els_HLT_Mu17_Ele8();
	if (els_HLT_Mu17_Ele8_TrailingLeg_branch != 0) els_HLT_Mu17_Ele8_TrailingLeg();
	if (els_HLT_Mu8_Ele17_branch != 0) els_HLT_Mu8_Ele17();
	if (els_id2012ext_loose_branch != 0) els_id2012ext_loose();
	if (els_id2012ext_medium_branch != 0) els_id2012ext_medium();
	if (els_id2012ext_tight_branch != 0) els_id2012ext_tight();
	if (els_id2012ext_veto_branch != 0) els_id2012ext_veto();
	if (els_id2012_loose_branch != 0) els_id2012_loose();
	if (els_id2012_medium_branch != 0) els_id2012_medium();
	if (els_id2012_tight_branch != 0) els_id2012_tight();
	if (els_id2012_veto_branch != 0) els_id2012_veto();
	if (hlt_prescales_branch != 0) hlt_prescales();
	if (mus_HLT_IsoMu24_eta2p1_branch != 0) mus_HLT_IsoMu24_eta2p1();
	if (mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_branch != 0) mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1();
	if (mus_HLT_Mu17_Ele8_branch != 0) mus_HLT_Mu17_Ele8();
	if (mus_HLT_Mu17_Ele8_LeadingLeg_branch != 0) mus_HLT_Mu17_Ele8_LeadingLeg();
	if (mus_HLT_Mu17_Mu8_branch != 0) mus_HLT_Mu17_Mu8();
	if (mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_branch != 0) mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen();
	if (mus_HLT_Mu17_Mu8_LeadingLeg_branch != 0) mus_HLT_Mu17_Mu8_LeadingLeg();
	if (mus_HLT_Mu17_Mu8_TrailingLeg_branch != 0) mus_HLT_Mu17_Mu8_TrailingLeg();
	if (mus_HLT_Mu17_TkMu8_branch != 0) mus_HLT_Mu17_TkMu8();
	if (mus_HLT_Mu17_TkMu8_LeadingLeg_branch != 0) mus_HLT_Mu17_TkMu8_LeadingLeg();
	if (mus_HLT_Mu17_TkMu8_TrailingLeg_branch != 0) mus_HLT_Mu17_TkMu8_TrailingLeg();
	if (mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_branch != 0) mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered();
	if (mus_HLT_Mu8_Ele17_branch != 0) mus_HLT_Mu8_Ele17();
	if (mus_HLT_Mu8_Ele17_TrailingLeg_branch != 0) mus_HLT_Mu8_Ele17_TrailingLeg();
}

	TBits &hlt_bits()
	{
		if (not hlt_bits_isLoaded) {
			if (hlt_bits_branch != 0) {
				hlt_bits_branch->GetEntry(index);
			} else { }
			hlt_bits_isLoaded = true;
		}
		return hlt_bits_;
	}
	vector<TString> &evt_CMS2tag()
	{
		if (not evt_CMS2tag_isLoaded) {
			if (evt_CMS2tag_branch != 0) {
				evt_CMS2tag_branch->GetEntry(index);
			} else { }
			evt_CMS2tag_isLoaded = true;
		}
		return evt_CMS2tag_;
	}
	vector<TString> &evt_dataset()
	{
		if (not evt_dataset_isLoaded) {
			if (evt_dataset_branch != 0) {
				evt_dataset_branch->GetEntry(index);
			} else { }
			evt_dataset_isLoaded = true;
		}
		return evt_dataset_;
	}
	vector<TString> &hlt_trigNames()
	{
		if (not hlt_trigNames_isLoaded) {
			if (hlt_trigNames_branch != 0) {
				hlt_trigNames_branch->GetEntry(index);
			} else { }
			hlt_trigNames_isLoaded = true;
		}
		return hlt_trigNames_;
	}
	bool &hcalnoise_HasBadRBXTS4TS5()
	{
		if (not hcalnoise_HasBadRBXTS4TS5_isLoaded) {
			if (hcalnoise_HasBadRBXTS4TS5_branch != 0) {
				hcalnoise_HasBadRBXTS4TS5_branch->GetEntry(index);
			} else { }
			hcalnoise_HasBadRBXTS4TS5_isLoaded = true;
		}
		return hcalnoise_HasBadRBXTS4TS5_;
	}
	bool &evt_cscTightHaloId()
	{
		if (not evt_cscTightHaloId_isLoaded) {
			if (evt_cscTightHaloId_branch != 0) {
				evt_cscTightHaloId_branch->GetEntry(index);
			} else { }
			evt_cscTightHaloId_isLoaded = true;
		}
		return evt_cscTightHaloId_;
	}
	bool &evt_hbheFilter()
	{
		if (not evt_hbheFilter_isLoaded) {
			if (evt_hbheFilter_branch != 0) {
				evt_hbheFilter_branch->GetEntry(index);
			} else { }
			evt_hbheFilter_isLoaded = true;
		}
		return evt_hbheFilter_;
	}
	bool &filt_cscBeamHalo()
	{
		if (not filt_cscBeamHalo_isLoaded) {
			if (filt_cscBeamHalo_branch != 0) {
				filt_cscBeamHalo_branch->GetEntry(index);
			} else { }
			filt_cscBeamHalo_isLoaded = true;
		}
		return filt_cscBeamHalo_;
	}
	bool &filt_ecalLaser()
	{
		if (not filt_ecalLaser_isLoaded) {
			if (filt_ecalLaser_branch != 0) {
				filt_ecalLaser_branch->GetEntry(index);
			} else { }
			filt_ecalLaser_isLoaded = true;
		}
		return filt_ecalLaser_;
	}
	bool &filt_ecalTP()
	{
		if (not filt_ecalTP_isLoaded) {
			if (filt_ecalTP_branch != 0) {
				filt_ecalTP_branch->GetEntry(index);
			} else { }
			filt_ecalTP_isLoaded = true;
		}
		return filt_ecalTP_;
	}
	bool &filt_eeBadSc()
	{
		if (not filt_eeBadSc_isLoaded) {
			if (filt_eeBadSc_branch != 0) {
				filt_eeBadSc_branch->GetEntry(index);
			} else { }
			filt_eeBadSc_isLoaded = true;
		}
		return filt_eeBadSc_;
	}
	bool &filt_goodVertices()
	{
		if (not filt_goodVertices_isLoaded) {
			if (filt_goodVertices_branch != 0) {
				filt_goodVertices_branch->GetEntry(index);
			} else { }
			filt_goodVertices_isLoaded = true;
		}
		return filt_goodVertices_;
	}
	bool &filt_hbheNoise()
	{
		if (not filt_hbheNoise_isLoaded) {
			if (filt_hbheNoise_branch != 0) {
				filt_hbheNoise_branch->GetEntry(index);
			} else { }
			filt_hbheNoise_isLoaded = true;
		}
		return filt_hbheNoise_;
	}
	bool &filt_hcalLaser()
	{
		if (not filt_hcalLaser_isLoaded) {
			if (filt_hcalLaser_branch != 0) {
				filt_hcalLaser_branch->GetEntry(index);
			} else { }
			filt_hcalLaser_isLoaded = true;
		}
		return filt_hcalLaser_;
	}
	bool &filt_metfilter()
	{
		if (not filt_metfilter_isLoaded) {
			if (filt_metfilter_branch != 0) {
				filt_metfilter_branch->GetEntry(index);
			} else { }
			filt_metfilter_isLoaded = true;
		}
		return filt_metfilter_;
	}
	bool &filt_trackingFailure()
	{
		if (not filt_trackingFailure_isLoaded) {
			if (filt_trackingFailure_branch != 0) {
				filt_trackingFailure_branch->GetEntry(index);
			} else { }
			filt_trackingFailure_isLoaded = true;
		}
		return filt_trackingFailure_;
	}
	bool &filt_trkPOGFilters()
	{
		if (not filt_trkPOGFilters_isLoaded) {
			if (filt_trkPOGFilters_branch != 0) {
				filt_trkPOGFilters_branch->GetEntry(index);
			} else { }
			filt_trkPOGFilters_isLoaded = true;
		}
		return filt_trkPOGFilters_;
	}
	bool &filt_trkPOG_logErrorTooManyClusters()
	{
		if (not filt_trkPOG_logErrorTooManyClusters_isLoaded) {
			if (filt_trkPOG_logErrorTooManyClusters_branch != 0) {
				filt_trkPOG_logErrorTooManyClusters_branch->GetEntry(index);
			} else { }
			filt_trkPOG_logErrorTooManyClusters_isLoaded = true;
		}
		return filt_trkPOG_logErrorTooManyClusters_;
	}
	bool &filt_trkPOG_manystripclus53X()
	{
		if (not filt_trkPOG_manystripclus53X_isLoaded) {
			if (filt_trkPOG_manystripclus53X_branch != 0) {
				filt_trkPOG_manystripclus53X_branch->GetEntry(index);
			} else { }
			filt_trkPOG_manystripclus53X_isLoaded = true;
		}
		return filt_trkPOG_manystripclus53X_;
	}
	bool &filt_trkPOG_toomanystripclus53X()
	{
		if (not filt_trkPOG_toomanystripclus53X_isLoaded) {
			if (filt_trkPOG_toomanystripclus53X_branch != 0) {
				filt_trkPOG_toomanystripclus53X_branch->GetEntry(index);
			} else { }
			filt_trkPOG_toomanystripclus53X_isLoaded = true;
		}
		return filt_trkPOG_toomanystripclus53X_;
	}
	vector<bool> &els_conv_vtx_flag()
	{
		if (not els_conv_vtx_flag_isLoaded) {
			if (els_conv_vtx_flag_branch != 0) {
				els_conv_vtx_flag_branch->GetEntry(index);
			} else { }
			els_conv_vtx_flag_isLoaded = true;
		}
		return els_conv_vtx_flag_;
	}
	vector<bool> &els_isGsfCtfScPixChargeConsistent()
	{
		if (not els_isGsfCtfScPixChargeConsistent_isLoaded) {
			if (els_isGsfCtfScPixChargeConsistent_branch != 0) {
				els_isGsfCtfScPixChargeConsistent_branch->GetEntry(index);
			} else { }
			els_isGsfCtfScPixChargeConsistent_isLoaded = true;
		}
		return els_isGsfCtfScPixChargeConsistent_;
	}
	vector<bool> &els_passingMvaPreselection()
	{
		if (not els_passingMvaPreselection_isLoaded) {
			if (els_passingMvaPreselection_branch != 0) {
				els_passingMvaPreselection_branch->GetEntry(index);
			} else { }
			els_passingMvaPreselection_isLoaded = true;
		}
		return els_passingMvaPreselection_;
	}
	vector<bool> &els_passingPflowPreselection()
	{
		if (not els_passingPflowPreselection_isLoaded) {
			if (els_passingPflowPreselection_branch != 0) {
				els_passingPflowPreselection_branch->GetEntry(index);
			} else { }
			els_passingPflowPreselection_isLoaded = true;
		}
		return els_passingPflowPreselection_;
	}
	vector<bool> &mus_isRPCMuon()
	{
		if (not mus_isRPCMuon_isLoaded) {
			if (mus_isRPCMuon_branch != 0) {
				mus_isRPCMuon_branch->GetEntry(index);
			} else { }
			mus_isRPCMuon_isLoaded = true;
		}
		return mus_isRPCMuon_;
	}
	vector<bool> &mus_tightMatch()
	{
		if (not mus_tightMatch_isLoaded) {
			if (mus_tightMatch_branch != 0) {
				mus_tightMatch_branch->GetEntry(index);
			} else { }
			mus_tightMatch_isLoaded = true;
		}
		return mus_tightMatch_;
	}
	vector<bool> &mus_updatedSta()
	{
		if (not mus_updatedSta_isLoaded) {
			if (mus_updatedSta_branch != 0) {
				mus_updatedSta_branch->GetEntry(index);
			} else { }
			mus_updatedSta_isLoaded = true;
		}
		return mus_updatedSta_;
	}
	vector<bool> &photons_haspixelSeed()
	{
		if (not photons_haspixelSeed_isLoaded) {
			if (photons_haspixelSeed_branch != 0) {
				photons_haspixelSeed_branch->GetEntry(index);
			} else { }
			photons_haspixelSeed_isLoaded = true;
		}
		return photons_haspixelSeed_;
	}
	float &evt_bs_Xwidth()
	{
		if (not evt_bs_Xwidth_isLoaded) {
			if (evt_bs_Xwidth_branch != 0) {
				evt_bs_Xwidth_branch->GetEntry(index);
			} else { }
			evt_bs_Xwidth_isLoaded = true;
		}
		return evt_bs_Xwidth_;
	}
	float &evt_bs_XwidthErr()
	{
		if (not evt_bs_XwidthErr_isLoaded) {
			if (evt_bs_XwidthErr_branch != 0) {
				evt_bs_XwidthErr_branch->GetEntry(index);
			} else { }
			evt_bs_XwidthErr_isLoaded = true;
		}
		return evt_bs_XwidthErr_;
	}
	float &evt_bs_Ywidth()
	{
		if (not evt_bs_Ywidth_isLoaded) {
			if (evt_bs_Ywidth_branch != 0) {
				evt_bs_Ywidth_branch->GetEntry(index);
			} else { }
			evt_bs_Ywidth_isLoaded = true;
		}
		return evt_bs_Ywidth_;
	}
	float &evt_bs_YwidthErr()
	{
		if (not evt_bs_YwidthErr_isLoaded) {
			if (evt_bs_YwidthErr_branch != 0) {
				evt_bs_YwidthErr_branch->GetEntry(index);
			} else { }
			evt_bs_YwidthErr_isLoaded = true;
		}
		return evt_bs_YwidthErr_;
	}
	float &evt_bs_dxdz()
	{
		if (not evt_bs_dxdz_isLoaded) {
			if (evt_bs_dxdz_branch != 0) {
				evt_bs_dxdz_branch->GetEntry(index);
			} else { }
			evt_bs_dxdz_isLoaded = true;
		}
		return evt_bs_dxdz_;
	}
	float &evt_bs_dxdzErr()
	{
		if (not evt_bs_dxdzErr_isLoaded) {
			if (evt_bs_dxdzErr_branch != 0) {
				evt_bs_dxdzErr_branch->GetEntry(index);
			} else { }
			evt_bs_dxdzErr_isLoaded = true;
		}
		return evt_bs_dxdzErr_;
	}
	float &evt_bs_dydz()
	{
		if (not evt_bs_dydz_isLoaded) {
			if (evt_bs_dydz_branch != 0) {
				evt_bs_dydz_branch->GetEntry(index);
			} else { }
			evt_bs_dydz_isLoaded = true;
		}
		return evt_bs_dydz_;
	}
	float &evt_bs_dydzErr()
	{
		if (not evt_bs_dydzErr_isLoaded) {
			if (evt_bs_dydzErr_branch != 0) {
				evt_bs_dydzErr_branch->GetEntry(index);
			} else { }
			evt_bs_dydzErr_isLoaded = true;
		}
		return evt_bs_dydzErr_;
	}
	float &evt_bs_sigmaZ()
	{
		if (not evt_bs_sigmaZ_isLoaded) {
			if (evt_bs_sigmaZ_branch != 0) {
				evt_bs_sigmaZ_branch->GetEntry(index);
			} else { }
			evt_bs_sigmaZ_isLoaded = true;
		}
		return evt_bs_sigmaZ_;
	}
	float &evt_bs_sigmaZErr()
	{
		if (not evt_bs_sigmaZErr_isLoaded) {
			if (evt_bs_sigmaZErr_branch != 0) {
				evt_bs_sigmaZErr_branch->GetEntry(index);
			} else { }
			evt_bs_sigmaZErr_isLoaded = true;
		}
		return evt_bs_sigmaZErr_;
	}
	float &evt_bs_xErr()
	{
		if (not evt_bs_xErr_isLoaded) {
			if (evt_bs_xErr_branch != 0) {
				evt_bs_xErr_branch->GetEntry(index);
			} else { }
			evt_bs_xErr_isLoaded = true;
		}
		return evt_bs_xErr_;
	}
	float &evt_bs_yErr()
	{
		if (not evt_bs_yErr_isLoaded) {
			if (evt_bs_yErr_branch != 0) {
				evt_bs_yErr_branch->GetEntry(index);
			} else { }
			evt_bs_yErr_isLoaded = true;
		}
		return evt_bs_yErr_;
	}
	float &evt_bs_zErr()
	{
		if (not evt_bs_zErr_isLoaded) {
			if (evt_bs_zErr_branch != 0) {
				evt_bs_zErr_branch->GetEntry(index);
			} else { }
			evt_bs_zErr_isLoaded = true;
		}
		return evt_bs_zErr_;
	}
	float &evt_bField()
	{
		if (not evt_bField_isLoaded) {
			if (evt_bField_branch != 0) {
				evt_bField_branch->GetEntry(index);
			} else { }
			evt_bField_isLoaded = true;
		}
		return evt_bField_;
	}
	float &evt_fixgrid_all_rho()
	{
		if (not evt_fixgrid_all_rho_isLoaded) {
			if (evt_fixgrid_all_rho_branch != 0) {
				evt_fixgrid_all_rho_branch->GetEntry(index);
			} else { }
			evt_fixgrid_all_rho_isLoaded = true;
		}
		return evt_fixgrid_all_rho_;
	}
	float &evt_fixgridfastjet_allcalo_rho()
	{
		if (not evt_fixgridfastjet_allcalo_rho_isLoaded) {
			if (evt_fixgridfastjet_allcalo_rho_branch != 0) {
				evt_fixgridfastjet_allcalo_rho_branch->GetEntry(index);
			} else { }
			evt_fixgridfastjet_allcalo_rho_isLoaded = true;
		}
		return evt_fixgridfastjet_allcalo_rho_;
	}
	float &evt_fixgridfastjet_all_rho()
	{
		if (not evt_fixgridfastjet_all_rho_isLoaded) {
			if (evt_fixgridfastjet_all_rho_branch != 0) {
				evt_fixgridfastjet_all_rho_branch->GetEntry(index);
			} else { }
			evt_fixgridfastjet_all_rho_isLoaded = true;
		}
		return evt_fixgridfastjet_all_rho_;
	}
	float &evt_fixgridfastjet_centralcalo_rho()
	{
		if (not evt_fixgridfastjet_centralcalo_rho_isLoaded) {
			if (evt_fixgridfastjet_centralcalo_rho_branch != 0) {
				evt_fixgridfastjet_centralcalo_rho_branch->GetEntry(index);
			} else { }
			evt_fixgridfastjet_centralcalo_rho_isLoaded = true;
		}
		return evt_fixgridfastjet_centralcalo_rho_;
	}
	float &evt_fixgridfastjet_centralchargedpileup_rho()
	{
		if (not evt_fixgridfastjet_centralchargedpileup_rho_isLoaded) {
			if (evt_fixgridfastjet_centralchargedpileup_rho_branch != 0) {
				evt_fixgridfastjet_centralchargedpileup_rho_branch->GetEntry(index);
			} else { }
			evt_fixgridfastjet_centralchargedpileup_rho_isLoaded = true;
		}
		return evt_fixgridfastjet_centralchargedpileup_rho_;
	}
	float &evt_fixgridfastjet_centralneutral_rho()
	{
		if (not evt_fixgridfastjet_centralneutral_rho_isLoaded) {
			if (evt_fixgridfastjet_centralneutral_rho_branch != 0) {
				evt_fixgridfastjet_centralneutral_rho_branch->GetEntry(index);
			} else { }
			evt_fixgridfastjet_centralneutral_rho_isLoaded = true;
		}
		return evt_fixgridfastjet_centralneutral_rho_;
	}
	float &evt_kfactor()
	{
		if (not evt_kfactor_isLoaded) {
			if (evt_kfactor_branch != 0) {
				evt_kfactor_branch->GetEntry(index);
			} else { }
			evt_kfactor_isLoaded = true;
		}
		return evt_kfactor_;
	}
	float &evt_scale1fb()
	{
		if (not evt_scale1fb_isLoaded) {
			if (evt_scale1fb_branch != 0) {
				evt_scale1fb_branch->GetEntry(index);
			} else { }
			evt_scale1fb_isLoaded = true;
		}
		return evt_scale1fb_;
	}
	float &evt_xsec_excl()
	{
		if (not evt_xsec_excl_isLoaded) {
			if (evt_xsec_excl_branch != 0) {
				evt_xsec_excl_branch->GetEntry(index);
			} else { }
			evt_xsec_excl_isLoaded = true;
		}
		return evt_xsec_excl_;
	}
	float &evt_xsec_incl()
	{
		if (not evt_xsec_incl_isLoaded) {
			if (evt_xsec_incl_branch != 0) {
				evt_xsec_incl_branch->GetEntry(index);
			} else { }
			evt_xsec_incl_isLoaded = true;
		}
		return evt_xsec_incl_;
	}
	float &gen_met()
	{
		if (not gen_met_isLoaded) {
			if (gen_met_branch != 0) {
				gen_met_branch->GetEntry(index);
			} else { }
			gen_met_isLoaded = true;
		}
		return gen_met_;
	}
	float &gen_metPhi()
	{
		if (not gen_metPhi_isLoaded) {
			if (gen_metPhi_branch != 0) {
				gen_metPhi_branch->GetEntry(index);
			} else { }
			gen_metPhi_isLoaded = true;
		}
		return gen_metPhi_;
	}
	float &genps_alphaQCD()
	{
		if (not genps_alphaQCD_isLoaded) {
			if (genps_alphaQCD_branch != 0) {
				genps_alphaQCD_branch->GetEntry(index);
			} else { }
			genps_alphaQCD_isLoaded = true;
		}
		return genps_alphaQCD_;
	}
	float &genps_pthat()
	{
		if (not genps_pthat_isLoaded) {
			if (genps_pthat_branch != 0) {
				genps_pthat_branch->GetEntry(index);
			} else { }
			genps_pthat_isLoaded = true;
		}
		return genps_pthat_;
	}
	float &genps_qScale()
	{
		if (not genps_qScale_isLoaded) {
			if (genps_qScale_branch != 0) {
				genps_qScale_branch->GetEntry(index);
			} else { }
			genps_qScale_isLoaded = true;
		}
		return genps_qScale_;
	}
	float &genps_weight()
	{
		if (not genps_weight_isLoaded) {
			if (genps_weight_branch != 0) {
				genps_weight_branch->GetEntry(index);
			} else { }
			genps_weight_isLoaded = true;
		}
		return genps_weight_;
	}
	float &gen_sumEt()
	{
		if (not gen_sumEt_isLoaded) {
			if (gen_sumEt_branch != 0) {
				gen_sumEt_branch->GetEntry(index);
			} else { }
			gen_sumEt_isLoaded = true;
		}
		return gen_sumEt_;
	}
	float &hcalnoise_GetRecHitEnergy()
	{
		if (not hcalnoise_GetRecHitEnergy_isLoaded) {
			if (hcalnoise_GetRecHitEnergy_branch != 0) {
				hcalnoise_GetRecHitEnergy_branch->GetEntry(index);
			} else { }
			hcalnoise_GetRecHitEnergy_isLoaded = true;
		}
		return hcalnoise_GetRecHitEnergy_;
	}
	float &hcalnoise_GetRecHitEnergy15()
	{
		if (not hcalnoise_GetRecHitEnergy15_isLoaded) {
			if (hcalnoise_GetRecHitEnergy15_branch != 0) {
				hcalnoise_GetRecHitEnergy15_branch->GetEntry(index);
			} else { }
			hcalnoise_GetRecHitEnergy15_isLoaded = true;
		}
		return hcalnoise_GetRecHitEnergy15_;
	}
	float &hcalnoise_GetTotalCalibCharge()
	{
		if (not hcalnoise_GetTotalCalibCharge_isLoaded) {
			if (hcalnoise_GetTotalCalibCharge_branch != 0) {
				hcalnoise_GetTotalCalibCharge_branch->GetEntry(index);
			} else { }
			hcalnoise_GetTotalCalibCharge_isLoaded = true;
		}
		return hcalnoise_GetTotalCalibCharge_;
	}
	float &hcalnoise_TS4TS5NoiseSumE()
	{
		if (not hcalnoise_TS4TS5NoiseSumE_isLoaded) {
			if (hcalnoise_TS4TS5NoiseSumE_branch != 0) {
				hcalnoise_TS4TS5NoiseSumE_branch->GetEntry(index);
			} else { }
			hcalnoise_TS4TS5NoiseSumE_isLoaded = true;
		}
		return hcalnoise_TS4TS5NoiseSumE_;
	}
	float &hcalnoise_TS4TS5NoiseSumEt()
	{
		if (not hcalnoise_TS4TS5NoiseSumEt_isLoaded) {
			if (hcalnoise_TS4TS5NoiseSumEt_branch != 0) {
				hcalnoise_TS4TS5NoiseSumEt_branch->GetEntry(index);
			} else { }
			hcalnoise_TS4TS5NoiseSumEt_isLoaded = true;
		}
		return hcalnoise_TS4TS5NoiseSumEt_;
	}
	float &hcalnoise_eventChargeFraction()
	{
		if (not hcalnoise_eventChargeFraction_isLoaded) {
			if (hcalnoise_eventChargeFraction_branch != 0) {
				hcalnoise_eventChargeFraction_branch->GetEntry(index);
			} else { }
			hcalnoise_eventChargeFraction_isLoaded = true;
		}
		return hcalnoise_eventChargeFraction_;
	}
	float &hcalnoise_eventEMEnergy()
	{
		if (not hcalnoise_eventEMEnergy_isLoaded) {
			if (hcalnoise_eventEMEnergy_branch != 0) {
				hcalnoise_eventEMEnergy_branch->GetEntry(index);
			} else { }
			hcalnoise_eventEMEnergy_isLoaded = true;
		}
		return hcalnoise_eventEMEnergy_;
	}
	float &hcalnoise_eventEMFraction()
	{
		if (not hcalnoise_eventEMFraction_isLoaded) {
			if (hcalnoise_eventEMFraction_branch != 0) {
				hcalnoise_eventEMFraction_branch->GetEntry(index);
			} else { }
			hcalnoise_eventEMFraction_isLoaded = true;
		}
		return hcalnoise_eventEMFraction_;
	}
	float &hcalnoise_eventHadEnergy()
	{
		if (not hcalnoise_eventHadEnergy_isLoaded) {
			if (hcalnoise_eventHadEnergy_branch != 0) {
				hcalnoise_eventHadEnergy_branch->GetEntry(index);
			} else { }
			hcalnoise_eventHadEnergy_isLoaded = true;
		}
		return hcalnoise_eventHadEnergy_;
	}
	float &hcalnoise_eventTrackEnergy()
	{
		if (not hcalnoise_eventTrackEnergy_isLoaded) {
			if (hcalnoise_eventTrackEnergy_branch != 0) {
				hcalnoise_eventTrackEnergy_branch->GetEntry(index);
			} else { }
			hcalnoise_eventTrackEnergy_isLoaded = true;
		}
		return hcalnoise_eventTrackEnergy_;
	}
	float &hcalnoise_flatNoiseSumE()
	{
		if (not hcalnoise_flatNoiseSumE_isLoaded) {
			if (hcalnoise_flatNoiseSumE_branch != 0) {
				hcalnoise_flatNoiseSumE_branch->GetEntry(index);
			} else { }
			hcalnoise_flatNoiseSumE_isLoaded = true;
		}
		return hcalnoise_flatNoiseSumE_;
	}
	float &hcalnoise_flatNoiseSumEt()
	{
		if (not hcalnoise_flatNoiseSumEt_isLoaded) {
			if (hcalnoise_flatNoiseSumEt_branch != 0) {
				hcalnoise_flatNoiseSumEt_branch->GetEntry(index);
			} else { }
			hcalnoise_flatNoiseSumEt_isLoaded = true;
		}
		return hcalnoise_flatNoiseSumEt_;
	}
	float &hcalnoise_isolatedNoiseSumE()
	{
		if (not hcalnoise_isolatedNoiseSumE_isLoaded) {
			if (hcalnoise_isolatedNoiseSumE_branch != 0) {
				hcalnoise_isolatedNoiseSumE_branch->GetEntry(index);
			} else { }
			hcalnoise_isolatedNoiseSumE_isLoaded = true;
		}
		return hcalnoise_isolatedNoiseSumE_;
	}
	float &hcalnoise_isolatedNoiseSumEt()
	{
		if (not hcalnoise_isolatedNoiseSumEt_isLoaded) {
			if (hcalnoise_isolatedNoiseSumEt_branch != 0) {
				hcalnoise_isolatedNoiseSumEt_branch->GetEntry(index);
			} else { }
			hcalnoise_isolatedNoiseSumEt_isLoaded = true;
		}
		return hcalnoise_isolatedNoiseSumEt_;
	}
	float &hcalnoise_max10GeVHitTime()
	{
		if (not hcalnoise_max10GeVHitTime_isLoaded) {
			if (hcalnoise_max10GeVHitTime_branch != 0) {
				hcalnoise_max10GeVHitTime_branch->GetEntry(index);
			} else { }
			hcalnoise_max10GeVHitTime_isLoaded = true;
		}
		return hcalnoise_max10GeVHitTime_;
	}
	float &hcalnoise_max25GeVHitTime()
	{
		if (not hcalnoise_max25GeVHitTime_isLoaded) {
			if (hcalnoise_max25GeVHitTime_branch != 0) {
				hcalnoise_max25GeVHitTime_branch->GetEntry(index);
			} else { }
			hcalnoise_max25GeVHitTime_isLoaded = true;
		}
		return hcalnoise_max25GeVHitTime_;
	}
	float &hcalnoise_maxE10TS()
	{
		if (not hcalnoise_maxE10TS_isLoaded) {
			if (hcalnoise_maxE10TS_branch != 0) {
				hcalnoise_maxE10TS_branch->GetEntry(index);
			} else { }
			hcalnoise_maxE10TS_isLoaded = true;
		}
		return hcalnoise_maxE10TS_;
	}
	float &hcalnoise_maxE2Over10TS()
	{
		if (not hcalnoise_maxE2Over10TS_isLoaded) {
			if (hcalnoise_maxE2Over10TS_branch != 0) {
				hcalnoise_maxE2Over10TS_branch->GetEntry(index);
			} else { }
			hcalnoise_maxE2Over10TS_isLoaded = true;
		}
		return hcalnoise_maxE2Over10TS_;
	}
	float &hcalnoise_maxE2TS()
	{
		if (not hcalnoise_maxE2TS_isLoaded) {
			if (hcalnoise_maxE2TS_branch != 0) {
				hcalnoise_maxE2TS_branch->GetEntry(index);
			} else { }
			hcalnoise_maxE2TS_isLoaded = true;
		}
		return hcalnoise_maxE2TS_;
	}
	float &hcalnoise_min10GeVHitTime()
	{
		if (not hcalnoise_min10GeVHitTime_isLoaded) {
			if (hcalnoise_min10GeVHitTime_branch != 0) {
				hcalnoise_min10GeVHitTime_branch->GetEntry(index);
			} else { }
			hcalnoise_min10GeVHitTime_isLoaded = true;
		}
		return hcalnoise_min10GeVHitTime_;
	}
	float &hcalnoise_min25GeVHitTime()
	{
		if (not hcalnoise_min25GeVHitTime_isLoaded) {
			if (hcalnoise_min25GeVHitTime_branch != 0) {
				hcalnoise_min25GeVHitTime_branch->GetEntry(index);
			} else { }
			hcalnoise_min25GeVHitTime_isLoaded = true;
		}
		return hcalnoise_min25GeVHitTime_;
	}
	float &hcalnoise_minE10TS()
	{
		if (not hcalnoise_minE10TS_isLoaded) {
			if (hcalnoise_minE10TS_branch != 0) {
				hcalnoise_minE10TS_branch->GetEntry(index);
			} else { }
			hcalnoise_minE10TS_isLoaded = true;
		}
		return hcalnoise_minE10TS_;
	}
	float &hcalnoise_minE2Over10TS()
	{
		if (not hcalnoise_minE2Over10TS_isLoaded) {
			if (hcalnoise_minE2Over10TS_branch != 0) {
				hcalnoise_minE2Over10TS_branch->GetEntry(index);
			} else { }
			hcalnoise_minE2Over10TS_isLoaded = true;
		}
		return hcalnoise_minE2Over10TS_;
	}
	float &hcalnoise_minE2TS()
	{
		if (not hcalnoise_minE2TS_isLoaded) {
			if (hcalnoise_minE2TS_branch != 0) {
				hcalnoise_minE2TS_branch->GetEntry(index);
			} else { }
			hcalnoise_minE2TS_isLoaded = true;
		}
		return hcalnoise_minE2TS_;
	}
	float &hcalnoise_minHPDEMF()
	{
		if (not hcalnoise_minHPDEMF_isLoaded) {
			if (hcalnoise_minHPDEMF_branch != 0) {
				hcalnoise_minHPDEMF_branch->GetEntry(index);
			} else { }
			hcalnoise_minHPDEMF_isLoaded = true;
		}
		return hcalnoise_minHPDEMF_;
	}
	float &hcalnoise_minRBXEMF()
	{
		if (not hcalnoise_minRBXEMF_isLoaded) {
			if (hcalnoise_minRBXEMF_branch != 0) {
				hcalnoise_minRBXEMF_branch->GetEntry(index);
			} else { }
			hcalnoise_minRBXEMF_isLoaded = true;
		}
		return hcalnoise_minRBXEMF_;
	}
	float &hcalnoise_rms10GeVHitTime()
	{
		if (not hcalnoise_rms10GeVHitTime_isLoaded) {
			if (hcalnoise_rms10GeVHitTime_branch != 0) {
				hcalnoise_rms10GeVHitTime_branch->GetEntry(index);
			} else { }
			hcalnoise_rms10GeVHitTime_isLoaded = true;
		}
		return hcalnoise_rms10GeVHitTime_;
	}
	float &hcalnoise_rms25GeVHitTime()
	{
		if (not hcalnoise_rms25GeVHitTime_isLoaded) {
			if (hcalnoise_rms25GeVHitTime_branch != 0) {
				hcalnoise_rms25GeVHitTime_branch->GetEntry(index);
			} else { }
			hcalnoise_rms25GeVHitTime_isLoaded = true;
		}
		return hcalnoise_rms25GeVHitTime_;
	}
	float &hcalnoise_spikeNoiseSumE()
	{
		if (not hcalnoise_spikeNoiseSumE_isLoaded) {
			if (hcalnoise_spikeNoiseSumE_branch != 0) {
				hcalnoise_spikeNoiseSumE_branch->GetEntry(index);
			} else { }
			hcalnoise_spikeNoiseSumE_isLoaded = true;
		}
		return hcalnoise_spikeNoiseSumE_;
	}
	float &hcalnoise_spikeNoiseSumEt()
	{
		if (not hcalnoise_spikeNoiseSumEt_isLoaded) {
			if (hcalnoise_spikeNoiseSumEt_branch != 0) {
				hcalnoise_spikeNoiseSumEt_branch->GetEntry(index);
			} else { }
			hcalnoise_spikeNoiseSumEt_isLoaded = true;
		}
		return hcalnoise_spikeNoiseSumEt_;
	}
	float &hcalnoise_triangleNoiseSumE()
	{
		if (not hcalnoise_triangleNoiseSumE_isLoaded) {
			if (hcalnoise_triangleNoiseSumE_branch != 0) {
				hcalnoise_triangleNoiseSumE_branch->GetEntry(index);
			} else { }
			hcalnoise_triangleNoiseSumE_isLoaded = true;
		}
		return hcalnoise_triangleNoiseSumE_;
	}
	float &hcalnoise_triangleNoiseSumEt()
	{
		if (not hcalnoise_triangleNoiseSumEt_isLoaded) {
			if (hcalnoise_triangleNoiseSumEt_branch != 0) {
				hcalnoise_triangleNoiseSumEt_branch->GetEntry(index);
			} else { }
			hcalnoise_triangleNoiseSumEt_isLoaded = true;
		}
		return hcalnoise_triangleNoiseSumEt_;
	}
	float &pdfinfo_pdf1()
	{
		if (not pdfinfo_pdf1_isLoaded) {
			if (pdfinfo_pdf1_branch != 0) {
				pdfinfo_pdf1_branch->GetEntry(index);
			} else { }
			pdfinfo_pdf1_isLoaded = true;
		}
		return pdfinfo_pdf1_;
	}
	float &pdfinfo_pdf2()
	{
		if (not pdfinfo_pdf2_isLoaded) {
			if (pdfinfo_pdf2_branch != 0) {
				pdfinfo_pdf2_branch->GetEntry(index);
			} else { }
			pdfinfo_pdf2_isLoaded = true;
		}
		return pdfinfo_pdf2_;
	}
	float &pdfinfo_scale()
	{
		if (not pdfinfo_scale_isLoaded) {
			if (pdfinfo_scale_branch != 0) {
				pdfinfo_scale_branch->GetEntry(index);
			} else { }
			pdfinfo_scale_isLoaded = true;
		}
		return pdfinfo_scale_;
	}
	float &pdfinfo_x1()
	{
		if (not pdfinfo_x1_isLoaded) {
			if (pdfinfo_x1_branch != 0) {
				pdfinfo_x1_branch->GetEntry(index);
			} else { }
			pdfinfo_x1_isLoaded = true;
		}
		return pdfinfo_x1_;
	}
	float &pdfinfo_x2()
	{
		if (not pdfinfo_x2_isLoaded) {
			if (pdfinfo_x2_branch != 0) {
				pdfinfo_x2_branch->GetEntry(index);
			} else { }
			pdfinfo_x2_isLoaded = true;
		}
		return pdfinfo_x2_;
	}
	float &evt_fixgrid_rho_all()
	{
		if (not evt_fixgrid_rho_all_isLoaded) {
			if (evt_fixgrid_rho_all_branch != 0) {
				evt_fixgrid_rho_all_branch->GetEntry(index);
			} else { }
			evt_fixgrid_rho_all_isLoaded = true;
		}
		return evt_fixgrid_rho_all_;
	}
	float &evt_fixgrid_rho_ctr()
	{
		if (not evt_fixgrid_rho_ctr_isLoaded) {
			if (evt_fixgrid_rho_ctr_branch != 0) {
				evt_fixgrid_rho_ctr_branch->GetEntry(index);
			} else { }
			evt_fixgrid_rho_ctr_isLoaded = true;
		}
		return evt_fixgrid_rho_ctr_;
	}
	float &evt_fixgrid_rho_fwd()
	{
		if (not evt_fixgrid_rho_fwd_isLoaded) {
			if (evt_fixgrid_rho_fwd_branch != 0) {
				evt_fixgrid_rho_fwd_branch->GetEntry(index);
			} else { }
			evt_fixgrid_rho_fwd_isLoaded = true;
		}
		return evt_fixgrid_rho_fwd_;
	}
	float &evt_pfmet()
	{
		if (not evt_pfmet_isLoaded) {
			if (evt_pfmet_branch != 0) {
				evt_pfmet_branch->GetEntry(index);
			} else { }
			evt_pfmet_isLoaded = true;
		}
		return evt_pfmet_;
	}
	float &evt_pfmetPhi()
	{
		if (not evt_pfmetPhi_isLoaded) {
			if (evt_pfmetPhi_branch != 0) {
				evt_pfmetPhi_branch->GetEntry(index);
			} else { }
			evt_pfmetPhi_isLoaded = true;
		}
		return evt_pfmetPhi_;
	}
	float &evt_pfmetSig()
	{
		if (not evt_pfmetSig_isLoaded) {
			if (evt_pfmetSig_branch != 0) {
				evt_pfmetSig_branch->GetEntry(index);
			} else { }
			evt_pfmetSig_isLoaded = true;
		}
		return evt_pfmetSig_;
	}
	float &evt_pfmetSignificance()
	{
		if (not evt_pfmetSignificance_isLoaded) {
			if (evt_pfmetSignificance_branch != 0) {
				evt_pfmetSignificance_branch->GetEntry(index);
			} else { }
			evt_pfmetSignificance_isLoaded = true;
		}
		return evt_pfmetSignificance_;
	}
	float &evt_pfsumet()
	{
		if (not evt_pfsumet_isLoaded) {
			if (evt_pfsumet_branch != 0) {
				evt_pfsumet_branch->GetEntry(index);
			} else { }
			evt_pfsumet_isLoaded = true;
		}
		return evt_pfsumet_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &evt_bsp4()
	{
		if (not evt_bsp4_isLoaded) {
			if (evt_bsp4_branch != 0) {
				evt_bsp4_branch->GetEntry(index);
			} else { }
			evt_bsp4_isLoaded = true;
		}
		return evt_bsp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_motherp4()
	{
		if (not els_mc_motherp4_isLoaded) {
			if (els_mc_motherp4_branch != 0) {
				els_mc_motherp4_branch->GetEntry(index);
			} else { }
			els_mc_motherp4_isLoaded = true;
		}
		return els_mc_motherp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_p4()
	{
		if (not els_mc_p4_isLoaded) {
			if (els_mc_p4_branch != 0) {
				els_mc_p4_branch->GetEntry(index);
			} else { }
			els_mc_p4_isLoaded = true;
		}
		return els_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_mc_motherp4()
	{
		if (not mus_mc_motherp4_isLoaded) {
			if (mus_mc_motherp4_branch != 0) {
				mus_mc_motherp4_branch->GetEntry(index);
			} else { }
			mus_mc_motherp4_isLoaded = true;
		}
		return mus_mc_motherp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_mc_p4()
	{
		if (not mus_mc_p4_isLoaded) {
			if (mus_mc_p4_branch != 0) {
				mus_mc_p4_branch->GetEntry(index);
			} else { }
			mus_mc_p4_isLoaded = true;
		}
		return mus_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_gp_p4()
	{
		if (not pfjets_mc_gp_p4_isLoaded) {
			if (pfjets_mc_gp_p4_branch != 0) {
				pfjets_mc_gp_p4_branch->GetEntry(index);
			} else { }
			pfjets_mc_gp_p4_isLoaded = true;
		}
		return pfjets_mc_gp_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_motherp4()
	{
		if (not pfjets_mc_motherp4_isLoaded) {
			if (pfjets_mc_motherp4_branch != 0) {
				pfjets_mc_motherp4_branch->GetEntry(index);
			} else { }
			pfjets_mc_motherp4_isLoaded = true;
		}
		return pfjets_mc_motherp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_p4()
	{
		if (not pfjets_mc_p4_isLoaded) {
			if (pfjets_mc_p4_branch != 0) {
				pfjets_mc_p4_branch->GetEntry(index);
			} else { }
			pfjets_mc_p4_isLoaded = true;
		}
		return pfjets_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_mc_motherp4()
	{
		if (not photons_mc_motherp4_isLoaded) {
			if (photons_mc_motherp4_branch != 0) {
				photons_mc_motherp4_branch->GetEntry(index);
			} else { }
			photons_mc_motherp4_isLoaded = true;
		}
		return photons_mc_motherp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_mc_p4()
	{
		if (not photons_mc_p4_isLoaded) {
			if (photons_mc_p4_branch != 0) {
				photons_mc_p4_branch->GetEntry(index);
			} else { }
			photons_mc_p4_isLoaded = true;
		}
		return photons_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_conv_pos_p4()
	{
		if (not els_conv_pos_p4_isLoaded) {
			if (els_conv_pos_p4_branch != 0) {
				els_conv_pos_p4_branch->GetEntry(index);
			} else { }
			els_conv_pos_p4_isLoaded = true;
		}
		return els_conv_pos_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_inner_position()
	{
		if (not els_inner_position_isLoaded) {
			if (els_inner_position_branch != 0) {
				els_inner_position_branch->GetEntry(index);
			} else { }
			els_inner_position_isLoaded = true;
		}
		return els_inner_position_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_patMatch_p4()
	{
		if (not els_mc_patMatch_p4_isLoaded) {
			if (els_mc_patMatch_p4_branch != 0) {
				els_mc_patMatch_p4_branch->GetEntry(index);
			} else { }
			els_mc_patMatch_p4_isLoaded = true;
		}
		return els_mc_patMatch_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_outer_position()
	{
		if (not els_outer_position_isLoaded) {
			if (els_outer_position_branch != 0) {
				els_outer_position_branch->GetEntry(index);
			} else { }
			els_outer_position_isLoaded = true;
		}
		return els_outer_position_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4()
	{
		if (not els_p4_isLoaded) {
			if (els_p4_branch != 0) {
				els_p4_branch->GetEntry(index);
			} else { }
			els_p4_isLoaded = true;
		}
		return els_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4In()
	{
		if (not els_p4In_isLoaded) {
			if (els_p4In_branch != 0) {
				els_p4In_branch->GetEntry(index);
			} else { }
			els_p4In_isLoaded = true;
		}
		return els_p4In_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4Out()
	{
		if (not els_p4Out_isLoaded) {
			if (els_p4Out_branch != 0) {
				els_p4Out_branch->GetEntry(index);
			} else { }
			els_p4Out_isLoaded = true;
		}
		return els_p4Out_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_p4()
	{
		if (not els_trk_p4_isLoaded) {
			if (els_trk_p4_branch != 0) {
				els_trk_p4_branch->GetEntry(index);
			} else { }
			els_trk_p4_isLoaded = true;
		}
		return els_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_vertex_p4()
	{
		if (not els_trk_vertex_p4_isLoaded) {
			if (els_trk_vertex_p4_branch != 0) {
				els_trk_vertex_p4_branch->GetEntry(index);
			} else { }
			els_trk_vertex_p4_isLoaded = true;
		}
		return els_trk_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_vertex_p4()
	{
		if (not els_vertex_p4_isLoaded) {
			if (els_vertex_p4_branch != 0) {
				els_vertex_p4_branch->GetEntry(index);
			} else { }
			els_vertex_p4_isLoaded = true;
		}
		return els_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genjets_p4NoMuNoNu()
	{
		if (not genjets_p4NoMuNoNu_isLoaded) {
			if (genjets_p4NoMuNoNu_branch != 0) {
				genjets_p4NoMuNoNu_branch->GetEntry(index);
			} else { }
			genjets_p4NoMuNoNu_isLoaded = true;
		}
		return genjets_p4NoMuNoNu_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_p4()
	{
		if (not genps_p4_isLoaded) {
			if (genps_p4_branch != 0) {
				genps_p4_branch->GetEntry(index);
			} else { }
			genps_p4_isLoaded = true;
		}
		return genps_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_prod_vtx()
	{
		if (not genps_prod_vtx_isLoaded) {
			if (genps_prod_vtx_branch != 0) {
				genps_prod_vtx_branch->GetEntry(index);
			} else { }
			genps_prod_vtx_isLoaded = true;
		}
		return genps_prod_vtx_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_p4()
	{
		if (not mus_gfit_p4_isLoaded) {
			if (mus_gfit_p4_branch != 0) {
				mus_gfit_p4_branch->GetEntry(index);
			} else { }
			mus_gfit_p4_isLoaded = true;
		}
		return mus_gfit_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_vertex_p4()
	{
		if (not mus_gfit_vertex_p4_isLoaded) {
			if (mus_gfit_vertex_p4_branch != 0) {
				mus_gfit_vertex_p4_branch->GetEntry(index);
			} else { }
			mus_gfit_vertex_p4_isLoaded = true;
		}
		return mus_gfit_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_mc_patMatch_p4()
	{
		if (not mus_mc_patMatch_p4_isLoaded) {
			if (mus_mc_patMatch_p4_branch != 0) {
				mus_mc_patMatch_p4_branch->GetEntry(index);
			} else { }
			mus_mc_patMatch_p4_isLoaded = true;
		}
		return mus_mc_patMatch_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_p4()
	{
		if (not mus_p4_isLoaded) {
			if (mus_p4_branch != 0) {
				mus_p4_branch->GetEntry(index);
			} else { }
			mus_p4_isLoaded = true;
		}
		return mus_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_pfp4()
	{
		if (not mus_pfp4_isLoaded) {
			if (mus_pfp4_branch != 0) {
				mus_pfp4_branch->GetEntry(index);
			} else { }
			mus_pfp4_isLoaded = true;
		}
		return mus_pfp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_sta_p4()
	{
		if (not mus_sta_p4_isLoaded) {
			if (mus_sta_p4_branch != 0) {
				mus_sta_p4_branch->GetEntry(index);
			} else { }
			mus_sta_p4_isLoaded = true;
		}
		return mus_sta_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_sta_vertex_p4()
	{
		if (not mus_sta_vertex_p4_isLoaded) {
			if (mus_sta_vertex_p4_branch != 0) {
				mus_sta_vertex_p4_branch->GetEntry(index);
			} else { }
			mus_sta_vertex_p4_isLoaded = true;
		}
		return mus_sta_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_trk_p4()
	{
		if (not mus_trk_p4_isLoaded) {
			if (mus_trk_p4_branch != 0) {
				mus_trk_p4_branch->GetEntry(index);
			} else { }
			mus_trk_p4_isLoaded = true;
		}
		return mus_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_vertex_p4()
	{
		if (not mus_vertex_p4_isLoaded) {
			if (mus_vertex_p4_branch != 0) {
				mus_vertex_p4_branch->GetEntry(index);
			} else { }
			mus_vertex_p4_isLoaded = true;
		}
		return mus_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfcands_p4()
	{
		if (not pfcands_p4_isLoaded) {
			if (pfcands_p4_branch != 0) {
				pfcands_p4_branch->GetEntry(index);
			} else { }
			pfcands_p4_isLoaded = true;
		}
		return pfcands_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_p4()
	{
		if (not pfjets_p4_isLoaded) {
			if (pfjets_p4_branch != 0) {
				pfjets_p4_branch->GetEntry(index);
			} else { }
			pfjets_p4_isLoaded = true;
		}
		return pfjets_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &taus_pf_lead_chargecand_p4()
	{
		if (not taus_pf_lead_chargecand_p4_isLoaded) {
			if (taus_pf_lead_chargecand_p4_branch != 0) {
				taus_pf_lead_chargecand_p4_branch->GetEntry(index);
			} else { }
			taus_pf_lead_chargecand_p4_isLoaded = true;
		}
		return taus_pf_lead_chargecand_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &taus_pf_lead_neutrcand_p4()
	{
		if (not taus_pf_lead_neutrcand_p4_isLoaded) {
			if (taus_pf_lead_neutrcand_p4_branch != 0) {
				taus_pf_lead_neutrcand_p4_branch->GetEntry(index);
			} else { }
			taus_pf_lead_neutrcand_p4_isLoaded = true;
		}
		return taus_pf_lead_neutrcand_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &taus_pf_p4()
	{
		if (not taus_pf_p4_isLoaded) {
			if (taus_pf_p4_branch != 0) {
				taus_pf_p4_branch->GetEntry(index);
			} else { }
			taus_pf_p4_isLoaded = true;
		}
		return taus_pf_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_p4()
	{
		if (not photons_p4_isLoaded) {
			if (photons_p4_branch != 0) {
				photons_p4_branch->GetEntry(index);
			} else { }
			photons_p4_isLoaded = true;
		}
		return photons_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &vtxs_position()
	{
		if (not vtxs_position_isLoaded) {
			if (vtxs_position_branch != 0) {
				vtxs_position_branch->GetEntry(index);
			} else { }
			vtxs_position_isLoaded = true;
		}
		return vtxs_position_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &genps_lepdaughter_p4()
	{
		if (not genps_lepdaughter_p4_isLoaded) {
			if (genps_lepdaughter_p4_branch != 0) {
				genps_lepdaughter_p4_branch->GetEntry(index);
			} else { }
			genps_lepdaughter_p4_isLoaded = true;
		}
		return genps_lepdaughter_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hlt_trigObjs_p4()
	{
		if (not hlt_trigObjs_p4_isLoaded) {
			if (hlt_trigObjs_p4_branch != 0) {
				hlt_trigObjs_p4_branch->GetEntry(index);
			} else { }
			hlt_trigObjs_p4_isLoaded = true;
		}
		return hlt_trigObjs_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &taus_pf_isocands_p4()
	{
		if (not taus_pf_isocands_p4_isLoaded) {
			if (taus_pf_isocands_p4_branch != 0) {
				taus_pf_isocands_p4_branch->GetEntry(index);
			} else { }
			taus_pf_isocands_p4_isLoaded = true;
		}
		return taus_pf_isocands_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &taus_pf_signalcands_p4()
	{
		if (not taus_pf_signalcands_p4_isLoaded) {
			if (taus_pf_signalcands_p4_branch != 0) {
				taus_pf_signalcands_p4_branch->GetEntry(index);
			} else { }
			taus_pf_signalcands_p4_isLoaded = true;
		}
		return taus_pf_signalcands_p4_;
	}
	vector<float> &evt_bs_covMatrix()
	{
		if (not evt_bs_covMatrix_isLoaded) {
			if (evt_bs_covMatrix_branch != 0) {
				evt_bs_covMatrix_branch->GetEntry(index);
			} else { }
			evt_bs_covMatrix_isLoaded = true;
		}
		return evt_bs_covMatrix_;
	}
	vector<float> &els_mc3dr()
	{
		if (not els_mc3dr_isLoaded) {
			if (els_mc3dr_branch != 0) {
				els_mc3dr_branch->GetEntry(index);
			} else { }
			els_mc3dr_isLoaded = true;
		}
		return els_mc3dr_;
	}
	vector<float> &els_mcdr()
	{
		if (not els_mcdr_isLoaded) {
			if (els_mcdr_branch != 0) {
				els_mcdr_branch->GetEntry(index);
			} else { }
			els_mcdr_isLoaded = true;
		}
		return els_mcdr_;
	}
	vector<float> &mus_mc3dr()
	{
		if (not mus_mc3dr_isLoaded) {
			if (mus_mc3dr_branch != 0) {
				mus_mc3dr_branch->GetEntry(index);
			} else { }
			mus_mc3dr_isLoaded = true;
		}
		return mus_mc3dr_;
	}
	vector<float> &mus_mcdr()
	{
		if (not mus_mcdr_isLoaded) {
			if (mus_mcdr_branch != 0) {
				mus_mcdr_branch->GetEntry(index);
			} else { }
			mus_mcdr_isLoaded = true;
		}
		return mus_mcdr_;
	}
	vector<float> &pfjets_mc3dr()
	{
		if (not pfjets_mc3dr_isLoaded) {
			if (pfjets_mc3dr_branch != 0) {
				pfjets_mc3dr_branch->GetEntry(index);
			} else { }
			pfjets_mc3dr_isLoaded = true;
		}
		return pfjets_mc3dr_;
	}
	vector<float> &pfjets_mcdr()
	{
		if (not pfjets_mcdr_isLoaded) {
			if (pfjets_mcdr_branch != 0) {
				pfjets_mcdr_branch->GetEntry(index);
			} else { }
			pfjets_mcdr_isLoaded = true;
		}
		return pfjets_mcdr_;
	}
	vector<float> &pfjets_mc_emEnergy()
	{
		if (not pfjets_mc_emEnergy_isLoaded) {
			if (pfjets_mc_emEnergy_branch != 0) {
				pfjets_mc_emEnergy_branch->GetEntry(index);
			} else { }
			pfjets_mc_emEnergy_isLoaded = true;
		}
		return pfjets_mc_emEnergy_;
	}
	vector<float> &pfjets_mc_gpdr()
	{
		if (not pfjets_mc_gpdr_isLoaded) {
			if (pfjets_mc_gpdr_branch != 0) {
				pfjets_mc_gpdr_branch->GetEntry(index);
			} else { }
			pfjets_mc_gpdr_isLoaded = true;
		}
		return pfjets_mc_gpdr_;
	}
	vector<float> &pfjets_mc_hadEnergy()
	{
		if (not pfjets_mc_hadEnergy_isLoaded) {
			if (pfjets_mc_hadEnergy_branch != 0) {
				pfjets_mc_hadEnergy_branch->GetEntry(index);
			} else { }
			pfjets_mc_hadEnergy_isLoaded = true;
		}
		return pfjets_mc_hadEnergy_;
	}
	vector<float> &pfjets_mc_invEnergy()
	{
		if (not pfjets_mc_invEnergy_isLoaded) {
			if (pfjets_mc_invEnergy_branch != 0) {
				pfjets_mc_invEnergy_branch->GetEntry(index);
			} else { }
			pfjets_mc_invEnergy_isLoaded = true;
		}
		return pfjets_mc_invEnergy_;
	}
	vector<float> &pfjets_mc_otherEnergy()
	{
		if (not pfjets_mc_otherEnergy_isLoaded) {
			if (pfjets_mc_otherEnergy_branch != 0) {
				pfjets_mc_otherEnergy_branch->GetEntry(index);
			} else { }
			pfjets_mc_otherEnergy_isLoaded = true;
		}
		return pfjets_mc_otherEnergy_;
	}
	vector<float> &photons_mc3dr()
	{
		if (not photons_mc3dr_isLoaded) {
			if (photons_mc3dr_branch != 0) {
				photons_mc3dr_branch->GetEntry(index);
			} else { }
			photons_mc3dr_isLoaded = true;
		}
		return photons_mc3dr_;
	}
	vector<float> &photons_mcdr()
	{
		if (not photons_mcdr_isLoaded) {
			if (photons_mcdr_branch != 0) {
				photons_mcdr_branch->GetEntry(index);
			} else { }
			photons_mcdr_isLoaded = true;
		}
		return photons_mcdr_;
	}
	vector<float> &els_chi2()
	{
		if (not els_chi2_isLoaded) {
			if (els_chi2_branch != 0) {
				els_chi2_branch->GetEntry(index);
			} else { }
			els_chi2_isLoaded = true;
		}
		return els_chi2_;
	}
	vector<float> &els_ckf_chi2()
	{
		if (not els_ckf_chi2_isLoaded) {
			if (els_ckf_chi2_branch != 0) {
				els_ckf_chi2_branch->GetEntry(index);
			} else { }
			els_ckf_chi2_isLoaded = true;
		}
		return els_ckf_chi2_;
	}
	vector<float> &els_ckf_ndof()
	{
		if (not els_ckf_ndof_isLoaded) {
			if (els_ckf_ndof_branch != 0) {
				els_ckf_ndof_branch->GetEntry(index);
			} else { }
			els_ckf_ndof_isLoaded = true;
		}
		return els_ckf_ndof_;
	}
	vector<float> &els_conv_dcot()
	{
		if (not els_conv_dcot_isLoaded) {
			if (els_conv_dcot_branch != 0) {
				els_conv_dcot_branch->GetEntry(index);
			} else { }
			els_conv_dcot_isLoaded = true;
		}
		return els_conv_dcot_;
	}
	vector<float> &els_conv_dist()
	{
		if (not els_conv_dist_isLoaded) {
			if (els_conv_dist_branch != 0) {
				els_conv_dist_branch->GetEntry(index);
			} else { }
			els_conv_dist_isLoaded = true;
		}
		return els_conv_dist_;
	}
	vector<float> &els_conv_old_dcot()
	{
		if (not els_conv_old_dcot_isLoaded) {
			if (els_conv_old_dcot_branch != 0) {
				els_conv_old_dcot_branch->GetEntry(index);
			} else { }
			els_conv_old_dcot_isLoaded = true;
		}
		return els_conv_old_dcot_;
	}
	vector<float> &els_conv_old_dist()
	{
		if (not els_conv_old_dist_isLoaded) {
			if (els_conv_old_dist_branch != 0) {
				els_conv_old_dist_branch->GetEntry(index);
			} else { }
			els_conv_old_dist_isLoaded = true;
		}
		return els_conv_old_dist_;
	}
	vector<float> &els_conv_old_radius()
	{
		if (not els_conv_old_radius_isLoaded) {
			if (els_conv_old_radius_branch != 0) {
				els_conv_old_radius_branch->GetEntry(index);
			} else { }
			els_conv_old_radius_isLoaded = true;
		}
		return els_conv_old_radius_;
	}
	vector<float> &els_conv_radius()
	{
		if (not els_conv_radius_isLoaded) {
			if (els_conv_radius_branch != 0) {
				els_conv_radius_branch->GetEntry(index);
			} else { }
			els_conv_radius_isLoaded = true;
		}
		return els_conv_radius_;
	}
	vector<float> &els_d0()
	{
		if (not els_d0_isLoaded) {
			if (els_d0_branch != 0) {
				els_d0_branch->GetEntry(index);
			} else { }
			els_d0_isLoaded = true;
		}
		return els_d0_;
	}
	vector<float> &els_d0Err()
	{
		if (not els_d0Err_isLoaded) {
			if (els_d0Err_branch != 0) {
				els_d0Err_branch->GetEntry(index);
			} else { }
			els_d0Err_isLoaded = true;
		}
		return els_d0Err_;
	}
	vector<float> &els_d0corr()
	{
		if (not els_d0corr_isLoaded) {
			if (els_d0corr_branch != 0) {
				els_d0corr_branch->GetEntry(index);
			} else { }
			els_d0corr_isLoaded = true;
		}
		return els_d0corr_;
	}
	vector<float> &els_d0corrPhi()
	{
		if (not els_d0corrPhi_isLoaded) {
			if (els_d0corrPhi_branch != 0) {
				els_d0corrPhi_branch->GetEntry(index);
			} else { }
			els_d0corrPhi_isLoaded = true;
		}
		return els_d0corrPhi_;
	}
	vector<float> &els_d0phiCov()
	{
		if (not els_d0phiCov_isLoaded) {
			if (els_d0phiCov_branch != 0) {
				els_d0phiCov_branch->GetEntry(index);
			} else { }
			els_d0phiCov_isLoaded = true;
		}
		return els_d0phiCov_;
	}
	vector<float> &els_dEtaIn()
	{
		if (not els_dEtaIn_isLoaded) {
			if (els_dEtaIn_branch != 0) {
				els_dEtaIn_branch->GetEntry(index);
			} else { }
			els_dEtaIn_isLoaded = true;
		}
		return els_dEtaIn_;
	}
	vector<float> &els_dEtaOut()
	{
		if (not els_dEtaOut_isLoaded) {
			if (els_dEtaOut_branch != 0) {
				els_dEtaOut_branch->GetEntry(index);
			} else { }
			els_dEtaOut_isLoaded = true;
		}
		return els_dEtaOut_;
	}
	vector<float> &els_dPhiIn()
	{
		if (not els_dPhiIn_isLoaded) {
			if (els_dPhiIn_branch != 0) {
				els_dPhiIn_branch->GetEntry(index);
			} else { }
			els_dPhiIn_isLoaded = true;
		}
		return els_dPhiIn_;
	}
	vector<float> &els_dPhiInPhiOut()
	{
		if (not els_dPhiInPhiOut_isLoaded) {
			if (els_dPhiInPhiOut_branch != 0) {
				els_dPhiInPhiOut_branch->GetEntry(index);
			} else { }
			els_dPhiInPhiOut_isLoaded = true;
		}
		return els_dPhiInPhiOut_;
	}
	vector<float> &els_dPhiOut()
	{
		if (not els_dPhiOut_isLoaded) {
			if (els_dPhiOut_branch != 0) {
				els_dPhiOut_branch->GetEntry(index);
			} else { }
			els_dPhiOut_isLoaded = true;
		}
		return els_dPhiOut_;
	}
	vector<float> &els_deltaEtaEleClusterTrackAtCalo()
	{
		if (not els_deltaEtaEleClusterTrackAtCalo_isLoaded) {
			if (els_deltaEtaEleClusterTrackAtCalo_branch != 0) {
				els_deltaEtaEleClusterTrackAtCalo_branch->GetEntry(index);
			} else { }
			els_deltaEtaEleClusterTrackAtCalo_isLoaded = true;
		}
		return els_deltaEtaEleClusterTrackAtCalo_;
	}
	vector<float> &els_deltaPhiEleClusterTrackAtCalo()
	{
		if (not els_deltaPhiEleClusterTrackAtCalo_isLoaded) {
			if (els_deltaPhiEleClusterTrackAtCalo_branch != 0) {
				els_deltaPhiEleClusterTrackAtCalo_branch->GetEntry(index);
			} else { }
			els_deltaPhiEleClusterTrackAtCalo_isLoaded = true;
		}
		return els_deltaPhiEleClusterTrackAtCalo_;
	}
	vector<float> &els_dxyPV()
	{
		if (not els_dxyPV_isLoaded) {
			if (els_dxyPV_branch != 0) {
				els_dxyPV_branch->GetEntry(index);
			} else { }
			els_dxyPV_isLoaded = true;
		}
		return els_dxyPV_;
	}
	vector<float> &els_dzPV()
	{
		if (not els_dzPV_isLoaded) {
			if (els_dzPV_branch != 0) {
				els_dzPV_branch->GetEntry(index);
			} else { }
			els_dzPV_isLoaded = true;
		}
		return els_dzPV_;
	}
	vector<float> &els_e1x5()
	{
		if (not els_e1x5_isLoaded) {
			if (els_e1x5_branch != 0) {
				els_e1x5_branch->GetEntry(index);
			} else { }
			els_e1x5_isLoaded = true;
		}
		return els_e1x5_;
	}
	vector<float> &els_e1x5_full5x5()
	{
		if (not els_e1x5_full5x5_isLoaded) {
			if (els_e1x5_full5x5_branch != 0) {
				els_e1x5_full5x5_branch->GetEntry(index);
			} else { }
			els_e1x5_full5x5_isLoaded = true;
		}
		return els_e1x5_full5x5_;
	}
	vector<float> &els_e2x5Max()
	{
		if (not els_e2x5Max_isLoaded) {
			if (els_e2x5Max_branch != 0) {
				els_e2x5Max_branch->GetEntry(index);
			} else { }
			els_e2x5Max_isLoaded = true;
		}
		return els_e2x5Max_;
	}
	vector<float> &els_e2x5Max_full5x5()
	{
		if (not els_e2x5Max_full5x5_isLoaded) {
			if (els_e2x5Max_full5x5_branch != 0) {
				els_e2x5Max_full5x5_branch->GetEntry(index);
			} else { }
			els_e2x5Max_full5x5_isLoaded = true;
		}
		return els_e2x5Max_full5x5_;
	}
	vector<float> &els_e3x3()
	{
		if (not els_e3x3_isLoaded) {
			if (els_e3x3_branch != 0) {
				els_e3x3_branch->GetEntry(index);
			} else { }
			els_e3x3_isLoaded = true;
		}
		return els_e3x3_;
	}
	vector<float> &els_e5x5()
	{
		if (not els_e5x5_isLoaded) {
			if (els_e5x5_branch != 0) {
				els_e5x5_branch->GetEntry(index);
			} else { }
			els_e5x5_isLoaded = true;
		}
		return els_e5x5_;
	}
	vector<float> &els_e5x5_full5x5()
	{
		if (not els_e5x5_full5x5_isLoaded) {
			if (els_e5x5_full5x5_branch != 0) {
				els_e5x5_full5x5_branch->GetEntry(index);
			} else { }
			els_e5x5_full5x5_isLoaded = true;
		}
		return els_e5x5_full5x5_;
	}
	vector<float> &els_eMax()
	{
		if (not els_eMax_isLoaded) {
			if (els_eMax_branch != 0) {
				els_eMax_branch->GetEntry(index);
			} else { }
			els_eMax_isLoaded = true;
		}
		return els_eMax_;
	}
	vector<float> &els_eOverPIn()
	{
		if (not els_eOverPIn_isLoaded) {
			if (els_eOverPIn_branch != 0) {
				els_eOverPIn_branch->GetEntry(index);
			} else { }
			els_eOverPIn_isLoaded = true;
		}
		return els_eOverPIn_;
	}
	vector<float> &els_eOverPOut()
	{
		if (not els_eOverPOut_isLoaded) {
			if (els_eOverPOut_branch != 0) {
				els_eOverPOut_branch->GetEntry(index);
			} else { }
			els_eOverPOut_isLoaded = true;
		}
		return els_eOverPOut_;
	}
	vector<float> &els_eSC()
	{
		if (not els_eSC_isLoaded) {
			if (els_eSC_branch != 0) {
				els_eSC_branch->GetEntry(index);
			} else { }
			els_eSC_isLoaded = true;
		}
		return els_eSC_;
	}
	vector<float> &els_eSCPresh()
	{
		if (not els_eSCPresh_isLoaded) {
			if (els_eSCPresh_branch != 0) {
				els_eSCPresh_branch->GetEntry(index);
			} else { }
			els_eSCPresh_isLoaded = true;
		}
		return els_eSCPresh_;
	}
	vector<float> &els_eSCRaw()
	{
		if (not els_eSCRaw_isLoaded) {
			if (els_eSCRaw_branch != 0) {
				els_eSCRaw_branch->GetEntry(index);
			} else { }
			els_eSCRaw_isLoaded = true;
		}
		return els_eSCRaw_;
	}
	vector<float> &els_eSeed()
	{
		if (not els_eSeed_isLoaded) {
			if (els_eSeed_branch != 0) {
				els_eSeed_branch->GetEntry(index);
			} else { }
			els_eSeed_isLoaded = true;
		}
		return els_eSeed_;
	}
	vector<float> &els_eSeedOverPIn()
	{
		if (not els_eSeedOverPIn_isLoaded) {
			if (els_eSeedOverPIn_branch != 0) {
				els_eSeedOverPIn_branch->GetEntry(index);
			} else { }
			els_eSeedOverPIn_isLoaded = true;
		}
		return els_eSeedOverPIn_;
	}
	vector<float> &els_eSeedOverPOut()
	{
		if (not els_eSeedOverPOut_isLoaded) {
			if (els_eSeedOverPOut_branch != 0) {
				els_eSeedOverPOut_branch->GetEntry(index);
			} else { }
			els_eSeedOverPOut_isLoaded = true;
		}
		return els_eSeedOverPOut_;
	}
	vector<float> &els_ecalEnergy()
	{
		if (not els_ecalEnergy_isLoaded) {
			if (els_ecalEnergy_branch != 0) {
				els_ecalEnergy_branch->GetEntry(index);
			} else { }
			els_ecalEnergy_isLoaded = true;
		}
		return els_ecalEnergy_;
	}
	vector<float> &els_ecalEnergyError()
	{
		if (not els_ecalEnergyError_isLoaded) {
			if (els_ecalEnergyError_branch != 0) {
				els_ecalEnergyError_branch->GetEntry(index);
			} else { }
			els_ecalEnergyError_isLoaded = true;
		}
		return els_ecalEnergyError_;
	}
	vector<float> &els_ecalIso()
	{
		if (not els_ecalIso_isLoaded) {
			if (els_ecalIso_branch != 0) {
				els_ecalIso_branch->GetEntry(index);
			} else { }
			els_ecalIso_isLoaded = true;
		}
		return els_ecalIso_;
	}
	vector<float> &els_ecalIso04()
	{
		if (not els_ecalIso04_isLoaded) {
			if (els_ecalIso04_branch != 0) {
				els_ecalIso04_branch->GetEntry(index);
			} else { }
			els_ecalIso04_isLoaded = true;
		}
		return els_ecalIso04_;
	}
	vector<float> &els_etaErr()
	{
		if (not els_etaErr_isLoaded) {
			if (els_etaErr_branch != 0) {
				els_etaErr_branch->GetEntry(index);
			} else { }
			els_etaErr_isLoaded = true;
		}
		return els_etaErr_;
	}
	vector<float> &els_etaSC()
	{
		if (not els_etaSC_isLoaded) {
			if (els_etaSC_branch != 0) {
				els_etaSC_branch->GetEntry(index);
			} else { }
			els_etaSC_isLoaded = true;
		}
		return els_etaSC_;
	}
	vector<float> &els_etaSCwidth()
	{
		if (not els_etaSCwidth_isLoaded) {
			if (els_etaSCwidth_branch != 0) {
				els_etaSCwidth_branch->GetEntry(index);
			} else { }
			els_etaSCwidth_isLoaded = true;
		}
		return els_etaSCwidth_;
	}
	vector<float> &els_fbrem()
	{
		if (not els_fbrem_isLoaded) {
			if (els_fbrem_branch != 0) {
				els_fbrem_branch->GetEntry(index);
			} else { }
			els_fbrem_isLoaded = true;
		}
		return els_fbrem_;
	}
	vector<float> &els_hOverE()
	{
		if (not els_hOverE_isLoaded) {
			if (els_hOverE_branch != 0) {
				els_hOverE_branch->GetEntry(index);
			} else { }
			els_hOverE_isLoaded = true;
		}
		return els_hOverE_;
	}
	vector<float> &els_hOverEBC()
	{
		if (not els_hOverEBC_isLoaded) {
			if (els_hOverEBC_branch != 0) {
				els_hOverEBC_branch->GetEntry(index);
			} else { }
			els_hOverEBC_isLoaded = true;
		}
		return els_hOverEBC_;
	}
	vector<float> &els_hcalDepth1OverEcal()
	{
		if (not els_hcalDepth1OverEcal_isLoaded) {
			if (els_hcalDepth1OverEcal_branch != 0) {
				els_hcalDepth1OverEcal_branch->GetEntry(index);
			} else { }
			els_hcalDepth1OverEcal_isLoaded = true;
		}
		return els_hcalDepth1OverEcal_;
	}
	vector<float> &els_hcalDepth1TowerSumEt()
	{
		if (not els_hcalDepth1TowerSumEt_isLoaded) {
			if (els_hcalDepth1TowerSumEt_branch != 0) {
				els_hcalDepth1TowerSumEt_branch->GetEntry(index);
			} else { }
			els_hcalDepth1TowerSumEt_isLoaded = true;
		}
		return els_hcalDepth1TowerSumEt_;
	}
	vector<float> &els_hcalDepth1TowerSumEt04()
	{
		if (not els_hcalDepth1TowerSumEt04_isLoaded) {
			if (els_hcalDepth1TowerSumEt04_branch != 0) {
				els_hcalDepth1TowerSumEt04_branch->GetEntry(index);
			} else { }
			els_hcalDepth1TowerSumEt04_isLoaded = true;
		}
		return els_hcalDepth1TowerSumEt04_;
	}
	vector<float> &els_hcalDepth2OverEcal()
	{
		if (not els_hcalDepth2OverEcal_isLoaded) {
			if (els_hcalDepth2OverEcal_branch != 0) {
				els_hcalDepth2OverEcal_branch->GetEntry(index);
			} else { }
			els_hcalDepth2OverEcal_isLoaded = true;
		}
		return els_hcalDepth2OverEcal_;
	}
	vector<float> &els_hcalDepth2TowerSumEt()
	{
		if (not els_hcalDepth2TowerSumEt_isLoaded) {
			if (els_hcalDepth2TowerSumEt_branch != 0) {
				els_hcalDepth2TowerSumEt_branch->GetEntry(index);
			} else { }
			els_hcalDepth2TowerSumEt_isLoaded = true;
		}
		return els_hcalDepth2TowerSumEt_;
	}
	vector<float> &els_hcalDepth2TowerSumEt04()
	{
		if (not els_hcalDepth2TowerSumEt04_isLoaded) {
			if (els_hcalDepth2TowerSumEt04_branch != 0) {
				els_hcalDepth2TowerSumEt04_branch->GetEntry(index);
			} else { }
			els_hcalDepth2TowerSumEt04_isLoaded = true;
		}
		return els_hcalDepth2TowerSumEt04_;
	}
	vector<float> &els_hcalIso()
	{
		if (not els_hcalIso_isLoaded) {
			if (els_hcalIso_branch != 0) {
				els_hcalIso_branch->GetEntry(index);
			} else { }
			els_hcalIso_isLoaded = true;
		}
		return els_hcalIso_;
	}
	vector<float> &els_hcalIso04()
	{
		if (not els_hcalIso04_isLoaded) {
			if (els_hcalIso04_branch != 0) {
				els_hcalIso04_branch->GetEntry(index);
			} else { }
			els_hcalIso04_isLoaded = true;
		}
		return els_hcalIso04_;
	}
	vector<float> &els_ip3d()
	{
		if (not els_ip3d_isLoaded) {
			if (els_ip3d_branch != 0) {
				els_ip3d_branch->GetEntry(index);
			} else { }
			els_ip3d_isLoaded = true;
		}
		return els_ip3d_;
	}
	vector<float> &els_ip3derr()
	{
		if (not els_ip3derr_isLoaded) {
			if (els_ip3derr_branch != 0) {
				els_ip3derr_branch->GetEntry(index);
			} else { }
			els_ip3derr_isLoaded = true;
		}
		return els_ip3derr_;
	}
	vector<float> &els_iso03_pf()
	{
		if (not els_iso03_pf_isLoaded) {
			if (els_iso03_pf_branch != 0) {
				els_iso03_pf_branch->GetEntry(index);
			} else { }
			els_iso03_pf_isLoaded = true;
		}
		return els_iso03_pf_;
	}
	vector<float> &els_iso03_pf2012_ch()
	{
		if (not els_iso03_pf2012_ch_isLoaded) {
			if (els_iso03_pf2012_ch_branch != 0) {
				els_iso03_pf2012_ch_branch->GetEntry(index);
			} else { }
			els_iso03_pf2012_ch_isLoaded = true;
		}
		return els_iso03_pf2012_ch_;
	}
	vector<float> &els_iso03_pf2012_em()
	{
		if (not els_iso03_pf2012_em_isLoaded) {
			if (els_iso03_pf2012_em_branch != 0) {
				els_iso03_pf2012_em_branch->GetEntry(index);
			} else { }
			els_iso03_pf2012_em_isLoaded = true;
		}
		return els_iso03_pf2012_em_;
	}
	vector<float> &els_iso03_pf2012ext_ch()
	{
		if (not els_iso03_pf2012ext_ch_isLoaded) {
			if (els_iso03_pf2012ext_ch_branch != 0) {
				els_iso03_pf2012ext_ch_branch->GetEntry(index);
			} else { }
			els_iso03_pf2012ext_ch_isLoaded = true;
		}
		return els_iso03_pf2012ext_ch_;
	}
	vector<float> &els_iso03_pf2012ext_em()
	{
		if (not els_iso03_pf2012ext_em_isLoaded) {
			if (els_iso03_pf2012ext_em_branch != 0) {
				els_iso03_pf2012ext_em_branch->GetEntry(index);
			} else { }
			els_iso03_pf2012ext_em_isLoaded = true;
		}
		return els_iso03_pf2012ext_em_;
	}
	vector<float> &els_iso03_pf2012ext_nh()
	{
		if (not els_iso03_pf2012ext_nh_isLoaded) {
			if (els_iso03_pf2012ext_nh_branch != 0) {
				els_iso03_pf2012ext_nh_branch->GetEntry(index);
			} else { }
			els_iso03_pf2012ext_nh_isLoaded = true;
		}
		return els_iso03_pf2012ext_nh_;
	}
	vector<float> &els_iso03_pf2012_nh()
	{
		if (not els_iso03_pf2012_nh_isLoaded) {
			if (els_iso03_pf2012_nh_branch != 0) {
				els_iso03_pf2012_nh_branch->GetEntry(index);
			} else { }
			els_iso03_pf2012_nh_isLoaded = true;
		}
		return els_iso03_pf2012_nh_;
	}
	vector<float> &els_iso03_pf_ch()
	{
		if (not els_iso03_pf_ch_isLoaded) {
			if (els_iso03_pf_ch_branch != 0) {
				els_iso03_pf_ch_branch->GetEntry(index);
			} else { }
			els_iso03_pf_ch_isLoaded = true;
		}
		return els_iso03_pf_ch_;
	}
	vector<float> &els_iso03_pf_gamma05()
	{
		if (not els_iso03_pf_gamma05_isLoaded) {
			if (els_iso03_pf_gamma05_branch != 0) {
				els_iso03_pf_gamma05_branch->GetEntry(index);
			} else { }
			els_iso03_pf_gamma05_isLoaded = true;
		}
		return els_iso03_pf_gamma05_;
	}
	vector<float> &els_iso03_pf_nhad05()
	{
		if (not els_iso03_pf_nhad05_isLoaded) {
			if (els_iso03_pf_nhad05_branch != 0) {
				els_iso03_pf_nhad05_branch->GetEntry(index);
			} else { }
			els_iso03_pf_nhad05_isLoaded = true;
		}
		return els_iso03_pf_nhad05_;
	}
	vector<float> &els_iso04_pf()
	{
		if (not els_iso04_pf_isLoaded) {
			if (els_iso04_pf_branch != 0) {
				els_iso04_pf_branch->GetEntry(index);
			} else { }
			els_iso04_pf_isLoaded = true;
		}
		return els_iso04_pf_;
	}
	vector<float> &els_iso04_pf2012_ch()
	{
		if (not els_iso04_pf2012_ch_isLoaded) {
			if (els_iso04_pf2012_ch_branch != 0) {
				els_iso04_pf2012_ch_branch->GetEntry(index);
			} else { }
			els_iso04_pf2012_ch_isLoaded = true;
		}
		return els_iso04_pf2012_ch_;
	}
	vector<float> &els_iso04_pf2012_em()
	{
		if (not els_iso04_pf2012_em_isLoaded) {
			if (els_iso04_pf2012_em_branch != 0) {
				els_iso04_pf2012_em_branch->GetEntry(index);
			} else { }
			els_iso04_pf2012_em_isLoaded = true;
		}
		return els_iso04_pf2012_em_;
	}
	vector<float> &els_iso04_pf2012ext_ch()
	{
		if (not els_iso04_pf2012ext_ch_isLoaded) {
			if (els_iso04_pf2012ext_ch_branch != 0) {
				els_iso04_pf2012ext_ch_branch->GetEntry(index);
			} else { }
			els_iso04_pf2012ext_ch_isLoaded = true;
		}
		return els_iso04_pf2012ext_ch_;
	}
	vector<float> &els_iso04_pf2012ext_em()
	{
		if (not els_iso04_pf2012ext_em_isLoaded) {
			if (els_iso04_pf2012ext_em_branch != 0) {
				els_iso04_pf2012ext_em_branch->GetEntry(index);
			} else { }
			els_iso04_pf2012ext_em_isLoaded = true;
		}
		return els_iso04_pf2012ext_em_;
	}
	vector<float> &els_iso04_pf2012ext_nh()
	{
		if (not els_iso04_pf2012ext_nh_isLoaded) {
			if (els_iso04_pf2012ext_nh_branch != 0) {
				els_iso04_pf2012ext_nh_branch->GetEntry(index);
			} else { }
			els_iso04_pf2012ext_nh_isLoaded = true;
		}
		return els_iso04_pf2012ext_nh_;
	}
	vector<float> &els_iso04_pf2012_nh()
	{
		if (not els_iso04_pf2012_nh_isLoaded) {
			if (els_iso04_pf2012_nh_branch != 0) {
				els_iso04_pf2012_nh_branch->GetEntry(index);
			} else { }
			els_iso04_pf2012_nh_isLoaded = true;
		}
		return els_iso04_pf2012_nh_;
	}
	vector<float> &els_iso04_pf_ch()
	{
		if (not els_iso04_pf_ch_isLoaded) {
			if (els_iso04_pf_ch_branch != 0) {
				els_iso04_pf_ch_branch->GetEntry(index);
			} else { }
			els_iso04_pf_ch_isLoaded = true;
		}
		return els_iso04_pf_ch_;
	}
	vector<float> &els_iso04_pf_gamma05()
	{
		if (not els_iso04_pf_gamma05_isLoaded) {
			if (els_iso04_pf_gamma05_branch != 0) {
				els_iso04_pf_gamma05_branch->GetEntry(index);
			} else { }
			els_iso04_pf_gamma05_isLoaded = true;
		}
		return els_iso04_pf_gamma05_;
	}
	vector<float> &els_iso04_pf_nhad05()
	{
		if (not els_iso04_pf_nhad05_isLoaded) {
			if (els_iso04_pf_nhad05_branch != 0) {
				els_iso04_pf_nhad05_branch->GetEntry(index);
			} else { }
			els_iso04_pf_nhad05_isLoaded = true;
		}
		return els_iso04_pf_nhad05_;
	}
	vector<float> &els_layer1_charge()
	{
		if (not els_layer1_charge_isLoaded) {
			if (els_layer1_charge_branch != 0) {
				els_layer1_charge_branch->GetEntry(index);
			} else { }
			els_layer1_charge_isLoaded = true;
		}
		return els_layer1_charge_;
	}
	vector<float> &els_lh()
	{
		if (not els_lh_isLoaded) {
			if (els_lh_branch != 0) {
				els_lh_branch->GetEntry(index);
			} else { }
			els_lh_isLoaded = true;
		}
		return els_lh_;
	}
	vector<float> &els_mc_patMatch_dr()
	{
		if (not els_mc_patMatch_dr_isLoaded) {
			if (els_mc_patMatch_dr_branch != 0) {
				els_mc_patMatch_dr_branch->GetEntry(index);
			} else { }
			els_mc_patMatch_dr_isLoaded = true;
		}
		return els_mc_patMatch_dr_;
	}
	vector<float> &els_mva()
	{
		if (not els_mva_isLoaded) {
			if (els_mva_branch != 0) {
				els_mva_branch->GetEntry(index);
			} else { }
			els_mva_isLoaded = true;
		}
		return els_mva_;
	}
	vector<float> &els_ndof()
	{
		if (not els_ndof_isLoaded) {
			if (els_ndof_branch != 0) {
				els_ndof_branch->GetEntry(index);
			} else { }
			els_ndof_isLoaded = true;
		}
		return els_ndof_;
	}
	vector<float> &els_pfChargedHadronIso()
	{
		if (not els_pfChargedHadronIso_isLoaded) {
			if (els_pfChargedHadronIso_branch != 0) {
				els_pfChargedHadronIso_branch->GetEntry(index);
			} else { }
			els_pfChargedHadronIso_isLoaded = true;
		}
		return els_pfChargedHadronIso_;
	}
	vector<float> &els_pfNeutralHadronIso()
	{
		if (not els_pfNeutralHadronIso_isLoaded) {
			if (els_pfNeutralHadronIso_branch != 0) {
				els_pfNeutralHadronIso_branch->GetEntry(index);
			} else { }
			els_pfNeutralHadronIso_isLoaded = true;
		}
		return els_pfNeutralHadronIso_;
	}
	vector<float> &els_pfPUIso()
	{
		if (not els_pfPUIso_isLoaded) {
			if (els_pfPUIso_branch != 0) {
				els_pfPUIso_branch->GetEntry(index);
			} else { }
			els_pfPUIso_isLoaded = true;
		}
		return els_pfPUIso_;
	}
	vector<float> &els_pfPhotonIso()
	{
		if (not els_pfPhotonIso_isLoaded) {
			if (els_pfPhotonIso_branch != 0) {
				els_pfPhotonIso_branch->GetEntry(index);
			} else { }
			els_pfPhotonIso_isLoaded = true;
		}
		return els_pfPhotonIso_;
	}
	vector<float> &els_phiErr()
	{
		if (not els_phiErr_isLoaded) {
			if (els_phiErr_branch != 0) {
				els_phiErr_branch->GetEntry(index);
			} else { }
			els_phiErr_isLoaded = true;
		}
		return els_phiErr_;
	}
	vector<float> &els_phiSC()
	{
		if (not els_phiSC_isLoaded) {
			if (els_phiSC_branch != 0) {
				els_phiSC_branch->GetEntry(index);
			} else { }
			els_phiSC_isLoaded = true;
		}
		return els_phiSC_;
	}
	vector<float> &els_phiSCwidth()
	{
		if (not els_phiSCwidth_isLoaded) {
			if (els_phiSCwidth_branch != 0) {
				els_phiSCwidth_branch->GetEntry(index);
			} else { }
			els_phiSCwidth_isLoaded = true;
		}
		return els_phiSCwidth_;
	}
	vector<float> &els_ptErr()
	{
		if (not els_ptErr_isLoaded) {
			if (els_ptErr_branch != 0) {
				els_ptErr_branch->GetEntry(index);
			} else { }
			els_ptErr_isLoaded = true;
		}
		return els_ptErr_;
	}
	vector<float> &els_ptErrGsf()
	{
		if (not els_ptErrGsf_isLoaded) {
			if (els_ptErrGsf_branch != 0) {
				els_ptErrGsf_branch->GetEntry(index);
			} else { }
			els_ptErrGsf_isLoaded = true;
		}
		return els_ptErrGsf_;
	}
	vector<float> &els_r9()
	{
		if (not els_r9_isLoaded) {
			if (els_r9_branch != 0) {
				els_r9_branch->GetEntry(index);
			} else { }
			els_r9_isLoaded = true;
		}
		return els_r9_;
	}
	vector<float> &els_r9_full5x5()
	{
		if (not els_r9_full5x5_isLoaded) {
			if (els_r9_full5x5_branch != 0) {
				els_r9_full5x5_branch->GetEntry(index);
			} else { }
			els_r9_full5x5_isLoaded = true;
		}
		return els_r9_full5x5_;
	}
	vector<float> &els_sigmaEtaEta()
	{
		if (not els_sigmaEtaEta_isLoaded) {
			if (els_sigmaEtaEta_branch != 0) {
				els_sigmaEtaEta_branch->GetEntry(index);
			} else { }
			els_sigmaEtaEta_isLoaded = true;
		}
		return els_sigmaEtaEta_;
	}
	vector<float> &els_sigmaEtaEta_full5x5()
	{
		if (not els_sigmaEtaEta_full5x5_isLoaded) {
			if (els_sigmaEtaEta_full5x5_branch != 0) {
				els_sigmaEtaEta_full5x5_branch->GetEntry(index);
			} else { }
			els_sigmaEtaEta_full5x5_isLoaded = true;
		}
		return els_sigmaEtaEta_full5x5_;
	}
	vector<float> &els_sigmaIEtaIEta()
	{
		if (not els_sigmaIEtaIEta_isLoaded) {
			if (els_sigmaIEtaIEta_branch != 0) {
				els_sigmaIEtaIEta_branch->GetEntry(index);
			} else { }
			els_sigmaIEtaIEta_isLoaded = true;
		}
		return els_sigmaIEtaIEta_;
	}
	vector<float> &els_sigmaIEtaIEtaSC()
	{
		if (not els_sigmaIEtaIEtaSC_isLoaded) {
			if (els_sigmaIEtaIEtaSC_branch != 0) {
				els_sigmaIEtaIEtaSC_branch->GetEntry(index);
			} else { }
			els_sigmaIEtaIEtaSC_isLoaded = true;
		}
		return els_sigmaIEtaIEtaSC_;
	}
	vector<float> &els_sigmaIEtaIEta_full5x5()
	{
		if (not els_sigmaIEtaIEta_full5x5_isLoaded) {
			if (els_sigmaIEtaIEta_full5x5_branch != 0) {
				els_sigmaIEtaIEta_full5x5_branch->GetEntry(index);
			} else { }
			els_sigmaIEtaIEta_full5x5_isLoaded = true;
		}
		return els_sigmaIEtaIEta_full5x5_;
	}
	vector<float> &els_sigmaIEtaIPhi()
	{
		if (not els_sigmaIEtaIPhi_isLoaded) {
			if (els_sigmaIEtaIPhi_branch != 0) {
				els_sigmaIEtaIPhi_branch->GetEntry(index);
			} else { }
			els_sigmaIEtaIPhi_isLoaded = true;
		}
		return els_sigmaIEtaIPhi_;
	}
	vector<float> &els_sigmaIPhiIPhi()
	{
		if (not els_sigmaIPhiIPhi_isLoaded) {
			if (els_sigmaIPhiIPhi_branch != 0) {
				els_sigmaIPhiIPhi_branch->GetEntry(index);
			} else { }
			els_sigmaIPhiIPhi_isLoaded = true;
		}
		return els_sigmaIPhiIPhi_;
	}
	vector<float> &els_sigmaIPhiIPhiSC()
	{
		if (not els_sigmaIPhiIPhiSC_isLoaded) {
			if (els_sigmaIPhiIPhiSC_branch != 0) {
				els_sigmaIPhiIPhiSC_branch->GetEntry(index);
			} else { }
			els_sigmaIPhiIPhiSC_isLoaded = true;
		}
		return els_sigmaIPhiIPhiSC_;
	}
	vector<float> &els_sigmaIPhiIPhi_full5x5()
	{
		if (not els_sigmaIPhiIPhi_full5x5_isLoaded) {
			if (els_sigmaIPhiIPhi_full5x5_branch != 0) {
				els_sigmaIPhiIPhi_full5x5_branch->GetEntry(index);
			} else { }
			els_sigmaIPhiIPhi_full5x5_isLoaded = true;
		}
		return els_sigmaIPhiIPhi_full5x5_;
	}
	vector<float> &els_sigmaIphiIphi()
	{
		if (not els_sigmaIphiIphi_isLoaded) {
			if (els_sigmaIphiIphi_branch != 0) {
				els_sigmaIphiIphi_branch->GetEntry(index);
			} else { }
			els_sigmaIphiIphi_isLoaded = true;
		}
		return els_sigmaIphiIphi_;
	}
	vector<float> &els_sigmaPhiPhi()
	{
		if (not els_sigmaPhiPhi_isLoaded) {
			if (els_sigmaPhiPhi_branch != 0) {
				els_sigmaPhiPhi_branch->GetEntry(index);
			} else { }
			els_sigmaPhiPhi_isLoaded = true;
		}
		return els_sigmaPhiPhi_;
	}
	vector<float> &els_tkIso()
	{
		if (not els_tkIso_isLoaded) {
			if (els_tkIso_branch != 0) {
				els_tkIso_branch->GetEntry(index);
			} else { }
			els_tkIso_isLoaded = true;
		}
		return els_tkIso_;
	}
	vector<float> &els_tkIso04()
	{
		if (not els_tkIso04_isLoaded) {
			if (els_tkIso04_branch != 0) {
				els_tkIso04_branch->GetEntry(index);
			} else { }
			els_tkIso04_isLoaded = true;
		}
		return els_tkIso04_;
	}
	vector<float> &els_trackMomentumError()
	{
		if (not els_trackMomentumError_isLoaded) {
			if (els_trackMomentumError_branch != 0) {
				els_trackMomentumError_branch->GetEntry(index);
			} else { }
			els_trackMomentumError_isLoaded = true;
		}
		return els_trackMomentumError_;
	}
	vector<float> &els_trkdr()
	{
		if (not els_trkdr_isLoaded) {
			if (els_trkdr_branch != 0) {
				els_trkdr_branch->GetEntry(index);
			} else { }
			els_trkdr_isLoaded = true;
		}
		return els_trkdr_;
	}
	vector<float> &els_trkshFrac()
	{
		if (not els_trkshFrac_isLoaded) {
			if (els_trkshFrac_branch != 0) {
				els_trkshFrac_branch->GetEntry(index);
			} else { }
			els_trkshFrac_isLoaded = true;
		}
		return els_trkshFrac_;
	}
	vector<float> &els_z0()
	{
		if (not els_z0_isLoaded) {
			if (els_z0_branch != 0) {
				els_z0_branch->GetEntry(index);
			} else { }
			els_z0_isLoaded = true;
		}
		return els_z0_;
	}
	vector<float> &els_z0Err()
	{
		if (not els_z0Err_isLoaded) {
			if (els_z0Err_branch != 0) {
				els_z0Err_branch->GetEntry(index);
			} else { }
			els_z0Err_isLoaded = true;
		}
		return els_z0Err_;
	}
	vector<float> &els_z0corr()
	{
		if (not els_z0corr_isLoaded) {
			if (els_z0corr_branch != 0) {
				els_z0corr_branch->GetEntry(index);
			} else { }
			els_z0corr_isLoaded = true;
		}
		return els_z0corr_;
	}
	vector<float> &mus_backToBackCompat()
	{
		if (not mus_backToBackCompat_isLoaded) {
			if (mus_backToBackCompat_branch != 0) {
				mus_backToBackCompat_branch->GetEntry(index);
			} else { }
			mus_backToBackCompat_isLoaded = true;
		}
		return mus_backToBackCompat_;
	}
	vector<float> &mus_best_dxyPV()
	{
		if (not mus_best_dxyPV_isLoaded) {
			if (mus_best_dxyPV_branch != 0) {
				mus_best_dxyPV_branch->GetEntry(index);
			} else { }
			mus_best_dxyPV_isLoaded = true;
		}
		return mus_best_dxyPV_;
	}
	vector<float> &mus_best_dzPV()
	{
		if (not mus_best_dzPV_isLoaded) {
			if (mus_best_dzPV_branch != 0) {
				mus_best_dzPV_branch->GetEntry(index);
			} else { }
			mus_best_dzPV_isLoaded = true;
		}
		return mus_best_dzPV_;
	}
	vector<float> &mus_caloCompatibility()
	{
		if (not mus_caloCompatibility_isLoaded) {
			if (mus_caloCompatibility_branch != 0) {
				mus_caloCompatibility_branch->GetEntry(index);
			} else { }
			mus_caloCompatibility_isLoaded = true;
		}
		return mus_caloCompatibility_;
	}
	vector<float> &mus_chi2()
	{
		if (not mus_chi2_isLoaded) {
			if (mus_chi2_branch != 0) {
				mus_chi2_branch->GetEntry(index);
			} else { }
			mus_chi2_isLoaded = true;
		}
		return mus_chi2_;
	}
	vector<float> &mus_chi2LocalMomentum()
	{
		if (not mus_chi2LocalMomentum_isLoaded) {
			if (mus_chi2LocalMomentum_branch != 0) {
				mus_chi2LocalMomentum_branch->GetEntry(index);
			} else { }
			mus_chi2LocalMomentum_isLoaded = true;
		}
		return mus_chi2LocalMomentum_;
	}
	vector<float> &mus_chi2LocalPosition()
	{
		if (not mus_chi2LocalPosition_isLoaded) {
			if (mus_chi2LocalPosition_branch != 0) {
				mus_chi2LocalPosition_branch->GetEntry(index);
			} else { }
			mus_chi2LocalPosition_isLoaded = true;
		}
		return mus_chi2LocalPosition_;
	}
	vector<float> &mus_cosmicCompat()
	{
		if (not mus_cosmicCompat_isLoaded) {
			if (mus_cosmicCompat_branch != 0) {
				mus_cosmicCompat_branch->GetEntry(index);
			} else { }
			mus_cosmicCompat_isLoaded = true;
		}
		return mus_cosmicCompat_;
	}
	vector<float> &mus_d0()
	{
		if (not mus_d0_isLoaded) {
			if (mus_d0_branch != 0) {
				mus_d0_branch->GetEntry(index);
			} else { }
			mus_d0_isLoaded = true;
		}
		return mus_d0_;
	}
	vector<float> &mus_d0Err()
	{
		if (not mus_d0Err_isLoaded) {
			if (mus_d0Err_branch != 0) {
				mus_d0Err_branch->GetEntry(index);
			} else { }
			mus_d0Err_isLoaded = true;
		}
		return mus_d0Err_;
	}
	vector<float> &mus_d0corr()
	{
		if (not mus_d0corr_isLoaded) {
			if (mus_d0corr_branch != 0) {
				mus_d0corr_branch->GetEntry(index);
			} else { }
			mus_d0corr_isLoaded = true;
		}
		return mus_d0corr_;
	}
	vector<float> &mus_d0corrPhi()
	{
		if (not mus_d0corrPhi_isLoaded) {
			if (mus_d0corrPhi_branch != 0) {
				mus_d0corrPhi_branch->GetEntry(index);
			} else { }
			mus_d0corrPhi_isLoaded = true;
		}
		return mus_d0corrPhi_;
	}
	vector<float> &mus_d0phiCov()
	{
		if (not mus_d0phiCov_isLoaded) {
			if (mus_d0phiCov_branch != 0) {
				mus_d0phiCov_branch->GetEntry(index);
			} else { }
			mus_d0phiCov_isLoaded = true;
		}
		return mus_d0phiCov_;
	}
	vector<float> &mus_dxyPV()
	{
		if (not mus_dxyPV_isLoaded) {
			if (mus_dxyPV_branch != 0) {
				mus_dxyPV_branch->GetEntry(index);
			} else { }
			mus_dxyPV_isLoaded = true;
		}
		return mus_dxyPV_;
	}
	vector<float> &mus_dzPV()
	{
		if (not mus_dzPV_isLoaded) {
			if (mus_dzPV_branch != 0) {
				mus_dzPV_branch->GetEntry(index);
			} else { }
			mus_dzPV_isLoaded = true;
		}
		return mus_dzPV_;
	}
	vector<float> &mus_ecal_time()
	{
		if (not mus_ecal_time_isLoaded) {
			if (mus_ecal_time_branch != 0) {
				mus_ecal_time_branch->GetEntry(index);
			} else { }
			mus_ecal_time_isLoaded = true;
		}
		return mus_ecal_time_;
	}
	vector<float> &mus_ecal_timeError()
	{
		if (not mus_ecal_timeError_isLoaded) {
			if (mus_ecal_timeError_branch != 0) {
				mus_ecal_timeError_branch->GetEntry(index);
			} else { }
			mus_ecal_timeError_isLoaded = true;
		}
		return mus_ecal_timeError_;
	}
	vector<float> &mus_e_em()
	{
		if (not mus_e_em_isLoaded) {
			if (mus_e_em_branch != 0) {
				mus_e_em_branch->GetEntry(index);
			} else { }
			mus_e_em_isLoaded = true;
		}
		return mus_e_em_;
	}
	vector<float> &mus_e_emS9()
	{
		if (not mus_e_emS9_isLoaded) {
			if (mus_e_emS9_branch != 0) {
				mus_e_emS9_branch->GetEntry(index);
			} else { }
			mus_e_emS9_isLoaded = true;
		}
		return mus_e_emS9_;
	}
	vector<float> &mus_e_had()
	{
		if (not mus_e_had_isLoaded) {
			if (mus_e_had_branch != 0) {
				mus_e_had_branch->GetEntry(index);
			} else { }
			mus_e_had_isLoaded = true;
		}
		return mus_e_had_;
	}
	vector<float> &mus_e_hadS9()
	{
		if (not mus_e_hadS9_isLoaded) {
			if (mus_e_hadS9_branch != 0) {
				mus_e_hadS9_branch->GetEntry(index);
			} else { }
			mus_e_hadS9_isLoaded = true;
		}
		return mus_e_hadS9_;
	}
	vector<float> &mus_e_ho()
	{
		if (not mus_e_ho_isLoaded) {
			if (mus_e_ho_branch != 0) {
				mus_e_ho_branch->GetEntry(index);
			} else { }
			mus_e_ho_isLoaded = true;
		}
		return mus_e_ho_;
	}
	vector<float> &mus_e_hoS9()
	{
		if (not mus_e_hoS9_isLoaded) {
			if (mus_e_hoS9_branch != 0) {
				mus_e_hoS9_branch->GetEntry(index);
			} else { }
			mus_e_hoS9_isLoaded = true;
		}
		return mus_e_hoS9_;
	}
	vector<float> &mus_emMax()
	{
		if (not mus_emMax_isLoaded) {
			if (mus_emMax_branch != 0) {
				mus_emMax_branch->GetEntry(index);
			} else { }
			mus_emMax_isLoaded = true;
		}
		return mus_emMax_;
	}
	vector<float> &mus_emS25()
	{
		if (not mus_emS25_isLoaded) {
			if (mus_emS25_branch != 0) {
				mus_emS25_branch->GetEntry(index);
			} else { }
			mus_emS25_isLoaded = true;
		}
		return mus_emS25_;
	}
	vector<float> &mus_etaErr()
	{
		if (not mus_etaErr_isLoaded) {
			if (mus_etaErr_branch != 0) {
				mus_etaErr_branch->GetEntry(index);
			} else { }
			mus_etaErr_isLoaded = true;
		}
		return mus_etaErr_;
	}
	vector<float> &mus_gfit_chi2()
	{
		if (not mus_gfit_chi2_isLoaded) {
			if (mus_gfit_chi2_branch != 0) {
				mus_gfit_chi2_branch->GetEntry(index);
			} else { }
			mus_gfit_chi2_isLoaded = true;
		}
		return mus_gfit_chi2_;
	}
	vector<float> &mus_gfit_d0()
	{
		if (not mus_gfit_d0_isLoaded) {
			if (mus_gfit_d0_branch != 0) {
				mus_gfit_d0_branch->GetEntry(index);
			} else { }
			mus_gfit_d0_isLoaded = true;
		}
		return mus_gfit_d0_;
	}
	vector<float> &mus_gfit_d0Err()
	{
		if (not mus_gfit_d0Err_isLoaded) {
			if (mus_gfit_d0Err_branch != 0) {
				mus_gfit_d0Err_branch->GetEntry(index);
			} else { }
			mus_gfit_d0Err_isLoaded = true;
		}
		return mus_gfit_d0Err_;
	}
	vector<float> &mus_gfit_d0corr()
	{
		if (not mus_gfit_d0corr_isLoaded) {
			if (mus_gfit_d0corr_branch != 0) {
				mus_gfit_d0corr_branch->GetEntry(index);
			} else { }
			mus_gfit_d0corr_isLoaded = true;
		}
		return mus_gfit_d0corr_;
	}
	vector<float> &mus_gfit_d0corrPhi()
	{
		if (not mus_gfit_d0corrPhi_isLoaded) {
			if (mus_gfit_d0corrPhi_branch != 0) {
				mus_gfit_d0corrPhi_branch->GetEntry(index);
			} else { }
			mus_gfit_d0corrPhi_isLoaded = true;
		}
		return mus_gfit_d0corrPhi_;
	}
	vector<float> &mus_gfit_d0phiCov()
	{
		if (not mus_gfit_d0phiCov_isLoaded) {
			if (mus_gfit_d0phiCov_branch != 0) {
				mus_gfit_d0phiCov_branch->GetEntry(index);
			} else { }
			mus_gfit_d0phiCov_isLoaded = true;
		}
		return mus_gfit_d0phiCov_;
	}
	vector<float> &mus_gfit_ndof()
	{
		if (not mus_gfit_ndof_isLoaded) {
			if (mus_gfit_ndof_branch != 0) {
				mus_gfit_ndof_branch->GetEntry(index);
			} else { }
			mus_gfit_ndof_isLoaded = true;
		}
		return mus_gfit_ndof_;
	}
	vector<float> &mus_gfit_qoverp()
	{
		if (not mus_gfit_qoverp_isLoaded) {
			if (mus_gfit_qoverp_branch != 0) {
				mus_gfit_qoverp_branch->GetEntry(index);
			} else { }
			mus_gfit_qoverp_isLoaded = true;
		}
		return mus_gfit_qoverp_;
	}
	vector<float> &mus_gfit_qoverpError()
	{
		if (not mus_gfit_qoverpError_isLoaded) {
			if (mus_gfit_qoverpError_branch != 0) {
				mus_gfit_qoverpError_branch->GetEntry(index);
			} else { }
			mus_gfit_qoverpError_isLoaded = true;
		}
		return mus_gfit_qoverpError_;
	}
	vector<float> &mus_gfit_z0()
	{
		if (not mus_gfit_z0_isLoaded) {
			if (mus_gfit_z0_branch != 0) {
				mus_gfit_z0_branch->GetEntry(index);
			} else { }
			mus_gfit_z0_isLoaded = true;
		}
		return mus_gfit_z0_;
	}
	vector<float> &mus_gfit_z0Err()
	{
		if (not mus_gfit_z0Err_isLoaded) {
			if (mus_gfit_z0Err_branch != 0) {
				mus_gfit_z0Err_branch->GetEntry(index);
			} else { }
			mus_gfit_z0Err_isLoaded = true;
		}
		return mus_gfit_z0Err_;
	}
	vector<float> &mus_gfit_z0corr()
	{
		if (not mus_gfit_z0corr_isLoaded) {
			if (mus_gfit_z0corr_branch != 0) {
				mus_gfit_z0corr_branch->GetEntry(index);
			} else { }
			mus_gfit_z0corr_isLoaded = true;
		}
		return mus_gfit_z0corr_;
	}
	vector<float> &mus_glbKink()
	{
		if (not mus_glbKink_isLoaded) {
			if (mus_glbKink_branch != 0) {
				mus_glbKink_branch->GetEntry(index);
			} else { }
			mus_glbKink_isLoaded = true;
		}
		return mus_glbKink_;
	}
	vector<float> &mus_glbTrackProbability()
	{
		if (not mus_glbTrackProbability_isLoaded) {
			if (mus_glbTrackProbability_branch != 0) {
				mus_glbTrackProbability_branch->GetEntry(index);
			} else { }
			mus_glbTrackProbability_isLoaded = true;
		}
		return mus_glbTrackProbability_;
	}
	vector<float> &mus_globalDeltaEtaPhi()
	{
		if (not mus_globalDeltaEtaPhi_isLoaded) {
			if (mus_globalDeltaEtaPhi_branch != 0) {
				mus_globalDeltaEtaPhi_branch->GetEntry(index);
			} else { }
			mus_globalDeltaEtaPhi_isLoaded = true;
		}
		return mus_globalDeltaEtaPhi_;
	}
	vector<float> &mus_hadMax()
	{
		if (not mus_hadMax_isLoaded) {
			if (mus_hadMax_branch != 0) {
				mus_hadMax_branch->GetEntry(index);
			} else { }
			mus_hadMax_isLoaded = true;
		}
		return mus_hadMax_;
	}
	vector<float> &mus_hcal_time()
	{
		if (not mus_hcal_time_isLoaded) {
			if (mus_hcal_time_branch != 0) {
				mus_hcal_time_branch->GetEntry(index);
			} else { }
			mus_hcal_time_isLoaded = true;
		}
		return mus_hcal_time_;
	}
	vector<float> &mus_hcal_timeError()
	{
		if (not mus_hcal_timeError_isLoaded) {
			if (mus_hcal_timeError_branch != 0) {
				mus_hcal_timeError_branch->GetEntry(index);
			} else { }
			mus_hcal_timeError_isLoaded = true;
		}
		return mus_hcal_timeError_;
	}
	vector<float> &mus_ip3d()
	{
		if (not mus_ip3d_isLoaded) {
			if (mus_ip3d_branch != 0) {
				mus_ip3d_branch->GetEntry(index);
			} else { }
			mus_ip3d_isLoaded = true;
		}
		return mus_ip3d_;
	}
	vector<float> &mus_ip3derr()
	{
		if (not mus_ip3derr_isLoaded) {
			if (mus_ip3derr_branch != 0) {
				mus_ip3derr_branch->GetEntry(index);
			} else { }
			mus_ip3derr_isLoaded = true;
		}
		return mus_ip3derr_;
	}
	vector<float> &mus_iso03_emEt()
	{
		if (not mus_iso03_emEt_isLoaded) {
			if (mus_iso03_emEt_branch != 0) {
				mus_iso03_emEt_branch->GetEntry(index);
			} else { }
			mus_iso03_emEt_isLoaded = true;
		}
		return mus_iso03_emEt_;
	}
	vector<float> &mus_iso03_hadEt()
	{
		if (not mus_iso03_hadEt_isLoaded) {
			if (mus_iso03_hadEt_branch != 0) {
				mus_iso03_hadEt_branch->GetEntry(index);
			} else { }
			mus_iso03_hadEt_isLoaded = true;
		}
		return mus_iso03_hadEt_;
	}
	vector<float> &mus_iso03_hoEt()
	{
		if (not mus_iso03_hoEt_isLoaded) {
			if (mus_iso03_hoEt_branch != 0) {
				mus_iso03_hoEt_branch->GetEntry(index);
			} else { }
			mus_iso03_hoEt_isLoaded = true;
		}
		return mus_iso03_hoEt_;
	}
	vector<float> &mus_iso03_sumPt()
	{
		if (not mus_iso03_sumPt_isLoaded) {
			if (mus_iso03_sumPt_branch != 0) {
				mus_iso03_sumPt_branch->GetEntry(index);
			} else { }
			mus_iso03_sumPt_isLoaded = true;
		}
		return mus_iso03_sumPt_;
	}
	vector<float> &mus_iso05_emEt()
	{
		if (not mus_iso05_emEt_isLoaded) {
			if (mus_iso05_emEt_branch != 0) {
				mus_iso05_emEt_branch->GetEntry(index);
			} else { }
			mus_iso05_emEt_isLoaded = true;
		}
		return mus_iso05_emEt_;
	}
	vector<float> &mus_iso05_hadEt()
	{
		if (not mus_iso05_hadEt_isLoaded) {
			if (mus_iso05_hadEt_branch != 0) {
				mus_iso05_hadEt_branch->GetEntry(index);
			} else { }
			mus_iso05_hadEt_isLoaded = true;
		}
		return mus_iso05_hadEt_;
	}
	vector<float> &mus_iso05_hoEt()
	{
		if (not mus_iso05_hoEt_isLoaded) {
			if (mus_iso05_hoEt_branch != 0) {
				mus_iso05_hoEt_branch->GetEntry(index);
			} else { }
			mus_iso05_hoEt_isLoaded = true;
		}
		return mus_iso05_hoEt_;
	}
	vector<float> &mus_iso05_sumPt()
	{
		if (not mus_iso05_sumPt_isLoaded) {
			if (mus_iso05_sumPt_branch != 0) {
				mus_iso05_sumPt_branch->GetEntry(index);
			} else { }
			mus_iso05_sumPt_isLoaded = true;
		}
		return mus_iso05_sumPt_;
	}
	vector<float> &mus_isoMeanDRR03_pf_ChargedHadronPt()
	{
		if (not mus_isoMeanDRR03_pf_ChargedHadronPt_isLoaded) {
			if (mus_isoMeanDRR03_pf_ChargedHadronPt_branch != 0) {
				mus_isoMeanDRR03_pf_ChargedHadronPt_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR03_pf_ChargedHadronPt_isLoaded = true;
		}
		return mus_isoMeanDRR03_pf_ChargedHadronPt_;
	}
	vector<float> &mus_isoMeanDRR03_pf_ChargedParticlePt()
	{
		if (not mus_isoMeanDRR03_pf_ChargedParticlePt_isLoaded) {
			if (mus_isoMeanDRR03_pf_ChargedParticlePt_branch != 0) {
				mus_isoMeanDRR03_pf_ChargedParticlePt_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR03_pf_ChargedParticlePt_isLoaded = true;
		}
		return mus_isoMeanDRR03_pf_ChargedParticlePt_;
	}
	vector<float> &mus_isoMeanDRR03_pf_NeutralHadronEt()
	{
		if (not mus_isoMeanDRR03_pf_NeutralHadronEt_isLoaded) {
			if (mus_isoMeanDRR03_pf_NeutralHadronEt_branch != 0) {
				mus_isoMeanDRR03_pf_NeutralHadronEt_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR03_pf_NeutralHadronEt_isLoaded = true;
		}
		return mus_isoMeanDRR03_pf_NeutralHadronEt_;
	}
	vector<float> &mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold()
	{
		if (not mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_isLoaded) {
			if (mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_branch != 0) {
				mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_isLoaded = true;
		}
		return mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold_;
	}
	vector<float> &mus_isoMeanDRR03_pf_PUPt()
	{
		if (not mus_isoMeanDRR03_pf_PUPt_isLoaded) {
			if (mus_isoMeanDRR03_pf_PUPt_branch != 0) {
				mus_isoMeanDRR03_pf_PUPt_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR03_pf_PUPt_isLoaded = true;
		}
		return mus_isoMeanDRR03_pf_PUPt_;
	}
	vector<float> &mus_isoMeanDRR03_pf_PhotonEt()
	{
		if (not mus_isoMeanDRR03_pf_PhotonEt_isLoaded) {
			if (mus_isoMeanDRR03_pf_PhotonEt_branch != 0) {
				mus_isoMeanDRR03_pf_PhotonEt_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR03_pf_PhotonEt_isLoaded = true;
		}
		return mus_isoMeanDRR03_pf_PhotonEt_;
	}
	vector<float> &mus_isoMeanDRR03_pf_PhotonEtHighThreshold()
	{
		if (not mus_isoMeanDRR03_pf_PhotonEtHighThreshold_isLoaded) {
			if (mus_isoMeanDRR03_pf_PhotonEtHighThreshold_branch != 0) {
				mus_isoMeanDRR03_pf_PhotonEtHighThreshold_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR03_pf_PhotonEtHighThreshold_isLoaded = true;
		}
		return mus_isoMeanDRR03_pf_PhotonEtHighThreshold_;
	}
	vector<float> &mus_isoMeanDRR04_pf_ChargedHadronPt()
	{
		if (not mus_isoMeanDRR04_pf_ChargedHadronPt_isLoaded) {
			if (mus_isoMeanDRR04_pf_ChargedHadronPt_branch != 0) {
				mus_isoMeanDRR04_pf_ChargedHadronPt_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR04_pf_ChargedHadronPt_isLoaded = true;
		}
		return mus_isoMeanDRR04_pf_ChargedHadronPt_;
	}
	vector<float> &mus_isoMeanDRR04_pf_ChargedParticlePt()
	{
		if (not mus_isoMeanDRR04_pf_ChargedParticlePt_isLoaded) {
			if (mus_isoMeanDRR04_pf_ChargedParticlePt_branch != 0) {
				mus_isoMeanDRR04_pf_ChargedParticlePt_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR04_pf_ChargedParticlePt_isLoaded = true;
		}
		return mus_isoMeanDRR04_pf_ChargedParticlePt_;
	}
	vector<float> &mus_isoMeanDRR04_pf_NeutralHadronEt()
	{
		if (not mus_isoMeanDRR04_pf_NeutralHadronEt_isLoaded) {
			if (mus_isoMeanDRR04_pf_NeutralHadronEt_branch != 0) {
				mus_isoMeanDRR04_pf_NeutralHadronEt_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR04_pf_NeutralHadronEt_isLoaded = true;
		}
		return mus_isoMeanDRR04_pf_NeutralHadronEt_;
	}
	vector<float> &mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold()
	{
		if (not mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_isLoaded) {
			if (mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_branch != 0) {
				mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_isLoaded = true;
		}
		return mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold_;
	}
	vector<float> &mus_isoMeanDRR04_pf_PUPt()
	{
		if (not mus_isoMeanDRR04_pf_PUPt_isLoaded) {
			if (mus_isoMeanDRR04_pf_PUPt_branch != 0) {
				mus_isoMeanDRR04_pf_PUPt_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR04_pf_PUPt_isLoaded = true;
		}
		return mus_isoMeanDRR04_pf_PUPt_;
	}
	vector<float> &mus_isoMeanDRR04_pf_PhotonEt()
	{
		if (not mus_isoMeanDRR04_pf_PhotonEt_isLoaded) {
			if (mus_isoMeanDRR04_pf_PhotonEt_branch != 0) {
				mus_isoMeanDRR04_pf_PhotonEt_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR04_pf_PhotonEt_isLoaded = true;
		}
		return mus_isoMeanDRR04_pf_PhotonEt_;
	}
	vector<float> &mus_isoMeanDRR04_pf_PhotonEtHighThreshold()
	{
		if (not mus_isoMeanDRR04_pf_PhotonEtHighThreshold_isLoaded) {
			if (mus_isoMeanDRR04_pf_PhotonEtHighThreshold_branch != 0) {
				mus_isoMeanDRR04_pf_PhotonEtHighThreshold_branch->GetEntry(index);
			} else { }
			mus_isoMeanDRR04_pf_PhotonEtHighThreshold_isLoaded = true;
		}
		return mus_isoMeanDRR04_pf_PhotonEtHighThreshold_;
	}
	vector<float> &mus_isoR03_pf_ChargedHadronPt()
	{
		if (not mus_isoR03_pf_ChargedHadronPt_isLoaded) {
			if (mus_isoR03_pf_ChargedHadronPt_branch != 0) {
				mus_isoR03_pf_ChargedHadronPt_branch->GetEntry(index);
			} else { }
			mus_isoR03_pf_ChargedHadronPt_isLoaded = true;
		}
		return mus_isoR03_pf_ChargedHadronPt_;
	}
	vector<float> &mus_isoR03_pf_ChargedParticlePt()
	{
		if (not mus_isoR03_pf_ChargedParticlePt_isLoaded) {
			if (mus_isoR03_pf_ChargedParticlePt_branch != 0) {
				mus_isoR03_pf_ChargedParticlePt_branch->GetEntry(index);
			} else { }
			mus_isoR03_pf_ChargedParticlePt_isLoaded = true;
		}
		return mus_isoR03_pf_ChargedParticlePt_;
	}
	vector<float> &mus_isoR03_pf_NeutralHadronEt()
	{
		if (not mus_isoR03_pf_NeutralHadronEt_isLoaded) {
			if (mus_isoR03_pf_NeutralHadronEt_branch != 0) {
				mus_isoR03_pf_NeutralHadronEt_branch->GetEntry(index);
			} else { }
			mus_isoR03_pf_NeutralHadronEt_isLoaded = true;
		}
		return mus_isoR03_pf_NeutralHadronEt_;
	}
	vector<float> &mus_isoR03_pf_NeutralHadronEtHighThreshold()
	{
		if (not mus_isoR03_pf_NeutralHadronEtHighThreshold_isLoaded) {
			if (mus_isoR03_pf_NeutralHadronEtHighThreshold_branch != 0) {
				mus_isoR03_pf_NeutralHadronEtHighThreshold_branch->GetEntry(index);
			} else { }
			mus_isoR03_pf_NeutralHadronEtHighThreshold_isLoaded = true;
		}
		return mus_isoR03_pf_NeutralHadronEtHighThreshold_;
	}
	vector<float> &mus_isoR03_pf_PUPt()
	{
		if (not mus_isoR03_pf_PUPt_isLoaded) {
			if (mus_isoR03_pf_PUPt_branch != 0) {
				mus_isoR03_pf_PUPt_branch->GetEntry(index);
			} else { }
			mus_isoR03_pf_PUPt_isLoaded = true;
		}
		return mus_isoR03_pf_PUPt_;
	}
	vector<float> &mus_isoR03_pf_PhotonEt()
	{
		if (not mus_isoR03_pf_PhotonEt_isLoaded) {
			if (mus_isoR03_pf_PhotonEt_branch != 0) {
				mus_isoR03_pf_PhotonEt_branch->GetEntry(index);
			} else { }
			mus_isoR03_pf_PhotonEt_isLoaded = true;
		}
		return mus_isoR03_pf_PhotonEt_;
	}
	vector<float> &mus_isoR03_pf_PhotonEtHighThreshold()
	{
		if (not mus_isoR03_pf_PhotonEtHighThreshold_isLoaded) {
			if (mus_isoR03_pf_PhotonEtHighThreshold_branch != 0) {
				mus_isoR03_pf_PhotonEtHighThreshold_branch->GetEntry(index);
			} else { }
			mus_isoR03_pf_PhotonEtHighThreshold_isLoaded = true;
		}
		return mus_isoR03_pf_PhotonEtHighThreshold_;
	}
	vector<float> &mus_isoR04_pf_ChargedHadronPt()
	{
		if (not mus_isoR04_pf_ChargedHadronPt_isLoaded) {
			if (mus_isoR04_pf_ChargedHadronPt_branch != 0) {
				mus_isoR04_pf_ChargedHadronPt_branch->GetEntry(index);
			} else { }
			mus_isoR04_pf_ChargedHadronPt_isLoaded = true;
		}
		return mus_isoR04_pf_ChargedHadronPt_;
	}
	vector<float> &mus_isoR04_pf_ChargedParticlePt()
	{
		if (not mus_isoR04_pf_ChargedParticlePt_isLoaded) {
			if (mus_isoR04_pf_ChargedParticlePt_branch != 0) {
				mus_isoR04_pf_ChargedParticlePt_branch->GetEntry(index);
			} else { }
			mus_isoR04_pf_ChargedParticlePt_isLoaded = true;
		}
		return mus_isoR04_pf_ChargedParticlePt_;
	}
	vector<float> &mus_isoR04_pf_NeutralHadronEt()
	{
		if (not mus_isoR04_pf_NeutralHadronEt_isLoaded) {
			if (mus_isoR04_pf_NeutralHadronEt_branch != 0) {
				mus_isoR04_pf_NeutralHadronEt_branch->GetEntry(index);
			} else { }
			mus_isoR04_pf_NeutralHadronEt_isLoaded = true;
		}
		return mus_isoR04_pf_NeutralHadronEt_;
	}
	vector<float> &mus_isoR04_pf_NeutralHadronEtHighThreshold()
	{
		if (not mus_isoR04_pf_NeutralHadronEtHighThreshold_isLoaded) {
			if (mus_isoR04_pf_NeutralHadronEtHighThreshold_branch != 0) {
				mus_isoR04_pf_NeutralHadronEtHighThreshold_branch->GetEntry(index);
			} else { }
			mus_isoR04_pf_NeutralHadronEtHighThreshold_isLoaded = true;
		}
		return mus_isoR04_pf_NeutralHadronEtHighThreshold_;
	}
	vector<float> &mus_isoR04_pf_PUPt()
	{
		if (not mus_isoR04_pf_PUPt_isLoaded) {
			if (mus_isoR04_pf_PUPt_branch != 0) {
				mus_isoR04_pf_PUPt_branch->GetEntry(index);
			} else { }
			mus_isoR04_pf_PUPt_isLoaded = true;
		}
		return mus_isoR04_pf_PUPt_;
	}
	vector<float> &mus_isoR04_pf_PhotonEt()
	{
		if (not mus_isoR04_pf_PhotonEt_isLoaded) {
			if (mus_isoR04_pf_PhotonEt_branch != 0) {
				mus_isoR04_pf_PhotonEt_branch->GetEntry(index);
			} else { }
			mus_isoR04_pf_PhotonEt_isLoaded = true;
		}
		return mus_isoR04_pf_PhotonEt_;
	}
	vector<float> &mus_isoR04_pf_PhotonEtHighThreshold()
	{
		if (not mus_isoR04_pf_PhotonEtHighThreshold_isLoaded) {
			if (mus_isoR04_pf_PhotonEtHighThreshold_branch != 0) {
				mus_isoR04_pf_PhotonEtHighThreshold_branch->GetEntry(index);
			} else { }
			mus_isoR04_pf_PhotonEtHighThreshold_isLoaded = true;
		}
		return mus_isoR04_pf_PhotonEtHighThreshold_;
	}
	vector<float> &mus_isoSumDRR03_pf_ChargedHadronPt()
	{
		if (not mus_isoSumDRR03_pf_ChargedHadronPt_isLoaded) {
			if (mus_isoSumDRR03_pf_ChargedHadronPt_branch != 0) {
				mus_isoSumDRR03_pf_ChargedHadronPt_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR03_pf_ChargedHadronPt_isLoaded = true;
		}
		return mus_isoSumDRR03_pf_ChargedHadronPt_;
	}
	vector<float> &mus_isoSumDRR03_pf_ChargedParticlePt()
	{
		if (not mus_isoSumDRR03_pf_ChargedParticlePt_isLoaded) {
			if (mus_isoSumDRR03_pf_ChargedParticlePt_branch != 0) {
				mus_isoSumDRR03_pf_ChargedParticlePt_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR03_pf_ChargedParticlePt_isLoaded = true;
		}
		return mus_isoSumDRR03_pf_ChargedParticlePt_;
	}
	vector<float> &mus_isoSumDRR03_pf_NeutralHadronEt()
	{
		if (not mus_isoSumDRR03_pf_NeutralHadronEt_isLoaded) {
			if (mus_isoSumDRR03_pf_NeutralHadronEt_branch != 0) {
				mus_isoSumDRR03_pf_NeutralHadronEt_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR03_pf_NeutralHadronEt_isLoaded = true;
		}
		return mus_isoSumDRR03_pf_NeutralHadronEt_;
	}
	vector<float> &mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold()
	{
		if (not mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_isLoaded) {
			if (mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_branch != 0) {
				mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_isLoaded = true;
		}
		return mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold_;
	}
	vector<float> &mus_isoSumDRR03_pf_PUPt()
	{
		if (not mus_isoSumDRR03_pf_PUPt_isLoaded) {
			if (mus_isoSumDRR03_pf_PUPt_branch != 0) {
				mus_isoSumDRR03_pf_PUPt_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR03_pf_PUPt_isLoaded = true;
		}
		return mus_isoSumDRR03_pf_PUPt_;
	}
	vector<float> &mus_isoSumDRR03_pf_PhotonEt()
	{
		if (not mus_isoSumDRR03_pf_PhotonEt_isLoaded) {
			if (mus_isoSumDRR03_pf_PhotonEt_branch != 0) {
				mus_isoSumDRR03_pf_PhotonEt_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR03_pf_PhotonEt_isLoaded = true;
		}
		return mus_isoSumDRR03_pf_PhotonEt_;
	}
	vector<float> &mus_isoSumDRR03_pf_PhotonEtHighThreshold()
	{
		if (not mus_isoSumDRR03_pf_PhotonEtHighThreshold_isLoaded) {
			if (mus_isoSumDRR03_pf_PhotonEtHighThreshold_branch != 0) {
				mus_isoSumDRR03_pf_PhotonEtHighThreshold_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR03_pf_PhotonEtHighThreshold_isLoaded = true;
		}
		return mus_isoSumDRR03_pf_PhotonEtHighThreshold_;
	}
	vector<float> &mus_isoSumDRR04_pf_ChargedHadronPt()
	{
		if (not mus_isoSumDRR04_pf_ChargedHadronPt_isLoaded) {
			if (mus_isoSumDRR04_pf_ChargedHadronPt_branch != 0) {
				mus_isoSumDRR04_pf_ChargedHadronPt_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR04_pf_ChargedHadronPt_isLoaded = true;
		}
		return mus_isoSumDRR04_pf_ChargedHadronPt_;
	}
	vector<float> &mus_isoSumDRR04_pf_ChargedParticlePt()
	{
		if (not mus_isoSumDRR04_pf_ChargedParticlePt_isLoaded) {
			if (mus_isoSumDRR04_pf_ChargedParticlePt_branch != 0) {
				mus_isoSumDRR04_pf_ChargedParticlePt_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR04_pf_ChargedParticlePt_isLoaded = true;
		}
		return mus_isoSumDRR04_pf_ChargedParticlePt_;
	}
	vector<float> &mus_isoSumDRR04_pf_NeutralHadronEt()
	{
		if (not mus_isoSumDRR04_pf_NeutralHadronEt_isLoaded) {
			if (mus_isoSumDRR04_pf_NeutralHadronEt_branch != 0) {
				mus_isoSumDRR04_pf_NeutralHadronEt_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR04_pf_NeutralHadronEt_isLoaded = true;
		}
		return mus_isoSumDRR04_pf_NeutralHadronEt_;
	}
	vector<float> &mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold()
	{
		if (not mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_isLoaded) {
			if (mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_branch != 0) {
				mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_isLoaded = true;
		}
		return mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold_;
	}
	vector<float> &mus_isoSumDRR04_pf_PUPt()
	{
		if (not mus_isoSumDRR04_pf_PUPt_isLoaded) {
			if (mus_isoSumDRR04_pf_PUPt_branch != 0) {
				mus_isoSumDRR04_pf_PUPt_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR04_pf_PUPt_isLoaded = true;
		}
		return mus_isoSumDRR04_pf_PUPt_;
	}
	vector<float> &mus_isoSumDRR04_pf_PhotonEt()
	{
		if (not mus_isoSumDRR04_pf_PhotonEt_isLoaded) {
			if (mus_isoSumDRR04_pf_PhotonEt_branch != 0) {
				mus_isoSumDRR04_pf_PhotonEt_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR04_pf_PhotonEt_isLoaded = true;
		}
		return mus_isoSumDRR04_pf_PhotonEt_;
	}
	vector<float> &mus_isoSumDRR04_pf_PhotonEtHighThreshold()
	{
		if (not mus_isoSumDRR04_pf_PhotonEtHighThreshold_isLoaded) {
			if (mus_isoSumDRR04_pf_PhotonEtHighThreshold_branch != 0) {
				mus_isoSumDRR04_pf_PhotonEtHighThreshold_branch->GetEntry(index);
			} else { }
			mus_isoSumDRR04_pf_PhotonEtHighThreshold_isLoaded = true;
		}
		return mus_isoSumDRR04_pf_PhotonEtHighThreshold_;
	}
	vector<float> &mus_iso_ecalvetoDep()
	{
		if (not mus_iso_ecalvetoDep_isLoaded) {
			if (mus_iso_ecalvetoDep_branch != 0) {
				mus_iso_ecalvetoDep_branch->GetEntry(index);
			} else { }
			mus_iso_ecalvetoDep_isLoaded = true;
		}
		return mus_iso_ecalvetoDep_;
	}
	vector<float> &mus_iso_hcalvetoDep()
	{
		if (not mus_iso_hcalvetoDep_isLoaded) {
			if (mus_iso_hcalvetoDep_branch != 0) {
				mus_iso_hcalvetoDep_branch->GetEntry(index);
			} else { }
			mus_iso_hcalvetoDep_isLoaded = true;
		}
		return mus_iso_hcalvetoDep_;
	}
	vector<float> &mus_iso_hovetoDep()
	{
		if (not mus_iso_hovetoDep_isLoaded) {
			if (mus_iso_hovetoDep_branch != 0) {
				mus_iso_hovetoDep_branch->GetEntry(index);
			} else { }
			mus_iso_hovetoDep_isLoaded = true;
		}
		return mus_iso_hovetoDep_;
	}
	vector<float> &mus_iso_trckvetoDep()
	{
		if (not mus_iso_trckvetoDep_isLoaded) {
			if (mus_iso_trckvetoDep_branch != 0) {
				mus_iso_trckvetoDep_branch->GetEntry(index);
			} else { }
			mus_iso_trckvetoDep_isLoaded = true;
		}
		return mus_iso_trckvetoDep_;
	}
	vector<float> &mus_localDistance()
	{
		if (not mus_localDistance_isLoaded) {
			if (mus_localDistance_branch != 0) {
				mus_localDistance_branch->GetEntry(index);
			} else { }
			mus_localDistance_isLoaded = true;
		}
		return mus_localDistance_;
	}
	vector<float> &mus_mc_patMatch_dr()
	{
		if (not mus_mc_patMatch_dr_isLoaded) {
			if (mus_mc_patMatch_dr_branch != 0) {
				mus_mc_patMatch_dr_branch->GetEntry(index);
			} else { }
			mus_mc_patMatch_dr_isLoaded = true;
		}
		return mus_mc_patMatch_dr_;
	}
	vector<float> &mus_ndof()
	{
		if (not mus_ndof_isLoaded) {
			if (mus_ndof_branch != 0) {
				mus_ndof_branch->GetEntry(index);
			} else { }
			mus_ndof_isLoaded = true;
		}
		return mus_ndof_;
	}
	vector<float> &mus_overlapCompat()
	{
		if (not mus_overlapCompat_isLoaded) {
			if (mus_overlapCompat_branch != 0) {
				mus_overlapCompat_branch->GetEntry(index);
			} else { }
			mus_overlapCompat_isLoaded = true;
		}
		return mus_overlapCompat_;
	}
	vector<float> &mus_phiErr()
	{
		if (not mus_phiErr_isLoaded) {
			if (mus_phiErr_branch != 0) {
				mus_phiErr_branch->GetEntry(index);
			} else { }
			mus_phiErr_isLoaded = true;
		}
		return mus_phiErr_;
	}
	vector<float> &mus_ptErr()
	{
		if (not mus_ptErr_isLoaded) {
			if (mus_ptErr_branch != 0) {
				mus_ptErr_branch->GetEntry(index);
			} else { }
			mus_ptErr_isLoaded = true;
		}
		return mus_ptErr_;
	}
	vector<float> &mus_qoverp()
	{
		if (not mus_qoverp_isLoaded) {
			if (mus_qoverp_branch != 0) {
				mus_qoverp_branch->GetEntry(index);
			} else { }
			mus_qoverp_isLoaded = true;
		}
		return mus_qoverp_;
	}
	vector<float> &mus_qoverpError()
	{
		if (not mus_qoverpError_isLoaded) {
			if (mus_qoverpError_branch != 0) {
				mus_qoverpError_branch->GetEntry(index);
			} else { }
			mus_qoverpError_isLoaded = true;
		}
		return mus_qoverpError_;
	}
	vector<float> &mus_segmCompatibility()
	{
		if (not mus_segmCompatibility_isLoaded) {
			if (mus_segmCompatibility_branch != 0) {
				mus_segmCompatibility_branch->GetEntry(index);
			} else { }
			mus_segmCompatibility_isLoaded = true;
		}
		return mus_segmCompatibility_;
	}
	vector<float> &mus_staRelChi2()
	{
		if (not mus_staRelChi2_isLoaded) {
			if (mus_staRelChi2_branch != 0) {
				mus_staRelChi2_branch->GetEntry(index);
			} else { }
			mus_staRelChi2_isLoaded = true;
		}
		return mus_staRelChi2_;
	}
	vector<float> &mus_sta_chi2()
	{
		if (not mus_sta_chi2_isLoaded) {
			if (mus_sta_chi2_branch != 0) {
				mus_sta_chi2_branch->GetEntry(index);
			} else { }
			mus_sta_chi2_isLoaded = true;
		}
		return mus_sta_chi2_;
	}
	vector<float> &mus_sta_d0()
	{
		if (not mus_sta_d0_isLoaded) {
			if (mus_sta_d0_branch != 0) {
				mus_sta_d0_branch->GetEntry(index);
			} else { }
			mus_sta_d0_isLoaded = true;
		}
		return mus_sta_d0_;
	}
	vector<float> &mus_sta_d0Err()
	{
		if (not mus_sta_d0Err_isLoaded) {
			if (mus_sta_d0Err_branch != 0) {
				mus_sta_d0Err_branch->GetEntry(index);
			} else { }
			mus_sta_d0Err_isLoaded = true;
		}
		return mus_sta_d0Err_;
	}
	vector<float> &mus_sta_d0corr()
	{
		if (not mus_sta_d0corr_isLoaded) {
			if (mus_sta_d0corr_branch != 0) {
				mus_sta_d0corr_branch->GetEntry(index);
			} else { }
			mus_sta_d0corr_isLoaded = true;
		}
		return mus_sta_d0corr_;
	}
	vector<float> &mus_sta_d0corrPhi()
	{
		if (not mus_sta_d0corrPhi_isLoaded) {
			if (mus_sta_d0corrPhi_branch != 0) {
				mus_sta_d0corrPhi_branch->GetEntry(index);
			} else { }
			mus_sta_d0corrPhi_isLoaded = true;
		}
		return mus_sta_d0corrPhi_;
	}
	vector<float> &mus_sta_d0phiCov()
	{
		if (not mus_sta_d0phiCov_isLoaded) {
			if (mus_sta_d0phiCov_branch != 0) {
				mus_sta_d0phiCov_branch->GetEntry(index);
			} else { }
			mus_sta_d0phiCov_isLoaded = true;
		}
		return mus_sta_d0phiCov_;
	}
	vector<float> &mus_sta_ndof()
	{
		if (not mus_sta_ndof_isLoaded) {
			if (mus_sta_ndof_branch != 0) {
				mus_sta_ndof_branch->GetEntry(index);
			} else { }
			mus_sta_ndof_isLoaded = true;
		}
		return mus_sta_ndof_;
	}
	vector<float> &mus_sta_qoverp()
	{
		if (not mus_sta_qoverp_isLoaded) {
			if (mus_sta_qoverp_branch != 0) {
				mus_sta_qoverp_branch->GetEntry(index);
			} else { }
			mus_sta_qoverp_isLoaded = true;
		}
		return mus_sta_qoverp_;
	}
	vector<float> &mus_sta_qoverpError()
	{
		if (not mus_sta_qoverpError_isLoaded) {
			if (mus_sta_qoverpError_branch != 0) {
				mus_sta_qoverpError_branch->GetEntry(index);
			} else { }
			mus_sta_qoverpError_isLoaded = true;
		}
		return mus_sta_qoverpError_;
	}
	vector<float> &mus_sta_z0()
	{
		if (not mus_sta_z0_isLoaded) {
			if (mus_sta_z0_branch != 0) {
				mus_sta_z0_branch->GetEntry(index);
			} else { }
			mus_sta_z0_isLoaded = true;
		}
		return mus_sta_z0_;
	}
	vector<float> &mus_sta_z0Err()
	{
		if (not mus_sta_z0Err_isLoaded) {
			if (mus_sta_z0Err_branch != 0) {
				mus_sta_z0Err_branch->GetEntry(index);
			} else { }
			mus_sta_z0Err_isLoaded = true;
		}
		return mus_sta_z0Err_;
	}
	vector<float> &mus_sta_z0corr()
	{
		if (not mus_sta_z0corr_isLoaded) {
			if (mus_sta_z0corr_branch != 0) {
				mus_sta_z0corr_branch->GetEntry(index);
			} else { }
			mus_sta_z0corr_isLoaded = true;
		}
		return mus_sta_z0corr_;
	}
	vector<float> &mus_timeAtIpInOut()
	{
		if (not mus_timeAtIpInOut_isLoaded) {
			if (mus_timeAtIpInOut_branch != 0) {
				mus_timeAtIpInOut_branch->GetEntry(index);
			} else { }
			mus_timeAtIpInOut_isLoaded = true;
		}
		return mus_timeAtIpInOut_;
	}
	vector<float> &mus_timeAtIpInOutErr()
	{
		if (not mus_timeAtIpInOutErr_isLoaded) {
			if (mus_timeAtIpInOutErr_branch != 0) {
				mus_timeAtIpInOutErr_branch->GetEntry(index);
			} else { }
			mus_timeAtIpInOutErr_isLoaded = true;
		}
		return mus_timeAtIpInOutErr_;
	}
	vector<float> &mus_timeAtIpOutIn()
	{
		if (not mus_timeAtIpOutIn_isLoaded) {
			if (mus_timeAtIpOutIn_branch != 0) {
				mus_timeAtIpOutIn_branch->GetEntry(index);
			} else { }
			mus_timeAtIpOutIn_isLoaded = true;
		}
		return mus_timeAtIpOutIn_;
	}
	vector<float> &mus_timeAtIpOutInErr()
	{
		if (not mus_timeAtIpOutInErr_isLoaded) {
			if (mus_timeAtIpOutInErr_branch != 0) {
				mus_timeAtIpOutInErr_branch->GetEntry(index);
			} else { }
			mus_timeAtIpOutInErr_isLoaded = true;
		}
		return mus_timeAtIpOutInErr_;
	}
	vector<float> &mus_timeCompat()
	{
		if (not mus_timeCompat_isLoaded) {
			if (mus_timeCompat_branch != 0) {
				mus_timeCompat_branch->GetEntry(index);
			} else { }
			mus_timeCompat_isLoaded = true;
		}
		return mus_timeCompat_;
	}
	vector<float> &mus_tower()
	{
		if (not mus_tower_isLoaded) {
			if (mus_tower_branch != 0) {
				mus_tower_branch->GetEntry(index);
			} else { }
			mus_tower_isLoaded = true;
		}
		return mus_tower_;
	}
	vector<float> &mus_towerS9()
	{
		if (not mus_towerS9_isLoaded) {
			if (mus_towerS9_branch != 0) {
				mus_towerS9_branch->GetEntry(index);
			} else { }
			mus_towerS9_isLoaded = true;
		}
		return mus_towerS9_;
	}
	vector<float> &mus_trkKink()
	{
		if (not mus_trkKink_isLoaded) {
			if (mus_trkKink_branch != 0) {
				mus_trkKink_branch->GetEntry(index);
			} else { }
			mus_trkKink_isLoaded = true;
		}
		return mus_trkKink_;
	}
	vector<float> &mus_trkRelChi2()
	{
		if (not mus_trkRelChi2_isLoaded) {
			if (mus_trkRelChi2_branch != 0) {
				mus_trkRelChi2_branch->GetEntry(index);
			} else { }
			mus_trkRelChi2_isLoaded = true;
		}
		return mus_trkRelChi2_;
	}
	vector<float> &mus_vertexCompat()
	{
		if (not mus_vertexCompat_isLoaded) {
			if (mus_vertexCompat_branch != 0) {
				mus_vertexCompat_branch->GetEntry(index);
			} else { }
			mus_vertexCompat_isLoaded = true;
		}
		return mus_vertexCompat_;
	}
	vector<float> &mus_vertexphi()
	{
		if (not mus_vertexphi_isLoaded) {
			if (mus_vertexphi_branch != 0) {
				mus_vertexphi_branch->GetEntry(index);
			} else { }
			mus_vertexphi_isLoaded = true;
		}
		return mus_vertexphi_;
	}
	vector<float> &mus_z0()
	{
		if (not mus_z0_isLoaded) {
			if (mus_z0_branch != 0) {
				mus_z0_branch->GetEntry(index);
			} else { }
			mus_z0_isLoaded = true;
		}
		return mus_z0_;
	}
	vector<float> &mus_z0Err()
	{
		if (not mus_z0Err_isLoaded) {
			if (mus_z0Err_branch != 0) {
				mus_z0Err_branch->GetEntry(index);
			} else { }
			mus_z0Err_isLoaded = true;
		}
		return mus_z0Err_;
	}
	vector<float> &mus_z0corr()
	{
		if (not mus_z0corr_isLoaded) {
			if (mus_z0corr_branch != 0) {
				mus_z0corr_branch->GetEntry(index);
			} else { }
			mus_z0corr_isLoaded = true;
		}
		return mus_z0corr_;
	}
	vector<float> &pfcands_dz()
	{
		if (not pfcands_dz_isLoaded) {
			if (pfcands_dz_branch != 0) {
				pfcands_dz_branch->GetEntry(index);
			} else { }
			pfcands_dz_isLoaded = true;
		}
		return pfcands_dz_;
	}
	vector<float> &pfjets_area()
	{
		if (not pfjets_area_isLoaded) {
			if (pfjets_area_branch != 0) {
				pfjets_area_branch->GetEntry(index);
			} else { }
			pfjets_area_isLoaded = true;
		}
		return pfjets_area_;
	}
	vector<float> &pfjets_chargedEmE()
	{
		if (not pfjets_chargedEmE_isLoaded) {
			if (pfjets_chargedEmE_branch != 0) {
				pfjets_chargedEmE_branch->GetEntry(index);
			} else { }
			pfjets_chargedEmE_isLoaded = true;
		}
		return pfjets_chargedEmE_;
	}
	vector<float> &pfjets_chargedHadronE()
	{
		if (not pfjets_chargedHadronE_isLoaded) {
			if (pfjets_chargedHadronE_branch != 0) {
				pfjets_chargedHadronE_branch->GetEntry(index);
			} else { }
			pfjets_chargedHadronE_isLoaded = true;
		}
		return pfjets_chargedHadronE_;
	}
	vector<float> &pfjets_combinedSecondaryVertexBJetTag()
	{
		if (not pfjets_combinedSecondaryVertexBJetTag_isLoaded) {
			if (pfjets_combinedSecondaryVertexBJetTag_branch != 0) {
				pfjets_combinedSecondaryVertexBJetTag_branch->GetEntry(index);
			} else { }
			pfjets_combinedSecondaryVertexBJetTag_isLoaded = true;
		}
		return pfjets_combinedSecondaryVertexBJetTag_;
	}
	vector<float> &pfjets_cor()
	{
		if (not pfjets_cor_isLoaded) {
			if (pfjets_cor_branch != 0) {
				pfjets_cor_branch->GetEntry(index);
			} else { }
			pfjets_cor_isLoaded = true;
		}
		return pfjets_cor_;
	}
	vector<float> &pfjets_corL1Fast()
	{
		if (not pfjets_corL1Fast_isLoaded) {
			if (pfjets_corL1Fast_branch != 0) {
				pfjets_corL1Fast_branch->GetEntry(index);
			} else { }
			pfjets_corL1Fast_isLoaded = true;
		}
		return pfjets_corL1Fast_;
	}
	vector<float> &pfjets_corL1FastL2L3()
	{
		if (not pfjets_corL1FastL2L3_isLoaded) {
			if (pfjets_corL1FastL2L3_branch != 0) {
				pfjets_corL1FastL2L3_branch->GetEntry(index);
			} else { }
			pfjets_corL1FastL2L3_isLoaded = true;
		}
		return pfjets_corL1FastL2L3_;
	}
	vector<float> &pfjets_corL1FastL2L3residual()
	{
		if (not pfjets_corL1FastL2L3residual_isLoaded) {
			if (pfjets_corL1FastL2L3residual_branch != 0) {
				pfjets_corL1FastL2L3residual_branch->GetEntry(index);
			} else { }
			pfjets_corL1FastL2L3residual_isLoaded = true;
		}
		return pfjets_corL1FastL2L3residual_;
	}
	vector<float> &pfjets_corL1L2L3()
	{
		if (not pfjets_corL1L2L3_isLoaded) {
			if (pfjets_corL1L2L3_branch != 0) {
				pfjets_corL1L2L3_branch->GetEntry(index);
			} else { }
			pfjets_corL1L2L3_isLoaded = true;
		}
		return pfjets_corL1L2L3_;
	}
	vector<float> &pfjets_electronE()
	{
		if (not pfjets_electronE_isLoaded) {
			if (pfjets_electronE_branch != 0) {
				pfjets_electronE_branch->GetEntry(index);
			} else { }
			pfjets_electronE_isLoaded = true;
		}
		return pfjets_electronE_;
	}
	vector<float> &pfjets_hfEmE()
	{
		if (not pfjets_hfEmE_isLoaded) {
			if (pfjets_hfEmE_branch != 0) {
				pfjets_hfEmE_branch->GetEntry(index);
			} else { }
			pfjets_hfEmE_isLoaded = true;
		}
		return pfjets_hfEmE_;
	}
	vector<float> &pfjets_hfHadronE()
	{
		if (not pfjets_hfHadronE_isLoaded) {
			if (pfjets_hfHadronE_branch != 0) {
				pfjets_hfHadronE_branch->GetEntry(index);
			} else { }
			pfjets_hfHadronE_isLoaded = true;
		}
		return pfjets_hfHadronE_;
	}
	vector<float> &pfjets_jetBProbabilityBJetTag()
	{
		if (not pfjets_jetBProbabilityBJetTag_isLoaded) {
			if (pfjets_jetBProbabilityBJetTag_branch != 0) {
				pfjets_jetBProbabilityBJetTag_branch->GetEntry(index);
			} else { }
			pfjets_jetBProbabilityBJetTag_isLoaded = true;
		}
		return pfjets_jetBProbabilityBJetTag_;
	}
	vector<float> &pfjets_jetProbabilityBJetTag()
	{
		if (not pfjets_jetProbabilityBJetTag_isLoaded) {
			if (pfjets_jetProbabilityBJetTag_branch != 0) {
				pfjets_jetProbabilityBJetTag_branch->GetEntry(index);
			} else { }
			pfjets_jetProbabilityBJetTag_isLoaded = true;
		}
		return pfjets_jetProbabilityBJetTag_;
	}
	vector<float> &pfjets_muonE()
	{
		if (not pfjets_muonE_isLoaded) {
			if (pfjets_muonE_branch != 0) {
				pfjets_muonE_branch->GetEntry(index);
			} else { }
			pfjets_muonE_isLoaded = true;
		}
		return pfjets_muonE_;
	}
	vector<float> &pfjets_neutralEmE()
	{
		if (not pfjets_neutralEmE_isLoaded) {
			if (pfjets_neutralEmE_branch != 0) {
				pfjets_neutralEmE_branch->GetEntry(index);
			} else { }
			pfjets_neutralEmE_isLoaded = true;
		}
		return pfjets_neutralEmE_;
	}
	vector<float> &pfjets_neutralHadronE()
	{
		if (not pfjets_neutralHadronE_isLoaded) {
			if (pfjets_neutralHadronE_branch != 0) {
				pfjets_neutralHadronE_branch->GetEntry(index);
			} else { }
			pfjets_neutralHadronE_isLoaded = true;
		}
		return pfjets_neutralHadronE_;
	}
	vector<float> &pfjets_photonE()
	{
		if (not pfjets_photonE_isLoaded) {
			if (pfjets_photonE_branch != 0) {
				pfjets_photonE_branch->GetEntry(index);
			} else { }
			pfjets_photonE_isLoaded = true;
		}
		return pfjets_photonE_;
	}
	vector<float> &pfjets_pileupJetId()
	{
		if (not pfjets_pileupJetId_isLoaded) {
			if (pfjets_pileupJetId_branch != 0) {
				pfjets_pileupJetId_branch->GetEntry(index);
			} else { }
			pfjets_pileupJetId_isLoaded = true;
		}
		return pfjets_pileupJetId_;
	}
	vector<float> &pfjets_simpleSecondaryVertexHighEffBJetTag()
	{
		if (not pfjets_simpleSecondaryVertexHighEffBJetTag_isLoaded) {
			if (pfjets_simpleSecondaryVertexHighEffBJetTag_branch != 0) {
				pfjets_simpleSecondaryVertexHighEffBJetTag_branch->GetEntry(index);
			} else { }
			pfjets_simpleSecondaryVertexHighEffBJetTag_isLoaded = true;
		}
		return pfjets_simpleSecondaryVertexHighEffBJetTag_;
	}
	vector<float> &pfjets_simpleSecondaryVertexHighPurBJetTags()
	{
		if (not pfjets_simpleSecondaryVertexHighPurBJetTags_isLoaded) {
			if (pfjets_simpleSecondaryVertexHighPurBJetTags_branch != 0) {
				pfjets_simpleSecondaryVertexHighPurBJetTags_branch->GetEntry(index);
			} else { }
			pfjets_simpleSecondaryVertexHighPurBJetTags_isLoaded = true;
		}
		return pfjets_simpleSecondaryVertexHighPurBJetTags_;
	}
	vector<float> &pfjets_trackCountingHighEffBJetTag()
	{
		if (not pfjets_trackCountingHighEffBJetTag_isLoaded) {
			if (pfjets_trackCountingHighEffBJetTag_branch != 0) {
				pfjets_trackCountingHighEffBJetTag_branch->GetEntry(index);
			} else { }
			pfjets_trackCountingHighEffBJetTag_isLoaded = true;
		}
		return pfjets_trackCountingHighEffBJetTag_;
	}
	vector<float> &pfjets_trackCountingHighPurBJetTag()
	{
		if (not pfjets_trackCountingHighPurBJetTag_isLoaded) {
			if (pfjets_trackCountingHighPurBJetTag_branch != 0) {
				pfjets_trackCountingHighPurBJetTag_branch->GetEntry(index);
			} else { }
			pfjets_trackCountingHighPurBJetTag_isLoaded = true;
		}
		return pfjets_trackCountingHighPurBJetTag_;
	}
	vector<float> &taus_pf_againstElectronDeadECAL()
	{
		if (not taus_pf_againstElectronDeadECAL_isLoaded) {
			if (taus_pf_againstElectronDeadECAL_branch != 0) {
				taus_pf_againstElectronDeadECAL_branch->GetEntry(index);
			} else { }
			taus_pf_againstElectronDeadECAL_isLoaded = true;
		}
		return taus_pf_againstElectronDeadECAL_;
	}
	vector<float> &taus_pf_againstElectronLoose()
	{
		if (not taus_pf_againstElectronLoose_isLoaded) {
			if (taus_pf_againstElectronLoose_branch != 0) {
				taus_pf_againstElectronLoose_branch->GetEntry(index);
			} else { }
			taus_pf_againstElectronLoose_isLoaded = true;
		}
		return taus_pf_againstElectronLoose_;
	}
	vector<float> &taus_pf_againstElectronMedium()
	{
		if (not taus_pf_againstElectronMedium_isLoaded) {
			if (taus_pf_againstElectronMedium_branch != 0) {
				taus_pf_againstElectronMedium_branch->GetEntry(index);
			} else { }
			taus_pf_againstElectronMedium_isLoaded = true;
		}
		return taus_pf_againstElectronMedium_;
	}
	vector<float> &taus_pf_againstElectronTight()
	{
		if (not taus_pf_againstElectronTight_isLoaded) {
			if (taus_pf_againstElectronTight_branch != 0) {
				taus_pf_againstElectronTight_branch->GetEntry(index);
			} else { }
			taus_pf_againstElectronTight_isLoaded = true;
		}
		return taus_pf_againstElectronTight_;
	}
	vector<float> &taus_pf_againstMuonLoose()
	{
		if (not taus_pf_againstMuonLoose_isLoaded) {
			if (taus_pf_againstMuonLoose_branch != 0) {
				taus_pf_againstMuonLoose_branch->GetEntry(index);
			} else { }
			taus_pf_againstMuonLoose_isLoaded = true;
		}
		return taus_pf_againstMuonLoose_;
	}
	vector<float> &taus_pf_againstMuonLoose2()
	{
		if (not taus_pf_againstMuonLoose2_isLoaded) {
			if (taus_pf_againstMuonLoose2_branch != 0) {
				taus_pf_againstMuonLoose2_branch->GetEntry(index);
			} else { }
			taus_pf_againstMuonLoose2_isLoaded = true;
		}
		return taus_pf_againstMuonLoose2_;
	}
	vector<float> &taus_pf_againstMuonLoose3()
	{
		if (not taus_pf_againstMuonLoose3_isLoaded) {
			if (taus_pf_againstMuonLoose3_branch != 0) {
				taus_pf_againstMuonLoose3_branch->GetEntry(index);
			} else { }
			taus_pf_againstMuonLoose3_isLoaded = true;
		}
		return taus_pf_againstMuonLoose3_;
	}
	vector<float> &taus_pf_againstMuonMedium()
	{
		if (not taus_pf_againstMuonMedium_isLoaded) {
			if (taus_pf_againstMuonMedium_branch != 0) {
				taus_pf_againstMuonMedium_branch->GetEntry(index);
			} else { }
			taus_pf_againstMuonMedium_isLoaded = true;
		}
		return taus_pf_againstMuonMedium_;
	}
	vector<float> &taus_pf_againstMuonMedium2()
	{
		if (not taus_pf_againstMuonMedium2_isLoaded) {
			if (taus_pf_againstMuonMedium2_branch != 0) {
				taus_pf_againstMuonMedium2_branch->GetEntry(index);
			} else { }
			taus_pf_againstMuonMedium2_isLoaded = true;
		}
		return taus_pf_againstMuonMedium2_;
	}
	vector<float> &taus_pf_againstMuonTight()
	{
		if (not taus_pf_againstMuonTight_isLoaded) {
			if (taus_pf_againstMuonTight_branch != 0) {
				taus_pf_againstMuonTight_branch->GetEntry(index);
			} else { }
			taus_pf_againstMuonTight_isLoaded = true;
		}
		return taus_pf_againstMuonTight_;
	}
	vector<float> &taus_pf_againstMuonTight2()
	{
		if (not taus_pf_againstMuonTight2_isLoaded) {
			if (taus_pf_againstMuonTight2_branch != 0) {
				taus_pf_againstMuonTight2_branch->GetEntry(index);
			} else { }
			taus_pf_againstMuonTight2_isLoaded = true;
		}
		return taus_pf_againstMuonTight2_;
	}
	vector<float> &taus_pf_againstMuonTight3()
	{
		if (not taus_pf_againstMuonTight3_isLoaded) {
			if (taus_pf_againstMuonTight3_branch != 0) {
				taus_pf_againstMuonTight3_branch->GetEntry(index);
			} else { }
			taus_pf_againstMuonTight3_isLoaded = true;
		}
		return taus_pf_againstMuonTight3_;
	}
	vector<float> &taus_pf_byCombinedIsolationDeltaBetaCorrRaw()
	{
		if (not taus_pf_byCombinedIsolationDeltaBetaCorrRaw_isLoaded) {
			if (taus_pf_byCombinedIsolationDeltaBetaCorrRaw_branch != 0) {
				taus_pf_byCombinedIsolationDeltaBetaCorrRaw_branch->GetEntry(index);
			} else { }
			taus_pf_byCombinedIsolationDeltaBetaCorrRaw_isLoaded = true;
		}
		return taus_pf_byCombinedIsolationDeltaBetaCorrRaw_;
	}
	vector<float> &taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits()
	{
		if (not taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_isLoaded) {
			if (taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_branch != 0) {
				taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_branch->GetEntry(index);
			} else { }
			taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_isLoaded = true;
		}
		return taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits_;
	}
	vector<float> &taus_pf_byDecayModeFinding()
	{
		if (not taus_pf_byDecayModeFinding_isLoaded) {
			if (taus_pf_byDecayModeFinding_branch != 0) {
				taus_pf_byDecayModeFinding_branch->GetEntry(index);
			} else { }
			taus_pf_byDecayModeFinding_isLoaded = true;
		}
		return taus_pf_byDecayModeFinding_;
	}
	vector<float> &taus_pf_byLooseCombinedIsolationDeltaBetaCorr()
	{
		if (not taus_pf_byLooseCombinedIsolationDeltaBetaCorr_isLoaded) {
			if (taus_pf_byLooseCombinedIsolationDeltaBetaCorr_branch != 0) {
				taus_pf_byLooseCombinedIsolationDeltaBetaCorr_branch->GetEntry(index);
			} else { }
			taus_pf_byLooseCombinedIsolationDeltaBetaCorr_isLoaded = true;
		}
		return taus_pf_byLooseCombinedIsolationDeltaBetaCorr_;
	}
	vector<float> &taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits()
	{
		if (not taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_isLoaded) {
			if (taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_branch != 0) {
				taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_branch->GetEntry(index);
			} else { }
			taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_isLoaded = true;
		}
		return taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits_;
	}
	vector<float> &taus_pf_byMediumCombinedIsolationDeltaBetaCorr()
	{
		if (not taus_pf_byMediumCombinedIsolationDeltaBetaCorr_isLoaded) {
			if (taus_pf_byMediumCombinedIsolationDeltaBetaCorr_branch != 0) {
				taus_pf_byMediumCombinedIsolationDeltaBetaCorr_branch->GetEntry(index);
			} else { }
			taus_pf_byMediumCombinedIsolationDeltaBetaCorr_isLoaded = true;
		}
		return taus_pf_byMediumCombinedIsolationDeltaBetaCorr_;
	}
	vector<float> &taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits()
	{
		if (not taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_isLoaded) {
			if (taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_branch != 0) {
				taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_branch->GetEntry(index);
			} else { }
			taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_isLoaded = true;
		}
		return taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits_;
	}
	vector<float> &taus_pf_byTightCombinedIsolationDeltaBetaCorr()
	{
		if (not taus_pf_byTightCombinedIsolationDeltaBetaCorr_isLoaded) {
			if (taus_pf_byTightCombinedIsolationDeltaBetaCorr_branch != 0) {
				taus_pf_byTightCombinedIsolationDeltaBetaCorr_branch->GetEntry(index);
			} else { }
			taus_pf_byTightCombinedIsolationDeltaBetaCorr_isLoaded = true;
		}
		return taus_pf_byTightCombinedIsolationDeltaBetaCorr_;
	}
	vector<float> &taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits()
	{
		if (not taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_isLoaded) {
			if (taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_branch != 0) {
				taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_branch->GetEntry(index);
			} else { }
			taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_isLoaded = true;
		}
		return taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits_;
	}
	vector<float> &taus_pf_byVLooseCombinedIsolationDeltaBetaCorr()
	{
		if (not taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_isLoaded) {
			if (taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_branch != 0) {
				taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_branch->GetEntry(index);
			} else { }
			taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_isLoaded = true;
		}
		return taus_pf_byVLooseCombinedIsolationDeltaBetaCorr_;
	}
	vector<float> &photons_e1x5()
	{
		if (not photons_e1x5_isLoaded) {
			if (photons_e1x5_branch != 0) {
				photons_e1x5_branch->GetEntry(index);
			} else { }
			photons_e1x5_isLoaded = true;
		}
		return photons_e1x5_;
	}
	vector<float> &photons_e2x5Max()
	{
		if (not photons_e2x5Max_isLoaded) {
			if (photons_e2x5Max_branch != 0) {
				photons_e2x5Max_branch->GetEntry(index);
			} else { }
			photons_e2x5Max_isLoaded = true;
		}
		return photons_e2x5Max_;
	}
	vector<float> &photons_e3x3()
	{
		if (not photons_e3x3_isLoaded) {
			if (photons_e3x3_branch != 0) {
				photons_e3x3_branch->GetEntry(index);
			} else { }
			photons_e3x3_isLoaded = true;
		}
		return photons_e3x3_;
	}
	vector<float> &photons_e5x5()
	{
		if (not photons_e5x5_isLoaded) {
			if (photons_e5x5_branch != 0) {
				photons_e5x5_branch->GetEntry(index);
			} else { }
			photons_e5x5_isLoaded = true;
		}
		return photons_e5x5_;
	}
	vector<float> &photons_ecalIso03()
	{
		if (not photons_ecalIso03_isLoaded) {
			if (photons_ecalIso03_branch != 0) {
				photons_ecalIso03_branch->GetEntry(index);
			} else { }
			photons_ecalIso03_isLoaded = true;
		}
		return photons_ecalIso03_;
	}
	vector<float> &photons_ecalIso04()
	{
		if (not photons_ecalIso04_isLoaded) {
			if (photons_ecalIso04_branch != 0) {
				photons_ecalIso04_branch->GetEntry(index);
			} else { }
			photons_ecalIso04_isLoaded = true;
		}
		return photons_ecalIso04_;
	}
	vector<float> &photons_full3x3_e3x3()
	{
		if (not photons_full3x3_e3x3_isLoaded) {
			if (photons_full3x3_e3x3_branch != 0) {
				photons_full3x3_e3x3_branch->GetEntry(index);
			} else { }
			photons_full3x3_e3x3_isLoaded = true;
		}
		return photons_full3x3_e3x3_;
	}
	vector<float> &photons_full5x5_e1x5()
	{
		if (not photons_full5x5_e1x5_isLoaded) {
			if (photons_full5x5_e1x5_branch != 0) {
				photons_full5x5_e1x5_branch->GetEntry(index);
			} else { }
			photons_full5x5_e1x5_isLoaded = true;
		}
		return photons_full5x5_e1x5_;
	}
	vector<float> &photons_full5x5_e2x5Max()
	{
		if (not photons_full5x5_e2x5Max_isLoaded) {
			if (photons_full5x5_e2x5Max_branch != 0) {
				photons_full5x5_e2x5Max_branch->GetEntry(index);
			} else { }
			photons_full5x5_e2x5Max_isLoaded = true;
		}
		return photons_full5x5_e2x5Max_;
	}
	vector<float> &photons_full5x5_e5x5()
	{
		if (not photons_full5x5_e5x5_isLoaded) {
			if (photons_full5x5_e5x5_branch != 0) {
				photons_full5x5_e5x5_branch->GetEntry(index);
			} else { }
			photons_full5x5_e5x5_isLoaded = true;
		}
		return photons_full5x5_e5x5_;
	}
	vector<float> &photons_full5x5_hOverE()
	{
		if (not photons_full5x5_hOverE_isLoaded) {
			if (photons_full5x5_hOverE_branch != 0) {
				photons_full5x5_hOverE_branch->GetEntry(index);
			} else { }
			photons_full5x5_hOverE_isLoaded = true;
		}
		return photons_full5x5_hOverE_;
	}
	vector<float> &photons_full5x5_sigmaEtaEta()
	{
		if (not photons_full5x5_sigmaEtaEta_isLoaded) {
			if (photons_full5x5_sigmaEtaEta_branch != 0) {
				photons_full5x5_sigmaEtaEta_branch->GetEntry(index);
			} else { }
			photons_full5x5_sigmaEtaEta_isLoaded = true;
		}
		return photons_full5x5_sigmaEtaEta_;
	}
	vector<float> &photons_full5x5_sigmaIEtaIEta()
	{
		if (not photons_full5x5_sigmaIEtaIEta_isLoaded) {
			if (photons_full5x5_sigmaIEtaIEta_branch != 0) {
				photons_full5x5_sigmaIEtaIEta_branch->GetEntry(index);
			} else { }
			photons_full5x5_sigmaIEtaIEta_isLoaded = true;
		}
		return photons_full5x5_sigmaIEtaIEta_;
	}
	vector<float> &photons_hOverE()
	{
		if (not photons_hOverE_isLoaded) {
			if (photons_hOverE_branch != 0) {
				photons_hOverE_branch->GetEntry(index);
			} else { }
			photons_hOverE_isLoaded = true;
		}
		return photons_hOverE_;
	}
	vector<float> &photons_hcalDepth1TowerSumEtBcConeDR03()
	{
		if (not photons_hcalDepth1TowerSumEtBcConeDR03_isLoaded) {
			if (photons_hcalDepth1TowerSumEtBcConeDR03_branch != 0) {
				photons_hcalDepth1TowerSumEtBcConeDR03_branch->GetEntry(index);
			} else { }
			photons_hcalDepth1TowerSumEtBcConeDR03_isLoaded = true;
		}
		return photons_hcalDepth1TowerSumEtBcConeDR03_;
	}
	vector<float> &photons_hcalDepth1TowerSumEtBcConeDR04()
	{
		if (not photons_hcalDepth1TowerSumEtBcConeDR04_isLoaded) {
			if (photons_hcalDepth1TowerSumEtBcConeDR04_branch != 0) {
				photons_hcalDepth1TowerSumEtBcConeDR04_branch->GetEntry(index);
			} else { }
			photons_hcalDepth1TowerSumEtBcConeDR04_isLoaded = true;
		}
		return photons_hcalDepth1TowerSumEtBcConeDR04_;
	}
	vector<float> &photons_hcalDepth2TowerSumEtBcConeDR03()
	{
		if (not photons_hcalDepth2TowerSumEtBcConeDR03_isLoaded) {
			if (photons_hcalDepth2TowerSumEtBcConeDR03_branch != 0) {
				photons_hcalDepth2TowerSumEtBcConeDR03_branch->GetEntry(index);
			} else { }
			photons_hcalDepth2TowerSumEtBcConeDR03_isLoaded = true;
		}
		return photons_hcalDepth2TowerSumEtBcConeDR03_;
	}
	vector<float> &photons_hcalDepth2TowerSumEtBcConeDR04()
	{
		if (not photons_hcalDepth2TowerSumEtBcConeDR04_isLoaded) {
			if (photons_hcalDepth2TowerSumEtBcConeDR04_branch != 0) {
				photons_hcalDepth2TowerSumEtBcConeDR04_branch->GetEntry(index);
			} else { }
			photons_hcalDepth2TowerSumEtBcConeDR04_isLoaded = true;
		}
		return photons_hcalDepth2TowerSumEtBcConeDR04_;
	}
	vector<float> &photons_hcalIso03()
	{
		if (not photons_hcalIso03_isLoaded) {
			if (photons_hcalIso03_branch != 0) {
				photons_hcalIso03_branch->GetEntry(index);
			} else { }
			photons_hcalIso03_isLoaded = true;
		}
		return photons_hcalIso03_;
	}
	vector<float> &photons_hcalIso04()
	{
		if (not photons_hcalIso04_isLoaded) {
			if (photons_hcalIso04_branch != 0) {
				photons_hcalIso04_branch->GetEntry(index);
			} else { }
			photons_hcalIso04_isLoaded = true;
		}
		return photons_hcalIso04_;
	}
	vector<float> &photons_hcalTowerSumEtBcConeDR03()
	{
		if (not photons_hcalTowerSumEtBcConeDR03_isLoaded) {
			if (photons_hcalTowerSumEtBcConeDR03_branch != 0) {
				photons_hcalTowerSumEtBcConeDR03_branch->GetEntry(index);
			} else { }
			photons_hcalTowerSumEtBcConeDR03_isLoaded = true;
		}
		return photons_hcalTowerSumEtBcConeDR03_;
	}
	vector<float> &photons_hcalTowerSumEtBcConeDR04()
	{
		if (not photons_hcalTowerSumEtBcConeDR04_isLoaded) {
			if (photons_hcalTowerSumEtBcConeDR04_branch != 0) {
				photons_hcalTowerSumEtBcConeDR04_branch->GetEntry(index);
			} else { }
			photons_hcalTowerSumEtBcConeDR04_isLoaded = true;
		}
		return photons_hcalTowerSumEtBcConeDR04_;
	}
	vector<float> &photons_ntkIsoHollow03()
	{
		if (not photons_ntkIsoHollow03_isLoaded) {
			if (photons_ntkIsoHollow03_branch != 0) {
				photons_ntkIsoHollow03_branch->GetEntry(index);
			} else { }
			photons_ntkIsoHollow03_isLoaded = true;
		}
		return photons_ntkIsoHollow03_;
	}
	vector<float> &photons_ntkIsoHollow04()
	{
		if (not photons_ntkIsoHollow04_isLoaded) {
			if (photons_ntkIsoHollow04_branch != 0) {
				photons_ntkIsoHollow04_branch->GetEntry(index);
			} else { }
			photons_ntkIsoHollow04_isLoaded = true;
		}
		return photons_ntkIsoHollow04_;
	}
	vector<float> &photons_ntkIsoSolid03()
	{
		if (not photons_ntkIsoSolid03_isLoaded) {
			if (photons_ntkIsoSolid03_branch != 0) {
				photons_ntkIsoSolid03_branch->GetEntry(index);
			} else { }
			photons_ntkIsoSolid03_isLoaded = true;
		}
		return photons_ntkIsoSolid03_;
	}
	vector<float> &photons_ntkIsoSolid04()
	{
		if (not photons_ntkIsoSolid04_isLoaded) {
			if (photons_ntkIsoSolid04_branch != 0) {
				photons_ntkIsoSolid04_branch->GetEntry(index);
			} else { }
			photons_ntkIsoSolid04_isLoaded = true;
		}
		return photons_ntkIsoSolid04_;
	}
	vector<float> &photons_sigmaEtaEta()
	{
		if (not photons_sigmaEtaEta_isLoaded) {
			if (photons_sigmaEtaEta_branch != 0) {
				photons_sigmaEtaEta_branch->GetEntry(index);
			} else { }
			photons_sigmaEtaEta_isLoaded = true;
		}
		return photons_sigmaEtaEta_;
	}
	vector<float> &photons_sigmaIEtaIEta()
	{
		if (not photons_sigmaIEtaIEta_isLoaded) {
			if (photons_sigmaIEtaIEta_branch != 0) {
				photons_sigmaIEtaIEta_branch->GetEntry(index);
			} else { }
			photons_sigmaIEtaIEta_isLoaded = true;
		}
		return photons_sigmaIEtaIEta_;
	}
	vector<float> &photons_tkIsoHollow03()
	{
		if (not photons_tkIsoHollow03_isLoaded) {
			if (photons_tkIsoHollow03_branch != 0) {
				photons_tkIsoHollow03_branch->GetEntry(index);
			} else { }
			photons_tkIsoHollow03_isLoaded = true;
		}
		return photons_tkIsoHollow03_;
	}
	vector<float> &photons_tkIsoHollow04()
	{
		if (not photons_tkIsoHollow04_isLoaded) {
			if (photons_tkIsoHollow04_branch != 0) {
				photons_tkIsoHollow04_branch->GetEntry(index);
			} else { }
			photons_tkIsoHollow04_isLoaded = true;
		}
		return photons_tkIsoHollow04_;
	}
	vector<float> &photons_tkIsoSolid03()
	{
		if (not photons_tkIsoSolid03_isLoaded) {
			if (photons_tkIsoSolid03_branch != 0) {
				photons_tkIsoSolid03_branch->GetEntry(index);
			} else { }
			photons_tkIsoSolid03_isLoaded = true;
		}
		return photons_tkIsoSolid03_;
	}
	vector<float> &photons_tkIsoSolid04()
	{
		if (not photons_tkIsoSolid04_isLoaded) {
			if (photons_tkIsoSolid04_branch != 0) {
				photons_tkIsoSolid04_branch->GetEntry(index);
			} else { }
			photons_tkIsoSolid04_isLoaded = true;
		}
		return photons_tkIsoSolid04_;
	}
	vector<float> &puInfo_trueNumInteractions()
	{
		if (not puInfo_trueNumInteractions_isLoaded) {
			if (puInfo_trueNumInteractions_branch != 0) {
				puInfo_trueNumInteractions_branch->GetEntry(index);
			} else { }
			puInfo_trueNumInteractions_isLoaded = true;
		}
		return puInfo_trueNumInteractions_;
	}
	vector<float> &convs_chi2()
	{
		if (not convs_chi2_isLoaded) {
			if (convs_chi2_branch != 0) {
				convs_chi2_branch->GetEntry(index);
			} else { }
			convs_chi2_isLoaded = true;
		}
		return convs_chi2_;
	}
	vector<float> &convs_dl()
	{
		if (not convs_dl_isLoaded) {
			if (convs_dl_branch != 0) {
				convs_dl_branch->GetEntry(index);
			} else { }
			convs_dl_isLoaded = true;
		}
		return convs_dl_;
	}
	vector<float> &convs_ndof()
	{
		if (not convs_ndof_isLoaded) {
			if (convs_ndof_branch != 0) {
				convs_ndof_branch->GetEntry(index);
			} else { }
			convs_ndof_isLoaded = true;
		}
		return convs_ndof_;
	}
	vector<float> &vtxs_chi2()
	{
		if (not vtxs_chi2_isLoaded) {
			if (vtxs_chi2_branch != 0) {
				vtxs_chi2_branch->GetEntry(index);
			} else { }
			vtxs_chi2_isLoaded = true;
		}
		return vtxs_chi2_;
	}
	vector<float> &vtxs_ndof()
	{
		if (not vtxs_ndof_isLoaded) {
			if (vtxs_ndof_branch != 0) {
				vtxs_ndof_branch->GetEntry(index);
			} else { }
			vtxs_ndof_isLoaded = true;
		}
		return vtxs_ndof_;
	}
	vector<float> &vtxs_sumpt()
	{
		if (not vtxs_sumpt_isLoaded) {
			if (vtxs_sumpt_branch != 0) {
				vtxs_sumpt_branch->GetEntry(index);
			} else { }
			vtxs_sumpt_isLoaded = true;
		}
		return vtxs_sumpt_;
	}
	vector<float> &vtxs_xError()
	{
		if (not vtxs_xError_isLoaded) {
			if (vtxs_xError_branch != 0) {
				vtxs_xError_branch->GetEntry(index);
			} else { }
			vtxs_xError_isLoaded = true;
		}
		return vtxs_xError_;
	}
	vector<float> &vtxs_yError()
	{
		if (not vtxs_yError_isLoaded) {
			if (vtxs_yError_branch != 0) {
				vtxs_yError_branch->GetEntry(index);
			} else { }
			vtxs_yError_isLoaded = true;
		}
		return vtxs_yError_;
	}
	vector<float> &vtxs_zError()
	{
		if (not vtxs_zError_isLoaded) {
			if (vtxs_zError_branch != 0) {
				vtxs_zError_branch->GetEntry(index);
			} else { }
			vtxs_zError_isLoaded = true;
		}
		return vtxs_zError_;
	}
	vector<vector<float> > &els_convs_dcot()
	{
		if (not els_convs_dcot_isLoaded) {
			if (els_convs_dcot_branch != 0) {
				els_convs_dcot_branch->GetEntry(index);
			} else { }
			els_convs_dcot_isLoaded = true;
		}
		return els_convs_dcot_;
	}
	vector<vector<float> > &els_convs_dist()
	{
		if (not els_convs_dist_isLoaded) {
			if (els_convs_dist_branch != 0) {
				els_convs_dist_branch->GetEntry(index);
			} else { }
			els_convs_dist_isLoaded = true;
		}
		return els_convs_dist_;
	}
	vector<vector<float> > &els_convs_radius()
	{
		if (not els_convs_radius_isLoaded) {
			if (els_convs_radius_branch != 0) {
				els_convs_radius_branch->GetEntry(index);
			} else { }
			els_convs_radius_isLoaded = true;
		}
		return els_convs_radius_;
	}
	vector<vector<float> > &mus_stationShowerDeltaR()
	{
		if (not mus_stationShowerDeltaR_isLoaded) {
			if (mus_stationShowerDeltaR_branch != 0) {
				mus_stationShowerDeltaR_branch->GetEntry(index);
			} else { }
			mus_stationShowerDeltaR_isLoaded = true;
		}
		return mus_stationShowerDeltaR_;
	}
	vector<vector<float> > &mus_stationShowerSizeT()
	{
		if (not mus_stationShowerSizeT_isLoaded) {
			if (mus_stationShowerSizeT_branch != 0) {
				mus_stationShowerSizeT_branch->GetEntry(index);
			} else { }
			mus_stationShowerSizeT_isLoaded = true;
		}
		return mus_stationShowerSizeT_;
	}
	vector<vector<float> > &puInfo_instLumi()
	{
		if (not puInfo_instLumi_isLoaded) {
			if (puInfo_instLumi_branch != 0) {
				puInfo_instLumi_branch->GetEntry(index);
			} else { }
			puInfo_instLumi_isLoaded = true;
		}
		return puInfo_instLumi_;
	}
	vector<vector<float> > &vtxs_covMatrix()
	{
		if (not vtxs_covMatrix_isLoaded) {
			if (vtxs_covMatrix_branch != 0) {
				vtxs_covMatrix_branch->GetEntry(index);
			} else { }
			vtxs_covMatrix_isLoaded = true;
		}
		return vtxs_covMatrix_;
	}
	int &evt_bsType()
	{
		if (not evt_bsType_isLoaded) {
			if (evt_bsType_branch != 0) {
				evt_bsType_branch->GetEntry(index);
			} else { }
			evt_bsType_isLoaded = true;
		}
		return evt_bsType_;
	}
	int &evt_bunchCrossing()
	{
		if (not evt_bunchCrossing_isLoaded) {
			if (evt_bunchCrossing_branch != 0) {
				evt_bunchCrossing_branch->GetEntry(index);
			} else { }
			evt_bunchCrossing_isLoaded = true;
		}
		return evt_bunchCrossing_;
	}
	int &evt_experimentType()
	{
		if (not evt_experimentType_isLoaded) {
			if (evt_experimentType_branch != 0) {
				evt_experimentType_branch->GetEntry(index);
			} else { }
			evt_experimentType_isLoaded = true;
		}
		return evt_experimentType_;
	}
	int &evt_isRealData()
	{
		if (not evt_isRealData_isLoaded) {
			if (evt_isRealData_branch != 0) {
				evt_isRealData_branch->GetEntry(index);
			} else { }
			evt_isRealData_isLoaded = true;
		}
		return evt_isRealData_;
	}
	int &evt_orbitNumber()
	{
		if (not evt_orbitNumber_isLoaded) {
			if (evt_orbitNumber_branch != 0) {
				evt_orbitNumber_branch->GetEntry(index);
			} else { }
			evt_orbitNumber_isLoaded = true;
		}
		return evt_orbitNumber_;
	}
	int &evt_storeNumber()
	{
		if (not evt_storeNumber_isLoaded) {
			if (evt_storeNumber_branch != 0) {
				evt_storeNumber_branch->GetEntry(index);
			} else { }
			evt_storeNumber_isLoaded = true;
		}
		return evt_storeNumber_;
	}
	int &hcalnoise_GetRecHitCount()
	{
		if (not hcalnoise_GetRecHitCount_isLoaded) {
			if (hcalnoise_GetRecHitCount_branch != 0) {
				hcalnoise_GetRecHitCount_branch->GetEntry(index);
			} else { }
			hcalnoise_GetRecHitCount_isLoaded = true;
		}
		return hcalnoise_GetRecHitCount_;
	}
	int &hcalnoise_GetRecHitCount15()
	{
		if (not hcalnoise_GetRecHitCount15_isLoaded) {
			if (hcalnoise_GetRecHitCount15_branch != 0) {
				hcalnoise_GetRecHitCount15_branch->GetEntry(index);
			} else { }
			hcalnoise_GetRecHitCount15_isLoaded = true;
		}
		return hcalnoise_GetRecHitCount15_;
	}
	int &hcalnoise_maxHPDHits()
	{
		if (not hcalnoise_maxHPDHits_isLoaded) {
			if (hcalnoise_maxHPDHits_branch != 0) {
				hcalnoise_maxHPDHits_branch->GetEntry(index);
			} else { }
			hcalnoise_maxHPDHits_isLoaded = true;
		}
		return hcalnoise_maxHPDHits_;
	}
	int &hcalnoise_maxHPDNoOtherHits()
	{
		if (not hcalnoise_maxHPDNoOtherHits_isLoaded) {
			if (hcalnoise_maxHPDNoOtherHits_branch != 0) {
				hcalnoise_maxHPDNoOtherHits_branch->GetEntry(index);
			} else { }
			hcalnoise_maxHPDNoOtherHits_isLoaded = true;
		}
		return hcalnoise_maxHPDNoOtherHits_;
	}
	int &hcalnoise_maxRBXHits()
	{
		if (not hcalnoise_maxRBXHits_isLoaded) {
			if (hcalnoise_maxRBXHits_branch != 0) {
				hcalnoise_maxRBXHits_branch->GetEntry(index);
			} else { }
			hcalnoise_maxRBXHits_isLoaded = true;
		}
		return hcalnoise_maxRBXHits_;
	}
	int &hcalnoise_maxZeros()
	{
		if (not hcalnoise_maxZeros_isLoaded) {
			if (hcalnoise_maxZeros_branch != 0) {
				hcalnoise_maxZeros_branch->GetEntry(index);
			} else { }
			hcalnoise_maxZeros_isLoaded = true;
		}
		return hcalnoise_maxZeros_;
	}
	int &hcalnoise_noiseFilterStatus()
	{
		if (not hcalnoise_noiseFilterStatus_isLoaded) {
			if (hcalnoise_noiseFilterStatus_branch != 0) {
				hcalnoise_noiseFilterStatus_branch->GetEntry(index);
			} else { }
			hcalnoise_noiseFilterStatus_isLoaded = true;
		}
		return hcalnoise_noiseFilterStatus_;
	}
	int &hcalnoise_noiseType()
	{
		if (not hcalnoise_noiseType_isLoaded) {
			if (hcalnoise_noiseType_branch != 0) {
				hcalnoise_noiseType_branch->GetEntry(index);
			} else { }
			hcalnoise_noiseType_isLoaded = true;
		}
		return hcalnoise_noiseType_;
	}
	int &hcalnoise_num10GeVHits()
	{
		if (not hcalnoise_num10GeVHits_isLoaded) {
			if (hcalnoise_num10GeVHits_branch != 0) {
				hcalnoise_num10GeVHits_branch->GetEntry(index);
			} else { }
			hcalnoise_num10GeVHits_isLoaded = true;
		}
		return hcalnoise_num10GeVHits_;
	}
	int &hcalnoise_num25GeVHits()
	{
		if (not hcalnoise_num25GeVHits_isLoaded) {
			if (hcalnoise_num25GeVHits_branch != 0) {
				hcalnoise_num25GeVHits_branch->GetEntry(index);
			} else { }
			hcalnoise_num25GeVHits_isLoaded = true;
		}
		return hcalnoise_num25GeVHits_;
	}
	int &hcalnoise_numFlatNoiseChannels()
	{
		if (not hcalnoise_numFlatNoiseChannels_isLoaded) {
			if (hcalnoise_numFlatNoiseChannels_branch != 0) {
				hcalnoise_numFlatNoiseChannels_branch->GetEntry(index);
			} else { }
			hcalnoise_numFlatNoiseChannels_isLoaded = true;
		}
		return hcalnoise_numFlatNoiseChannels_;
	}
	int &hcalnoise_numIsolatedNoiseChannels()
	{
		if (not hcalnoise_numIsolatedNoiseChannels_isLoaded) {
			if (hcalnoise_numIsolatedNoiseChannels_branch != 0) {
				hcalnoise_numIsolatedNoiseChannels_branch->GetEntry(index);
			} else { }
			hcalnoise_numIsolatedNoiseChannels_isLoaded = true;
		}
		return hcalnoise_numIsolatedNoiseChannels_;
	}
	int &hcalnoise_numProblematicRBXs()
	{
		if (not hcalnoise_numProblematicRBXs_isLoaded) {
			if (hcalnoise_numProblematicRBXs_branch != 0) {
				hcalnoise_numProblematicRBXs_branch->GetEntry(index);
			} else { }
			hcalnoise_numProblematicRBXs_isLoaded = true;
		}
		return hcalnoise_numProblematicRBXs_;
	}
	int &hcalnoise_numSpikeNoiseChannels()
	{
		if (not hcalnoise_numSpikeNoiseChannels_isLoaded) {
			if (hcalnoise_numSpikeNoiseChannels_branch != 0) {
				hcalnoise_numSpikeNoiseChannels_branch->GetEntry(index);
			} else { }
			hcalnoise_numSpikeNoiseChannels_isLoaded = true;
		}
		return hcalnoise_numSpikeNoiseChannels_;
	}
	int &hcalnoise_numTS4TS5NoiseChannels()
	{
		if (not hcalnoise_numTS4TS5NoiseChannels_isLoaded) {
			if (hcalnoise_numTS4TS5NoiseChannels_branch != 0) {
				hcalnoise_numTS4TS5NoiseChannels_branch->GetEntry(index);
			} else { }
			hcalnoise_numTS4TS5NoiseChannels_isLoaded = true;
		}
		return hcalnoise_numTS4TS5NoiseChannels_;
	}
	int &hcalnoise_numTriangleNoiseChannels()
	{
		if (not hcalnoise_numTriangleNoiseChannels_isLoaded) {
			if (hcalnoise_numTriangleNoiseChannels_branch != 0) {
				hcalnoise_numTriangleNoiseChannels_branch->GetEntry(index);
			} else { }
			hcalnoise_numTriangleNoiseChannels_isLoaded = true;
		}
		return hcalnoise_numTriangleNoiseChannels_;
	}
	int &hcalnoise_passHighLevelNoiseFilter()
	{
		if (not hcalnoise_passHighLevelNoiseFilter_isLoaded) {
			if (hcalnoise_passHighLevelNoiseFilter_branch != 0) {
				hcalnoise_passHighLevelNoiseFilter_branch->GetEntry(index);
			} else { }
			hcalnoise_passHighLevelNoiseFilter_isLoaded = true;
		}
		return hcalnoise_passHighLevelNoiseFilter_;
	}
	int &hcalnoise_passLooseNoiseFilter()
	{
		if (not hcalnoise_passLooseNoiseFilter_isLoaded) {
			if (hcalnoise_passLooseNoiseFilter_branch != 0) {
				hcalnoise_passLooseNoiseFilter_branch->GetEntry(index);
			} else { }
			hcalnoise_passLooseNoiseFilter_isLoaded = true;
		}
		return hcalnoise_passLooseNoiseFilter_;
	}
	int &hcalnoise_passTightNoiseFilter()
	{
		if (not hcalnoise_passTightNoiseFilter_isLoaded) {
			if (hcalnoise_passTightNoiseFilter_branch != 0) {
				hcalnoise_passTightNoiseFilter_branch->GetEntry(index);
			} else { }
			hcalnoise_passTightNoiseFilter_isLoaded = true;
		}
		return hcalnoise_passTightNoiseFilter_;
	}
	int &pdfinfo_id1()
	{
		if (not pdfinfo_id1_isLoaded) {
			if (pdfinfo_id1_branch != 0) {
				pdfinfo_id1_branch->GetEntry(index);
			} else { }
			pdfinfo_id1_isLoaded = true;
		}
		return pdfinfo_id1_;
	}
	int &pdfinfo_id2()
	{
		if (not pdfinfo_id2_isLoaded) {
			if (pdfinfo_id2_branch != 0) {
				pdfinfo_id2_branch->GetEntry(index);
			} else { }
			pdfinfo_id2_isLoaded = true;
		}
		return pdfinfo_id2_;
	}
	vector<int> &els_mc3_id()
	{
		if (not els_mc3_id_isLoaded) {
			if (els_mc3_id_branch != 0) {
				els_mc3_id_branch->GetEntry(index);
			} else { }
			els_mc3_id_isLoaded = true;
		}
		return els_mc3_id_;
	}
	vector<int> &els_mc3idx()
	{
		if (not els_mc3idx_isLoaded) {
			if (els_mc3idx_branch != 0) {
				els_mc3idx_branch->GetEntry(index);
			} else { }
			els_mc3idx_isLoaded = true;
		}
		return els_mc3idx_;
	}
	vector<int> &els_mc3_motherid()
	{
		if (not els_mc3_motherid_isLoaded) {
			if (els_mc3_motherid_branch != 0) {
				els_mc3_motherid_branch->GetEntry(index);
			} else { }
			els_mc3_motherid_isLoaded = true;
		}
		return els_mc3_motherid_;
	}
	vector<int> &els_mc3_motheridx()
	{
		if (not els_mc3_motheridx_isLoaded) {
			if (els_mc3_motheridx_branch != 0) {
				els_mc3_motheridx_branch->GetEntry(index);
			} else { }
			els_mc3_motheridx_isLoaded = true;
		}
		return els_mc3_motheridx_;
	}
	vector<int> &els_mc_id()
	{
		if (not els_mc_id_isLoaded) {
			if (els_mc_id_branch != 0) {
				els_mc_id_branch->GetEntry(index);
			} else { }
			els_mc_id_isLoaded = true;
		}
		return els_mc_id_;
	}
	vector<int> &els_mcidx()
	{
		if (not els_mcidx_isLoaded) {
			if (els_mcidx_branch != 0) {
				els_mcidx_branch->GetEntry(index);
			} else { }
			els_mcidx_isLoaded = true;
		}
		return els_mcidx_;
	}
	vector<int> &els_mc_motherid()
	{
		if (not els_mc_motherid_isLoaded) {
			if (els_mc_motherid_branch != 0) {
				els_mc_motherid_branch->GetEntry(index);
			} else { }
			els_mc_motherid_isLoaded = true;
		}
		return els_mc_motherid_;
	}
	vector<int> &mus_mc3_id()
	{
		if (not mus_mc3_id_isLoaded) {
			if (mus_mc3_id_branch != 0) {
				mus_mc3_id_branch->GetEntry(index);
			} else { }
			mus_mc3_id_isLoaded = true;
		}
		return mus_mc3_id_;
	}
	vector<int> &mus_mc3idx()
	{
		if (not mus_mc3idx_isLoaded) {
			if (mus_mc3idx_branch != 0) {
				mus_mc3idx_branch->GetEntry(index);
			} else { }
			mus_mc3idx_isLoaded = true;
		}
		return mus_mc3idx_;
	}
	vector<int> &mus_mc3_motherid()
	{
		if (not mus_mc3_motherid_isLoaded) {
			if (mus_mc3_motherid_branch != 0) {
				mus_mc3_motherid_branch->GetEntry(index);
			} else { }
			mus_mc3_motherid_isLoaded = true;
		}
		return mus_mc3_motherid_;
	}
	vector<int> &mus_mc3_motheridx()
	{
		if (not mus_mc3_motheridx_isLoaded) {
			if (mus_mc3_motheridx_branch != 0) {
				mus_mc3_motheridx_branch->GetEntry(index);
			} else { }
			mus_mc3_motheridx_isLoaded = true;
		}
		return mus_mc3_motheridx_;
	}
	vector<int> &mus_mc_id()
	{
		if (not mus_mc_id_isLoaded) {
			if (mus_mc_id_branch != 0) {
				mus_mc_id_branch->GetEntry(index);
			} else { }
			mus_mc_id_isLoaded = true;
		}
		return mus_mc_id_;
	}
	vector<int> &mus_mcidx()
	{
		if (not mus_mcidx_isLoaded) {
			if (mus_mcidx_branch != 0) {
				mus_mcidx_branch->GetEntry(index);
			} else { }
			mus_mcidx_isLoaded = true;
		}
		return mus_mcidx_;
	}
	vector<int> &mus_mc_motherid()
	{
		if (not mus_mc_motherid_isLoaded) {
			if (mus_mc_motherid_branch != 0) {
				mus_mc_motherid_branch->GetEntry(index);
			} else { }
			mus_mc_motherid_isLoaded = true;
		}
		return mus_mc_motherid_;
	}
	vector<int> &pfjets_mc3_id()
	{
		if (not pfjets_mc3_id_isLoaded) {
			if (pfjets_mc3_id_branch != 0) {
				pfjets_mc3_id_branch->GetEntry(index);
			} else { }
			pfjets_mc3_id_isLoaded = true;
		}
		return pfjets_mc3_id_;
	}
	vector<int> &pfjets_mc3idx()
	{
		if (not pfjets_mc3idx_isLoaded) {
			if (pfjets_mc3idx_branch != 0) {
				pfjets_mc3idx_branch->GetEntry(index);
			} else { }
			pfjets_mc3idx_isLoaded = true;
		}
		return pfjets_mc3idx_;
	}
	vector<int> &pfjets_mc_gpidx()
	{
		if (not pfjets_mc_gpidx_isLoaded) {
			if (pfjets_mc_gpidx_branch != 0) {
				pfjets_mc_gpidx_branch->GetEntry(index);
			} else { }
			pfjets_mc_gpidx_isLoaded = true;
		}
		return pfjets_mc_gpidx_;
	}
	vector<int> &pfjets_mc_id()
	{
		if (not pfjets_mc_id_isLoaded) {
			if (pfjets_mc_id_branch != 0) {
				pfjets_mc_id_branch->GetEntry(index);
			} else { }
			pfjets_mc_id_isLoaded = true;
		}
		return pfjets_mc_id_;
	}
	vector<int> &pfjets_mcidx()
	{
		if (not pfjets_mcidx_isLoaded) {
			if (pfjets_mcidx_branch != 0) {
				pfjets_mcidx_branch->GetEntry(index);
			} else { }
			pfjets_mcidx_isLoaded = true;
		}
		return pfjets_mcidx_;
	}
	vector<int> &pfjets_mc_motherid()
	{
		if (not pfjets_mc_motherid_isLoaded) {
			if (pfjets_mc_motherid_branch != 0) {
				pfjets_mc_motherid_branch->GetEntry(index);
			} else { }
			pfjets_mc_motherid_isLoaded = true;
		}
		return pfjets_mc_motherid_;
	}
	vector<int> &photons_mc3_id()
	{
		if (not photons_mc3_id_isLoaded) {
			if (photons_mc3_id_branch != 0) {
				photons_mc3_id_branch->GetEntry(index);
			} else { }
			photons_mc3_id_isLoaded = true;
		}
		return photons_mc3_id_;
	}
	vector<int> &photons_mc3idx()
	{
		if (not photons_mc3idx_isLoaded) {
			if (photons_mc3idx_branch != 0) {
				photons_mc3idx_branch->GetEntry(index);
			} else { }
			photons_mc3idx_isLoaded = true;
		}
		return photons_mc3idx_;
	}
	vector<int> &photons_mc3_motherid()
	{
		if (not photons_mc3_motherid_isLoaded) {
			if (photons_mc3_motherid_branch != 0) {
				photons_mc3_motherid_branch->GetEntry(index);
			} else { }
			photons_mc3_motherid_isLoaded = true;
		}
		return photons_mc3_motherid_;
	}
	vector<int> &photons_mc3_motheridx()
	{
		if (not photons_mc3_motheridx_isLoaded) {
			if (photons_mc3_motheridx_branch != 0) {
				photons_mc3_motheridx_branch->GetEntry(index);
			} else { }
			photons_mc3_motheridx_isLoaded = true;
		}
		return photons_mc3_motheridx_;
	}
	vector<int> &photons_mc_id()
	{
		if (not photons_mc_id_isLoaded) {
			if (photons_mc_id_branch != 0) {
				photons_mc_id_branch->GetEntry(index);
			} else { }
			photons_mc_id_isLoaded = true;
		}
		return photons_mc_id_;
	}
	vector<int> &photons_mcidx()
	{
		if (not photons_mcidx_isLoaded) {
			if (photons_mcidx_branch != 0) {
				photons_mcidx_branch->GetEntry(index);
			} else { }
			photons_mcidx_isLoaded = true;
		}
		return photons_mcidx_;
	}
	vector<int> &photons_mc_motherid()
	{
		if (not photons_mc_motherid_isLoaded) {
			if (photons_mc_motherid_branch != 0) {
				photons_mc_motherid_branch->GetEntry(index);
			} else { }
			photons_mc_motherid_isLoaded = true;
		}
		return photons_mc_motherid_;
	}
	vector<int> &els_category()
	{
		if (not els_category_isLoaded) {
			if (els_category_branch != 0) {
				els_category_branch->GetEntry(index);
			} else { }
			els_category_isLoaded = true;
		}
		return els_category_;
	}
	vector<int> &els_charge()
	{
		if (not els_charge_isLoaded) {
			if (els_charge_branch != 0) {
				els_charge_branch->GetEntry(index);
			} else { }
			els_charge_isLoaded = true;
		}
		return els_charge_;
	}
	vector<int> &els_ckf_charge()
	{
		if (not els_ckf_charge_isLoaded) {
			if (els_ckf_charge_branch != 0) {
				els_ckf_charge_branch->GetEntry(index);
			} else { }
			els_ckf_charge_isLoaded = true;
		}
		return els_ckf_charge_;
	}
	vector<int> &els_ckf_laywithmeas()
	{
		if (not els_ckf_laywithmeas_isLoaded) {
			if (els_ckf_laywithmeas_branch != 0) {
				els_ckf_laywithmeas_branch->GetEntry(index);
			} else { }
			els_ckf_laywithmeas_isLoaded = true;
		}
		return els_ckf_laywithmeas_;
	}
	vector<int> &els_class()
	{
		if (not els_class_isLoaded) {
			if (els_class_branch != 0) {
				els_class_branch->GetEntry(index);
			} else { }
			els_class_isLoaded = true;
		}
		return els_class_;
	}
	vector<int> &els_conv_delMissHits()
	{
		if (not els_conv_delMissHits_isLoaded) {
			if (els_conv_delMissHits_branch != 0) {
				els_conv_delMissHits_branch->GetEntry(index);
			} else { }
			els_conv_delMissHits_isLoaded = true;
		}
		return els_conv_delMissHits_;
	}
	vector<int> &els_conv_flag()
	{
		if (not els_conv_flag_isLoaded) {
			if (els_conv_flag_branch != 0) {
				els_conv_flag_branch->GetEntry(index);
			} else { }
			els_conv_flag_isLoaded = true;
		}
		return els_conv_flag_;
	}
	vector<int> &els_conv_gsftkidx()
	{
		if (not els_conv_gsftkidx_isLoaded) {
			if (els_conv_gsftkidx_branch != 0) {
				els_conv_gsftkidx_branch->GetEntry(index);
			} else { }
			els_conv_gsftkidx_isLoaded = true;
		}
		return els_conv_gsftkidx_;
	}
	vector<int> &els_conv_old_delMissHits()
	{
		if (not els_conv_old_delMissHits_isLoaded) {
			if (els_conv_old_delMissHits_branch != 0) {
				els_conv_old_delMissHits_branch->GetEntry(index);
			} else { }
			els_conv_old_delMissHits_isLoaded = true;
		}
		return els_conv_old_delMissHits_;
	}
	vector<int> &els_conv_old_flag()
	{
		if (not els_conv_old_flag_isLoaded) {
			if (els_conv_old_flag_branch != 0) {
				els_conv_old_flag_branch->GetEntry(index);
			} else { }
			els_conv_old_flag_isLoaded = true;
		}
		return els_conv_old_flag_;
	}
	vector<int> &els_conv_old_gsftkidx()
	{
		if (not els_conv_old_gsftkidx_isLoaded) {
			if (els_conv_old_gsftkidx_branch != 0) {
				els_conv_old_gsftkidx_branch->GetEntry(index);
			} else { }
			els_conv_old_gsftkidx_isLoaded = true;
		}
		return els_conv_old_gsftkidx_;
	}
	vector<int> &els_conv_old_tkidx()
	{
		if (not els_conv_old_tkidx_isLoaded) {
			if (els_conv_old_tkidx_branch != 0) {
				els_conv_old_tkidx_branch->GetEntry(index);
			} else { }
			els_conv_old_tkidx_isLoaded = true;
		}
		return els_conv_old_tkidx_;
	}
	vector<int> &els_conv_tkidx()
	{
		if (not els_conv_tkidx_isLoaded) {
			if (els_conv_tkidx_branch != 0) {
				els_conv_tkidx_branch->GetEntry(index);
			} else { }
			els_conv_tkidx_isLoaded = true;
		}
		return els_conv_tkidx_;
	}
	vector<int> &els_exp_innerlayers()
	{
		if (not els_exp_innerlayers_isLoaded) {
			if (els_exp_innerlayers_branch != 0) {
				els_exp_innerlayers_branch->GetEntry(index);
			} else { }
			els_exp_innerlayers_isLoaded = true;
		}
		return els_exp_innerlayers_;
	}
	vector<int> &els_exp_outerlayers()
	{
		if (not els_exp_outerlayers_isLoaded) {
			if (els_exp_outerlayers_branch != 0) {
				els_exp_outerlayers_branch->GetEntry(index);
			} else { }
			els_exp_outerlayers_isLoaded = true;
		}
		return els_exp_outerlayers_;
	}
	vector<int> &els_fiduciality()
	{
		if (not els_fiduciality_isLoaded) {
			if (els_fiduciality_branch != 0) {
				els_fiduciality_branch->GetEntry(index);
			} else { }
			els_fiduciality_isLoaded = true;
		}
		return els_fiduciality_;
	}
	vector<int> &els_gsftrkidx()
	{
		if (not els_gsftrkidx_isLoaded) {
			if (els_gsftrkidx_branch != 0) {
				els_gsftrkidx_branch->GetEntry(index);
			} else { }
			els_gsftrkidx_isLoaded = true;
		}
		return els_gsftrkidx_;
	}
	vector<int> &els_layer1_det()
	{
		if (not els_layer1_det_isLoaded) {
			if (els_layer1_det_branch != 0) {
				els_layer1_det_branch->GetEntry(index);
			} else { }
			els_layer1_det_isLoaded = true;
		}
		return els_layer1_det_;
	}
	vector<int> &els_layer1_layer()
	{
		if (not els_layer1_layer_isLoaded) {
			if (els_layer1_layer_branch != 0) {
				els_layer1_layer_branch->GetEntry(index);
			} else { }
			els_layer1_layer_isLoaded = true;
		}
		return els_layer1_layer_;
	}
	vector<int> &els_layer1_sizerphi()
	{
		if (not els_layer1_sizerphi_isLoaded) {
			if (els_layer1_sizerphi_branch != 0) {
				els_layer1_sizerphi_branch->GetEntry(index);
			} else { }
			els_layer1_sizerphi_isLoaded = true;
		}
		return els_layer1_sizerphi_;
	}
	vector<int> &els_layer1_sizerz()
	{
		if (not els_layer1_sizerz_isLoaded) {
			if (els_layer1_sizerz_branch != 0) {
				els_layer1_sizerz_branch->GetEntry(index);
			} else { }
			els_layer1_sizerz_isLoaded = true;
		}
		return els_layer1_sizerz_;
	}
	vector<int> &els_lostHits()
	{
		if (not els_lostHits_isLoaded) {
			if (els_lostHits_branch != 0) {
				els_lostHits_branch->GetEntry(index);
			} else { }
			els_lostHits_isLoaded = true;
		}
		return els_lostHits_;
	}
	vector<int> &els_lost_pixelhits()
	{
		if (not els_lost_pixelhits_isLoaded) {
			if (els_lost_pixelhits_branch != 0) {
				els_lost_pixelhits_branch->GetEntry(index);
			} else { }
			els_lost_pixelhits_isLoaded = true;
		}
		return els_lost_pixelhits_;
	}
	vector<int> &els_mc_patMatch_id()
	{
		if (not els_mc_patMatch_id_isLoaded) {
			if (els_mc_patMatch_id_branch != 0) {
				els_mc_patMatch_id_branch->GetEntry(index);
			} else { }
			els_mc_patMatch_id_isLoaded = true;
		}
		return els_mc_patMatch_id_;
	}
	vector<int> &els_nSeed()
	{
		if (not els_nSeed_isLoaded) {
			if (els_nSeed_branch != 0) {
				els_nSeed_branch->GetEntry(index);
			} else { }
			els_nSeed_isLoaded = true;
		}
		return els_nSeed_;
	}
	vector<int> &els_nlayers()
	{
		if (not els_nlayers_isLoaded) {
			if (els_nlayers_branch != 0) {
				els_nlayers_branch->GetEntry(index);
			} else { }
			els_nlayers_isLoaded = true;
		}
		return els_nlayers_;
	}
	vector<int> &els_nlayers3D()
	{
		if (not els_nlayers3D_isLoaded) {
			if (els_nlayers3D_branch != 0) {
				els_nlayers3D_branch->GetEntry(index);
			} else { }
			els_nlayers3D_isLoaded = true;
		}
		return els_nlayers3D_;
	}
	vector<int> &els_nlayersLost()
	{
		if (not els_nlayersLost_isLoaded) {
			if (els_nlayersLost_branch != 0) {
				els_nlayersLost_branch->GetEntry(index);
			} else { }
			els_nlayersLost_isLoaded = true;
		}
		return els_nlayersLost_;
	}
	vector<int> &els_sccharge()
	{
		if (not els_sccharge_isLoaded) {
			if (els_sccharge_branch != 0) {
				els_sccharge_branch->GetEntry(index);
			} else { }
			els_sccharge_isLoaded = true;
		}
		return els_sccharge_;
	}
	vector<int> &els_scindex()
	{
		if (not els_scindex_isLoaded) {
			if (els_scindex_branch != 0) {
				els_scindex_branch->GetEntry(index);
			} else { }
			els_scindex_isLoaded = true;
		}
		return els_scindex_;
	}
	vector<int> &els_trk_charge()
	{
		if (not els_trk_charge_isLoaded) {
			if (els_trk_charge_branch != 0) {
				els_trk_charge_branch->GetEntry(index);
			} else { }
			els_trk_charge_isLoaded = true;
		}
		return els_trk_charge_;
	}
	vector<int> &els_trkidx()
	{
		if (not els_trkidx_isLoaded) {
			if (els_trkidx_branch != 0) {
				els_trkidx_branch->GetEntry(index);
			} else { }
			els_trkidx_isLoaded = true;
		}
		return els_trkidx_;
	}
	vector<int> &els_type()
	{
		if (not els_type_isLoaded) {
			if (els_type_branch != 0) {
				els_type_branch->GetEntry(index);
			} else { }
			els_type_isLoaded = true;
		}
		return els_type_;
	}
	vector<int> &els_validHits()
	{
		if (not els_validHits_isLoaded) {
			if (els_validHits_branch != 0) {
				els_validHits_branch->GetEntry(index);
			} else { }
			els_validHits_isLoaded = true;
		}
		return els_validHits_;
	}
	vector<int> &els_valid_pixelhits()
	{
		if (not els_valid_pixelhits_isLoaded) {
			if (els_valid_pixelhits_branch != 0) {
				els_valid_pixelhits_branch->GetEntry(index);
			} else { }
			els_valid_pixelhits_isLoaded = true;
		}
		return els_valid_pixelhits_;
	}
	vector<int> &genps_id()
	{
		if (not genps_id_isLoaded) {
			if (genps_id_branch != 0) {
				genps_id_branch->GetEntry(index);
			} else { }
			genps_id_isLoaded = true;
		}
		return genps_id_;
	}
	vector<int> &genps_id_mother()
	{
		if (not genps_id_mother_isLoaded) {
			if (genps_id_mother_branch != 0) {
				genps_id_mother_branch->GetEntry(index);
			} else { }
			genps_id_mother_isLoaded = true;
		}
		return genps_id_mother_;
	}
	vector<int> &genps_status()
	{
		if (not genps_status_isLoaded) {
			if (genps_status_branch != 0) {
				genps_status_branch->GetEntry(index);
			} else { }
			genps_status_isLoaded = true;
		}
		return genps_status_;
	}
	vector<int> &mus_algo()
	{
		if (not mus_algo_isLoaded) {
			if (mus_algo_branch != 0) {
				mus_algo_branch->GetEntry(index);
			} else { }
			mus_algo_isLoaded = true;
		}
		return mus_algo_;
	}
	vector<int> &mus_charge()
	{
		if (not mus_charge_isLoaded) {
			if (mus_charge_branch != 0) {
				mus_charge_branch->GetEntry(index);
			} else { }
			mus_charge_isLoaded = true;
		}
		return mus_charge_;
	}
	vector<int> &mus_ecal_rawId()
	{
		if (not mus_ecal_rawId_isLoaded) {
			if (mus_ecal_rawId_branch != 0) {
				mus_ecal_rawId_branch->GetEntry(index);
			} else { }
			mus_ecal_rawId_isLoaded = true;
		}
		return mus_ecal_rawId_;
	}
	vector<int> &mus_exp_innerlayers()
	{
		if (not mus_exp_innerlayers_isLoaded) {
			if (mus_exp_innerlayers_branch != 0) {
				mus_exp_innerlayers_branch->GetEntry(index);
			} else { }
			mus_exp_innerlayers_isLoaded = true;
		}
		return mus_exp_innerlayers_;
	}
	vector<int> &mus_exp_outerlayers()
	{
		if (not mus_exp_outerlayers_isLoaded) {
			if (mus_exp_outerlayers_branch != 0) {
				mus_exp_outerlayers_branch->GetEntry(index);
			} else { }
			mus_exp_outerlayers_isLoaded = true;
		}
		return mus_exp_outerlayers_;
	}
	vector<int> &mus_gfit_algo()
	{
		if (not mus_gfit_algo_isLoaded) {
			if (mus_gfit_algo_branch != 0) {
				mus_gfit_algo_branch->GetEntry(index);
			} else { }
			mus_gfit_algo_isLoaded = true;
		}
		return mus_gfit_algo_;
	}
	vector<int> &mus_gfit_etaErr()
	{
		if (not mus_gfit_etaErr_isLoaded) {
			if (mus_gfit_etaErr_branch != 0) {
				mus_gfit_etaErr_branch->GetEntry(index);
			} else { }
			mus_gfit_etaErr_isLoaded = true;
		}
		return mus_gfit_etaErr_;
	}
	vector<int> &mus_gfit_exp_innerlayer()
	{
		if (not mus_gfit_exp_innerlayer_isLoaded) {
			if (mus_gfit_exp_innerlayer_branch != 0) {
				mus_gfit_exp_innerlayer_branch->GetEntry(index);
			} else { }
			mus_gfit_exp_innerlayer_isLoaded = true;
		}
		return mus_gfit_exp_innerlayer_;
	}
	vector<int> &mus_gfit_exp_outerlayer()
	{
		if (not mus_gfit_exp_outerlayer_isLoaded) {
			if (mus_gfit_exp_outerlayer_branch != 0) {
				mus_gfit_exp_outerlayer_branch->GetEntry(index);
			} else { }
			mus_gfit_exp_outerlayer_isLoaded = true;
		}
		return mus_gfit_exp_outerlayer_;
	}
	vector<int> &mus_gfit_lostHits()
	{
		if (not mus_gfit_lostHits_isLoaded) {
			if (mus_gfit_lostHits_branch != 0) {
				mus_gfit_lostHits_branch->GetEntry(index);
			} else { }
			mus_gfit_lostHits_isLoaded = true;
		}
		return mus_gfit_lostHits_;
	}
	vector<int> &mus_gfit_lostPixelHits()
	{
		if (not mus_gfit_lostPixelHits_isLoaded) {
			if (mus_gfit_lostPixelHits_branch != 0) {
				mus_gfit_lostPixelHits_branch->GetEntry(index);
			} else { }
			mus_gfit_lostPixelHits_isLoaded = true;
		}
		return mus_gfit_lostPixelHits_;
	}
	vector<int> &mus_gfit_nlayers()
	{
		if (not mus_gfit_nlayers_isLoaded) {
			if (mus_gfit_nlayers_branch != 0) {
				mus_gfit_nlayers_branch->GetEntry(index);
			} else { }
			mus_gfit_nlayers_isLoaded = true;
		}
		return mus_gfit_nlayers_;
	}
	vector<int> &mus_gfit_nlayers3D()
	{
		if (not mus_gfit_nlayers3D_isLoaded) {
			if (mus_gfit_nlayers3D_branch != 0) {
				mus_gfit_nlayers3D_branch->GetEntry(index);
			} else { }
			mus_gfit_nlayers3D_isLoaded = true;
		}
		return mus_gfit_nlayers3D_;
	}
	vector<int> &mus_gfit_nlayersLost()
	{
		if (not mus_gfit_nlayersLost_isLoaded) {
			if (mus_gfit_nlayersLost_branch != 0) {
				mus_gfit_nlayersLost_branch->GetEntry(index);
			} else { }
			mus_gfit_nlayersLost_isLoaded = true;
		}
		return mus_gfit_nlayersLost_;
	}
	vector<int> &mus_gfit_phiErr()
	{
		if (not mus_gfit_phiErr_isLoaded) {
			if (mus_gfit_phiErr_branch != 0) {
				mus_gfit_phiErr_branch->GetEntry(index);
			} else { }
			mus_gfit_phiErr_isLoaded = true;
		}
		return mus_gfit_phiErr_;
	}
	vector<int> &mus_gfit_ptErr()
	{
		if (not mus_gfit_ptErr_isLoaded) {
			if (mus_gfit_ptErr_branch != 0) {
				mus_gfit_ptErr_branch->GetEntry(index);
			} else { }
			mus_gfit_ptErr_isLoaded = true;
		}
		return mus_gfit_ptErr_;
	}
	vector<int> &mus_gfit_qualityMask()
	{
		if (not mus_gfit_qualityMask_isLoaded) {
			if (mus_gfit_qualityMask_branch != 0) {
				mus_gfit_qualityMask_branch->GetEntry(index);
			} else { }
			mus_gfit_qualityMask_isLoaded = true;
		}
		return mus_gfit_qualityMask_;
	}
	vector<int> &mus_gfit_trk_charge()
	{
		if (not mus_gfit_trk_charge_isLoaded) {
			if (mus_gfit_trk_charge_branch != 0) {
				mus_gfit_trk_charge_branch->GetEntry(index);
			} else { }
			mus_gfit_trk_charge_isLoaded = true;
		}
		return mus_gfit_trk_charge_;
	}
	vector<int> &mus_gfit_validHits()
	{
		if (not mus_gfit_validHits_isLoaded) {
			if (mus_gfit_validHits_branch != 0) {
				mus_gfit_validHits_branch->GetEntry(index);
			} else { }
			mus_gfit_validHits_isLoaded = true;
		}
		return mus_gfit_validHits_;
	}
	vector<int> &mus_gfit_validPixelHits()
	{
		if (not mus_gfit_validPixelHits_isLoaded) {
			if (mus_gfit_validPixelHits_branch != 0) {
				mus_gfit_validPixelHits_branch->GetEntry(index);
			} else { }
			mus_gfit_validPixelHits_isLoaded = true;
		}
		return mus_gfit_validPixelHits_;
	}
	vector<int> &mus_gfit_validSTAHits()
	{
		if (not mus_gfit_validSTAHits_isLoaded) {
			if (mus_gfit_validSTAHits_branch != 0) {
				mus_gfit_validSTAHits_branch->GetEntry(index);
			} else { }
			mus_gfit_validSTAHits_isLoaded = true;
		}
		return mus_gfit_validSTAHits_;
	}
	vector<int> &mus_gfit_validSiHits()
	{
		if (not mus_gfit_validSiHits_isLoaded) {
			if (mus_gfit_validSiHits_branch != 0) {
				mus_gfit_validSiHits_branch->GetEntry(index);
			} else { }
			mus_gfit_validSiHits_isLoaded = true;
		}
		return mus_gfit_validSiHits_;
	}
	vector<int> &mus_goodmask()
	{
		if (not mus_goodmask_isLoaded) {
			if (mus_goodmask_branch != 0) {
				mus_goodmask_branch->GetEntry(index);
			} else { }
			mus_goodmask_isLoaded = true;
		}
		return mus_goodmask_;
	}
	vector<int> &mus_hcal_rawId()
	{
		if (not mus_hcal_rawId_isLoaded) {
			if (mus_hcal_rawId_branch != 0) {
				mus_hcal_rawId_branch->GetEntry(index);
			} else { }
			mus_hcal_rawId_isLoaded = true;
		}
		return mus_hcal_rawId_;
	}
	vector<int> &mus_iso03_njets()
	{
		if (not mus_iso03_njets_isLoaded) {
			if (mus_iso03_njets_branch != 0) {
				mus_iso03_njets_branch->GetEntry(index);
			} else { }
			mus_iso03_njets_isLoaded = true;
		}
		return mus_iso03_njets_;
	}
	vector<int> &mus_iso03_ntrk()
	{
		if (not mus_iso03_ntrk_isLoaded) {
			if (mus_iso03_ntrk_branch != 0) {
				mus_iso03_ntrk_branch->GetEntry(index);
			} else { }
			mus_iso03_ntrk_isLoaded = true;
		}
		return mus_iso03_ntrk_;
	}
	vector<int> &mus_iso05_ntrk()
	{
		if (not mus_iso05_ntrk_isLoaded) {
			if (mus_iso05_ntrk_branch != 0) {
				mus_iso05_ntrk_branch->GetEntry(index);
			} else { }
			mus_iso05_ntrk_isLoaded = true;
		}
		return mus_iso05_ntrk_;
	}
	vector<int> &mus_lostHits()
	{
		if (not mus_lostHits_isLoaded) {
			if (mus_lostHits_branch != 0) {
				mus_lostHits_branch->GetEntry(index);
			} else { }
			mus_lostHits_isLoaded = true;
		}
		return mus_lostHits_;
	}
	vector<int> &mus_lostPixelHits()
	{
		if (not mus_lostPixelHits_isLoaded) {
			if (mus_lostPixelHits_branch != 0) {
				mus_lostPixelHits_branch->GetEntry(index);
			} else { }
			mus_lostPixelHits_isLoaded = true;
		}
		return mus_lostPixelHits_;
	}
	vector<int> &mus_mc_patMatch_id()
	{
		if (not mus_mc_patMatch_id_isLoaded) {
			if (mus_mc_patMatch_id_branch != 0) {
				mus_mc_patMatch_id_branch->GetEntry(index);
			} else { }
			mus_mc_patMatch_id_isLoaded = true;
		}
		return mus_mc_patMatch_id_;
	}
	vector<int> &mus_muonBestTrackType()
	{
		if (not mus_muonBestTrackType_isLoaded) {
			if (mus_muonBestTrackType_branch != 0) {
				mus_muonBestTrackType_branch->GetEntry(index);
			} else { }
			mus_muonBestTrackType_isLoaded = true;
		}
		return mus_muonBestTrackType_;
	}
	vector<int> &mus_nOverlaps()
	{
		if (not mus_nOverlaps_isLoaded) {
			if (mus_nOverlaps_branch != 0) {
				mus_nOverlaps_branch->GetEntry(index);
			} else { }
			mus_nOverlaps_isLoaded = true;
		}
		return mus_nOverlaps_;
	}
	vector<int> &mus_nlayers()
	{
		if (not mus_nlayers_isLoaded) {
			if (mus_nlayers_branch != 0) {
				mus_nlayers_branch->GetEntry(index);
			} else { }
			mus_nlayers_isLoaded = true;
		}
		return mus_nlayers_;
	}
	vector<int> &mus_nlayers3D()
	{
		if (not mus_nlayers3D_isLoaded) {
			if (mus_nlayers3D_branch != 0) {
				mus_nlayers3D_branch->GetEntry(index);
			} else { }
			mus_nlayers3D_isLoaded = true;
		}
		return mus_nlayers3D_;
	}
	vector<int> &mus_nlayersLost()
	{
		if (not mus_nlayersLost_isLoaded) {
			if (mus_nlayersLost_branch != 0) {
				mus_nlayersLost_branch->GetEntry(index);
			} else { }
			mus_nlayersLost_isLoaded = true;
		}
		return mus_nlayersLost_;
	}
	vector<int> &mus_nmatches()
	{
		if (not mus_nmatches_isLoaded) {
			if (mus_nmatches_branch != 0) {
				mus_nmatches_branch->GetEntry(index);
			} else { }
			mus_nmatches_isLoaded = true;
		}
		return mus_nmatches_;
	}
	vector<int> &mus_numberOfMatchedStations()
	{
		if (not mus_numberOfMatchedStations_isLoaded) {
			if (mus_numberOfMatchedStations_branch != 0) {
				mus_numberOfMatchedStations_branch->GetEntry(index);
			} else { }
			mus_numberOfMatchedStations_isLoaded = true;
		}
		return mus_numberOfMatchedStations_;
	}
	vector<int> &mus_overlap0()
	{
		if (not mus_overlap0_isLoaded) {
			if (mus_overlap0_branch != 0) {
				mus_overlap0_branch->GetEntry(index);
			} else { }
			mus_overlap0_isLoaded = true;
		}
		return mus_overlap0_;
	}
	vector<int> &mus_overlap1()
	{
		if (not mus_overlap1_isLoaded) {
			if (mus_overlap1_branch != 0) {
				mus_overlap1_branch->GetEntry(index);
			} else { }
			mus_overlap1_isLoaded = true;
		}
		return mus_overlap1_;
	}
	vector<int> &mus_pfcharge()
	{
		if (not mus_pfcharge_isLoaded) {
			if (mus_pfcharge_branch != 0) {
				mus_pfcharge_branch->GetEntry(index);
			} else { }
			mus_pfcharge_isLoaded = true;
		}
		return mus_pfcharge_;
	}
	vector<int> &mus_pfparticleId()
	{
		if (not mus_pfparticleId_isLoaded) {
			if (mus_pfparticleId_branch != 0) {
				mus_pfparticleId_branch->GetEntry(index);
			} else { }
			mus_pfparticleId_isLoaded = true;
		}
		return mus_pfparticleId_;
	}
	vector<int> &mus_pid_PFMuon()
	{
		if (not mus_pid_PFMuon_isLoaded) {
			if (mus_pid_PFMuon_branch != 0) {
				mus_pid_PFMuon_branch->GetEntry(index);
			} else { }
			mus_pid_PFMuon_isLoaded = true;
		}
		return mus_pid_PFMuon_;
	}
	vector<int> &mus_pid_TM2DCompatibilityLoose()
	{
		if (not mus_pid_TM2DCompatibilityLoose_isLoaded) {
			if (mus_pid_TM2DCompatibilityLoose_branch != 0) {
				mus_pid_TM2DCompatibilityLoose_branch->GetEntry(index);
			} else { }
			mus_pid_TM2DCompatibilityLoose_isLoaded = true;
		}
		return mus_pid_TM2DCompatibilityLoose_;
	}
	vector<int> &mus_pid_TM2DCompatibilityTight()
	{
		if (not mus_pid_TM2DCompatibilityTight_isLoaded) {
			if (mus_pid_TM2DCompatibilityTight_branch != 0) {
				mus_pid_TM2DCompatibilityTight_branch->GetEntry(index);
			} else { }
			mus_pid_TM2DCompatibilityTight_isLoaded = true;
		}
		return mus_pid_TM2DCompatibilityTight_;
	}
	vector<int> &mus_pid_TMLastStationLoose()
	{
		if (not mus_pid_TMLastStationLoose_isLoaded) {
			if (mus_pid_TMLastStationLoose_branch != 0) {
				mus_pid_TMLastStationLoose_branch->GetEntry(index);
			} else { }
			mus_pid_TMLastStationLoose_isLoaded = true;
		}
		return mus_pid_TMLastStationLoose_;
	}
	vector<int> &mus_pid_TMLastStationTight()
	{
		if (not mus_pid_TMLastStationTight_isLoaded) {
			if (mus_pid_TMLastStationTight_branch != 0) {
				mus_pid_TMLastStationTight_branch->GetEntry(index);
			} else { }
			mus_pid_TMLastStationTight_isLoaded = true;
		}
		return mus_pid_TMLastStationTight_;
	}
	vector<int> &mus_qualityMask()
	{
		if (not mus_qualityMask_isLoaded) {
			if (mus_qualityMask_branch != 0) {
				mus_qualityMask_branch->GetEntry(index);
			} else { }
			mus_qualityMask_isLoaded = true;
		}
		return mus_qualityMask_;
	}
	vector<int> &mus_sta_algo()
	{
		if (not mus_sta_algo_isLoaded) {
			if (mus_sta_algo_branch != 0) {
				mus_sta_algo_branch->GetEntry(index);
			} else { }
			mus_sta_algo_isLoaded = true;
		}
		return mus_sta_algo_;
	}
	vector<int> &mus_sta_etaErr()
	{
		if (not mus_sta_etaErr_isLoaded) {
			if (mus_sta_etaErr_branch != 0) {
				mus_sta_etaErr_branch->GetEntry(index);
			} else { }
			mus_sta_etaErr_isLoaded = true;
		}
		return mus_sta_etaErr_;
	}
	vector<int> &mus_sta_exp_innerlayer()
	{
		if (not mus_sta_exp_innerlayer_isLoaded) {
			if (mus_sta_exp_innerlayer_branch != 0) {
				mus_sta_exp_innerlayer_branch->GetEntry(index);
			} else { }
			mus_sta_exp_innerlayer_isLoaded = true;
		}
		return mus_sta_exp_innerlayer_;
	}
	vector<int> &mus_sta_exp_outerlayer()
	{
		if (not mus_sta_exp_outerlayer_isLoaded) {
			if (mus_sta_exp_outerlayer_branch != 0) {
				mus_sta_exp_outerlayer_branch->GetEntry(index);
			} else { }
			mus_sta_exp_outerlayer_isLoaded = true;
		}
		return mus_sta_exp_outerlayer_;
	}
	vector<int> &mus_sta_lostHits()
	{
		if (not mus_sta_lostHits_isLoaded) {
			if (mus_sta_lostHits_branch != 0) {
				mus_sta_lostHits_branch->GetEntry(index);
			} else { }
			mus_sta_lostHits_isLoaded = true;
		}
		return mus_sta_lostHits_;
	}
	vector<int> &mus_sta_lostPixelHits()
	{
		if (not mus_sta_lostPixelHits_isLoaded) {
			if (mus_sta_lostPixelHits_branch != 0) {
				mus_sta_lostPixelHits_branch->GetEntry(index);
			} else { }
			mus_sta_lostPixelHits_isLoaded = true;
		}
		return mus_sta_lostPixelHits_;
	}
	vector<int> &mus_sta_nlayers()
	{
		if (not mus_sta_nlayers_isLoaded) {
			if (mus_sta_nlayers_branch != 0) {
				mus_sta_nlayers_branch->GetEntry(index);
			} else { }
			mus_sta_nlayers_isLoaded = true;
		}
		return mus_sta_nlayers_;
	}
	vector<int> &mus_sta_nlayers3D()
	{
		if (not mus_sta_nlayers3D_isLoaded) {
			if (mus_sta_nlayers3D_branch != 0) {
				mus_sta_nlayers3D_branch->GetEntry(index);
			} else { }
			mus_sta_nlayers3D_isLoaded = true;
		}
		return mus_sta_nlayers3D_;
	}
	vector<int> &mus_sta_nlayersLost()
	{
		if (not mus_sta_nlayersLost_isLoaded) {
			if (mus_sta_nlayersLost_branch != 0) {
				mus_sta_nlayersLost_branch->GetEntry(index);
			} else { }
			mus_sta_nlayersLost_isLoaded = true;
		}
		return mus_sta_nlayersLost_;
	}
	vector<int> &mus_sta_phiErr()
	{
		if (not mus_sta_phiErr_isLoaded) {
			if (mus_sta_phiErr_branch != 0) {
				mus_sta_phiErr_branch->GetEntry(index);
			} else { }
			mus_sta_phiErr_isLoaded = true;
		}
		return mus_sta_phiErr_;
	}
	vector<int> &mus_sta_ptErr()
	{
		if (not mus_sta_ptErr_isLoaded) {
			if (mus_sta_ptErr_branch != 0) {
				mus_sta_ptErr_branch->GetEntry(index);
			} else { }
			mus_sta_ptErr_isLoaded = true;
		}
		return mus_sta_ptErr_;
	}
	vector<int> &mus_sta_qualityMask()
	{
		if (not mus_sta_qualityMask_isLoaded) {
			if (mus_sta_qualityMask_branch != 0) {
				mus_sta_qualityMask_branch->GetEntry(index);
			} else { }
			mus_sta_qualityMask_isLoaded = true;
		}
		return mus_sta_qualityMask_;
	}
	vector<int> &mus_sta_trk_charge()
	{
		if (not mus_sta_trk_charge_isLoaded) {
			if (mus_sta_trk_charge_branch != 0) {
				mus_sta_trk_charge_branch->GetEntry(index);
			} else { }
			mus_sta_trk_charge_isLoaded = true;
		}
		return mus_sta_trk_charge_;
	}
	vector<int> &mus_sta_validHits()
	{
		if (not mus_sta_validHits_isLoaded) {
			if (mus_sta_validHits_branch != 0) {
				mus_sta_validHits_branch->GetEntry(index);
			} else { }
			mus_sta_validHits_isLoaded = true;
		}
		return mus_sta_validHits_;
	}
	vector<int> &mus_sta_validPixelHits()
	{
		if (not mus_sta_validPixelHits_isLoaded) {
			if (mus_sta_validPixelHits_branch != 0) {
				mus_sta_validPixelHits_branch->GetEntry(index);
			} else { }
			mus_sta_validPixelHits_isLoaded = true;
		}
		return mus_sta_validPixelHits_;
	}
	vector<int> &mus_timeDirection()
	{
		if (not mus_timeDirection_isLoaded) {
			if (mus_timeDirection_branch != 0) {
				mus_timeDirection_branch->GetEntry(index);
			} else { }
			mus_timeDirection_isLoaded = true;
		}
		return mus_timeDirection_;
	}
	vector<int> &mus_timeNumStationsUsed()
	{
		if (not mus_timeNumStationsUsed_isLoaded) {
			if (mus_timeNumStationsUsed_branch != 0) {
				mus_timeNumStationsUsed_branch->GetEntry(index);
			} else { }
			mus_timeNumStationsUsed_isLoaded = true;
		}
		return mus_timeNumStationsUsed_;
	}
	vector<int> &mus_trk_charge()
	{
		if (not mus_trk_charge_isLoaded) {
			if (mus_trk_charge_branch != 0) {
				mus_trk_charge_branch->GetEntry(index);
			} else { }
			mus_trk_charge_isLoaded = true;
		}
		return mus_trk_charge_;
	}
	vector<int> &mus_trkidx()
	{
		if (not mus_trkidx_isLoaded) {
			if (mus_trkidx_branch != 0) {
				mus_trkidx_branch->GetEntry(index);
			} else { }
			mus_trkidx_isLoaded = true;
		}
		return mus_trkidx_;
	}
	vector<int> &mus_type()
	{
		if (not mus_type_isLoaded) {
			if (mus_type_branch != 0) {
				mus_type_branch->GetEntry(index);
			} else { }
			mus_type_isLoaded = true;
		}
		return mus_type_;
	}
	vector<int> &mus_validHits()
	{
		if (not mus_validHits_isLoaded) {
			if (mus_validHits_branch != 0) {
				mus_validHits_branch->GetEntry(index);
			} else { }
			mus_validHits_isLoaded = true;
		}
		return mus_validHits_;
	}
	vector<int> &mus_validPixelHits()
	{
		if (not mus_validPixelHits_isLoaded) {
			if (mus_validPixelHits_branch != 0) {
				mus_validPixelHits_branch->GetEntry(index);
			} else { }
			mus_validPixelHits_isLoaded = true;
		}
		return mus_validPixelHits_;
	}
	vector<int> &pfcands_charge()
	{
		if (not pfcands_charge_isLoaded) {
			if (pfcands_charge_branch != 0) {
				pfcands_charge_branch->GetEntry(index);
			} else { }
			pfcands_charge_isLoaded = true;
		}
		return pfcands_charge_;
	}
	vector<int> &pfcands_particleId()
	{
		if (not pfcands_particleId_isLoaded) {
			if (pfcands_particleId_branch != 0) {
				pfcands_particleId_branch->GetEntry(index);
			} else { }
			pfcands_particleId_isLoaded = true;
		}
		return pfcands_particleId_;
	}
	vector<int> &pfjets_chargedHadronMultiplicity()
	{
		if (not pfjets_chargedHadronMultiplicity_isLoaded) {
			if (pfjets_chargedHadronMultiplicity_branch != 0) {
				pfjets_chargedHadronMultiplicity_branch->GetEntry(index);
			} else { }
			pfjets_chargedHadronMultiplicity_isLoaded = true;
		}
		return pfjets_chargedHadronMultiplicity_;
	}
	vector<int> &pfjets_chargedMultiplicity()
	{
		if (not pfjets_chargedMultiplicity_isLoaded) {
			if (pfjets_chargedMultiplicity_branch != 0) {
				pfjets_chargedMultiplicity_branch->GetEntry(index);
			} else { }
			pfjets_chargedMultiplicity_isLoaded = true;
		}
		return pfjets_chargedMultiplicity_;
	}
	vector<int> &pfjets_electronMultiplicity()
	{
		if (not pfjets_electronMultiplicity_isLoaded) {
			if (pfjets_electronMultiplicity_branch != 0) {
				pfjets_electronMultiplicity_branch->GetEntry(index);
			} else { }
			pfjets_electronMultiplicity_isLoaded = true;
		}
		return pfjets_electronMultiplicity_;
	}
	vector<int> &pfjets_hfEmMultiplicity()
	{
		if (not pfjets_hfEmMultiplicity_isLoaded) {
			if (pfjets_hfEmMultiplicity_branch != 0) {
				pfjets_hfEmMultiplicity_branch->GetEntry(index);
			} else { }
			pfjets_hfEmMultiplicity_isLoaded = true;
		}
		return pfjets_hfEmMultiplicity_;
	}
	vector<int> &pfjets_hfHadronMultiplicity()
	{
		if (not pfjets_hfHadronMultiplicity_isLoaded) {
			if (pfjets_hfHadronMultiplicity_branch != 0) {
				pfjets_hfHadronMultiplicity_branch->GetEntry(index);
			} else { }
			pfjets_hfHadronMultiplicity_isLoaded = true;
		}
		return pfjets_hfHadronMultiplicity_;
	}
	vector<int> &pfjets_muonMultiplicity()
	{
		if (not pfjets_muonMultiplicity_isLoaded) {
			if (pfjets_muonMultiplicity_branch != 0) {
				pfjets_muonMultiplicity_branch->GetEntry(index);
			} else { }
			pfjets_muonMultiplicity_isLoaded = true;
		}
		return pfjets_muonMultiplicity_;
	}
	vector<int> &pfjets_neutralHadronMultiplicity()
	{
		if (not pfjets_neutralHadronMultiplicity_isLoaded) {
			if (pfjets_neutralHadronMultiplicity_branch != 0) {
				pfjets_neutralHadronMultiplicity_branch->GetEntry(index);
			} else { }
			pfjets_neutralHadronMultiplicity_isLoaded = true;
		}
		return pfjets_neutralHadronMultiplicity_;
	}
	vector<int> &pfjets_neutralMultiplicity()
	{
		if (not pfjets_neutralMultiplicity_isLoaded) {
			if (pfjets_neutralMultiplicity_branch != 0) {
				pfjets_neutralMultiplicity_branch->GetEntry(index);
			} else { }
			pfjets_neutralMultiplicity_isLoaded = true;
		}
		return pfjets_neutralMultiplicity_;
	}
	vector<int> &pfjets_partonFlavour()
	{
		if (not pfjets_partonFlavour_isLoaded) {
			if (pfjets_partonFlavour_branch != 0) {
				pfjets_partonFlavour_branch->GetEntry(index);
			} else { }
			pfjets_partonFlavour_isLoaded = true;
		}
		return pfjets_partonFlavour_;
	}
	vector<int> &pfjets_photonMultiplicity()
	{
		if (not pfjets_photonMultiplicity_isLoaded) {
			if (pfjets_photonMultiplicity_branch != 0) {
				pfjets_photonMultiplicity_branch->GetEntry(index);
			} else { }
			pfjets_photonMultiplicity_isLoaded = true;
		}
		return pfjets_photonMultiplicity_;
	}
	vector<int> &taus_pf_charge()
	{
		if (not taus_pf_charge_isLoaded) {
			if (taus_pf_charge_branch != 0) {
				taus_pf_charge_branch->GetEntry(index);
			} else { }
			taus_pf_charge_isLoaded = true;
		}
		return taus_pf_charge_;
	}
	vector<int> &photons_fiduciality()
	{
		if (not photons_fiduciality_isLoaded) {
			if (photons_fiduciality_branch != 0) {
				photons_fiduciality_branch->GetEntry(index);
			} else { }
			photons_fiduciality_isLoaded = true;
		}
		return photons_fiduciality_;
	}
	vector<int> &puInfo_bunchCrossing()
	{
		if (not puInfo_bunchCrossing_isLoaded) {
			if (puInfo_bunchCrossing_branch != 0) {
				puInfo_bunchCrossing_branch->GetEntry(index);
			} else { }
			puInfo_bunchCrossing_isLoaded = true;
		}
		return puInfo_bunchCrossing_;
	}
	vector<int> &puInfo_nPUvertices()
	{
		if (not puInfo_nPUvertices_isLoaded) {
			if (puInfo_nPUvertices_branch != 0) {
				puInfo_nPUvertices_branch->GetEntry(index);
			} else { }
			puInfo_nPUvertices_isLoaded = true;
		}
		return puInfo_nPUvertices_;
	}
	vector<int> &convs_algo()
	{
		if (not convs_algo_isLoaded) {
			if (convs_algo_branch != 0) {
				convs_algo_branch->GetEntry(index);
			} else { }
			convs_algo_isLoaded = true;
		}
		return convs_algo_;
	}
	vector<int> &convs_isConverted()
	{
		if (not convs_isConverted_isLoaded) {
			if (convs_isConverted_branch != 0) {
				convs_isConverted_branch->GetEntry(index);
			} else { }
			convs_isConverted_isLoaded = true;
		}
		return convs_isConverted_;
	}
	vector<int> &convs_quality()
	{
		if (not convs_quality_isLoaded) {
			if (convs_quality_branch != 0) {
				convs_quality_branch->GetEntry(index);
			} else { }
			convs_quality_isLoaded = true;
		}
		return convs_quality_;
	}
	vector<int> &vtxs_isFake()
	{
		if (not vtxs_isFake_isLoaded) {
			if (vtxs_isFake_branch != 0) {
				vtxs_isFake_branch->GetEntry(index);
			} else { }
			vtxs_isFake_isLoaded = true;
		}
		return vtxs_isFake_;
	}
	vector<int> &vtxs_isValid()
	{
		if (not vtxs_isValid_isLoaded) {
			if (vtxs_isValid_branch != 0) {
				vtxs_isValid_branch->GetEntry(index);
			} else { }
			vtxs_isValid_isLoaded = true;
		}
		return vtxs_isValid_;
	}
	vector<int> &vtxs_tracksSize()
	{
		if (not vtxs_tracksSize_isLoaded) {
			if (vtxs_tracksSize_branch != 0) {
				vtxs_tracksSize_branch->GetEntry(index);
			} else { }
			vtxs_tracksSize_isLoaded = true;
		}
		return vtxs_tracksSize_;
	}
	vector<vector<int> > &els_convs_delMissHits()
	{
		if (not els_convs_delMissHits_isLoaded) {
			if (els_convs_delMissHits_branch != 0) {
				els_convs_delMissHits_branch->GetEntry(index);
			} else { }
			els_convs_delMissHits_isLoaded = true;
		}
		return els_convs_delMissHits_;
	}
	vector<vector<int> > &els_convs_flag()
	{
		if (not els_convs_flag_isLoaded) {
			if (els_convs_flag_branch != 0) {
				els_convs_flag_branch->GetEntry(index);
			} else { }
			els_convs_flag_isLoaded = true;
		}
		return els_convs_flag_;
	}
	vector<vector<int> > &els_convs_gsftkidx()
	{
		if (not els_convs_gsftkidx_isLoaded) {
			if (els_convs_gsftkidx_branch != 0) {
				els_convs_gsftkidx_branch->GetEntry(index);
			} else { }
			els_convs_gsftkidx_isLoaded = true;
		}
		return els_convs_gsftkidx_;
	}
	vector<vector<int> > &els_convs_tkidx()
	{
		if (not els_convs_tkidx_isLoaded) {
			if (els_convs_tkidx_branch != 0) {
				els_convs_tkidx_branch->GetEntry(index);
			} else { }
			els_convs_tkidx_isLoaded = true;
		}
		return els_convs_tkidx_;
	}
	vector<vector<int> > &els_PFCand_idx()
	{
		if (not els_PFCand_idx_isLoaded) {
			if (els_PFCand_idx_branch != 0) {
				els_PFCand_idx_branch->GetEntry(index);
			} else { }
			els_PFCand_idx_isLoaded = true;
		}
		return els_PFCand_idx_;
	}
	vector<vector<int> > &genps_lepdaughter_id()
	{
		if (not genps_lepdaughter_id_isLoaded) {
			if (genps_lepdaughter_id_branch != 0) {
				genps_lepdaughter_id_branch->GetEntry(index);
			} else { }
			genps_lepdaughter_id_isLoaded = true;
		}
		return genps_lepdaughter_id_;
	}
	vector<vector<int> > &genps_lepdaughter_idx()
	{
		if (not genps_lepdaughter_idx_isLoaded) {
			if (genps_lepdaughter_idx_branch != 0) {
				genps_lepdaughter_idx_branch->GetEntry(index);
			} else { }
			genps_lepdaughter_idx_isLoaded = true;
		}
		return genps_lepdaughter_idx_;
	}
	vector<vector<int> > &hlt_trigObjs_id()
	{
		if (not hlt_trigObjs_id_isLoaded) {
			if (hlt_trigObjs_id_branch != 0) {
				hlt_trigObjs_id_branch->GetEntry(index);
			} else { }
			hlt_trigObjs_id_isLoaded = true;
		}
		return hlt_trigObjs_id_;
	}
	vector<vector<int> > &mus_nStationCorrelatedHits()
	{
		if (not mus_nStationCorrelatedHits_isLoaded) {
			if (mus_nStationCorrelatedHits_branch != 0) {
				mus_nStationCorrelatedHits_branch->GetEntry(index);
			} else { }
			mus_nStationCorrelatedHits_isLoaded = true;
		}
		return mus_nStationCorrelatedHits_;
	}
	vector<vector<int> > &mus_nStationHits()
	{
		if (not mus_nStationHits_isLoaded) {
			if (mus_nStationHits_branch != 0) {
				mus_nStationHits_branch->GetEntry(index);
			} else { }
			mus_nStationHits_isLoaded = true;
		}
		return mus_nStationHits_;
	}
	vector<vector<int> > &pfjets_pfcandIndicies()
	{
		if (not pfjets_pfcandIndicies_isLoaded) {
			if (pfjets_pfcandIndicies_branch != 0) {
				pfjets_pfcandIndicies_branch->GetEntry(index);
			} else { }
			pfjets_pfcandIndicies_isLoaded = true;
		}
		return pfjets_pfcandIndicies_;
	}
	vector<vector<int> > &convs_nHitsBeforeVtx()
	{
		if (not convs_nHitsBeforeVtx_isLoaded) {
			if (convs_nHitsBeforeVtx_branch != 0) {
				convs_nHitsBeforeVtx_branch->GetEntry(index);
			} else { }
			convs_nHitsBeforeVtx_isLoaded = true;
		}
		return convs_nHitsBeforeVtx_;
	}
	vector<vector<int> > &convs_tkalgo()
	{
		if (not convs_tkalgo_isLoaded) {
			if (convs_tkalgo_branch != 0) {
				convs_tkalgo_branch->GetEntry(index);
			} else { }
			convs_tkalgo_isLoaded = true;
		}
		return convs_tkalgo_;
	}
	vector<vector<int> > &convs_tkidx()
	{
		if (not convs_tkidx_isLoaded) {
			if (convs_tkidx_branch != 0) {
				convs_tkidx_branch->GetEntry(index);
			} else { }
			convs_tkidx_isLoaded = true;
		}
		return convs_tkidx_;
	}
	unsigned int &els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version()
	{
		if (not els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_isLoaded) {
			if (els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch != 0) {
				els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_branch->GetEntry(index);
			} else { }
			els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_isLoaded = true;
		}
		return els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version_;
	}
	unsigned int &els_HLT_Ele17_Ele8_LeadingLeg_version()
	{
		if (not els_HLT_Ele17_Ele8_LeadingLeg_version_isLoaded) {
			if (els_HLT_Ele17_Ele8_LeadingLeg_version_branch != 0) {
				els_HLT_Ele17_Ele8_LeadingLeg_version_branch->GetEntry(index);
			} else { }
			els_HLT_Ele17_Ele8_LeadingLeg_version_isLoaded = true;
		}
		return els_HLT_Ele17_Ele8_LeadingLeg_version_;
	}
	unsigned int &els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version()
	{
		if (not els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_isLoaded) {
			if (els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch != 0) {
				els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_branch->GetEntry(index);
			} else { }
			els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_isLoaded = true;
		}
		return els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version_;
	}
	unsigned int &els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version()
	{
		if (not els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_isLoaded) {
			if (els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch != 0) {
				els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_branch->GetEntry(index);
			} else { }
			els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_isLoaded = true;
		}
		return els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version_;
	}
	unsigned int &els_HLT_Ele17_Ele8_TrailingLeg_version()
	{
		if (not els_HLT_Ele17_Ele8_TrailingLeg_version_isLoaded) {
			if (els_HLT_Ele17_Ele8_TrailingLeg_version_branch != 0) {
				els_HLT_Ele17_Ele8_TrailingLeg_version_branch->GetEntry(index);
			} else { }
			els_HLT_Ele17_Ele8_TrailingLeg_version_isLoaded = true;
		}
		return els_HLT_Ele17_Ele8_TrailingLeg_version_;
	}
	unsigned int &els_HLT_Ele17_Ele8_version()
	{
		if (not els_HLT_Ele17_Ele8_version_isLoaded) {
			if (els_HLT_Ele17_Ele8_version_branch != 0) {
				els_HLT_Ele17_Ele8_version_branch->GetEntry(index);
			} else { }
			els_HLT_Ele17_Ele8_version_isLoaded = true;
		}
		return els_HLT_Ele17_Ele8_version_;
	}
	unsigned int &els_HLT_Ele20_SC4_Mass50_LeadingLeg_version()
	{
		if (not els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_isLoaded) {
			if (els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch != 0) {
				els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_branch->GetEntry(index);
			} else { }
			els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_isLoaded = true;
		}
		return els_HLT_Ele20_SC4_Mass50_LeadingLeg_version_;
	}
	unsigned int &els_HLT_Ele20_SC4_Mass50_TrailingLeg_version()
	{
		if (not els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_isLoaded) {
			if (els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch != 0) {
				els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_branch->GetEntry(index);
			} else { }
			els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_isLoaded = true;
		}
		return els_HLT_Ele20_SC4_Mass50_TrailingLeg_version_;
	}
	unsigned int &els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version()
	{
		if (not els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_isLoaded) {
			if (els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch != 0) {
				els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_branch->GetEntry(index);
			} else { }
			els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_isLoaded = true;
		}
		return els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version_;
	}
	unsigned int &els_HLT_Ele27_WP80_version()
	{
		if (not els_HLT_Ele27_WP80_version_isLoaded) {
			if (els_HLT_Ele27_WP80_version_branch != 0) {
				els_HLT_Ele27_WP80_version_branch->GetEntry(index);
			} else { }
			els_HLT_Ele27_WP80_version_isLoaded = true;
		}
		return els_HLT_Ele27_WP80_version_;
	}
	unsigned int &els_HLT_Ele32_SC17_Mass50_LeadingLeg_version()
	{
		if (not els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_isLoaded) {
			if (els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch != 0) {
				els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_branch->GetEntry(index);
			} else { }
			els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_isLoaded = true;
		}
		return els_HLT_Ele32_SC17_Mass50_LeadingLeg_version_;
	}
	unsigned int &els_HLT_Ele32_SC17_Mass50_TrailingLeg_version()
	{
		if (not els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_isLoaded) {
			if (els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch != 0) {
				els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_branch->GetEntry(index);
			} else { }
			els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_isLoaded = true;
		}
		return els_HLT_Ele32_SC17_Mass50_TrailingLeg_version_;
	}
	unsigned int &els_HLT_Mu17_Ele8_TrailingLeg_version()
	{
		if (not els_HLT_Mu17_Ele8_TrailingLeg_version_isLoaded) {
			if (els_HLT_Mu17_Ele8_TrailingLeg_version_branch != 0) {
				els_HLT_Mu17_Ele8_TrailingLeg_version_branch->GetEntry(index);
			} else { }
			els_HLT_Mu17_Ele8_TrailingLeg_version_isLoaded = true;
		}
		return els_HLT_Mu17_Ele8_TrailingLeg_version_;
	}
	unsigned int &els_HLT_Mu17_Ele8_version()
	{
		if (not els_HLT_Mu17_Ele8_version_isLoaded) {
			if (els_HLT_Mu17_Ele8_version_branch != 0) {
				els_HLT_Mu17_Ele8_version_branch->GetEntry(index);
			} else { }
			els_HLT_Mu17_Ele8_version_isLoaded = true;
		}
		return els_HLT_Mu17_Ele8_version_;
	}
	unsigned int &els_HLT_Mu8_Ele17_version()
	{
		if (not els_HLT_Mu8_Ele17_version_isLoaded) {
			if (els_HLT_Mu8_Ele17_version_branch != 0) {
				els_HLT_Mu8_Ele17_version_branch->GetEntry(index);
			} else { }
			els_HLT_Mu8_Ele17_version_isLoaded = true;
		}
		return els_HLT_Mu8_Ele17_version_;
	}
	unsigned int &evt_nels()
	{
		if (not evt_nels_isLoaded) {
			if (evt_nels_branch != 0) {
				evt_nels_branch->GetEntry(index);
			} else { }
			evt_nels_isLoaded = true;
		}
		return evt_nels_;
	}
	unsigned int &evt_detectorStatus()
	{
		if (not evt_detectorStatus_isLoaded) {
			if (evt_detectorStatus_branch != 0) {
				evt_detectorStatus_branch->GetEntry(index);
			} else { }
			evt_detectorStatus_isLoaded = true;
		}
		return evt_detectorStatus_;
	}
	unsigned int &evt_event()
	{
		if (not evt_event_isLoaded) {
			if (evt_event_branch != 0) {
				evt_event_branch->GetEntry(index);
			} else { }
			evt_event_isLoaded = true;
		}
		return evt_event_;
	}
	unsigned int &evt_lumiBlock()
	{
		if (not evt_lumiBlock_isLoaded) {
			if (evt_lumiBlock_branch != 0) {
				evt_lumiBlock_branch->GetEntry(index);
			} else { }
			evt_lumiBlock_isLoaded = true;
		}
		return evt_lumiBlock_;
	}
	unsigned int &evt_run()
	{
		if (not evt_run_isLoaded) {
			if (evt_run_branch != 0) {
				evt_run_branch->GetEntry(index);
			} else { }
			evt_run_isLoaded = true;
		}
		return evt_run_;
	}
	unsigned int &evt_ngenjetsNoMuNoNu()
	{
		if (not evt_ngenjetsNoMuNoNu_isLoaded) {
			if (evt_ngenjetsNoMuNoNu_branch != 0) {
				evt_ngenjetsNoMuNoNu_branch->GetEntry(index);
			} else { }
			evt_ngenjetsNoMuNoNu_isLoaded = true;
		}
		return evt_ngenjetsNoMuNoNu_;
	}
	unsigned int &genps_signalProcessID()
	{
		if (not genps_signalProcessID_isLoaded) {
			if (genps_signalProcessID_branch != 0) {
				genps_signalProcessID_branch->GetEntry(index);
			} else { }
			genps_signalProcessID_isLoaded = true;
		}
		return genps_signalProcessID_;
	}
	unsigned int &mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version()
	{
		if (not mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_isLoaded) {
			if (mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch != 0) {
				mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_branch->GetEntry(index);
			} else { }
			mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_isLoaded = true;
		}
		return mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version_;
	}
	unsigned int &mus_HLT_IsoMu24_eta2p1_version()
	{
		if (not mus_HLT_IsoMu24_eta2p1_version_isLoaded) {
			if (mus_HLT_IsoMu24_eta2p1_version_branch != 0) {
				mus_HLT_IsoMu24_eta2p1_version_branch->GetEntry(index);
			} else { }
			mus_HLT_IsoMu24_eta2p1_version_isLoaded = true;
		}
		return mus_HLT_IsoMu24_eta2p1_version_;
	}
	unsigned int &mus_HLT_Mu17_Ele8_LeadingLeg_version()
	{
		if (not mus_HLT_Mu17_Ele8_LeadingLeg_version_isLoaded) {
			if (mus_HLT_Mu17_Ele8_LeadingLeg_version_branch != 0) {
				mus_HLT_Mu17_Ele8_LeadingLeg_version_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_Ele8_LeadingLeg_version_isLoaded = true;
		}
		return mus_HLT_Mu17_Ele8_LeadingLeg_version_;
	}
	unsigned int &mus_HLT_Mu17_Ele8_version()
	{
		if (not mus_HLT_Mu17_Ele8_version_isLoaded) {
			if (mus_HLT_Mu17_Ele8_version_branch != 0) {
				mus_HLT_Mu17_Ele8_version_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_Ele8_version_isLoaded = true;
		}
		return mus_HLT_Mu17_Ele8_version_;
	}
	unsigned int &mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version()
	{
		if (not mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_isLoaded) {
			if (mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch != 0) {
				mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_isLoaded = true;
		}
		return mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version_;
	}
	unsigned int &mus_HLT_Mu17_Mu8_LeadingLeg_version()
	{
		if (not mus_HLT_Mu17_Mu8_LeadingLeg_version_isLoaded) {
			if (mus_HLT_Mu17_Mu8_LeadingLeg_version_branch != 0) {
				mus_HLT_Mu17_Mu8_LeadingLeg_version_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_Mu8_LeadingLeg_version_isLoaded = true;
		}
		return mus_HLT_Mu17_Mu8_LeadingLeg_version_;
	}
	unsigned int &mus_HLT_Mu17_Mu8_TrailingLeg_version()
	{
		if (not mus_HLT_Mu17_Mu8_TrailingLeg_version_isLoaded) {
			if (mus_HLT_Mu17_Mu8_TrailingLeg_version_branch != 0) {
				mus_HLT_Mu17_Mu8_TrailingLeg_version_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_Mu8_TrailingLeg_version_isLoaded = true;
		}
		return mus_HLT_Mu17_Mu8_TrailingLeg_version_;
	}
	unsigned int &mus_HLT_Mu17_Mu8_version()
	{
		if (not mus_HLT_Mu17_Mu8_version_isLoaded) {
			if (mus_HLT_Mu17_Mu8_version_branch != 0) {
				mus_HLT_Mu17_Mu8_version_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_Mu8_version_isLoaded = true;
		}
		return mus_HLT_Mu17_Mu8_version_;
	}
	unsigned int &mus_HLT_Mu17_TkMu8_LeadingLeg_version()
	{
		if (not mus_HLT_Mu17_TkMu8_LeadingLeg_version_isLoaded) {
			if (mus_HLT_Mu17_TkMu8_LeadingLeg_version_branch != 0) {
				mus_HLT_Mu17_TkMu8_LeadingLeg_version_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_TkMu8_LeadingLeg_version_isLoaded = true;
		}
		return mus_HLT_Mu17_TkMu8_LeadingLeg_version_;
	}
	unsigned int &mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version()
	{
		if (not mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_isLoaded) {
			if (mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch != 0) {
				mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_isLoaded = true;
		}
		return mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version_;
	}
	unsigned int &mus_HLT_Mu17_TkMu8_TrailingLeg_version()
	{
		if (not mus_HLT_Mu17_TkMu8_TrailingLeg_version_isLoaded) {
			if (mus_HLT_Mu17_TkMu8_TrailingLeg_version_branch != 0) {
				mus_HLT_Mu17_TkMu8_TrailingLeg_version_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_TkMu8_TrailingLeg_version_isLoaded = true;
		}
		return mus_HLT_Mu17_TkMu8_TrailingLeg_version_;
	}
	unsigned int &mus_HLT_Mu17_TkMu8_version()
	{
		if (not mus_HLT_Mu17_TkMu8_version_isLoaded) {
			if (mus_HLT_Mu17_TkMu8_version_branch != 0) {
				mus_HLT_Mu17_TkMu8_version_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_TkMu8_version_isLoaded = true;
		}
		return mus_HLT_Mu17_TkMu8_version_;
	}
	unsigned int &mus_HLT_Mu8_Ele17_TrailingLeg_version()
	{
		if (not mus_HLT_Mu8_Ele17_TrailingLeg_version_isLoaded) {
			if (mus_HLT_Mu8_Ele17_TrailingLeg_version_branch != 0) {
				mus_HLT_Mu8_Ele17_TrailingLeg_version_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu8_Ele17_TrailingLeg_version_isLoaded = true;
		}
		return mus_HLT_Mu8_Ele17_TrailingLeg_version_;
	}
	unsigned int &mus_HLT_Mu8_Ele17_version()
	{
		if (not mus_HLT_Mu8_Ele17_version_isLoaded) {
			if (mus_HLT_Mu8_Ele17_version_branch != 0) {
				mus_HLT_Mu8_Ele17_version_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu8_Ele17_version_isLoaded = true;
		}
		return mus_HLT_Mu8_Ele17_version_;
	}
	unsigned int &evt_nphotons()
	{
		if (not evt_nphotons_isLoaded) {
			if (evt_nphotons_branch != 0) {
				evt_nphotons_branch->GetEntry(index);
			} else { }
			evt_nphotons_isLoaded = true;
		}
		return evt_nphotons_;
	}
	unsigned int &evt_nvtxs()
	{
		if (not evt_nvtxs_isLoaded) {
			if (evt_nvtxs_branch != 0) {
				evt_nvtxs_branch->GetEntry(index);
			} else { }
			evt_nvtxs_isLoaded = true;
		}
		return evt_nvtxs_;
	}
	vector<unsigned int> &els_HLT_Ele17_Ele8()
	{
		if (not els_HLT_Ele17_Ele8_isLoaded) {
			if (els_HLT_Ele17_Ele8_branch != 0) {
				els_HLT_Ele17_Ele8_branch->GetEntry(index);
			} else { }
			els_HLT_Ele17_Ele8_isLoaded = true;
		}
		return els_HLT_Ele17_Ele8_;
	}
	vector<unsigned int> &els_HLT_Ele17_Ele8_L1sL1DoubleEG137()
	{
		if (not els_HLT_Ele17_Ele8_L1sL1DoubleEG137_isLoaded) {
			if (els_HLT_Ele17_Ele8_L1sL1DoubleEG137_branch != 0) {
				els_HLT_Ele17_Ele8_L1sL1DoubleEG137_branch->GetEntry(index);
			} else { }
			els_HLT_Ele17_Ele8_L1sL1DoubleEG137_isLoaded = true;
		}
		return els_HLT_Ele17_Ele8_L1sL1DoubleEG137_;
	}
	vector<unsigned int> &els_HLT_Ele17_Ele8_LeadingLeg()
	{
		if (not els_HLT_Ele17_Ele8_LeadingLeg_isLoaded) {
			if (els_HLT_Ele17_Ele8_LeadingLeg_branch != 0) {
				els_HLT_Ele17_Ele8_LeadingLeg_branch->GetEntry(index);
			} else { }
			els_HLT_Ele17_Ele8_LeadingLeg_isLoaded = true;
		}
		return els_HLT_Ele17_Ele8_LeadingLeg_;
	}
	vector<unsigned int> &els_HLT_Ele17_Ele8_Mass50_LeadingLeg()
	{
		if (not els_HLT_Ele17_Ele8_Mass50_LeadingLeg_isLoaded) {
			if (els_HLT_Ele17_Ele8_Mass50_LeadingLeg_branch != 0) {
				els_HLT_Ele17_Ele8_Mass50_LeadingLeg_branch->GetEntry(index);
			} else { }
			els_HLT_Ele17_Ele8_Mass50_LeadingLeg_isLoaded = true;
		}
		return els_HLT_Ele17_Ele8_Mass50_LeadingLeg_;
	}
	vector<unsigned int> &els_HLT_Ele17_Ele8_Mass50_TrailingLeg()
	{
		if (not els_HLT_Ele17_Ele8_Mass50_TrailingLeg_isLoaded) {
			if (els_HLT_Ele17_Ele8_Mass50_TrailingLeg_branch != 0) {
				els_HLT_Ele17_Ele8_Mass50_TrailingLeg_branch->GetEntry(index);
			} else { }
			els_HLT_Ele17_Ele8_Mass50_TrailingLeg_isLoaded = true;
		}
		return els_HLT_Ele17_Ele8_Mass50_TrailingLeg_;
	}
	vector<unsigned int> &els_HLT_Ele17_Ele8_TrailingLeg()
	{
		if (not els_HLT_Ele17_Ele8_TrailingLeg_isLoaded) {
			if (els_HLT_Ele17_Ele8_TrailingLeg_branch != 0) {
				els_HLT_Ele17_Ele8_TrailingLeg_branch->GetEntry(index);
			} else { }
			els_HLT_Ele17_Ele8_TrailingLeg_isLoaded = true;
		}
		return els_HLT_Ele17_Ele8_TrailingLeg_;
	}
	vector<unsigned int> &els_HLT_Ele20_SC4_Mass50_LeadingLeg()
	{
		if (not els_HLT_Ele20_SC4_Mass50_LeadingLeg_isLoaded) {
			if (els_HLT_Ele20_SC4_Mass50_LeadingLeg_branch != 0) {
				els_HLT_Ele20_SC4_Mass50_LeadingLeg_branch->GetEntry(index);
			} else { }
			els_HLT_Ele20_SC4_Mass50_LeadingLeg_isLoaded = true;
		}
		return els_HLT_Ele20_SC4_Mass50_LeadingLeg_;
	}
	vector<unsigned int> &els_HLT_Ele20_SC4_Mass50_TrailingLeg()
	{
		if (not els_HLT_Ele20_SC4_Mass50_TrailingLeg_isLoaded) {
			if (els_HLT_Ele20_SC4_Mass50_TrailingLeg_branch != 0) {
				els_HLT_Ele20_SC4_Mass50_TrailingLeg_branch->GetEntry(index);
			} else { }
			els_HLT_Ele20_SC4_Mass50_TrailingLeg_isLoaded = true;
		}
		return els_HLT_Ele20_SC4_Mass50_TrailingLeg_;
	}
	vector<unsigned int> &els_HLT_Ele27_WP80()
	{
		if (not els_HLT_Ele27_WP80_isLoaded) {
			if (els_HLT_Ele27_WP80_branch != 0) {
				els_HLT_Ele27_WP80_branch->GetEntry(index);
			} else { }
			els_HLT_Ele27_WP80_isLoaded = true;
		}
		return els_HLT_Ele27_WP80_;
	}
	vector<unsigned int> &els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22()
	{
		if (not els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_isLoaded) {
			if (els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_branch != 0) {
				els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_branch->GetEntry(index);
			} else { }
			els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_isLoaded = true;
		}
		return els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_;
	}
	vector<unsigned int> &els_HLT_Ele32_SC17_Mass50_LeadingLeg()
	{
		if (not els_HLT_Ele32_SC17_Mass50_LeadingLeg_isLoaded) {
			if (els_HLT_Ele32_SC17_Mass50_LeadingLeg_branch != 0) {
				els_HLT_Ele32_SC17_Mass50_LeadingLeg_branch->GetEntry(index);
			} else { }
			els_HLT_Ele32_SC17_Mass50_LeadingLeg_isLoaded = true;
		}
		return els_HLT_Ele32_SC17_Mass50_LeadingLeg_;
	}
	vector<unsigned int> &els_HLT_Ele32_SC17_Mass50_TrailingLeg()
	{
		if (not els_HLT_Ele32_SC17_Mass50_TrailingLeg_isLoaded) {
			if (els_HLT_Ele32_SC17_Mass50_TrailingLeg_branch != 0) {
				els_HLT_Ele32_SC17_Mass50_TrailingLeg_branch->GetEntry(index);
			} else { }
			els_HLT_Ele32_SC17_Mass50_TrailingLeg_isLoaded = true;
		}
		return els_HLT_Ele32_SC17_Mass50_TrailingLeg_;
	}
	vector<unsigned int> &els_HLT_Mu17_Ele8()
	{
		if (not els_HLT_Mu17_Ele8_isLoaded) {
			if (els_HLT_Mu17_Ele8_branch != 0) {
				els_HLT_Mu17_Ele8_branch->GetEntry(index);
			} else { }
			els_HLT_Mu17_Ele8_isLoaded = true;
		}
		return els_HLT_Mu17_Ele8_;
	}
	vector<unsigned int> &els_HLT_Mu17_Ele8_TrailingLeg()
	{
		if (not els_HLT_Mu17_Ele8_TrailingLeg_isLoaded) {
			if (els_HLT_Mu17_Ele8_TrailingLeg_branch != 0) {
				els_HLT_Mu17_Ele8_TrailingLeg_branch->GetEntry(index);
			} else { }
			els_HLT_Mu17_Ele8_TrailingLeg_isLoaded = true;
		}
		return els_HLT_Mu17_Ele8_TrailingLeg_;
	}
	vector<unsigned int> &els_HLT_Mu8_Ele17()
	{
		if (not els_HLT_Mu8_Ele17_isLoaded) {
			if (els_HLT_Mu8_Ele17_branch != 0) {
				els_HLT_Mu8_Ele17_branch->GetEntry(index);
			} else { }
			els_HLT_Mu8_Ele17_isLoaded = true;
		}
		return els_HLT_Mu8_Ele17_;
	}
	vector<unsigned int> &els_id2012ext_loose()
	{
		if (not els_id2012ext_loose_isLoaded) {
			if (els_id2012ext_loose_branch != 0) {
				els_id2012ext_loose_branch->GetEntry(index);
			} else { }
			els_id2012ext_loose_isLoaded = true;
		}
		return els_id2012ext_loose_;
	}
	vector<unsigned int> &els_id2012ext_medium()
	{
		if (not els_id2012ext_medium_isLoaded) {
			if (els_id2012ext_medium_branch != 0) {
				els_id2012ext_medium_branch->GetEntry(index);
			} else { }
			els_id2012ext_medium_isLoaded = true;
		}
		return els_id2012ext_medium_;
	}
	vector<unsigned int> &els_id2012ext_tight()
	{
		if (not els_id2012ext_tight_isLoaded) {
			if (els_id2012ext_tight_branch != 0) {
				els_id2012ext_tight_branch->GetEntry(index);
			} else { }
			els_id2012ext_tight_isLoaded = true;
		}
		return els_id2012ext_tight_;
	}
	vector<unsigned int> &els_id2012ext_veto()
	{
		if (not els_id2012ext_veto_isLoaded) {
			if (els_id2012ext_veto_branch != 0) {
				els_id2012ext_veto_branch->GetEntry(index);
			} else { }
			els_id2012ext_veto_isLoaded = true;
		}
		return els_id2012ext_veto_;
	}
	vector<unsigned int> &els_id2012_loose()
	{
		if (not els_id2012_loose_isLoaded) {
			if (els_id2012_loose_branch != 0) {
				els_id2012_loose_branch->GetEntry(index);
			} else { }
			els_id2012_loose_isLoaded = true;
		}
		return els_id2012_loose_;
	}
	vector<unsigned int> &els_id2012_medium()
	{
		if (not els_id2012_medium_isLoaded) {
			if (els_id2012_medium_branch != 0) {
				els_id2012_medium_branch->GetEntry(index);
			} else { }
			els_id2012_medium_isLoaded = true;
		}
		return els_id2012_medium_;
	}
	vector<unsigned int> &els_id2012_tight()
	{
		if (not els_id2012_tight_isLoaded) {
			if (els_id2012_tight_branch != 0) {
				els_id2012_tight_branch->GetEntry(index);
			} else { }
			els_id2012_tight_isLoaded = true;
		}
		return els_id2012_tight_;
	}
	vector<unsigned int> &els_id2012_veto()
	{
		if (not els_id2012_veto_isLoaded) {
			if (els_id2012_veto_branch != 0) {
				els_id2012_veto_branch->GetEntry(index);
			} else { }
			els_id2012_veto_isLoaded = true;
		}
		return els_id2012_veto_;
	}
	vector<unsigned int> &hlt_prescales()
	{
		if (not hlt_prescales_isLoaded) {
			if (hlt_prescales_branch != 0) {
				hlt_prescales_branch->GetEntry(index);
			} else { }
			hlt_prescales_isLoaded = true;
		}
		return hlt_prescales_;
	}
	vector<unsigned int> &mus_HLT_IsoMu24_eta2p1()
	{
		if (not mus_HLT_IsoMu24_eta2p1_isLoaded) {
			if (mus_HLT_IsoMu24_eta2p1_branch != 0) {
				mus_HLT_IsoMu24_eta2p1_branch->GetEntry(index);
			} else { }
			mus_HLT_IsoMu24_eta2p1_isLoaded = true;
		}
		return mus_HLT_IsoMu24_eta2p1_;
	}
	vector<unsigned int> &mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1()
	{
		if (not mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_isLoaded) {
			if (mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_branch != 0) {
				mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_branch->GetEntry(index);
			} else { }
			mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_isLoaded = true;
		}
		return mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_;
	}
	vector<unsigned int> &mus_HLT_Mu17_Ele8()
	{
		if (not mus_HLT_Mu17_Ele8_isLoaded) {
			if (mus_HLT_Mu17_Ele8_branch != 0) {
				mus_HLT_Mu17_Ele8_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_Ele8_isLoaded = true;
		}
		return mus_HLT_Mu17_Ele8_;
	}
	vector<unsigned int> &mus_HLT_Mu17_Ele8_LeadingLeg()
	{
		if (not mus_HLT_Mu17_Ele8_LeadingLeg_isLoaded) {
			if (mus_HLT_Mu17_Ele8_LeadingLeg_branch != 0) {
				mus_HLT_Mu17_Ele8_LeadingLeg_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_Ele8_LeadingLeg_isLoaded = true;
		}
		return mus_HLT_Mu17_Ele8_LeadingLeg_;
	}
	vector<unsigned int> &mus_HLT_Mu17_Mu8()
	{
		if (not mus_HLT_Mu17_Mu8_isLoaded) {
			if (mus_HLT_Mu17_Mu8_branch != 0) {
				mus_HLT_Mu17_Mu8_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_Mu8_isLoaded = true;
		}
		return mus_HLT_Mu17_Mu8_;
	}
	vector<unsigned int> &mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen()
	{
		if (not mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_isLoaded) {
			if (mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_branch != 0) {
				mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_isLoaded = true;
		}
		return mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_;
	}
	vector<unsigned int> &mus_HLT_Mu17_Mu8_LeadingLeg()
	{
		if (not mus_HLT_Mu17_Mu8_LeadingLeg_isLoaded) {
			if (mus_HLT_Mu17_Mu8_LeadingLeg_branch != 0) {
				mus_HLT_Mu17_Mu8_LeadingLeg_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_Mu8_LeadingLeg_isLoaded = true;
		}
		return mus_HLT_Mu17_Mu8_LeadingLeg_;
	}
	vector<unsigned int> &mus_HLT_Mu17_Mu8_TrailingLeg()
	{
		if (not mus_HLT_Mu17_Mu8_TrailingLeg_isLoaded) {
			if (mus_HLT_Mu17_Mu8_TrailingLeg_branch != 0) {
				mus_HLT_Mu17_Mu8_TrailingLeg_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_Mu8_TrailingLeg_isLoaded = true;
		}
		return mus_HLT_Mu17_Mu8_TrailingLeg_;
	}
	vector<unsigned int> &mus_HLT_Mu17_TkMu8()
	{
		if (not mus_HLT_Mu17_TkMu8_isLoaded) {
			if (mus_HLT_Mu17_TkMu8_branch != 0) {
				mus_HLT_Mu17_TkMu8_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_TkMu8_isLoaded = true;
		}
		return mus_HLT_Mu17_TkMu8_;
	}
	vector<unsigned int> &mus_HLT_Mu17_TkMu8_LeadingLeg()
	{
		if (not mus_HLT_Mu17_TkMu8_LeadingLeg_isLoaded) {
			if (mus_HLT_Mu17_TkMu8_LeadingLeg_branch != 0) {
				mus_HLT_Mu17_TkMu8_LeadingLeg_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_TkMu8_LeadingLeg_isLoaded = true;
		}
		return mus_HLT_Mu17_TkMu8_LeadingLeg_;
	}
	vector<unsigned int> &mus_HLT_Mu17_TkMu8_TrailingLeg()
	{
		if (not mus_HLT_Mu17_TkMu8_TrailingLeg_isLoaded) {
			if (mus_HLT_Mu17_TkMu8_TrailingLeg_branch != 0) {
				mus_HLT_Mu17_TkMu8_TrailingLeg_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_TkMu8_TrailingLeg_isLoaded = true;
		}
		return mus_HLT_Mu17_TkMu8_TrailingLeg_;
	}
	vector<unsigned int> &mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered()
	{
		if (not mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_isLoaded) {
			if (mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_branch != 0) {
				mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_isLoaded = true;
		}
		return mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_;
	}
	vector<unsigned int> &mus_HLT_Mu8_Ele17()
	{
		if (not mus_HLT_Mu8_Ele17_isLoaded) {
			if (mus_HLT_Mu8_Ele17_branch != 0) {
				mus_HLT_Mu8_Ele17_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu8_Ele17_isLoaded = true;
		}
		return mus_HLT_Mu8_Ele17_;
	}
	vector<unsigned int> &mus_HLT_Mu8_Ele17_TrailingLeg()
	{
		if (not mus_HLT_Mu8_Ele17_TrailingLeg_isLoaded) {
			if (mus_HLT_Mu8_Ele17_TrailingLeg_branch != 0) {
				mus_HLT_Mu8_Ele17_TrailingLeg_branch->GetEntry(index);
			} else { }
			mus_HLT_Mu8_Ele17_TrailingLeg_isLoaded = true;
		}
		return mus_HLT_Mu8_Ele17_TrailingLeg_;
	}
	bool passHLTTrigger(TString trigName) {
		int trigIndx;
		vector<TString>::const_iterator begin_it = hlt_trigNames().begin();
		vector<TString>::const_iterator end_it = hlt_trigNames().end();
		vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
		if(found_it != end_it)
			trigIndx = found_it - begin_it;
		else {
			cout << "Cannot find Trigger " << trigName << endl; 
			return 0;
		}

	return hlt_bits().TestBitNumber(trigIndx);
	}

  static void progress( int nEventsTotal, int nEventsChain ){
    int period = 1000;
    if(nEventsTotal%1000 == 0) {
      // xterm magic from L. Vacavant and A. Cerri
      if (isatty(1)) {
        if( ( nEventsChain - nEventsTotal ) > period ){
          float frac = (float)nEventsTotal/(nEventsChain*0.01);
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", frac);
          fflush(stdout);
        }
        else {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", 100.);
          cout << endl;
        }
      }
    }
  }
  
};

#ifndef __CINT__
extern CMS2 cms2;
#endif

namespace tas {
	TBits &hlt_bits();
	vector<TString> &evt_CMS2tag();
	vector<TString> &evt_dataset();
	vector<TString> &hlt_trigNames();
	bool &hcalnoise_HasBadRBXTS4TS5();
	bool &evt_cscTightHaloId();
	bool &evt_hbheFilter();
	bool &filt_cscBeamHalo();
	bool &filt_ecalLaser();
	bool &filt_ecalTP();
	bool &filt_eeBadSc();
	bool &filt_goodVertices();
	bool &filt_hbheNoise();
	bool &filt_hcalLaser();
	bool &filt_metfilter();
	bool &filt_trackingFailure();
	bool &filt_trkPOGFilters();
	bool &filt_trkPOG_logErrorTooManyClusters();
	bool &filt_trkPOG_manystripclus53X();
	bool &filt_trkPOG_toomanystripclus53X();
	vector<bool> &els_conv_vtx_flag();
	vector<bool> &els_isGsfCtfScPixChargeConsistent();
	vector<bool> &els_passingMvaPreselection();
	vector<bool> &els_passingPflowPreselection();
	vector<bool> &mus_isRPCMuon();
	vector<bool> &mus_tightMatch();
	vector<bool> &mus_updatedSta();
	vector<bool> &photons_haspixelSeed();
	float &evt_bs_Xwidth();
	float &evt_bs_XwidthErr();
	float &evt_bs_Ywidth();
	float &evt_bs_YwidthErr();
	float &evt_bs_dxdz();
	float &evt_bs_dxdzErr();
	float &evt_bs_dydz();
	float &evt_bs_dydzErr();
	float &evt_bs_sigmaZ();
	float &evt_bs_sigmaZErr();
	float &evt_bs_xErr();
	float &evt_bs_yErr();
	float &evt_bs_zErr();
	float &evt_bField();
	float &evt_fixgrid_all_rho();
	float &evt_fixgridfastjet_allcalo_rho();
	float &evt_fixgridfastjet_all_rho();
	float &evt_fixgridfastjet_centralcalo_rho();
	float &evt_fixgridfastjet_centralchargedpileup_rho();
	float &evt_fixgridfastjet_centralneutral_rho();
	float &evt_kfactor();
	float &evt_scale1fb();
	float &evt_xsec_excl();
	float &evt_xsec_incl();
	float &gen_met();
	float &gen_metPhi();
	float &genps_alphaQCD();
	float &genps_pthat();
	float &genps_qScale();
	float &genps_weight();
	float &gen_sumEt();
	float &hcalnoise_GetRecHitEnergy();
	float &hcalnoise_GetRecHitEnergy15();
	float &hcalnoise_GetTotalCalibCharge();
	float &hcalnoise_TS4TS5NoiseSumE();
	float &hcalnoise_TS4TS5NoiseSumEt();
	float &hcalnoise_eventChargeFraction();
	float &hcalnoise_eventEMEnergy();
	float &hcalnoise_eventEMFraction();
	float &hcalnoise_eventHadEnergy();
	float &hcalnoise_eventTrackEnergy();
	float &hcalnoise_flatNoiseSumE();
	float &hcalnoise_flatNoiseSumEt();
	float &hcalnoise_isolatedNoiseSumE();
	float &hcalnoise_isolatedNoiseSumEt();
	float &hcalnoise_max10GeVHitTime();
	float &hcalnoise_max25GeVHitTime();
	float &hcalnoise_maxE10TS();
	float &hcalnoise_maxE2Over10TS();
	float &hcalnoise_maxE2TS();
	float &hcalnoise_min10GeVHitTime();
	float &hcalnoise_min25GeVHitTime();
	float &hcalnoise_minE10TS();
	float &hcalnoise_minE2Over10TS();
	float &hcalnoise_minE2TS();
	float &hcalnoise_minHPDEMF();
	float &hcalnoise_minRBXEMF();
	float &hcalnoise_rms10GeVHitTime();
	float &hcalnoise_rms25GeVHitTime();
	float &hcalnoise_spikeNoiseSumE();
	float &hcalnoise_spikeNoiseSumEt();
	float &hcalnoise_triangleNoiseSumE();
	float &hcalnoise_triangleNoiseSumEt();
	float &pdfinfo_pdf1();
	float &pdfinfo_pdf2();
	float &pdfinfo_scale();
	float &pdfinfo_x1();
	float &pdfinfo_x2();
	float &evt_fixgrid_rho_all();
	float &evt_fixgrid_rho_ctr();
	float &evt_fixgrid_rho_fwd();
	float &evt_pfmet();
	float &evt_pfmetPhi();
	float &evt_pfmetSig();
	float &evt_pfmetSignificance();
	float &evt_pfsumet();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &evt_bsp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_motherp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_mc_motherp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_mc_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_gp_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_motherp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_mc_motherp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_mc_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_conv_pos_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_inner_position();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_patMatch_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_outer_position();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4In();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4Out();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_vertex_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_vertex_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genjets_p4NoMuNoNu();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_prod_vtx();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_vertex_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_mc_patMatch_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_pfp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_sta_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_sta_vertex_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_trk_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_vertex_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfcands_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &taus_pf_lead_chargecand_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &taus_pf_lead_neutrcand_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &taus_pf_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &vtxs_position();
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &genps_lepdaughter_p4();
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hlt_trigObjs_p4();
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &taus_pf_isocands_p4();
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &taus_pf_signalcands_p4();
	vector<float> &evt_bs_covMatrix();
	vector<float> &els_mc3dr();
	vector<float> &els_mcdr();
	vector<float> &mus_mc3dr();
	vector<float> &mus_mcdr();
	vector<float> &pfjets_mc3dr();
	vector<float> &pfjets_mcdr();
	vector<float> &pfjets_mc_emEnergy();
	vector<float> &pfjets_mc_gpdr();
	vector<float> &pfjets_mc_hadEnergy();
	vector<float> &pfjets_mc_invEnergy();
	vector<float> &pfjets_mc_otherEnergy();
	vector<float> &photons_mc3dr();
	vector<float> &photons_mcdr();
	vector<float> &els_chi2();
	vector<float> &els_ckf_chi2();
	vector<float> &els_ckf_ndof();
	vector<float> &els_conv_dcot();
	vector<float> &els_conv_dist();
	vector<float> &els_conv_old_dcot();
	vector<float> &els_conv_old_dist();
	vector<float> &els_conv_old_radius();
	vector<float> &els_conv_radius();
	vector<float> &els_d0();
	vector<float> &els_d0Err();
	vector<float> &els_d0corr();
	vector<float> &els_d0corrPhi();
	vector<float> &els_d0phiCov();
	vector<float> &els_dEtaIn();
	vector<float> &els_dEtaOut();
	vector<float> &els_dPhiIn();
	vector<float> &els_dPhiInPhiOut();
	vector<float> &els_dPhiOut();
	vector<float> &els_deltaEtaEleClusterTrackAtCalo();
	vector<float> &els_deltaPhiEleClusterTrackAtCalo();
	vector<float> &els_dxyPV();
	vector<float> &els_dzPV();
	vector<float> &els_e1x5();
	vector<float> &els_e1x5_full5x5();
	vector<float> &els_e2x5Max();
	vector<float> &els_e2x5Max_full5x5();
	vector<float> &els_e3x3();
	vector<float> &els_e5x5();
	vector<float> &els_e5x5_full5x5();
	vector<float> &els_eMax();
	vector<float> &els_eOverPIn();
	vector<float> &els_eOverPOut();
	vector<float> &els_eSC();
	vector<float> &els_eSCPresh();
	vector<float> &els_eSCRaw();
	vector<float> &els_eSeed();
	vector<float> &els_eSeedOverPIn();
	vector<float> &els_eSeedOverPOut();
	vector<float> &els_ecalEnergy();
	vector<float> &els_ecalEnergyError();
	vector<float> &els_ecalIso();
	vector<float> &els_ecalIso04();
	vector<float> &els_etaErr();
	vector<float> &els_etaSC();
	vector<float> &els_etaSCwidth();
	vector<float> &els_fbrem();
	vector<float> &els_hOverE();
	vector<float> &els_hOverEBC();
	vector<float> &els_hcalDepth1OverEcal();
	vector<float> &els_hcalDepth1TowerSumEt();
	vector<float> &els_hcalDepth1TowerSumEt04();
	vector<float> &els_hcalDepth2OverEcal();
	vector<float> &els_hcalDepth2TowerSumEt();
	vector<float> &els_hcalDepth2TowerSumEt04();
	vector<float> &els_hcalIso();
	vector<float> &els_hcalIso04();
	vector<float> &els_ip3d();
	vector<float> &els_ip3derr();
	vector<float> &els_iso03_pf();
	vector<float> &els_iso03_pf2012_ch();
	vector<float> &els_iso03_pf2012_em();
	vector<float> &els_iso03_pf2012ext_ch();
	vector<float> &els_iso03_pf2012ext_em();
	vector<float> &els_iso03_pf2012ext_nh();
	vector<float> &els_iso03_pf2012_nh();
	vector<float> &els_iso03_pf_ch();
	vector<float> &els_iso03_pf_gamma05();
	vector<float> &els_iso03_pf_nhad05();
	vector<float> &els_iso04_pf();
	vector<float> &els_iso04_pf2012_ch();
	vector<float> &els_iso04_pf2012_em();
	vector<float> &els_iso04_pf2012ext_ch();
	vector<float> &els_iso04_pf2012ext_em();
	vector<float> &els_iso04_pf2012ext_nh();
	vector<float> &els_iso04_pf2012_nh();
	vector<float> &els_iso04_pf_ch();
	vector<float> &els_iso04_pf_gamma05();
	vector<float> &els_iso04_pf_nhad05();
	vector<float> &els_layer1_charge();
	vector<float> &els_lh();
	vector<float> &els_mc_patMatch_dr();
	vector<float> &els_mva();
	vector<float> &els_ndof();
	vector<float> &els_pfChargedHadronIso();
	vector<float> &els_pfNeutralHadronIso();
	vector<float> &els_pfPUIso();
	vector<float> &els_pfPhotonIso();
	vector<float> &els_phiErr();
	vector<float> &els_phiSC();
	vector<float> &els_phiSCwidth();
	vector<float> &els_ptErr();
	vector<float> &els_ptErrGsf();
	vector<float> &els_r9();
	vector<float> &els_r9_full5x5();
	vector<float> &els_sigmaEtaEta();
	vector<float> &els_sigmaEtaEta_full5x5();
	vector<float> &els_sigmaIEtaIEta();
	vector<float> &els_sigmaIEtaIEtaSC();
	vector<float> &els_sigmaIEtaIEta_full5x5();
	vector<float> &els_sigmaIEtaIPhi();
	vector<float> &els_sigmaIPhiIPhi();
	vector<float> &els_sigmaIPhiIPhiSC();
	vector<float> &els_sigmaIPhiIPhi_full5x5();
	vector<float> &els_sigmaIphiIphi();
	vector<float> &els_sigmaPhiPhi();
	vector<float> &els_tkIso();
	vector<float> &els_tkIso04();
	vector<float> &els_trackMomentumError();
	vector<float> &els_trkdr();
	vector<float> &els_trkshFrac();
	vector<float> &els_z0();
	vector<float> &els_z0Err();
	vector<float> &els_z0corr();
	vector<float> &mus_backToBackCompat();
	vector<float> &mus_best_dxyPV();
	vector<float> &mus_best_dzPV();
	vector<float> &mus_caloCompatibility();
	vector<float> &mus_chi2();
	vector<float> &mus_chi2LocalMomentum();
	vector<float> &mus_chi2LocalPosition();
	vector<float> &mus_cosmicCompat();
	vector<float> &mus_d0();
	vector<float> &mus_d0Err();
	vector<float> &mus_d0corr();
	vector<float> &mus_d0corrPhi();
	vector<float> &mus_d0phiCov();
	vector<float> &mus_dxyPV();
	vector<float> &mus_dzPV();
	vector<float> &mus_ecal_time();
	vector<float> &mus_ecal_timeError();
	vector<float> &mus_e_em();
	vector<float> &mus_e_emS9();
	vector<float> &mus_e_had();
	vector<float> &mus_e_hadS9();
	vector<float> &mus_e_ho();
	vector<float> &mus_e_hoS9();
	vector<float> &mus_emMax();
	vector<float> &mus_emS25();
	vector<float> &mus_etaErr();
	vector<float> &mus_gfit_chi2();
	vector<float> &mus_gfit_d0();
	vector<float> &mus_gfit_d0Err();
	vector<float> &mus_gfit_d0corr();
	vector<float> &mus_gfit_d0corrPhi();
	vector<float> &mus_gfit_d0phiCov();
	vector<float> &mus_gfit_ndof();
	vector<float> &mus_gfit_qoverp();
	vector<float> &mus_gfit_qoverpError();
	vector<float> &mus_gfit_z0();
	vector<float> &mus_gfit_z0Err();
	vector<float> &mus_gfit_z0corr();
	vector<float> &mus_glbKink();
	vector<float> &mus_glbTrackProbability();
	vector<float> &mus_globalDeltaEtaPhi();
	vector<float> &mus_hadMax();
	vector<float> &mus_hcal_time();
	vector<float> &mus_hcal_timeError();
	vector<float> &mus_ip3d();
	vector<float> &mus_ip3derr();
	vector<float> &mus_iso03_emEt();
	vector<float> &mus_iso03_hadEt();
	vector<float> &mus_iso03_hoEt();
	vector<float> &mus_iso03_sumPt();
	vector<float> &mus_iso05_emEt();
	vector<float> &mus_iso05_hadEt();
	vector<float> &mus_iso05_hoEt();
	vector<float> &mus_iso05_sumPt();
	vector<float> &mus_isoMeanDRR03_pf_ChargedHadronPt();
	vector<float> &mus_isoMeanDRR03_pf_ChargedParticlePt();
	vector<float> &mus_isoMeanDRR03_pf_NeutralHadronEt();
	vector<float> &mus_isoMeanDRR03_pf_NeutralHadronEtHighThreshold();
	vector<float> &mus_isoMeanDRR03_pf_PUPt();
	vector<float> &mus_isoMeanDRR03_pf_PhotonEt();
	vector<float> &mus_isoMeanDRR03_pf_PhotonEtHighThreshold();
	vector<float> &mus_isoMeanDRR04_pf_ChargedHadronPt();
	vector<float> &mus_isoMeanDRR04_pf_ChargedParticlePt();
	vector<float> &mus_isoMeanDRR04_pf_NeutralHadronEt();
	vector<float> &mus_isoMeanDRR04_pf_NeutralHadronEtHighThreshold();
	vector<float> &mus_isoMeanDRR04_pf_PUPt();
	vector<float> &mus_isoMeanDRR04_pf_PhotonEt();
	vector<float> &mus_isoMeanDRR04_pf_PhotonEtHighThreshold();
	vector<float> &mus_isoR03_pf_ChargedHadronPt();
	vector<float> &mus_isoR03_pf_ChargedParticlePt();
	vector<float> &mus_isoR03_pf_NeutralHadronEt();
	vector<float> &mus_isoR03_pf_NeutralHadronEtHighThreshold();
	vector<float> &mus_isoR03_pf_PUPt();
	vector<float> &mus_isoR03_pf_PhotonEt();
	vector<float> &mus_isoR03_pf_PhotonEtHighThreshold();
	vector<float> &mus_isoR04_pf_ChargedHadronPt();
	vector<float> &mus_isoR04_pf_ChargedParticlePt();
	vector<float> &mus_isoR04_pf_NeutralHadronEt();
	vector<float> &mus_isoR04_pf_NeutralHadronEtHighThreshold();
	vector<float> &mus_isoR04_pf_PUPt();
	vector<float> &mus_isoR04_pf_PhotonEt();
	vector<float> &mus_isoR04_pf_PhotonEtHighThreshold();
	vector<float> &mus_isoSumDRR03_pf_ChargedHadronPt();
	vector<float> &mus_isoSumDRR03_pf_ChargedParticlePt();
	vector<float> &mus_isoSumDRR03_pf_NeutralHadronEt();
	vector<float> &mus_isoSumDRR03_pf_NeutralHadronEtHighThreshold();
	vector<float> &mus_isoSumDRR03_pf_PUPt();
	vector<float> &mus_isoSumDRR03_pf_PhotonEt();
	vector<float> &mus_isoSumDRR03_pf_PhotonEtHighThreshold();
	vector<float> &mus_isoSumDRR04_pf_ChargedHadronPt();
	vector<float> &mus_isoSumDRR04_pf_ChargedParticlePt();
	vector<float> &mus_isoSumDRR04_pf_NeutralHadronEt();
	vector<float> &mus_isoSumDRR04_pf_NeutralHadronEtHighThreshold();
	vector<float> &mus_isoSumDRR04_pf_PUPt();
	vector<float> &mus_isoSumDRR04_pf_PhotonEt();
	vector<float> &mus_isoSumDRR04_pf_PhotonEtHighThreshold();
	vector<float> &mus_iso_ecalvetoDep();
	vector<float> &mus_iso_hcalvetoDep();
	vector<float> &mus_iso_hovetoDep();
	vector<float> &mus_iso_trckvetoDep();
	vector<float> &mus_localDistance();
	vector<float> &mus_mc_patMatch_dr();
	vector<float> &mus_ndof();
	vector<float> &mus_overlapCompat();
	vector<float> &mus_phiErr();
	vector<float> &mus_ptErr();
	vector<float> &mus_qoverp();
	vector<float> &mus_qoverpError();
	vector<float> &mus_segmCompatibility();
	vector<float> &mus_staRelChi2();
	vector<float> &mus_sta_chi2();
	vector<float> &mus_sta_d0();
	vector<float> &mus_sta_d0Err();
	vector<float> &mus_sta_d0corr();
	vector<float> &mus_sta_d0corrPhi();
	vector<float> &mus_sta_d0phiCov();
	vector<float> &mus_sta_ndof();
	vector<float> &mus_sta_qoverp();
	vector<float> &mus_sta_qoverpError();
	vector<float> &mus_sta_z0();
	vector<float> &mus_sta_z0Err();
	vector<float> &mus_sta_z0corr();
	vector<float> &mus_timeAtIpInOut();
	vector<float> &mus_timeAtIpInOutErr();
	vector<float> &mus_timeAtIpOutIn();
	vector<float> &mus_timeAtIpOutInErr();
	vector<float> &mus_timeCompat();
	vector<float> &mus_tower();
	vector<float> &mus_towerS9();
	vector<float> &mus_trkKink();
	vector<float> &mus_trkRelChi2();
	vector<float> &mus_vertexCompat();
	vector<float> &mus_vertexphi();
	vector<float> &mus_z0();
	vector<float> &mus_z0Err();
	vector<float> &mus_z0corr();
	vector<float> &pfcands_dz();
	vector<float> &pfjets_area();
	vector<float> &pfjets_chargedEmE();
	vector<float> &pfjets_chargedHadronE();
	vector<float> &pfjets_combinedSecondaryVertexBJetTag();
	vector<float> &pfjets_cor();
	vector<float> &pfjets_corL1Fast();
	vector<float> &pfjets_corL1FastL2L3();
	vector<float> &pfjets_corL1FastL2L3residual();
	vector<float> &pfjets_corL1L2L3();
	vector<float> &pfjets_electronE();
	vector<float> &pfjets_hfEmE();
	vector<float> &pfjets_hfHadronE();
	vector<float> &pfjets_jetBProbabilityBJetTag();
	vector<float> &pfjets_jetProbabilityBJetTag();
	vector<float> &pfjets_muonE();
	vector<float> &pfjets_neutralEmE();
	vector<float> &pfjets_neutralHadronE();
	vector<float> &pfjets_photonE();
	vector<float> &pfjets_pileupJetId();
	vector<float> &pfjets_simpleSecondaryVertexHighEffBJetTag();
	vector<float> &pfjets_simpleSecondaryVertexHighPurBJetTags();
	vector<float> &pfjets_trackCountingHighEffBJetTag();
	vector<float> &pfjets_trackCountingHighPurBJetTag();
	vector<float> &taus_pf_againstElectronDeadECAL();
	vector<float> &taus_pf_againstElectronLoose();
	vector<float> &taus_pf_againstElectronMedium();
	vector<float> &taus_pf_againstElectronTight();
	vector<float> &taus_pf_againstMuonLoose();
	vector<float> &taus_pf_againstMuonLoose2();
	vector<float> &taus_pf_againstMuonLoose3();
	vector<float> &taus_pf_againstMuonMedium();
	vector<float> &taus_pf_againstMuonMedium2();
	vector<float> &taus_pf_againstMuonTight();
	vector<float> &taus_pf_againstMuonTight2();
	vector<float> &taus_pf_againstMuonTight3();
	vector<float> &taus_pf_byCombinedIsolationDeltaBetaCorrRaw();
	vector<float> &taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits();
	vector<float> &taus_pf_byDecayModeFinding();
	vector<float> &taus_pf_byLooseCombinedIsolationDeltaBetaCorr();
	vector<float> &taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits();
	vector<float> &taus_pf_byMediumCombinedIsolationDeltaBetaCorr();
	vector<float> &taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits();
	vector<float> &taus_pf_byTightCombinedIsolationDeltaBetaCorr();
	vector<float> &taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits();
	vector<float> &taus_pf_byVLooseCombinedIsolationDeltaBetaCorr();
	vector<float> &photons_e1x5();
	vector<float> &photons_e2x5Max();
	vector<float> &photons_e3x3();
	vector<float> &photons_e5x5();
	vector<float> &photons_ecalIso03();
	vector<float> &photons_ecalIso04();
	vector<float> &photons_full3x3_e3x3();
	vector<float> &photons_full5x5_e1x5();
	vector<float> &photons_full5x5_e2x5Max();
	vector<float> &photons_full5x5_e5x5();
	vector<float> &photons_full5x5_hOverE();
	vector<float> &photons_full5x5_sigmaEtaEta();
	vector<float> &photons_full5x5_sigmaIEtaIEta();
	vector<float> &photons_hOverE();
	vector<float> &photons_hcalDepth1TowerSumEtBcConeDR03();
	vector<float> &photons_hcalDepth1TowerSumEtBcConeDR04();
	vector<float> &photons_hcalDepth2TowerSumEtBcConeDR03();
	vector<float> &photons_hcalDepth2TowerSumEtBcConeDR04();
	vector<float> &photons_hcalIso03();
	vector<float> &photons_hcalIso04();
	vector<float> &photons_hcalTowerSumEtBcConeDR03();
	vector<float> &photons_hcalTowerSumEtBcConeDR04();
	vector<float> &photons_ntkIsoHollow03();
	vector<float> &photons_ntkIsoHollow04();
	vector<float> &photons_ntkIsoSolid03();
	vector<float> &photons_ntkIsoSolid04();
	vector<float> &photons_sigmaEtaEta();
	vector<float> &photons_sigmaIEtaIEta();
	vector<float> &photons_tkIsoHollow03();
	vector<float> &photons_tkIsoHollow04();
	vector<float> &photons_tkIsoSolid03();
	vector<float> &photons_tkIsoSolid04();
	vector<float> &puInfo_trueNumInteractions();
	vector<float> &convs_chi2();
	vector<float> &convs_dl();
	vector<float> &convs_ndof();
	vector<float> &vtxs_chi2();
	vector<float> &vtxs_ndof();
	vector<float> &vtxs_sumpt();
	vector<float> &vtxs_xError();
	vector<float> &vtxs_yError();
	vector<float> &vtxs_zError();
	vector<vector<float> > &els_convs_dcot();
	vector<vector<float> > &els_convs_dist();
	vector<vector<float> > &els_convs_radius();
	vector<vector<float> > &mus_stationShowerDeltaR();
	vector<vector<float> > &mus_stationShowerSizeT();
	vector<vector<float> > &puInfo_instLumi();
	vector<vector<float> > &vtxs_covMatrix();
	int &evt_bsType();
	int &evt_bunchCrossing();
	int &evt_experimentType();
	int &evt_isRealData();
	int &evt_orbitNumber();
	int &evt_storeNumber();
	int &hcalnoise_GetRecHitCount();
	int &hcalnoise_GetRecHitCount15();
	int &hcalnoise_maxHPDHits();
	int &hcalnoise_maxHPDNoOtherHits();
	int &hcalnoise_maxRBXHits();
	int &hcalnoise_maxZeros();
	int &hcalnoise_noiseFilterStatus();
	int &hcalnoise_noiseType();
	int &hcalnoise_num10GeVHits();
	int &hcalnoise_num25GeVHits();
	int &hcalnoise_numFlatNoiseChannels();
	int &hcalnoise_numIsolatedNoiseChannels();
	int &hcalnoise_numProblematicRBXs();
	int &hcalnoise_numSpikeNoiseChannels();
	int &hcalnoise_numTS4TS5NoiseChannels();
	int &hcalnoise_numTriangleNoiseChannels();
	int &hcalnoise_passHighLevelNoiseFilter();
	int &hcalnoise_passLooseNoiseFilter();
	int &hcalnoise_passTightNoiseFilter();
	int &pdfinfo_id1();
	int &pdfinfo_id2();
	vector<int> &els_mc3_id();
	vector<int> &els_mc3idx();
	vector<int> &els_mc3_motherid();
	vector<int> &els_mc3_motheridx();
	vector<int> &els_mc_id();
	vector<int> &els_mcidx();
	vector<int> &els_mc_motherid();
	vector<int> &mus_mc3_id();
	vector<int> &mus_mc3idx();
	vector<int> &mus_mc3_motherid();
	vector<int> &mus_mc3_motheridx();
	vector<int> &mus_mc_id();
	vector<int> &mus_mcidx();
	vector<int> &mus_mc_motherid();
	vector<int> &pfjets_mc3_id();
	vector<int> &pfjets_mc3idx();
	vector<int> &pfjets_mc_gpidx();
	vector<int> &pfjets_mc_id();
	vector<int> &pfjets_mcidx();
	vector<int> &pfjets_mc_motherid();
	vector<int> &photons_mc3_id();
	vector<int> &photons_mc3idx();
	vector<int> &photons_mc3_motherid();
	vector<int> &photons_mc3_motheridx();
	vector<int> &photons_mc_id();
	vector<int> &photons_mcidx();
	vector<int> &photons_mc_motherid();
	vector<int> &els_category();
	vector<int> &els_charge();
	vector<int> &els_ckf_charge();
	vector<int> &els_ckf_laywithmeas();
	vector<int> &els_class();
	vector<int> &els_conv_delMissHits();
	vector<int> &els_conv_flag();
	vector<int> &els_conv_gsftkidx();
	vector<int> &els_conv_old_delMissHits();
	vector<int> &els_conv_old_flag();
	vector<int> &els_conv_old_gsftkidx();
	vector<int> &els_conv_old_tkidx();
	vector<int> &els_conv_tkidx();
	vector<int> &els_exp_innerlayers();
	vector<int> &els_exp_outerlayers();
	vector<int> &els_fiduciality();
	vector<int> &els_gsftrkidx();
	vector<int> &els_layer1_det();
	vector<int> &els_layer1_layer();
	vector<int> &els_layer1_sizerphi();
	vector<int> &els_layer1_sizerz();
	vector<int> &els_lostHits();
	vector<int> &els_lost_pixelhits();
	vector<int> &els_mc_patMatch_id();
	vector<int> &els_nSeed();
	vector<int> &els_nlayers();
	vector<int> &els_nlayers3D();
	vector<int> &els_nlayersLost();
	vector<int> &els_sccharge();
	vector<int> &els_scindex();
	vector<int> &els_trk_charge();
	vector<int> &els_trkidx();
	vector<int> &els_type();
	vector<int> &els_validHits();
	vector<int> &els_valid_pixelhits();
	vector<int> &genps_id();
	vector<int> &genps_id_mother();
	vector<int> &genps_status();
	vector<int> &mus_algo();
	vector<int> &mus_charge();
	vector<int> &mus_ecal_rawId();
	vector<int> &mus_exp_innerlayers();
	vector<int> &mus_exp_outerlayers();
	vector<int> &mus_gfit_algo();
	vector<int> &mus_gfit_etaErr();
	vector<int> &mus_gfit_exp_innerlayer();
	vector<int> &mus_gfit_exp_outerlayer();
	vector<int> &mus_gfit_lostHits();
	vector<int> &mus_gfit_lostPixelHits();
	vector<int> &mus_gfit_nlayers();
	vector<int> &mus_gfit_nlayers3D();
	vector<int> &mus_gfit_nlayersLost();
	vector<int> &mus_gfit_phiErr();
	vector<int> &mus_gfit_ptErr();
	vector<int> &mus_gfit_qualityMask();
	vector<int> &mus_gfit_trk_charge();
	vector<int> &mus_gfit_validHits();
	vector<int> &mus_gfit_validPixelHits();
	vector<int> &mus_gfit_validSTAHits();
	vector<int> &mus_gfit_validSiHits();
	vector<int> &mus_goodmask();
	vector<int> &mus_hcal_rawId();
	vector<int> &mus_iso03_njets();
	vector<int> &mus_iso03_ntrk();
	vector<int> &mus_iso05_ntrk();
	vector<int> &mus_lostHits();
	vector<int> &mus_lostPixelHits();
	vector<int> &mus_mc_patMatch_id();
	vector<int> &mus_muonBestTrackType();
	vector<int> &mus_nOverlaps();
	vector<int> &mus_nlayers();
	vector<int> &mus_nlayers3D();
	vector<int> &mus_nlayersLost();
	vector<int> &mus_nmatches();
	vector<int> &mus_numberOfMatchedStations();
	vector<int> &mus_overlap0();
	vector<int> &mus_overlap1();
	vector<int> &mus_pfcharge();
	vector<int> &mus_pfparticleId();
	vector<int> &mus_pid_PFMuon();
	vector<int> &mus_pid_TM2DCompatibilityLoose();
	vector<int> &mus_pid_TM2DCompatibilityTight();
	vector<int> &mus_pid_TMLastStationLoose();
	vector<int> &mus_pid_TMLastStationTight();
	vector<int> &mus_qualityMask();
	vector<int> &mus_sta_algo();
	vector<int> &mus_sta_etaErr();
	vector<int> &mus_sta_exp_innerlayer();
	vector<int> &mus_sta_exp_outerlayer();
	vector<int> &mus_sta_lostHits();
	vector<int> &mus_sta_lostPixelHits();
	vector<int> &mus_sta_nlayers();
	vector<int> &mus_sta_nlayers3D();
	vector<int> &mus_sta_nlayersLost();
	vector<int> &mus_sta_phiErr();
	vector<int> &mus_sta_ptErr();
	vector<int> &mus_sta_qualityMask();
	vector<int> &mus_sta_trk_charge();
	vector<int> &mus_sta_validHits();
	vector<int> &mus_sta_validPixelHits();
	vector<int> &mus_timeDirection();
	vector<int> &mus_timeNumStationsUsed();
	vector<int> &mus_trk_charge();
	vector<int> &mus_trkidx();
	vector<int> &mus_type();
	vector<int> &mus_validHits();
	vector<int> &mus_validPixelHits();
	vector<int> &pfcands_charge();
	vector<int> &pfcands_particleId();
	vector<int> &pfjets_chargedHadronMultiplicity();
	vector<int> &pfjets_chargedMultiplicity();
	vector<int> &pfjets_electronMultiplicity();
	vector<int> &pfjets_hfEmMultiplicity();
	vector<int> &pfjets_hfHadronMultiplicity();
	vector<int> &pfjets_muonMultiplicity();
	vector<int> &pfjets_neutralHadronMultiplicity();
	vector<int> &pfjets_neutralMultiplicity();
	vector<int> &pfjets_partonFlavour();
	vector<int> &pfjets_photonMultiplicity();
	vector<int> &taus_pf_charge();
	vector<int> &photons_fiduciality();
	vector<int> &puInfo_bunchCrossing();
	vector<int> &puInfo_nPUvertices();
	vector<int> &convs_algo();
	vector<int> &convs_isConverted();
	vector<int> &convs_quality();
	vector<int> &vtxs_isFake();
	vector<int> &vtxs_isValid();
	vector<int> &vtxs_tracksSize();
	vector<vector<int> > &els_convs_delMissHits();
	vector<vector<int> > &els_convs_flag();
	vector<vector<int> > &els_convs_gsftkidx();
	vector<vector<int> > &els_convs_tkidx();
	vector<vector<int> > &els_PFCand_idx();
	vector<vector<int> > &genps_lepdaughter_id();
	vector<vector<int> > &genps_lepdaughter_idx();
	vector<vector<int> > &hlt_trigObjs_id();
	vector<vector<int> > &mus_nStationCorrelatedHits();
	vector<vector<int> > &mus_nStationHits();
	vector<vector<int> > &pfjets_pfcandIndicies();
	vector<vector<int> > &convs_nHitsBeforeVtx();
	vector<vector<int> > &convs_tkalgo();
	vector<vector<int> > &convs_tkidx();
	unsigned int &els_HLT_Ele17_Ele8_L1sL1DoubleEG137_version();
	unsigned int &els_HLT_Ele17_Ele8_LeadingLeg_version();
	unsigned int &els_HLT_Ele17_Ele8_Mass50_LeadingLeg_version();
	unsigned int &els_HLT_Ele17_Ele8_Mass50_TrailingLeg_version();
	unsigned int &els_HLT_Ele17_Ele8_TrailingLeg_version();
	unsigned int &els_HLT_Ele17_Ele8_version();
	unsigned int &els_HLT_Ele20_SC4_Mass50_LeadingLeg_version();
	unsigned int &els_HLT_Ele20_SC4_Mass50_TrailingLeg_version();
	unsigned int &els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22_version();
	unsigned int &els_HLT_Ele27_WP80_version();
	unsigned int &els_HLT_Ele32_SC17_Mass50_LeadingLeg_version();
	unsigned int &els_HLT_Ele32_SC17_Mass50_TrailingLeg_version();
	unsigned int &els_HLT_Mu17_Ele8_TrailingLeg_version();
	unsigned int &els_HLT_Mu17_Ele8_version();
	unsigned int &els_HLT_Mu8_Ele17_version();
	unsigned int &evt_nels();
	unsigned int &evt_detectorStatus();
	unsigned int &evt_event();
	unsigned int &evt_lumiBlock();
	unsigned int &evt_run();
	unsigned int &evt_ngenjetsNoMuNoNu();
	unsigned int &genps_signalProcessID();
	unsigned int &mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1_version();
	unsigned int &mus_HLT_IsoMu24_eta2p1_version();
	unsigned int &mus_HLT_Mu17_Ele8_LeadingLeg_version();
	unsigned int &mus_HLT_Mu17_Ele8_version();
	unsigned int &mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen_version();
	unsigned int &mus_HLT_Mu17_Mu8_LeadingLeg_version();
	unsigned int &mus_HLT_Mu17_Mu8_TrailingLeg_version();
	unsigned int &mus_HLT_Mu17_Mu8_version();
	unsigned int &mus_HLT_Mu17_TkMu8_LeadingLeg_version();
	unsigned int &mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered_version();
	unsigned int &mus_HLT_Mu17_TkMu8_TrailingLeg_version();
	unsigned int &mus_HLT_Mu17_TkMu8_version();
	unsigned int &mus_HLT_Mu8_Ele17_TrailingLeg_version();
	unsigned int &mus_HLT_Mu8_Ele17_version();
	unsigned int &evt_nphotons();
	unsigned int &evt_nvtxs();
	vector<unsigned int> &els_HLT_Ele17_Ele8();
	vector<unsigned int> &els_HLT_Ele17_Ele8_L1sL1DoubleEG137();
	vector<unsigned int> &els_HLT_Ele17_Ele8_LeadingLeg();
	vector<unsigned int> &els_HLT_Ele17_Ele8_Mass50_LeadingLeg();
	vector<unsigned int> &els_HLT_Ele17_Ele8_Mass50_TrailingLeg();
	vector<unsigned int> &els_HLT_Ele17_Ele8_TrailingLeg();
	vector<unsigned int> &els_HLT_Ele20_SC4_Mass50_LeadingLeg();
	vector<unsigned int> &els_HLT_Ele20_SC4_Mass50_TrailingLeg();
	vector<unsigned int> &els_HLT_Ele27_WP80();
	vector<unsigned int> &els_HLT_Ele27_WP80_L1sL1SingleEG20ORL1SingleEG22();
	vector<unsigned int> &els_HLT_Ele32_SC17_Mass50_LeadingLeg();
	vector<unsigned int> &els_HLT_Ele32_SC17_Mass50_TrailingLeg();
	vector<unsigned int> &els_HLT_Mu17_Ele8();
	vector<unsigned int> &els_HLT_Mu17_Ele8_TrailingLeg();
	vector<unsigned int> &els_HLT_Mu8_Ele17();
	vector<unsigned int> &els_id2012ext_loose();
	vector<unsigned int> &els_id2012ext_medium();
	vector<unsigned int> &els_id2012ext_tight();
	vector<unsigned int> &els_id2012ext_veto();
	vector<unsigned int> &els_id2012_loose();
	vector<unsigned int> &els_id2012_medium();
	vector<unsigned int> &els_id2012_tight();
	vector<unsigned int> &els_id2012_veto();
	vector<unsigned int> &hlt_prescales();
	vector<unsigned int> &mus_HLT_IsoMu24_eta2p1();
	vector<unsigned int> &mus_HLT_IsoMu24_eta2p1_L1sMu16Eta2p1();
	vector<unsigned int> &mus_HLT_Mu17_Ele8();
	vector<unsigned int> &mus_HLT_Mu17_Ele8_LeadingLeg();
	vector<unsigned int> &mus_HLT_Mu17_Mu8();
	vector<unsigned int> &mus_HLT_Mu17_Mu8_L1sL1DoubleMu10MuOpen();
	vector<unsigned int> &mus_HLT_Mu17_Mu8_LeadingLeg();
	vector<unsigned int> &mus_HLT_Mu17_Mu8_TrailingLeg();
	vector<unsigned int> &mus_HLT_Mu17_TkMu8();
	vector<unsigned int> &mus_HLT_Mu17_TkMu8_LeadingLeg();
	vector<unsigned int> &mus_HLT_Mu17_TkMu8_TrailingLeg();
	vector<unsigned int> &mus_HLT_Mu17_TkMu8_TrailingLegTrkFiltered();
	vector<unsigned int> &mus_HLT_Mu8_Ele17();
	vector<unsigned int> &mus_HLT_Mu8_Ele17_TrailingLeg();
	bool passHLTTrigger(TString trigName);
}
#endif
