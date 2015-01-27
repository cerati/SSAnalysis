#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TF1.h"
#include "Math/VectorUtil.h" 
#include "SSCORE/CMS2.h"
#include "SSCORE/selections.h"
#include "MT2AG.h"
//#include "structAG.h"
//#include "looper.h" 
//#include "fromCore.h"
//#include "structAG.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef vector<pair<const LorentzVector *, double> > jets_with_corr_t;

using namespace std;

//Enums
enum anal_type_t { HIGH_PT, LOW_PT, VLOW_PT };
enum hyp_type_t { EE, MM, EM, UNASSIGNED }; 
enum MT2_type_t { MT2W, MT2B }; 
enum signal_t { T1TTTT, T5TTTT, T5LNU, T5VV, T7BTW, T4TW, T6TTWW_X08, T6TTWW_X05, T6TTHH, NONE};
enum btag_type_t { CSVL, CSVM, CSVT };

//Structs
struct btagInfo_t { LorentzVector jetcor; int flavorPhys; int flavorAlgo; bool bjetFlags; };
struct btags_t { vector<LorentzVector> btags; vector <bool> isBtagged; };
struct hyp_result_t { int best_hyp; int hyp_class; };
struct particle_t { int id; LorentzVector p4; int idx; };
struct val_err_t { float value; float error; };

//Comparisons
//bool compare_btagInfo(const btagInfo_t a, const btagInfo_t b){ return a.jetcor.pt() > b.jetcor.pt(); }

//Declare functions in other files
//val_err_t gluino_cs(float mGluino);
//val_err_t stop_sbottom_cs(float mGluino);

//Classes
class babyMaker {

  public:
    babyMaker(bool debug = 0) {
      path = ".";
      verbose = debug;
      evt_cut = 13398;
    }
    void MakeBabyNtuple(const char* output_name);
    void InitBabyNtuple();
    void CloseBabyNtuple () { BabyFile->cd();BabyTree->Write();BabyFile->Close();}
    int ProcessBaby();
    hyp_result_t chooseBestHyp();
    int isGoodHyp(int iHyp, int analType = 2);
    pair <particle_t, int> getThirdLepton(int hyp);
    vector <particle_t> getGenPair();

  protected:
    TFile* BabyFile;
    TTree* BabyTree;

  private:

    //Switches
    TString path;
    bool verbose;
    unsigned int evt_cut;

    //MET
    float met;
    float metPhi;
 
    //Meta Variables
    int event;
    int lumi;
    int run;
    bool is_real_data;
    float scale1fb;     
    float xsec;         
    float kfactor;      
    TString filename;

    //Filters
    bool filt_csc;
    bool filt_hbhe;
    bool filt_hcallaser;
    bool filt_ecaltp;
    bool filt_trkfail;
    bool filt_eebadsc;

    //Gen MET 
    float gen_met;      
    float gen_met_phi;  

    //Jets
    int njets;
    float ht;
    vector <LorentzVector> jets;
    vector <float> jets_disc;

    //Hyp Class -- in this order
       //3 for num-num SS leptons
       //2 for num-den SS leptons
       //1 for den-den SS leptons
       //4 for num-num OS leptons
       //0 otherwise (not saved)
    int hyp_class;

    //Leptons
    LorentzVector lep1_p4;
    LorentzVector lep2_p4;
    LorentzVector dilep_p4;
    int lep1_id;
    int lep2_id;
    int lep1_idx;
    int lep2_idx;
    int hyp_type; 

    //Lepton Mother (old isFromW function)
    int lep1_motherID;
    int lep2_motherID;

    //Lepton Truth ID
      //From src/MatchUtilities/matchCandToGen function
      //delta-R matches to nearest particle other than neutrino or LSP
      //must be within 0.2
    int lep1_mc_id;
    int lep2_mc_id;

    //b-tags
    vector <LorentzVector> btags;
    vector <float> btags_disc;
    int nbtags;

    //Scale factors (from 8 TeV, outdated)
    float sf_dilepTrig_hpt; 
    float sf_dilepTrig_lpt; 
    float sf_dilepTrig_vlpt; 
    float sf_dilep_eff;

    //mT for both leptons, mt2
    float mt;
    float mt_l2;
    float mt2;

    //SUSY sparms
    float mLSP;
    float mGluino;
    float mSbottom;
    float mChargino;
 
    //Look at gen-level particles and choose favorite hypothesis
    //Limited usefulness as-is; probably not smart enough to reject non-prompt leptons
    int lep1_id_gen;
    int lep2_id_gen;
    LorentzVector lep1_p4_gen;
    LorentzVector lep2_p4_gen;

    //Third lepton -- quality is 2 for good, 1 for fakeable, 0 for veto
    int lep3_id;
    int lep3_idx;
    int lep3_quality;
    LorentzVector lep3_p4;

    //Isolation
    float lep1_iso;
    float lep2_iso;

    //Gen Leptons
    vector <LorentzVector> genps_p4;
    vector <int> genps_id;
    vector <int> genps_id_mother;
    vector <int> genps_status;
	vector <int> genps_id_grandma;

   //Leptons pass numerator ID
    bool lep1_passes_id;
    bool lep2_passes_id;

    //Imparct parameter
    float lep1_dxyPV;
    float lep2_dxyPV;
    float lep1_d0_err;
    float lep2_d0_err;
    float lep1_dZ;
    float lep2_dZ;
    float lep1_ip3d;
    float lep1_ip3d_err;
    float lep2_ip3d;
    float lep2_ip3d_err;

    //nGoodElectrons passing certain pT cuts
    int nGoodElectrons7;
    int nGoodElectrons10;
    int nGoodElectrons25;
    int nGoodMuons5;
    int nGoodMuons10;
    int nGoodMuons25;
};
