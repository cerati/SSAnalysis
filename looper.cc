// C++
#include <iostream>
#include <vector>

// ROOT
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TMath.h"

#include "./looper.h"

// CMS3
#include "./SSCORE/CMS2.h"
#include "./SSCORE/selections.h"

using namespace tas;
using namespace std;

struct TightLep {
  TightLep(int pdgid_, int idx_):pdgid(pdgid_),idx(idx_){}
  int charge() {return -1*pdgid/abs(pdgid);}
  int pdgId() {return pdgid;}
  float pt() {return abs(pdgid)==11 ? els_p4().at(idx).pt() : mus_p4().at(idx).pt();}
  LorentzVector p4() {return abs(pdgid)==11 ? els_p4().at(idx) : mus_p4().at(idx);}
  int mc3_id() {return abs(pdgid)==11 ? els_mc3_id().at(idx) : mus_mc3_id().at(idx);}
  int mc3_motherid() {return abs(pdgid)==11 ? els_mc3_motherid().at(idx) : mus_mc3_motherid().at(idx);}
  int mc_id() {
    return abs(pdgid)==11 ? els_mc_id().at(idx) : mus_mc_id().at(idx);
  }
  int mc_motherid() {return abs(pdgid)==11 ? els_mc_motherid().at(idx) : mus_mc_motherid().at(idx);}
private:
  int pdgid, idx;
};

struct DilepHyp {
  DilepHyp(TightLep lepone_, TightLep leptwo_):
    leadlep(lepone_),trailep(leptwo_) {
    if (lepone_.pt()<leptwo_.pt()) {
      leadlep = leptwo_;
      trailep = lepone_;
    }
  }
  int charge() {return leadlep.charge()+trailep.charge();}
  LorentzVector p4() {return leadlep.p4()+trailep.p4();}
  TightLep leadLep() {return leadlep;}
  TightLep traiLep() {return trailep;}
private:
  TightLep leadlep, trailep;
};

int looper::ScanChain( TChain* chain, const char* prefix, bool isData, int nEvents) {

  makebaby       = false;
  makehist       = true;
  maketext       = false;
  
  if (makebaby) MakeBabyNtuple( Form( "%s_baby.root", prefix ) );
  if (makehist) CreateOutputFile( Form( "%s_histos.root", prefix ) );

  // File Loop
  if( nEvents == -1 ) nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    // Get File Content
    TFile f( currentFile->GetTitle() );
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    // Event Loop
    unsigned int nEvents = tree->GetEntries();

    for(unsigned int event = 0; event < nEvents; ++event) {
    
      // Get Event Content
      cms2.GetEntry(event);
      ++nEventsTotal;

      // progress feedback to user
      if (nEventsTotal % 1000 == 0) {
	// xterm magic from L. Vacavant and A. Cerri
	if (isatty(1)) {
	  printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
		 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
	  fflush(stdout);
	}
      }

      if (makebaby) InitBabyNtuple();

      //-------- DO THE ANALYSIS HERE --------//

      //print info
      if (maketext) printEvent();

      //fill baby
      run_   = evt_run();
      ls_    = evt_lumiBlock();
      evt_   = evt_event();
      weight_ = isData ? 1. : evt_scale1fb();

      //start selection
      //event

      //jets, met, ht

      //leptons
      vector<TightLep> tightleps;
      for (unsigned int elidx=0;elidx<els_p4().size();++elidx) {
	//tight electron selection
	if (els_p4().at(elidx).pt()<20. || isTightElectron(elidx)==0 || threeChargeAgree(elidx)==0) continue;
	//if (els_p4().at(elidx).pt()<20. || isMediumElectron(elidx)==0 || threeChargeAgree(elidx)==0) continue;
	//if (els_p4().at(elidx).pt()<20. || isLooseElectron(elidx)==0 || threeChargeAgree(elidx)==0) continue;
	TightLep tightel(-1*els_charge().at(elidx)*11,elidx);
	tightleps.push_back(tightel);
      }
      for (unsigned int muidx=0;muidx<mus_p4().size();++muidx) {
	//tight muon selection
	if (mus_p4().at(muidx).pt()<20. || isTightMuon(muidx)==0 || muRelIso03(muidx)>0.15 ) continue;
	TightLep tightmu(-1*mus_charge().at(muidx)*13,muidx);
	tightleps.push_back(tightmu);
      }
      //cout << "muon size=" << mus_p4().size() << " electron size=" << els_p4().size() << endl;
      if (tightleps.size()!=2) continue;
      cout << "Tight leptons = " << tightleps.size() << " pdgids=" << tightleps[0].pdgId() << ", " << tightleps[1].pdgId() << endl;

      DilepHyp hyp(tightleps[0],tightleps[1]);
      cout << "Total charge = " << hyp.charge() << " m=" << hyp.p4().mass() << endl;
      if (hyp.p4().mass()<8) continue;

      //if (makebaby) FillBabyNtuple();

      //fill histos
      if (makehist) {
	//makeFillHisto1D<TH1F,int>("dummy","dummy",2,0,2,1);
	makeFillHisto1D<TH1F,int>("hyp_charge","hyp_charge",7,-3.5,3.5,hyp.charge());
	makeFillHisto1D<TH1F,int>("hyp_mll","hyp_mll",100,0,1000,hyp.p4().mass());
	makeFillHisto1D<TH1F,int>("hyp_ptll","hyp_ptll",100,0,1000,hyp.p4().pt());

	int type = -1;
	if (abs(hyp.traiLep().pdgId())==13 && abs(hyp.leadLep().pdgId())==13) type=0;
	if (abs(hyp.traiLep().pdgId())==13 && abs(hyp.leadLep().pdgId())==11) type=1;
	if (abs(hyp.traiLep().pdgId())==11 && abs(hyp.leadLep().pdgId())==13) type=2;
	if (abs(hyp.traiLep().pdgId())==11 && abs(hyp.leadLep().pdgId())==11) type=3;
	makeFillHisto1D<TH1F,int>("hyp_type","hyp_type",5,0,5,type);

	if (hyp.charge()!=0) {
	  //same sign
	  makeFillHisto1D<TH1F,int>("hyp_ss_mll","hyp_ss_mll",100,0,1000,hyp.p4().mass());
	  makeFillHisto1D<TH1F,int>("hyp_ss_ptll","hyp_ss_ptll",100,0,1000,hyp.p4().pt());
	  makeFillHisto1D<TH1F,int>("hyp_ss_type","hyp_ss_type",5,0,5,type);
	  //check fakes (mother not W)
	  if (abs(hyp.traiLep().mc_motherid())!=24) {
	    if (abs(hyp.traiLep().pdgId())==13) {
	      //trailing muon
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_mu_mc","hyp_ss_trail_mu_mc",11001,-5500.5,5500.5,hyp.traiLep().mc_id());
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_mu_mc_mother","hyp_ss_trail_mu_mc_mother",11001,-5500.5,5500.5,hyp.traiLep().mc_motherid());
	      //cout << evt_event() << " " << hyp.traiLep().pt() << " " << hyp.traiLep().p4().eta() << " " << hyp.traiLep().p4().phi() << endl; 
	      //makeFillHisto1D<TH1F,int>("hyp_ss_trail_mu_mc3","hyp_ss_trail_mu_mc3",11001,-5500.5,5500.5,hyp.traiLep().mc3_id());
	      //makeFillHisto1D<TH1F,int>("hyp_ss_trail_mu_mc3_mother","hyp_ss_trail_mu_mc3_mother",11001,-5500.5,5500.5,hyp.traiLep().mc3_motherid());
	    }
	    if (abs(hyp.traiLep().pdgId())==11) {
	      //trailing elec
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_el_mc","hyp_ss_trail_el_mc",11001,-5500.5,5500.5,hyp.traiLep().mc_id());
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_el_mc_mother","hyp_ss_trail_el_mc_mother",11001,-5500.5,5500.5,hyp.traiLep().mc_motherid());
	      //makeFillHisto1D<TH1F,int>("hyp_ss_trail_el_mc3","hyp_ss_trail_el_mc3",11001,-5500.5,5500.5,hyp.traiLep().mc3_id());
	      //makeFillHisto1D<TH1F,int>("hyp_ss_trail_el_mc3_mother","hyp_ss_trail_el_mc3_mother",11001,-5500.5,5500.5,hyp.traiLep().mc3_motherid());
	    }
	  }
	  if (abs(hyp.leadLep().mc_motherid())!=24) {
	    if (abs(hyp.leadLep().pdgId())==13) {
	      //leading muon
	      makeFillHisto1D<TH1F,int>("hyp_ss_lead_mu_mc","hyp_ss_lead_mu_mc",11001,-5500.5,5500.5,hyp.leadLep().mc_id());
	      makeFillHisto1D<TH1F,int>("hyp_ss_lead_mu_mc_mother","hyp_ss_lead_mu_mc_mother",11001,-5500.5,5500.5,hyp.leadLep().mc_motherid());
	      //cout << evt_event() << " " << hyp.leadLep().pt() << " " << hyp.leadLep().p4().eta() << " " << hyp.leadLep().p4().phi() << endl; 
	      //makeFillHisto1D<TH1F,int>("hyp_ss_lead_mu_mc3","hyp_ss_lead_mu_mc3",11001,-5500.5,5500.5,hyp.leadLep().mc3_id());
	      //makeFillHisto1D<TH1F,int>("hyp_ss_lead_mu_mc3_mother","hyp_ss_lead_mu_mc3_mother",11001,-5500.5,5500.5,hyp.leadLep().mc3_motherid());
	    }
	    if (abs(hyp.leadLep().pdgId())==11) {
	      //leading elec
	      makeFillHisto1D<TH1F,int>("hyp_ss_lead_el_mc","hyp_ss_lead_el_mc",11001,-5500.5,5500.5,hyp.leadLep().mc_id());
	      makeFillHisto1D<TH1F,int>("hyp_ss_lead_el_mc_mother","hyp_ss_lead_el_mc_mother",11001,-5500.5,5500.5,hyp.leadLep().mc_motherid());
	      //makeFillHisto1D<TH1F,int>("hyp_ss_lead_el_mc3","hyp_ss_lead_el_mc3",11001,-5500.5,5500.5,hyp.leadLep().mc3_id());
	      //makeFillHisto1D<TH1F,int>("hyp_ss_lead_el_mc3_mother","hyp_ss_lead_el_mc3_mother",11001,-5500.5,5500.5,hyp.leadLep().mc3_motherid());
	    }
	  }
	} else {
	  //opposite sign
	  makeFillHisto1D<TH1F,int>("hyp_os_mll","hyp_os_mll",100,0,1000,hyp.p4().mass());
	  makeFillHisto1D<TH1F,int>("hyp_os_ptll","hyp_os_ptll",100,0,1000,hyp.p4().pt());
	  makeFillHisto1D<TH1F,int>("hyp_os_type","hyp_os_type",5,0,5,type);
	}
      }//if (makehist)

    }

    delete tree;
    f.Close();
  }

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  if (makebaby) CloseBabyNtuple();
  if (makehist) SaveHistos();

  return 0;
}

void looper::printEvent(  ostream& ostr ){
  ostr << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl; 
}

// Book the baby ntuple
void looper::MakeBabyNtuple(const char *babyFilename) {
  babyFile_ = new TFile(Form("%s", babyFilename), "RECREATE");
  babyTree_ = new TTree("tree", "A Baby Ntuple");
  
  babyTree_->Branch("run", &run_, "run/I");
  babyTree_->Branch("ls", &ls_, "ls/I");
  babyTree_->Branch("evt", &evt_, "evt/I");
  babyTree_->Branch("weight", &weight_, "weight/F");

  babyTree_->Branch("p4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &p4_);
}

// Init the baby
void looper::InitBabyNtuple() { 

  run_ = -9999;
  ls_ = -9999;
  evt_ = -9999;
  weight_ = -9999.;

  babyTree_->SetBranchAddress("p4",          &p4_);
  *p4_ = LorentzVector();

}

void looper::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}
