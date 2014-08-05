// C++
#include <iostream>
#include <vector>
#include <unistd.h> //isatty

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

int looper::ScanChain( TChain* chain, const char* prefix, bool isData, int nEvents) {

  makebaby       = false;
  makehist       = true;
  maketext       = false;
  makeskim       = 0;
  
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

    //Skimmed output file - needs to be before cms2.Init(tree)
    TFile *skim_file = 0;
    TTree* skim_tree = 0;
    if (makeskim) {
      TString skim_file_name = TString(currentFile->GetTitle());
      skim_file_name.ReplaceAll(".root","_newskim.root");
      skim_file_name.ReplaceAll("/hadoop/cms/store/group/snt/csa14/MC_CMS3_V07-00-03/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/merged/","./");
      skim_file = new TFile(skim_file_name,"recreate");
      //locally
      skim_tree = (TTree*) tree->CloneTree(0, "fast");
      skim_tree->SetDirectory(skim_file);
      cms2.Init(skim_tree);
    }

    // Event Loop
    cms2.Init(tree);
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

      if (mus_dxyPV().size()!=mus_dzPV().size()) {
	cout << "run=" << run_ << " evt=" << evt_ << " mus_dxyPV().size()=" << mus_dxyPV().size() << " mus_dzPV().size()=" << mus_dzPV().size() << endl;
	continue;
      }

      //start selection
      //vertex
      //fixme: check leptons are from this vertex!
      if (firstGoodVertex()!=0) continue;

      //met
      float met = evt_pfmet();
      if (met<30.) continue;

      //fakable objects
      vector<Lep> fobs;
      for (unsigned int elidx=0;elidx<els_p4().size();++elidx) {
	//electron fo selection
	if (isFakableElectron(elidx)==0) continue;
	Lep foel(-1*els_charge().at(elidx)*11,elidx);
	fobs.push_back(foel);
      }
      for (unsigned int muidx=0;muidx<mus_p4().size();++muidx) {
	//muon fo selection
	if (isFakableMuon(muidx)==0) continue;
	Lep fomu(-1*mus_charge().at(muidx)*13,muidx);
	fobs.push_back(fomu);
      }

      if (fobs.size()!=2) continue;
      DilepHyp hyp(fobs[0],fobs[1]);
      if (hyp.p4().mass()<8) continue;

      unsigned int ac_base = analysisCategory(hyp.leadLep(),hyp.traiLep());

      //jets, ht, btags
      int njets = 0;
      int nbtag = 0;
      float ht=0;
      //fixme: should add corrections
      for (unsigned int pfjidx=0;pfjidx<pfjets_p4().size();++pfjidx) {
	if (fabs(pfjets_p4()[pfjidx].eta())>2.4) continue;
	if (isLoosePFJet(pfjidx)==false) continue;
	//add pu jet id pfjets_pileupJetId()>????
	bool isLep = false;
	for (unsigned int fo=0;fo<fobs.size();++fo) {
	  if (deltaR(pfjets_p4()[pfjidx],fobs[fo].p4())<0.4) {
	    isLep =true;
	    break;
	  }
	}
	if (isLep) continue;
	if (pfjets_p4()[pfjidx].pt()>40) {
	  njets++;
	  ht+=pfjets_p4()[pfjidx].pt();
	  if (pfjets_combinedSecondaryVertexBJetTag()[pfjidx]>0.679) nbtag++;
	}
      }

      passesBaselineCuts(njets, nbtag, met, ht, ac_base);
      if (ac_base==0) continue;
      int br = baselineRegion(nbtag);
      unsigned int ac_sig = ac_base;
      passesSignalRegionCuts(ht, ac_sig);
      int sr = ac_sig!=0 ? signalRegion(njets, nbtag, met, ht) : 0;

      //write skim here (only ss)
      if (makeskim) {
	if (hyp.charge()!=0) {
	  cms2.LoadAllBranches();
	  skim_file->cd(); 
	  skim_tree->Fill();
	}
      }

      //leptons
      vector<Lep> goodleps;
      for (unsigned int elidx=0;elidx<els_p4().size();++elidx) {
	//medium electron selection
	if (isGoodElectron(elidx)==0) continue;
	Lep goodel(-1*els_charge().at(elidx)*11,elidx);
	goodleps.push_back(goodel);
      }
      for (unsigned int muidx=0;muidx<mus_p4().size();++muidx) {
	//good muon selection
	if (isGoodMuon(muidx)==0) continue;
	Lep goodmu(-1*mus_charge().at(muidx)*13,muidx);
	goodleps.push_back(goodmu);
      }

      //if (goodleps.size()==0 || goodleps.size()>2) continue;

      //compute fake rate in ss ttbar (should be moved after which selection? ss at least?)
      if (hyp.charge()!=0) {
	for (unsigned int fo=0;fo<fobs.size();++fo) {
	  if (isFromW(fobs[fo])) continue;
	  int pdgid = abs(fobs[fo].pdgId());
	  //denominator
	  makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_den":"fr_el_den"),(pdgid==13?"fr_mu_den":"fr_el_den"),10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()));
	  //numerator
	  for (unsigned int gl=0;gl<goodleps.size();++gl) {
	    if (abs(goodleps[gl].pdgId())==pdgid && goodleps[gl].idx()==fobs[fo].idx()) {
	      //cout << "goodleps.size()=" << goodleps.size() << endl;
	      makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_num":"fr_el_num"),(pdgid==13?"fr_mu_num":"fr_el_num"),10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()));
	      break;
	    }
	  }
	}
      }

      if (goodleps.size()!=2) continue;

      //cout << "Good leptons = " << goodleps.size() << " pdgids=" << goodleps[0].pdgId() << ", " << goodleps[1].pdgId() << endl;
      //cout << "Total charge = " << hyp.charge() << " m=" << hyp.p4().mass() << endl;

      //if (makebaby) FillBabyNtuple();

      //fill histos
      if (makehist) {
	//makeFillHisto1D<TH1F,int>("dummy","dummy",2,0,2,1);
	makeFillHisto1D<TH1F,int>("hyp_charge","hyp_charge",7,-3.5,3.5,hyp.charge());
	makeFillHisto1D<TH1F,float>("hyp_mll","hyp_mll",100,0,1000,hyp.p4().mass());
	makeFillHisto1D<TH1F,float>("hyp_ptll","hyp_ptll",100,0,1000,hyp.p4().pt());
	makeFillHisto1D<TH1F,int>("hyp_njets","hyp_njets",20,0,20,njets);
	makeFillHisto1D<TH1F,float>("hyp_ht","hyp_ht",50,0,2000,ht);
	makeFillHisto1D<TH1F,float>("hyp_met","hyp_met",50,0,500,met);

	int type = -1;
	if (abs(hyp.traiLep().pdgId())==13 && abs(hyp.leadLep().pdgId())==13) type=0;
	if (abs(hyp.traiLep().pdgId())==13 && abs(hyp.leadLep().pdgId())==11) type=1;
	if (abs(hyp.traiLep().pdgId())==11 && abs(hyp.leadLep().pdgId())==13) type=2;
	if (abs(hyp.traiLep().pdgId())==11 && abs(hyp.leadLep().pdgId())==11) type=3;
	makeFillHisto1D<TH1F,int>("hyp_type","hyp_type",5,0,5,type);

	if (hyp.charge()!=0) {
	  //same sign
	  cout << endl << "NEW SS EVENT" << endl << endl;

 	  cout << "lead lep id=" << hyp.leadLep().pdgId() << " p4=" << hyp.leadLep().p4() << " mcid=" << hyp.leadLep().mc_id() << " mother_id=" << hyp.leadLep().mc_motherid() << endl;
	  cout << "trai lep id=" << hyp.traiLep().pdgId() << " p4=" << hyp.traiLep().p4() << " mcid=" << hyp.traiLep().mc_id() << " mother_id=" << hyp.traiLep().mc_motherid() //<< endl;
	       << " mcidx=" << hyp.traiLep().mcidx() << " mc3_motheridx=" << hyp.traiLep().mc3_motheridx() << " mc3_mother_id=" << hyp.traiLep().mc3_motherid() 
	       << " genps_id_mother()[hyp.traiLep().mc3_motheridx()]=" << genps_id_mother()[hyp.traiLep().mc3_motheridx()]
	       << endl;

	  makeFillHisto1D<TH1F,int>("hyp_ss_njets","hyp_ss_njets",20,0,20,njets);
	  makeFillHisto1D<TH1F,float>("hyp_ss_ht","hyp_ss_ht",50,0,2000,ht);
	  makeFillHisto1D<TH1F,float>("hyp_ss_met","hyp_ss_met",50,0,500,met);

	  makeFillHisto1D<TH1F,float>("hyp_ss_mll","hyp_ss_mll",100,0,1000,hyp.p4().mass());
	  makeFillHisto1D<TH1F,float>("hyp_ss_ptll","hyp_ss_ptll",100,0,1000,hyp.p4().pt());
	  makeFillHisto1D<TH1F,int>("hyp_ss_type","hyp_ss_type",5,0,5,type);

	  if (ac_base & 1<<HighPt) {
	    makeFillHisto1D<TH1F,int>("hyp_highpt_sr","hyp_highpt_sr",30,0,30,br);
	    makeFillHisto1D<TH1F,int>("hyp_highpt_sr","hyp_highpt_sr",30,0,30,sr);
	  }
	  if (ac_base & 1<<LowPt) {
	    makeFillHisto1D<TH1F,int>("hyp_lowpt_sr","hyp_lowpt_sr",30,0,30,br);
	    makeFillHisto1D<TH1F,int>("hyp_lowpt_sr","hyp_lowpt_sr",30,0,30,sr);
	  }
	  if (ac_base & 1<<VeryLowPt) {
	    makeFillHisto1D<TH1F,int>("hyp_verylowpt_sr","hyp_verylowpt_sr",30,0,30,br);
	    makeFillHisto1D<TH1F,int>("hyp_verylowpt_sr","hyp_verylowpt_sr",30,0,30,sr);
	  }
	  //check fakes (mother not W)
	  if ( !isFromW(hyp.traiLep()) ) {
	    if (abs(hyp.traiLep().pdgId())==13) {
	      //trailing muon
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_mu_mc","hyp_ss_trail_mu_mc",11001,-5500.5,5500.5,hyp.traiLep().mc_id());
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_mu_mc_mother","hyp_ss_trail_mu_mc_mother",11001,-5500.5,5500.5,hyp.traiLep().mc_motherid());
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_mu_mc3","hyp_ss_trail_mu_mc3",11001,-5500.5,5500.5,hyp.traiLep().mc3_id());
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_mu_mc3_mother","hyp_ss_trail_mu_mc3_mother",11001,-5500.5,5500.5,hyp.traiLep().mc3_motherid());
	    }
	    if (abs(hyp.traiLep().pdgId())==11) {
	      //trailing elec
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_el_mc","hyp_ss_trail_el_mc",11001,-5500.5,5500.5,hyp.traiLep().mc_id());
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_el_mc_mother","hyp_ss_trail_el_mc_mother",11001,-5500.5,5500.5,hyp.traiLep().mc_motherid());
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_el_mc3","hyp_ss_trail_el_mc3",11001,-5500.5,5500.5,hyp.traiLep().mc3_id());
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_el_mc3_mother","hyp_ss_trail_el_mc3_mother",11001,-5500.5,5500.5,hyp.traiLep().mc3_motherid());
	    }
	  }
	  if ( !isFromW(hyp.leadLep()) ) {
	    if (abs(hyp.leadLep().pdgId())==13) {
	      //leading muon
	      makeFillHisto1D<TH1F,int>("hyp_ss_lead_mu_mc","hyp_ss_lead_mu_mc",11001,-5500.5,5500.5,hyp.leadLep().mc_id());
	      makeFillHisto1D<TH1F,int>("hyp_ss_lead_mu_mc_mother","hyp_ss_lead_mu_mc_mother",11001,-5500.5,5500.5,hyp.leadLep().mc_motherid());
	      makeFillHisto1D<TH1F,int>("hyp_ss_lead_mu_mc3","hyp_ss_lead_mu_mc3",11001,-5500.5,5500.5,hyp.leadLep().mc3_id());
	      makeFillHisto1D<TH1F,int>("hyp_ss_lead_mu_mc3_mother","hyp_ss_lead_mu_mc3_mother",11001,-5500.5,5500.5,hyp.leadLep().mc3_motherid());
	    }
	    if (abs(hyp.leadLep().pdgId())==11) {
	      //leading elec
	      makeFillHisto1D<TH1F,int>("hyp_ss_lead_el_mc","hyp_ss_lead_el_mc",11001,-5500.5,5500.5,hyp.leadLep().mc_id());
	      makeFillHisto1D<TH1F,int>("hyp_ss_lead_el_mc_mother","hyp_ss_lead_el_mc_mother",11001,-5500.5,5500.5,hyp.leadLep().mc_motherid());
	      makeFillHisto1D<TH1F,int>("hyp_ss_lead_el_mc3","hyp_ss_lead_el_mc3",11001,-5500.5,5500.5,hyp.leadLep().mc3_id());
	      makeFillHisto1D<TH1F,int>("hyp_ss_lead_el_mc3_mother","hyp_ss_lead_el_mc3_mother",11001,-5500.5,5500.5,hyp.leadLep().mc3_motherid());
	    }
	  }

	  //dump genps
	  for (unsigned int gp_idx=0;gp_idx<genps_p4().size();++gp_idx) {
	    cout << "gp idx=" << gp_idx << " id=" << genps_id()[gp_idx] 
		 << " status=" << genps_status()[gp_idx]
		 << " id_mother=" << genps_id_mother()[gp_idx]
		 << " p4=" << genps_p4()[gp_idx]
		 << endl;    
	  }


	} else {
	  //opposite sign
	  makeFillHisto1D<TH1F,float>("hyp_os_mll","hyp_os_mll",100,0,1000,hyp.p4().mass());
	  makeFillHisto1D<TH1F,float>("hyp_os_ptll","hyp_os_ptll",100,0,1000,hyp.p4().pt());
	  makeFillHisto1D<TH1F,int>("hyp_os_type","hyp_os_type",5,0,5,type);
	}
      }//if (makehist)

    }

    if (makeskim) {
      skim_file->cd(); 
      skim_tree->Write(); 
      skim_file->Close();
      delete skim_file;
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
