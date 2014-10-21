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

//fixme: put WF and FSR in different categories
enum LeptonCategories { Prompt = 0, PromptWS = 1, PromptWF = 2, PromptFSR = 2, 
			FakeLightTrue = 3, FakeC = 4, FakeB = 5, FakeLightFake = 6, FakeHiPtGamma = 7, 
			FakeUnknown = 8, FakeLowPtGamma = 9, All9999 = 10,
			Other = 11, End = 12};

int looper::ScanChain( TChain* chain, TString prefix, TString postfix, bool isData, TString whatTest, int nEvents) {

  makebaby       = 0;
  makehist       = 1;
  maketext       = 0;

  makeQCDtest    = 0;
  makeDYtest     = 0;
  makeSSskim     = 0;
  makeQCDskim    = 0;
  if (whatTest=="QCDtest") makeQCDtest    = 1;
  if (whatTest=="DYtest" ) makeDYtest     = 1;
  if (whatTest=="SSskim" ) makeSSskim     = 1;
  if (whatTest=="QCDskim") makeQCDskim    = 1;

  bool debug = 0;  

  if (postfix!="") postfix = "_"+postfix;
  if (makebaby) MakeBabyNtuple( Form( "%s_baby%s.root", prefix.Data(), postfix.Data() ) );
  if (makehist) CreateOutputFile( Form( "%s_histos%s.root", prefix.Data(), postfix.Data() ) );

  TFile* fr_file=TFile::Open("fakeRates_qcd_pt-50to170.root");

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

    cout << "processing file: " << currentFile->GetTitle() << endl;

    //Skimmed output file - needs to be before cms2.Init(tree)
    TFile *skim_file = 0;
    TTree* skim_tree = 0;
    if (makeSSskim || makeQCDskim) {
      TString skim_file_name = TString(currentFile->GetTitle());
      if (makeSSskim) {
	skim_file_name.ReplaceAll(".root","_skimSS.root");
	skim_file_name.Remove(0,skim_file_name.Last('/'));
	skim_file_name.Prepend('.');
      } else if (makeQCDskim) {
	skim_file_name.ReplaceAll(".root","_skimQCD.root");
	skim_file_name.Remove(0,skim_file_name.First('Q'));
	skim_file_name.Remove(skim_file_name.First('T')-1,skim_file_name.Last('/')-skim_file_name.First('T')+1);
      }
      skim_file = new TFile(skim_file_name,"recreate");
      //locally
      skim_tree = (TTree*) tree->CloneTree(0, "fast");
      skim_tree->SetDirectory(skim_file);
      cms2.Init(skim_tree);
    }

    // Event Loop
    cms2.Init(tree);
    unsigned int nEvents = tree->GetEntries();
    bool newfile = true;

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

      float lumi = 5.;

      //fill baby
      run_   = evt_run();
      ls_    = evt_lumiBlock();
      evt_   = evt_event();
      weight_ = isData ? 1. : lumi*evt_scale1fb();
      if (prefix=="wj") weight_*=30;//fixme using only a subset of events
      if (prefix=="dy") weight_*=25;//fixme using only a subset of events

      if (newfile) {
	cout << "weight=" << weight_ << endl;
	newfile = false;
      }

      if (debug) cout << "file=" << currentFile->GetTitle() << " run=" << run_ << " evt=" << evt_ << endl;

      if (mus_dxyPV().size()!=mus_dzPV().size()) {
	cout << "run=" << run_ << " evt=" << evt_ << " mus_dxyPV().size()=" << mus_dxyPV().size() << " mus_dzPV().size()=" << mus_dzPV().size() << endl;
	continue;
      }


      //start selection
      //vertex
      //fixme: check leptons are from this vertex!
      if (firstGoodVertex()!=0) continue;

      //met
      if (debug) cout << "met" << endl;
      float met = evt_pfmet();
      if (met<30. && !makeQCDskim && !makeQCDtest && !makeDYtest) continue;

      //fakable objects
      if (debug) cout << "fobs" << endl;
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
      if (fobs.size()==0) continue;

      //write skim here (only qcd)
      if (makeQCDskim) {
	if (debug) cout << "qcd skim" << endl;
	if (fobs.size()!=0) {
	  cms2.LoadAllBranches();
	  skim_file->cd(); 
	  skim_tree->Fill();
	  continue;
	}
      }

      //leptons
      if (debug) cout << "goodleps" << endl;
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

      //jets, ht, btags
      if (debug) cout << "jets" << endl;
      int njets = 0;
      int nbtag = 0;
      float ht=0;
      //fixme: should add corrections
      float drcut = 0.4;
      if (makeQCDtest) drcut = 1.0;
      for (unsigned int pfjidx=0;pfjidx<pfjets_p4().size();++pfjidx) {
	if (fabs(pfjets_p4()[pfjidx].eta())>2.4) continue;
	if (isLoosePFJet(pfjidx)==false) continue;
	//add pu jet id pfjets_pileupJetId()>????
	bool isLep = false;
	//fixme: should be checked agains fo or good leptons?
	for (unsigned int fo=0;fo<fobs.size();++fo) {
	  if (deltaR(pfjets_p4()[pfjidx],fobs[fo].p4())<drcut) {
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

      if (makeQCDtest) {
	if (debug) cout << "qcdtest" << endl;
	if (fobs.size()>1)continue; 
	if (fobs[0].pt()<20) continue; 
	makeFillHisto1D<TH1F,int>("njets","njets",20,0,20,njets,weight_);
	for (unsigned int ipu=0;ipu<puInfo_bunchCrossing().size();ipu++)
	  if (puInfo_bunchCrossing()[ipu]==0) 
	    makeFillHisto1D<TH1F,float>("npu_true","npu_true",20,0,100,puInfo_trueNumInteractions()[ipu],weight_);
	if (njets==0) continue;
	float evtmt = mt(fobs[0].pt(),met,deltaPhi(fobs[0].p4().phi(),evt_pfmetPhi()));
	float genmt = mt(fobs[0].pt(),gen_met(),deltaPhi(fobs[0].p4().phi(),gen_metPhi()));
	float tkmet = trackerMET(0.2).met;
	float minmet = std::min(tkmet,met);
	makeFillHisto1D<TH1F,float>("fo_pt","fo_pt",20,0,100,fobs[0].pt(),weight_);
	makeFillHisto1D<TH1F,float>("fo_eta","fo_eta",10,-2.5,2.5,fobs[0].eta(),weight_);
	makeFillHisto1D<TH1F,float>("evt_met","evt_met",10,0,100,met,weight_);
	makeFillHisto1D<TH1F,float>("evt_tkmet","evt_tkmet",10,0,100,tkmet,weight_);
	makeFillHisto1D<TH1F,float>("evt_minmet","evt_minmet",10,0,100,minmet,weight_);
	makeFillHisto1D<TH1F,float>("evt_mt","evt_mt",20,0,200,evtmt,weight_);
	makeFillHisto2D<TH2F,float>("evt_mt_vs_met","evt_mt_vs_met",10,0,100,met,20,0,200,evtmt,weight_);
	makeFillHisto2D<TH2F,float>("evt_mt_vs_pt", "evt_mt_vs_pt",20,0,100,fobs[0].pt(),20,0,200,evtmt,weight_);
	makeFillHisto1D<TH1F,float>("gen_met","gen_met",10,0,100,gen_met(),weight_);
	makeFillHisto1D<TH1F,float>("gen_mt","gen_mt",20,0,200,genmt,weight_);
	makeFillHisto2D<TH2F,float>("gen_mt_vs_met","gen_mt_vs_met",10,0,100,gen_met(),20,0,200,genmt,weight_);
	if (met<20 && evtmt<20) makeFillHisto1D<TH1F,int>("pass_metmt","pass_metmt",2,0,2,1,weight_);
	else makeFillHisto1D<TH1F,int>("pass_metmt","pass_metmt",2,0,2,0,weight_);
 	if (minmet<20 && evtmt<20) makeFillHisto1D<TH1F,int>("pass_minmetmt","pass_minmetmt",2,0,2,1,weight_);
	else makeFillHisto1D<TH1F,int>("pass_minmetmt","pass_minmetmt",2,0,2,0,weight_);
	if (gen_met()<20 && genmt<20) makeFillHisto1D<TH1F,int>("pass_gen_metmt","pass_gen_metmt",2,0,2,1,weight_);
	else makeFillHisto1D<TH1F,int>("pass_gen_metmt","pass_gen_metmt",2,0,2,0,weight_);
	//compute fake rate
	//reco level selection
	if (met<20. && evtmt<20.) {
	  for (unsigned int fo=0;fo<fobs.size();++fo) {
	    //if (isFromW(fobs[fo])) continue;//fixme
	    //if (!isFromB(fobs[fo])) continue;//fixme
	    int pdgid = abs(fobs[fo].pdgId());
	    //denominator
	    makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_den":"fr_el_den"),(pdgid==13?"fr_mu_den":"fr_el_den"),
					10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),weight_);
	    makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_den_now":"fr_el_den_now"),(pdgid==13?"fr_mu_den_now":"fr_el_den_now"),
					10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),1.);
	    //numerator
	    for (unsigned int gl=0;gl<goodleps.size();++gl) {
	      if (abs(goodleps[gl].pdgId())==pdgid && goodleps[gl].idx()==fobs[fo].idx()) {
		//cout << "goodleps.size()=" << goodleps.size() << endl;
		makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_num":"fr_el_num"),(pdgid==13?"fr_mu_num":"fr_el_num"),
					    10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),weight_);
		makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_num_now":"fr_el_num_now"),(pdgid==13?"fr_mu_num_now":"fr_el_num_now"),
					    10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),1.);
		//std::cout << "scale1fb=" << evt_scale1fb() << endl;
		break;
	      }
	    }
	  }
	}
	//reco level selection (minmet)
	if (minmet<20. && evtmt<20.) {
	  for (unsigned int fo=0;fo<fobs.size();++fo) {
	    //if (isFromW(fobs[fo])) continue;//fixme
	    //if (!isFromB(fobs[fo])) continue;//fixme
	    int pdgid = abs(fobs[fo].pdgId());
	    //denominator
	    makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_den_min":"fr_el_den_min"),(pdgid==13?"fr_mu_den_min":"fr_el_den_min"),
					10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),weight_);
	    makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_den_min_now":"fr_el_den_min_now"),(pdgid==13?"fr_mu_den_min_now":"fr_el_den_min_now"),
					10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),1.);
	    //numerator
	    for (unsigned int gl=0;gl<goodleps.size();++gl) {
	      if (abs(goodleps[gl].pdgId())==pdgid && goodleps[gl].idx()==fobs[fo].idx()) {
		//cout << "goodleps.size()=" << goodleps.size() << endl;
		makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_num_min":"fr_el_num_min"),(pdgid==13?"fr_mu_num_min":"fr_el_num_min"),
					    10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),weight_);
		makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_num_min_now":"fr_el_num_min_now"),(pdgid==13?"fr_mu_num_min_now":"fr_el_num_min_now"),
					    10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),1.);
		//std::cout << "scale1fb=" << evt_scale1fb() << endl;
		break;
	      }
	    }
	  }
	}
	//gen level selection
	if (gen_met()<20. && genmt<20.) {
	  for (unsigned int fo=0;fo<fobs.size();++fo) {
	    //if (isFromW(fobs[fo])) continue;//fixme
	    //if (!isFromB(fobs[fo])) continue;//fixme
	    int pdgid = abs(fobs[fo].pdgId());
	    //denominator
	    makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_den_gen":"fr_el_den_gen"),(pdgid==13?"fr_mu_den_gen":"fr_el_den_gen"),
					10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),weight_);
	    makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_den_gen_now":"fr_el_den_gen_now"),(pdgid==13?"fr_mu_den_gen_now":"fr_el_den_gen_now"),
					10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),1.);
	    //numerator
	    for (unsigned int gl=0;gl<goodleps.size();++gl) {
	      if (abs(goodleps[gl].pdgId())==pdgid && goodleps[gl].idx()==fobs[fo].idx()) {
		//cout << "goodleps.size()=" << goodleps.size() << endl;
		makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_num_gen":"fr_el_num_gen"),(pdgid==13?"fr_mu_num_gen":"fr_el_num_gen"),
					    10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),weight_);
		makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_num_gen_now":"fr_el_num_gen_now"),(pdgid==13?"fr_mu_num_gen_now":"fr_el_num_gen_now"),
					    10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),1.);
		//std::cout << "scale1fb=" << evt_scale1fb() << endl;
		break;
	      }
	    }
	  }
	}
	continue;
      }

      if (makeDYtest) {
	if (debug) cout << "dytest" << endl;
	//check lepton selection efficiency
	for (unsigned int gp=0;gp<genps_id().size();++gp) {
	  int pdgid = abs(genps_id()[gp]);
	  if (pdgid!=13 && pdgid!=11) continue;
	  if (genps_id_mother()[gp]!=23) continue;
	  if (genps_status()[gp]!=1) continue;//is this needed?
	  if (genps_p4()[gp].eta()>2.5) continue;
	  if (genps_p4()[gp].pt()<5) continue;
	  if (pdgid==11 && genps_p4()[gp].pt()<10) continue;
	  //denominator
	  makeFillHisto2D<TH2F,float>((pdgid==13?"ef_mu_den":"ef_el_den"),(pdgid==13?"ef_mu_den":"ef_el_den"),10,0.,50.,genps_p4()[gp].pt(),5,0.,2.5,fabs(genps_p4()[gp].eta()),weight_);
	  //numerator
	  for (unsigned int gl=0;gl<goodleps.size();++gl) {
	    if (abs(goodleps[gl].pdgId())==pdgid && deltaR(goodleps[gl].p4(),genps_p4()[gp])<0.1) {
	      //cout << "goodleps.size()=" << goodleps.size() << endl;
	      makeFillHisto2D<TH2F,float>((pdgid==13?"ef_mu_num":"ef_el_num"),(pdgid==13?"ef_mu_num":"ef_el_num"),10,0.,50.,genps_p4()[gp].pt(),5,0.,2.5,fabs(genps_p4()[gp].eta()),weight_);
	      //charge flip
	      if (goodleps[gl].pdgId()==-genps_id()[gp]) makeFillHisto2D<TH2F,float>((pdgid==13?"flip_mu":"flip_el"),(pdgid==13?"flip_mu":"flip_el"),10,0.,50.,genps_p4()[gp].pt(),5,0.,2.5,fabs(genps_p4()[gp].eta()),weight_);
	      break;
	    }
	  }
	}
	continue;
      }

      if (fobs.size()!=2) continue;
      DilepHyp hyp(fobs[0],fobs[1]);
      if (hyp.p4().mass()<8) continue;

      unsigned int ac_base = analysisCategory(hyp.leadLep(),hyp.traiLep());

      passesBaselineCuts(njets, nbtag, met, ht, ac_base);
      if (ac_base==0) continue;
      int br = baselineRegion(nbtag);
      unsigned int ac_sig = ac_base;
      passesSignalRegionCuts(ht, ac_sig);
      int sr = ac_sig!=0 ? signalRegion(njets, nbtag, met, ht) : 0;

      //write skim here (only ss)
      if (debug) cout << "ss skim" << endl;
      if (makeSSskim) {
	if (hyp.charge()!=0) {
	  cms2.LoadAllBranches();
	  skim_file->cd(); 
	  skim_tree->Fill();
	}
      }

      //compute fake rate in ss ttbar
      if (hyp.charge()!=0) {
	for (unsigned int fo=0;fo<fobs.size();++fo) {
	  if (isFromW(fobs[fo])) continue;
	  int pdgid = abs(fobs[fo].pdgId());
	  //denominator
	  makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_den":"fr_el_den"),(pdgid==13?"fr_mu_den":"fr_el_den"),10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),weight_);
	  //numerator
	  bool isNumerator = false;
	  for (unsigned int gl=0;gl<goodleps.size();++gl) {
	    if (abs(goodleps[gl].pdgId())==pdgid && goodleps[gl].idx()==fobs[fo].idx()) {
	      //cout << "goodleps.size()=" << goodleps.size() << endl;
	      makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_num":"fr_el_num"),(pdgid==13?"fr_mu_num":"fr_el_num"),10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),weight_);
	      isNumerator = true;
	      break;
	    }
	  }
	  if (!isNumerator) {
	    TH2F* fr_h = (TH2F*) fr_file->Get((pdgid==13?"fr_mu_gen":"fr_el_gen"));
	    float maxPt=fr_h->GetXaxis()->GetBinUpEdge(fr_h->GetXaxis()->GetNbins())-0.01; 
	    float maxEta=fr_h->GetYaxis()->GetBinUpEdge(fr_h->GetYaxis()->GetNbins())-0.01;
	    float pt = fobs[fo].pt();
	    float eta = fabs(fobs[fo].eta());
	    if (pt>maxPt) pt=maxPt;
	    if (eta>maxEta) eta=maxEta;
	    float fr = fr_h->GetBinContent(fr_h->FindBin(pt,eta));
	    float fre = fr_h->GetBinError(fr_h->FindBin(pt,eta));
	    float frW = fr/(1.-fr);
	    std::cout << "fake rate=" << fr << " +/- " << fre << std::endl;
	    makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_close":"fr_el_close"),(pdgid==13?"fr_mu_close":"fr_el_close"),10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),weight_*frW);
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
	makeFillHisto1D<TH1F,int>("hyp_charge","hyp_charge",7,-3.5,3.5,hyp.charge(),weight_);
	makeFillHisto1D<TH1F,float>("hyp_mll","hyp_mll",100,0,1000,hyp.p4().mass(),weight_);
	makeFillHisto1D<TH1F,float>("hyp_ptll","hyp_ptll",100,0,1000,hyp.p4().pt(),weight_);
	makeFillHisto1D<TH1F,int>("hyp_njets","hyp_njets",20,0,20,njets,weight_);
	makeFillHisto1D<TH1F,int>("hyp_nbtag","hyp_nbtag",20,0,20,nbtag,weight_);
	makeFillHisto1D<TH1F,float>("hyp_ht","hyp_ht",50,0,2000,ht,weight_);
	makeFillHisto1D<TH1F,float>("hyp_met","hyp_met",50,0,500,met,weight_);

	int type = -1;
	if (abs(hyp.traiLep().pdgId())==13 && abs(hyp.leadLep().pdgId())==13) type=0;
	if (abs(hyp.traiLep().pdgId())==13 && abs(hyp.leadLep().pdgId())==11) type=1;
	if (abs(hyp.traiLep().pdgId())==11 && abs(hyp.leadLep().pdgId())==13) type=2;
	if (abs(hyp.traiLep().pdgId())==11 && abs(hyp.leadLep().pdgId())==11) type=3;
	makeFillHisto1D<TH1F,int>("hyp_type","hyp_type",5,0,5,type,weight_);

	TString ll = abs(hyp.leadLep().pdgId())==13 ? "mu" : "el";
	TString lt = abs(hyp.traiLep().pdgId())==13 ? "mu" : "el";

	if (hyp.charge()!=0) {
	  //same sign
	  cout << endl << "NEW SS EVENT" << endl << endl;

	  /*
 	  cout << "lead lep id=" << hyp.leadLep().pdgId() << " p4=" << hyp.leadLep().p4() 
	       << " mcid=" << hyp.leadLep().mc_id() << " mcp4=" << hyp.leadLep().mc_p4() << " mother_id=" << hyp.leadLep().mc_motherid()
	       << " mc3idx=" << hyp.leadLep().mc3idx() << " mc3_id=" << hyp.leadLep().mc3_id() 
	       << " mc3_motheridx=" << hyp.leadLep().mc3_motheridx() << " mc3_mother_id=" << hyp.leadLep().mc3_motherid()
	       << " genps_id_mother()[hyp.leadLep().mc3_motheridx()]=" << genps_id_mother()[hyp.leadLep().mc3_motheridx()]
	       << endl;
	  cout << "trai lep id=" << hyp.traiLep().pdgId() << " p4=" << hyp.traiLep().p4() 
	       << " mcid=" << hyp.traiLep().mc_id() << " mcp4=" << hyp.traiLep().mc_p4()  << " mother_id=" << hyp.traiLep().mc_motherid()
	       << " mc3idx=" << hyp.traiLep().mc3idx() << " mc3_id=" << hyp.traiLep().mc3_id() 
	       << " mc3_motheridx=" << hyp.traiLep().mc3_motheridx() << " mc3_mother_id=" << hyp.traiLep().mc3_motherid()
	       << " genps_id_mother()[hyp.traiLep().mc3_motheridx()]=" << genps_id_mother()[hyp.traiLep().mc3_motheridx()]
	       << endl;
	  */

	  makeFillHisto1D<TH1F,int>("hyp_ss_njets","hyp_ss_njets",20,0,20,njets,weight_);
	  makeFillHisto1D<TH1F,int>("hyp_ss_nbtag","hyp_ss_nbtag",20,0,20,nbtag,weight_);
	  makeFillHisto1D<TH1F,float>("hyp_ss_ht","hyp_ss_ht",50,0,2000,ht,weight_);
	  makeFillHisto1D<TH1F,float>("hyp_ss_met","hyp_ss_met",50,0,500,met,weight_);

	  makeFillHisto1D<TH1F,float>("hyp_ss_mll","hyp_ss_mll",100,0,1000,hyp.p4().mass(),weight_);
	  makeFillHisto1D<TH1F,float>("hyp_ss_ptll","hyp_ss_ptll",100,0,1000,hyp.p4().pt(),weight_);
	  makeFillHisto1D<TH1F,int>("hyp_ss_type","hyp_ss_type",5,0,5,type,weight_);

	  if (ac_base & 1<<HighPt) {
	    makeFillHisto1D<TH1F,int>("hyp_highpt_sr","hyp_highpt_sr",30,0,30,br,weight_);
	    makeFillHisto1D<TH1F,int>("hyp_highpt_sr","hyp_highpt_sr",30,0,30,sr,weight_);
	  }
	  if (ac_base & 1<<LowPt) {
	    makeFillHisto1D<TH1F,int>("hyp_lowpt_sr","hyp_lowpt_sr",30,0,30,br,weight_);
	    makeFillHisto1D<TH1F,int>("hyp_lowpt_sr","hyp_lowpt_sr",30,0,30,sr,weight_);
	  }
	  if (ac_base & 1<<VeryLowPt) {
	    makeFillHisto1D<TH1F,int>("hyp_verylowpt_sr","hyp_verylowpt_sr",30,0,30,br,weight_);
	    makeFillHisto1D<TH1F,int>("hyp_verylowpt_sr","hyp_verylowpt_sr",30,0,30,sr,weight_);
	  }

	  bool isLeadPrompt = 0;
	  bool isTrailPrompt = 0;
	  int leadType = -1;
	  int trailType = -1;

	  if ( !isFromW(hyp.traiLep()) ) {
	    //check fakes (mother not W)

	    makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_mc","hyp_ss_trail_"+lt+"_mc",11001,-5500.5,5500.5,hyp.traiLep().mc_id(),weight_);
	    makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_mc_mother","hyp_ss_trail_"+lt+"_mc_mother",11001,-5500.5,5500.5,hyp.traiLep().mc_motherid(),weight_);
	    makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_mc3","hyp_ss_trail_"+lt+"_mc3",11001,-5500.5,5500.5,hyp.traiLep().mc3_id(),weight_);
	    makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_mc3_mother","hyp_ss_trail_"+lt+"_mc3_mother",11001,-5500.5,5500.5,hyp.traiLep().mc3_motherid(),weight_);

	    if (isFromB(hyp.traiLep()) ) trailType=FakeB;
	    else if (isFromC(hyp.traiLep()) ) trailType=FakeC;
	    else if (isFromLight(hyp.traiLep())) trailType=FakeLightTrue;
	    else if (isFromLightFake(hyp.traiLep())) trailType=FakeLightFake;
	    else if (hyp.traiLep().mc_id()==22 && hyp.traiLep().mc_p4().pt()>hyp.traiLep().pt()/2.) {
	      trailType=FakeHiPtGamma;
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_hiptgamma_mc","hyp_ss_trail_"+lt+"_hiptgamma_mc",11001,-5500.5,5500.5,hyp.traiLep().mc_id(),weight_);
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_hiptgamma_mc_mother","hyp_ss_trail_"+lt+"_hiptgamma_mc_mother",11001,-5500.5,5500.5,hyp.traiLep().mc_motherid(),weight_);
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_hiptgamma_mc3","hyp_ss_trail_"+lt+"_hiptgamma_mc3",11001,-5500.5,5500.5,hyp.traiLep().mc3_id(),weight_);
	      makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_hiptgamma_mc3_mother","hyp_ss_trail_"+lt+"_hiptgamma_mc3_mother",11001,-5500.5,5500.5,hyp.traiLep().mc3_motherid(),weight_);
	      makeFillHisto1D<TH1F,float>("hyp_ss_trail_"+lt+"_hiptgamma_pt","hyp_ss_trail_"+lt+"_hiptgamma_pt",50,0,100,hyp.traiLep().pt(),weight_);
	      makeFillHisto1D<TH1F,float>("hyp_ss_trail_"+lt+"_hiptgamma_ptmc","hyp_ss_trail_"+lt+"_hiptgamma_ptmc",50,0,100,hyp.traiLep().mc_p4().pt(),weight_);
	      makeFillHisto1D<TH1F,float>("hyp_ss_trail_"+lt+"_hiptgamma_deltaptoverpt","hyp_ss_trail_"+lt+"_hiptgamma_deltaptoverpt",100,-0.5,0.5,(hyp.traiLep().pt()-hyp.traiLep().mc_p4().pt())/hyp.traiLep().mc_p4().pt(),weight_);
	    } else  {
	      trailType=FakeUnknown;
	      if (hyp.traiLep().mc_id()==22 && hyp.traiLep().mc_p4().pt()<1.) {
		makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_category","hyp_ss_trail_"+lt+"_category",End,0,End,FakeLowPtGamma,weight_);
		makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_lowptgamma_mc","hyp_ss_trail_"+lt+"_lowptgamma_mc",11001,-5500.5,5500.5,hyp.traiLep().mc_id(),weight_);
		makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_lowptgamma_mc_mother","hyp_ss_trail_"+lt+"_lowptgamma_mc_mother",11001,-5500.5,5500.5,hyp.traiLep().mc_motherid(),weight_);
		makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_lowptgamma_mc3","hyp_ss_trail_"+lt+"_lowptgamma_mc3",11001,-5500.5,5500.5,hyp.traiLep().mc3_id(),weight_);
		makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_lowptgamma_mc3_mother","hyp_ss_trail_"+lt+"_lowptgamma_mc3_mother",11001,-5500.5,5500.5,hyp.traiLep().mc3_motherid(),weight_);
		makeFillHisto1D<TH1F,float>("hyp_ss_trail_"+lt+"_category_lowptgamma_pt","hyp_ss_trail_"+lt+"_category_lowptgamma_pt",50,0,100,hyp.traiLep().pt(),weight_);
	      } else if (hyp.traiLep().mc_id()==-9999 && hyp.traiLep().mc_motherid()==-9999 && hyp.traiLep().mc3_id()==-9999 && hyp.traiLep().mc3_motherid()==-9999) {
		makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_category","hyp_ss_trail_"+lt+"_category",End,0,End,All9999,weight_);
		makeFillHisto1D<TH1F,float>("hyp_ss_trail_"+lt+"_category_all9999_pt","hyp_ss_trail_"+lt+"_category_all9999_pt",50,0,100,hyp.traiLep().pt(),weight_);
	      }
	      else {
		makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_category","hyp_ss_trail_"+lt+"_category",End,0,End,Other,weight_);
		makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_other_mc","hyp_ss_trail_"+lt+"_other_mc",11001,-5500.5,5500.5,hyp.traiLep().mc_id(),weight_);
		makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_other_mc_mother","hyp_ss_trail_"+lt+"_other_mc_mother",11001,-5500.5,5500.5,hyp.traiLep().mc_motherid(),weight_);
		makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_other_mc3","hyp_ss_trail_"+lt+"_other_mc3",11001,-5500.5,5500.5,hyp.traiLep().mc3_id(),weight_);
		makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_other_mc3_mother","hyp_ss_trail_"+lt+"_other_mc3_mother",11001,-5500.5,5500.5,hyp.traiLep().mc3_motherid(),weight_);
		makeFillHisto1D<TH1F,float>("hyp_ss_trail_"+lt+"_category_other_pt","hyp_ss_trail_"+lt+"_category_other_pt",50,0,100,hyp.traiLep().pt(),weight_);
		/*
		cout << "UNKNOWN FAKE LEPTON " << lt << endl;
		cout << "trai lep id=" << hyp.traiLep().pdgId() << " p4=" << hyp.traiLep().p4() 
		     << " mcid=" << hyp.traiLep().mc_id() << " mcp4=" << hyp.traiLep().mc_p4()  << " mother_id=" << hyp.traiLep().mc_motherid()
		     << " mc3idx=" << hyp.traiLep().mc3idx() << " mc3_id=" << hyp.traiLep().mc3_id() 
		     << " mc3_motheridx=" << hyp.traiLep().mc3_motheridx() << " mc3_mother_id=" << hyp.traiLep().mc3_motherid()
		     << " genps_id_mother()[hyp.traiLep().mc3_motheridx()]=" << genps_id_mother()[hyp.traiLep().mc3_motheridx()]
		     << endl;
		*/
	      }
	    }

	  } else {

	    isTrailPrompt = 1;
	    if ( hyp.traiLep().pdgId()==hyp.traiLep().mc_id() ) trailType=Prompt;
	    else if ( abs(hyp.traiLep().pdgId())==abs(hyp.traiLep().mc_id()) ) trailType=PromptWS;
	    else if ( abs(hyp.traiLep().pdgId())==11 ? abs(hyp.traiLep().mc_id())==13 : abs(hyp.traiLep().mc_id())==11 ) trailType=PromptWF;
	    else if ( hyp.traiLep().mc_id()==22) trailType=PromptFSR;
	    else cout << "UNKNOWN PROMPT LEPTON " << lt << endl;

	  }
	  if ( !isFromW(hyp.leadLep()) ) {

	    makeFillHisto1D<TH1F,int>("hyp_ss_lead_"+ll+"_mc","hyp_ss_lead_"+ll+"_mc",11001,-5500.5,5500.5,hyp.leadLep().mc_id(),weight_);
	    makeFillHisto1D<TH1F,int>("hyp_ss_lead_"+ll+"_mc_mother","hyp_ss_lead_"+ll+"_mc_mother",11001,-5500.5,5500.5,hyp.leadLep().mc_motherid(),weight_);
	    makeFillHisto1D<TH1F,int>("hyp_ss_lead_"+ll+"_mc3","hyp_ss_lead_"+ll+"_mc3",11001,-5500.5,5500.5,hyp.leadLep().mc3_id(),weight_);
	    makeFillHisto1D<TH1F,int>("hyp_ss_lead_"+ll+"_mc3_mother","hyp_ss_lead_"+ll+"_mc3_mother",11001,-5500.5,5500.5,hyp.leadLep().mc3_motherid(),weight_);

	    if (isFromB(hyp.leadLep()) ) leadType=FakeB;
	    else if (isFromC(hyp.leadLep()) ) leadType=FakeC;
	    else if (isFromLight(hyp.leadLep())) leadType=FakeLightTrue;
	    else if (isFromLightFake(hyp.leadLep())) leadType=FakeLightFake;
	    else if (hyp.leadLep().mc_id()==22 && hyp.leadLep().mc_p4().pt()>hyp.leadLep().pt()/2.) leadType=FakeHiPtGamma;
	    else {
	      leadType=FakeUnknown;
	      if (hyp.leadLep().mc_id()==22 && hyp.leadLep().mc_p4().pt()<1.) 
		makeFillHisto1D<TH1F,int>("hyp_ss_lead_"+ll+"_category","hyp_ss_lead_"+ll+"_category",End,0,End,FakeLowPtGamma,weight_);
	      else if (hyp.leadLep().mc_id()==-9999 && hyp.leadLep().mc_motherid()==-9999 && hyp.leadLep().mc3_id()==-9999 && hyp.leadLep().mc3_motherid()==-9999)
		makeFillHisto1D<TH1F,int>("hyp_ss_lead_"+ll+"_category","hyp_ss_lead_"+ll+"_category",End,0,End,All9999,weight_);
	      else {
		makeFillHisto1D<TH1F,int>("hyp_ss_lead_"+ll+"_category","hyp_ss_lead_"+ll+"_category",End,0,End,Other,weight_);
		/*
		cout << "UNKNOWN FAKE LEPTON " << ll << endl;
		cout << "lead lep id=" << hyp.leadLep().pdgId() << " p4=" << hyp.leadLep().p4() 
		     << " mcid=" << hyp.leadLep().mc_id() << " mcp4=" << hyp.leadLep().mc_p4() << " mother_id=" << hyp.leadLep().mc_motherid()
		     << " mc3idx=" << hyp.leadLep().mc3idx() << " mc3_id=" << hyp.leadLep().mc3_id() 
		     << " mc3_motheridx=" << hyp.leadLep().mc3_motheridx() << " mc3_mother_id=" << hyp.leadLep().mc3_motherid()
		     << " genps_id_mother()[hyp.leadLep().mc3_motheridx()]=" << genps_id_mother()[hyp.leadLep().mc3_motheridx()]
		     << endl;
		*/
	      }
	    }

	  } else {

	    isLeadPrompt = 1;
	    if ( hyp.leadLep().pdgId()==hyp.leadLep().mc_id() ) leadType=Prompt;
	    else if ( abs(hyp.leadLep().pdgId())==abs(hyp.leadLep().mc_id()) ) leadType=PromptWS;
	    else if ( abs(hyp.leadLep().pdgId())==11 ? abs(hyp.leadLep().mc_id())==13 : abs(hyp.leadLep().mc_id())==11 ) leadType=PromptWF;
	    else if ( hyp.leadLep().mc_id()==22) leadType=PromptFSR;
	    else cout << "UNKNOWN PROMPT LEPTON " << ll << endl;
	    
	  }

	  if (trailType>-1) makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_category","hyp_ss_trail_"+lt+"_category",End,0,End,trailType,weight_);
	  if (leadType>-1 ) makeFillHisto1D<TH1F,int>("hyp_ss_lead_"+ll+"_category","hyp_ss_lead_"+ll+"_category",End,0,End,leadType,weight_);

	  if (abs(hyp.traiLep().pdgId())==13 && hyp.traiLep().pt()>10 && trailType>-1) makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_ptGt10_category","hyp_ss_trail_"+lt+"_ptGt10_category",End,0,End,trailType,weight_);
	  else if (abs(hyp.traiLep().pdgId())==13 && hyp.traiLep().pt()<10 && trailType>-1) makeFillHisto1D<TH1F,int>("hyp_ss_trail_"+lt+"_ptLt10_category","hyp_ss_trail_"+lt+"_ptLt10_category",End,0,End,trailType,weight_);

	  int prompttype = -1;
	  if (isLeadPrompt && isTrailPrompt) prompttype=0;
	  if (isLeadPrompt && !isTrailPrompt) prompttype=1;
	  if (!isLeadPrompt && isTrailPrompt) prompttype=2;
	  if (!isLeadPrompt && !isTrailPrompt) prompttype=3;
	  makeFillHisto1D<TH1F,int>("hyp_prompttype","hyp_prompttype",5,0,5,prompttype,weight_);

	  if (isLeadPrompt && isTrailPrompt) {
	    int leadprompttype = -1;
	    if (leadType==0 && trailType==0) {
	      leadprompttype=0;
	      cout << "UNEXPECTED DOUBLE PROMPT" << endl;
	    } else if (leadType==1 || trailType==1) leadprompttype=1;
	    else if (leadType==2 || trailType==2) leadprompttype=1;
	    else  leadprompttype=3;
	    makeFillHisto1D<TH1F,int>("hyp_leadprompttype","hyp_leadprompttype",5,0,5,leadprompttype,weight_);
	  }

	  /*
	  //dump genps
	  for (unsigned int gp_idx=0;gp_idx<genps_p4().size();++gp_idx) {
	    cout << "gp idx=" << gp_idx << " id=" << genps_id()[gp_idx] 
		 << " status=" << genps_status()[gp_idx]
		 << " id_mother=" << genps_id_mother()[gp_idx]
		 << " p4=" << genps_p4()[gp_idx]
		 << endl;    
	  }
	  */

	} else {
	  //opposite sign
	  makeFillHisto1D<TH1F,float>("hyp_os_mll","hyp_os_mll",100,0,1000,hyp.p4().mass(),weight_);
	  makeFillHisto1D<TH1F,float>("hyp_os_ptll","hyp_os_ptll",100,0,1000,hyp.p4().pt(),weight_);
	  makeFillHisto1D<TH1F,int>("hyp_os_type","hyp_os_type",5,0,5,type,weight_);
	}
      }//if (makehist)

    }

    if (makeSSskim || makeQCDskim) {
      skim_file->cd(); 
      skim_tree->Write(); 
      skim_file->Close();
      delete skim_file;
    }
    delete tree;
    f.Close();
  }

  fr_file->Close();

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
