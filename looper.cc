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

bool ptsort (int i,int j) { return (genps_p4()[i].pt()>genps_p4()[j].pt()); }

bool lepsort (Lep i,Lep j) { 
  if ( abs(i.pdgId())==abs(j.pdgId()) ) return ( i.pt()>j.pt() ); //sort by pt if same flavor
  else return ( abs(i.pdgId())>abs(j.pdgId()) ); //prefer muons over electrons, but check that mu have pt>25//fixme, need to sync // && i.pt()>ptCutHigh
}

bool jetptsort (Jet i,Jet j) { return (i.pt()>j.pt()); }

//fixme: put WF and FSR in different categories
enum LeptonCategories { Prompt = 0, PromptWS = 1, PromptWF = 2, PromptFSR = 2, 
			FakeLightTrue = 3, FakeC = 4, FakeB = 5, FakeLightFake = 6, FakeHiPtGamma = 7, 
			FakeUnknown = 8, FakeLowPtGamma = 9, All9999 = 10,
			Other = 11, End = 12};

int looper::ScanChain( TChain* chain, TString prefix, TString postfix, bool isData, TString whatTest, int nEvents, vector<int> evtToDebug) {

  makebaby       = 0;
  makehist       = 1;
  maketext       = 0;

  bool makeQCDtest    = 0;
  bool makeDYtest     = 0;
  bool makeWZtest     = 0;
  bool makeSSskim     = 0;
  bool makeQCDskim    = 0;
  bool makeSyncTest   = 0;
  if (whatTest=="QCDtest") makeQCDtest    = 1;
  if (whatTest=="DYtest" ) makeDYtest     = 1;
  if (whatTest=="WZtest" ) makeWZtest     = 1;
  if (whatTest=="SSskim" ){makeSSskim     = 1; makehist = 0;}
  if (whatTest=="QCDskim") makeQCDskim    = 1;
  if (whatTest=="SyncTest")makeSyncTest   = 1;

  cout << "Processing " << prefix << " " << postfix << " " << whatTest << endl;

  bool debug = 0;  
  if (debug) cout << "running with flags: makeQCDtest=" << makeQCDtest << " makeDYtest=" << makeDYtest 
		  << " makeWZtest=" << makeWZtest << " makeSSskim=" << makeSSskim << " makeQCDskim=" << makeQCDskim << " makeSyncTest=" << makeSyncTest
		  << " makebaby=" << makebaby << " makehist=" << makehist << " maketext=" << maketext
		  << endl;

  if (postfix!="") postfix = "_"+postfix;
  if (makebaby) MakeBabyNtuple( Form( "%s_baby%s.root", prefix.Data(), postfix.Data() ) );
  if (makehist) CreateOutputFile( Form( "%s_histos%s.root", prefix.Data(), postfix.Data() ) );

  TFile* fr_file=0;
  if (!makeDYtest&&!makeQCDtest&&!makeSSskim&&!makeQCDskim) 
    fr_file=TFile::Open("fakeRates_qcd_pt-50to170.root");

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

    if (debug) cout << "processing file: " << currentFile->GetTitle() << endl;

    //Skimmed output file - needs to be before cms2.Init(tree)
    TFile *skim_file = 0;
    TTree* skim_tree = 0;
    if (makeSSskim || makeQCDskim) {
      TString skim_file_name = TString(currentFile->GetTitle());
      if (makeSSskim) {
	skim_file_name.ReplaceAll(".root","_skimSS_"+prefix+".root");
	skim_file_name.Remove(0,skim_file_name.Last('/'));
	skim_file_name.Prepend("./"+prefix+"_skimSS/");
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

      float lumi = 10.;

      if (evtToDebug.size()>0) {
	bool pass = false;
	for (unsigned int evt=0;evt<evtToDebug.size();evt++) {
	  if (evt_event()==abs(evtToDebug[evt])) {
	    pass = true;
	    break;
	  }
	}
	if (!pass) continue;
	std::cout << "event =" << evt_event() << std::endl;
	cout << "file=" << currentFile->GetTitle() << " run=" << run_ << " evt=" << evt_ << endl;      
	debug=1;
      }

      //fill baby
      run_   = evt_run();
      ls_    = evt_lumiBlock();
      evt_   = evt_event();
      weight_ = isData ? 1. : lumi*evt_scale1fb();
      //if (prefix=="wj") weight_*=30;//fixme using only a subset of events
      //if (prefix=="dy") weight_*=25;//fixme using only a subset of events
      if (makeSyncTest) weight_ = 1.;

      if (newfile) {
	if (debug) cout << "weight=" << weight_ << endl;
	newfile = false;
      }

      if (debug) cout << "file=" << currentFile->GetTitle() << " run=" << run_ << " evt=" << evt_ << endl;

      if (mus_dxyPV().size()!=mus_dzPV().size()) {
	cout << "run=" << run_ << " evt=" << evt_ << " mus_dxyPV().size()=" << mus_dxyPV().size() << " mus_dzPV().size()=" << mus_dzPV().size() << endl;
	continue;
      }

      //make sure there are genuine SS leptons for cut_flow_ss
      int qp=0,qn=0;
      int qpe=0,qne=0;
      int qpm=0,qnm=0;
      for (unsigned int gp=0;gp<genps_id().size();++gp) {
	int pdgid = abs(genps_id()[gp]);
	if (pdgid!=13 && pdgid!=11) continue;
	if (genps_id_mother()[gp]!=23 && abs(genps_id_mother()[gp])!=24) continue;
	if (genps_status()[gp]!=1) continue;//is this needed?
	if (fabs(genps_p4()[gp].eta())>2.4) continue;
	if (genps_p4()[gp].pt()<5) continue;
	if (pdgid==11 && genps_p4()[gp].pt()<10) continue;
	if (genps_charge()[gp]==1 ) qp++;
	if (genps_charge()[gp]==-1) qn++;
	if (pdgid==11 && genps_charge()[gp]==1 ) qpe++;
	if (pdgid==11 && genps_charge()[gp]==-1) qne++;
	if (pdgid==13 && genps_charge()[gp]==1 ) qpm++;
	if (pdgid==13 && genps_charge()[gp]==-1) qnm++;
      }
      bool isGenSS = false;
      bool isGenSSee = false;
      bool isGenSSmm = false;
      if (qp>1 || qn>1) isGenSS = true;
      if (qpe>1 || qne>1) isGenSSee = true;
      if (qpm>1 || qnm>1) isGenSSmm = true;

      if (makeSSskim&&isGenSS) {
	  cms2.LoadAllBranches();
	  skim_file->cd(); 
	  skim_tree->Fill();
	  continue;	
      }

      makeFillHisto1D<TH1F,int>("cut_flow","cut_flow",50,0,50,0,weight_);
      if (isGenSS) makeFillHisto1D<TH1F,int>("cut_flow_ss","cut_flow_ss",50,0,50,0,weight_);
      if (isGenSSee) makeFillHisto1D<TH1F,int>("cut_flow_ssee","cut_flow_ssee",50,0,50,0,weight_);
      if (isGenSSmm) makeFillHisto1D<TH1F,int>("cut_flow_ssmm","cut_flow_ssmm",50,0,50,0,weight_);

      //start selection
      //vertex
      //fixme: check leptons are from this vertex!
      if (firstGoodVertex()<0) {
        if (debug) cout << "does not pass firstGoodVertex" << endl;
        continue;
      }
      makeFillHisto1D<TH1F,int>("cut_flow","cut_flow",50,0,50,1,weight_);
      if (isGenSS) makeFillHisto1D<TH1F,int>("cut_flow_ss","cut_flow_ss",50,0,50,1,weight_);
      if (isGenSSee) makeFillHisto1D<TH1F,int>("cut_flow_ssee","cut_flow_ssee",50,0,50,1,weight_);
      if (isGenSSmm) makeFillHisto1D<TH1F,int>("cut_flow_ssmm","cut_flow_ssmm",50,0,50,1,weight_);

      //fakable objects
      if (debug) cout << "fobs" << endl;
      vector<Lep> fobs;
      for (unsigned int elidx=0;elidx<els_p4().size();++elidx) {
	//electron fo selection
	if (debug) cout << "el pt=" << els_p4()[elidx].pt() << " eta=" << els_p4()[elidx].eta() << " phi=" << els_p4()[elidx].phi() << " q=" << els_charge()[elidx] << endl;
	if (isFakableElectron(elidx)==0) continue;
	if (debug) cout << "pass FO selection" << endl;
	Lep foel(-1*els_charge().at(elidx)*11,elidx);
	fobs.push_back(foel);
      }
      for (unsigned int muidx=0;muidx<mus_p4().size();++muidx) {
	//muon fo selection
	if (debug) cout << "mu pt=" << mus_p4()[muidx].pt() << " eta=" << mus_p4()[muidx].eta() << " phi=" << mus_p4()[muidx].phi() << " q=" << mus_charge()[muidx]<< endl;
	if (isFakableMuon(muidx)==0) continue;
	if (debug) cout << "pass FO selection" << endl;
	Lep fomu(-1*mus_charge().at(muidx)*13,muidx);
	fobs.push_back(fomu);
      }
      if (fobs.size()==0 && !makeDYtest) continue;
      makeFillHisto1D<TH1F,int>("cut_flow","cut_flow",50,0,50,2,weight_);
      if (isGenSS) makeFillHisto1D<TH1F,int>("cut_flow_ss","cut_flow_ss",50,0,50,2,weight_);
      if (isGenSSee) makeFillHisto1D<TH1F,int>("cut_flow_ssee","cut_flow_ssee",50,0,50,2,weight_);
      if (isGenSSmm) makeFillHisto1D<TH1F,int>("cut_flow_ssmm","cut_flow_ssmm",50,0,50,2,weight_);

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
      for (unsigned int fo=0;fo<fobs.size();++fo) {
	if (abs(fobs[fo].pdgId())==13 && isGoodMuon(fobs[fo].idx())==0) continue;
	if (abs(fobs[fo].pdgId())==11 && isGoodElectron(fobs[fo].idx())==0) continue;
      	if (debug) cout << "good lep id=" << fobs[fo].pdgId() << " pt=" << fobs[fo].pt() 
			<< " eta=" << fobs[fo].eta() << " phi=" << fobs[fo].p4().phi() << " q=" << fobs[fo].charge()<< endl;
      	goodleps.push_back(fobs[fo]);
      }

      //veto leptons
      if (debug) cout << "vetoleps" << endl;
      vector<Lep> vetoleps;
      for (unsigned int elidx=0;elidx<els_p4().size();++elidx) {
	//medium electron selection
	if (isGoodVetoElectron(elidx)==0) continue;
	Lep vetoel(-1*els_charge().at(elidx)*11,elidx);
	vetoleps.push_back(vetoel);
      }
      for (unsigned int muidx=0;muidx<mus_p4().size();++muidx) {
	//veto muon selection
	if (isGoodVetoMuon(muidx)==0) continue;
	Lep vetomu(-1*mus_charge().at(muidx)*13,muidx);
	vetoleps.push_back(vetomu);
      }

      //jets, ht, btags
      if (debug) cout << "jets" << endl;
      int njets = 0;
      vector<Jet> jets;
      vector<Jet> alljets;
      vector<Jet> lepjets;
      int nbtag = 0;
      int nbtag_pt35 = 0;
      int nbtag_pt30 = 0;
      int nbtag_pt25 = 0;
      int nbtag_pt20 = 0;
      vector<Jet> btags;
      float ht=0;
      //fixme: should add corrections
      float drcut = 0.4;
      if (makeQCDtest) drcut = 1.0;
      for (unsigned int pfjidx=0;pfjidx<pfjets_p4().size();++pfjidx) {
        Jet jet(pfjidx);
        if (debug) cout << "jet pt=" << jet.pt() << " eta=" << jet.eta() << " phi=" << jet.phi() << endl;
      }
      for (unsigned int pfjidx=0;pfjidx<pfjets_p4().size();++pfjidx) {
        Jet jet(pfjidx);
        if (debug) cout << "jet pt=" << jet.pt() << " eta=" << jet.eta() << endl;
        if (fabs(jet.eta())>2.4) continue;
        if (debug) cout << "jet pass eta" << endl;
        //add pu jet id pfjets_pileupJetId()>????
        bool isLep = false;
        //fixme: should be checked agains fo or good leptons?
        if (makeQCDtest) {
          for (unsigned int fo=0;fo<fobs.size();++fo) {
            if (deltaR(jet.p4(),fobs[fo].p4())<drcut) {
              isLep =true;
              break;
            }
          }
        } else {
          //for (unsigned int gl=0;gl<goodleps.size();++gl) {
          for (unsigned int gl=0;gl<vetoleps.size();++gl) {//fixme not sure this is the best choice
	    Lep lep = vetoleps[gl];
	    if (lep.pt()<10) continue;
            if (deltaR(jet.p4(),lep.p4())<drcut) {
              isLep =true;
	      if (debug) cout << "jet overlapping with lepton" << endl;
              break;
            }
          }
          for (unsigned int fo=0;fo<fobs.size();++fo) {
            if (deltaR(jet.p4(),fobs[fo].p4())<0.7) {
	      lepjets.push_back(jet);
	      if (debug) cout << "found lepjet inerfering with fo" << endl;
              //break;
            }
          }
        }
        if (isLep) continue;
        if (debug) cout << "jet pass lepton clean" << endl;
        if (isLoosePFJet(pfjidx)==false) continue;
        if (debug) cout << "jet pass loose pf id" << endl;
        float jetpt = jet.pt();
        if (jetpt>40.) {
          if (debug) cout << "jet pass pT cut" << endl;
          njets++;
          jets.push_back(jet);
          ht+=jetpt;
          if (jet.isBtag()) {
            nbtag++;
            btags.push_back(jet);
          }
        }
	alljets.push_back(jet);
	if (jetpt>20. && jet.csv()>0.679) nbtag_pt20++;
	if (jetpt>25. && jet.csv()>0.679) nbtag_pt25++;
	if (jetpt>30. && jet.csv()>0.679) nbtag_pt30++;
	if (jetpt>35. && jet.csv()>0.679) nbtag_pt35++;
      }
      std::sort(jets.begin(),jets.end(),jetptsort);
      std::sort(btags.begin(),btags.end(),jetptsort);

      /*
      //add back to goodleps those fobs not passing iso but passing ptrel wrt lepjet
      for (unsigned int fo=0;fo<fobs.size();++fo) {
	if (fabs(fobs[fo].dxyPV())>0.01 || (abs(fobs[fo].pdgId())==13&&fabs(fobs[fo].dxyPV())>0.005)) continue;
	if (fabs(fobs[fo].dzPV())>0.1) continue;
	if (fobs[fo].relIso03()<0.1) continue;//ok, this is inverted here
	int lepjetidx = -1;
	float mindr = 0.7;
	for (unsigned int j=0;j<lepjets.size();++j) {
	  float dr = deltaR(lepjets[j].p4(),fobs[fo].p4());
	  if (dr<mindr) {
	    mindr = dr;
	    lepjetidx = j;
	  }
	} 
	if (lepjetidx>=0) {
	  float sinA = fabs(fobs[fo].p4().x()*lepjets[lepjetidx].p4().y()-fobs[fo].p4().y()*lepjets[lepjetidx].p4().x())/(sqrt(fobs[fo].p4().x()*fobs[fo].p4().x()+fobs[fo].p4().y()*fobs[fo].p4().y())*sqrt(lepjets[lepjetidx].p4().x()*lepjets[lepjetidx].p4().x()+lepjets[lepjetidx].p4().y()*lepjets[lepjetidx].p4().y()));//fixme fabs? 
	  float ptrel = fobs[fo].pt()*sinA;
	  if (ptrel>8.) goodleps.push_back(fobs[fo]);
	}
      }
      */

      //met
      if (debug) cout << "met" << endl;
      float met = evt_pfmet();
      if (met<30. && ht<500. && !makeQCDskim && !makeQCDtest && !makeDYtest) continue;
      makeFillHisto1D<TH1F,int>("cut_flow","cut_flow",50,0,50,3,weight_);
      if (isGenSS) makeFillHisto1D<TH1F,int>("cut_flow_ss","cut_flow_ss",50,0,50,3,weight_);
      if (isGenSSee) makeFillHisto1D<TH1F,int>("cut_flow_ssee","cut_flow_ssee",50,0,50,3,weight_);
      if (isGenSSmm) makeFillHisto1D<TH1F,int>("cut_flow_ssmm","cut_flow_ssmm",50,0,50,3,weight_);

      if (makeQCDtest) {
	if (debug) cout << "qcdtest" << endl;
	runQCDtest(fobs, goodleps, njets, met);
	continue;
      }

      if (makeDYtest) {
	if (debug) cout << "dytest" << endl;
	runDYtest(fobs, goodleps, njets, met);
	continue;
      }

      if (makeWZtest) {
	if (debug) cout << "WZtest" << endl;
	runWZtest(vetoleps, fobs, goodleps, njets, met);
	continue;
      }

      vector<Lep> hypleps;
      //select best hyp in case >3 good leps
      if (goodleps.size()>2) {
	vector<Lep> goodlepsp, goodlepsn;
	for (unsigned int gl=0;gl<goodleps.size();++gl) {
	  if (goodleps[gl].pdgId()>0) goodlepsn.push_back(goodleps[gl]);
	  else goodlepsp.push_back(goodleps[gl]);
	}
	//sort goodleps by muon and then by pt
	std::sort(goodlepsp.begin(),goodlepsp.end(),lepsort);
	std::sort(goodlepsn.begin(),goodlepsn.end(),lepsort);
	//take first SS hyp
	if (goodlepsp.size()<2) {
	  hypleps.push_back(goodlepsn[0]);
	  hypleps.push_back(goodlepsn[1]);
	} else if (goodlepsn.size()<2) {
	  hypleps.push_back(goodlepsp[0]);
	  hypleps.push_back(goodlepsp[1]);
	} else {
	  if ( abs(goodlepsn[0].pdgId()+goodlepsn[1].pdgId())>abs(goodlepsp[0].pdgId()+goodlepsp[1].pdgId()) ) {
	    hypleps.push_back(goodlepsn[0]);
	    hypleps.push_back(goodlepsn[1]);
	  } else if ( abs(goodlepsn[0].pdgId()+goodlepsn[1].pdgId())<abs(goodlepsp[0].pdgId()+goodlepsp[1].pdgId()) ) {
	    hypleps.push_back(goodlepsp[0]);
	    hypleps.push_back(goodlepsp[1]);
	  } else if ( (goodlepsn[0].pt()+goodlepsn[1].pt())>(goodlepsp[0].pt()+goodlepsp[1].pt()) ) {
	    hypleps.push_back(goodlepsn[0]);
	    hypleps.push_back(goodlepsn[1]);
	  } else {
	    hypleps.push_back(goodlepsp[0]);
	    hypleps.push_back(goodlepsp[1]);
	  }
	}
      } else if (goodleps.size()==2) {
	hypleps.push_back(goodleps[0]);
	hypleps.push_back(goodleps[1]);	
      }      
      if (debug) cout << "hypleps size=" << hypleps.size() << endl;

      if (fobs.size()<2) {
	if (debug) cout << "skip, fobs size=" << fobs.size() << " hypleps size=" << hypleps.size() << endl;
	if (isGenSS) {
	  for (unsigned int gp=0;gp<genps_id().size();++gp) {
	    int pdgid = abs(genps_id()[gp]);
	    if (pdgid!=13 && pdgid!=11) continue;
	    if (genps_id_mother()[gp]!=23 && abs(genps_id_mother()[gp])!=24) continue;
	    if (genps_status()[gp]!=1) continue;//is this needed?
	    if (fabs(genps_p4()[gp].eta())>2.4) continue;
	    if (genps_p4()[gp].pt()<5) continue;
	    if (pdgid==11 && genps_p4()[gp].pt()<10) continue;
	    if (debug) cout << "gen lep id=" << genps_id()[gp] << " eta=" << genps_p4()[gp].eta() << " pt=" << genps_p4()[gp].pt() << endl;	  
	    if (pdgid==11) {
	      for (unsigned int elidx=0;elidx<els_p4().size();++elidx) {
		if (fabs(ROOT::Math::VectorUtil::DeltaR(els_p4()[elidx],genps_p4()[gp]))>0.1) continue;
		if (debug) cout << "el pt=" << els_p4()[elidx].pt() << " eta=" << els_p4()[elidx].eta() << " q=" << els_charge()[elidx] << " iso=" << eleRelIso03(elidx)<< endl;
		int failmode = 0;
		if (isElectronFO(elidx)==0) failmode = isElectronFO_debug(elidx);
		if (threeChargeAgree(elidx)==0) {
		  failmode = 10;
		  bool wrongGsf = (genps_charge()[gp]!=els_charge()[elidx]);
		  bool wrongTrk = (genps_charge()[gp]!=els_trk_charge()[elidx]);
		  bool wrongScl = (genps_charge()[gp]!=els_sccharge()[elidx]);
		  makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode","el_fo_q3fail_mode",7,0,7,0,weight_);
		  if (wrongGsf) makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode","el_fo_q3fail_mode",7,0,7,1,weight_);
		  else makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode","el_fo_q3fail_mode",7,0,7,2,weight_);
		  if (wrongTrk) makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode","el_fo_q3fail_mode",7,0,7,3,weight_);
		  else makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode","el_fo_q3fail_mode",7,0,7,4,weight_);
		  if (wrongScl) makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode","el_fo_q3fail_mode",7,0,7,5,weight_);
		  else makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode","el_fo_q3fail_mode",7,0,7,6,weight_);
		  if (fabs(els_etaSC().at(elidx)) <= 1.479){
		    makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_barrel","el_fo_q3fail_mode_barrel",7,0,7,0,weight_);
		    if (wrongGsf) makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_barrel","el_fo_q3fail_mode_barrel",7,0,7,1,weight_);
		    else makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_barrel","el_fo_q3fail_mode_barrel",7,0,7,2,weight_);
		    if (wrongTrk) makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_barrel","el_fo_q3fail_mode_barrel",7,0,7,3,weight_);
		    else makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_barrel","el_fo_q3fail_mode_barrel",7,0,7,4,weight_);
		    if (wrongScl) makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_barrel","el_fo_q3fail_mode_barrel",7,0,7,5,weight_);
		    else makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_barrel","el_fo_q3fail_mode_barrel",7,0,7,6,weight_);
		  } else {
		    makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_endcap","el_fo_q3fail_mode_endcap",7,0,7,0,weight_);
		    if (wrongGsf) makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_endcap","el_fo_q3fail_mode_endcap",7,0,7,1,weight_);
		    else makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_endcap","el_fo_q3fail_mode_endcap",7,0,7,2,weight_);
		    if (wrongTrk) makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_endcap","el_fo_q3fail_mode_endcap",7,0,7,3,weight_);
		    else makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_endcap","el_fo_q3fail_mode_endcap",7,0,7,4,weight_);
		    if (wrongScl) makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_endcap","el_fo_q3fail_mode_endcap",7,0,7,5,weight_);
		    else makeFillHisto1D<TH1F,int>("el_fo_q3fail_mode_endcap","el_fo_q3fail_mode_endcap",7,0,7,6,weight_);
		  }
		}
		if (fabs(els_etaSC().at(elidx)) <= 1.479)
		  makeFillHisto1D<TH1F,float>("el_fo_fail_sietaieta_barrel","el_fo_fail_sietaieta_barrel",25,0.,0.025,els_sigmaIEtaIEta_full5x5().at(elidx),weight_);
		else
		  makeFillHisto1D<TH1F,float>("el_fo_fail_sietaieta_endcap","el_fo_fail_sietaieta_endcap",25,0.,0.050,els_sigmaIEtaIEta_full5x5().at(elidx),weight_);
		if (fabs(els_p4().at(elidx).eta())>2.4) failmode = 11;
		if (els_p4().at(elidx).pt()<10.) failmode = 12;
		if (isFakableElectron(elidx)==0) {
		  makeFillHisto1D<TH1F,int>("el_fo_fail_mode","el_fo_fail_mode",15,0,15,failmode,weight_);
		  if (fabs(els_etaSC().at(elidx)) <= 1.479) makeFillHisto1D<TH1F,int>("el_fo_fail_mode_barrel","el_fo_fail_mode_barrel",15,0,15,failmode,weight_);
		  else makeFillHisto1D<TH1F,int>("el_fo_fail_mode_endcap","el_fo_fail_mode_endcap",15,0,15,failmode,weight_);
		  continue;
		}
		if (debug) cout << "pass FO selection" << endl;
	      }
	    }
	  }
	  if (debug) {
	    for (unsigned int fo=0;fo<fobs.size();++fo) {
	      cout<< "fo lep id=" << fobs[fo].pdgId() << " eta=" << fobs[fo].eta() << " pt=" << fobs[fo].pt() << " relIso=" << fobs[fo].relIso03()<< endl;
	    }
	  }
	}
	continue;
      }
      makeFillHisto1D<TH1F,int>("cut_flow","cut_flow",50,0,50,4,weight_);
      if (isGenSS) makeFillHisto1D<TH1F,int>("cut_flow_ss","cut_flow_ss",50,0,50,4,weight_);
      if (isGenSSee) makeFillHisto1D<TH1F,int>("cut_flow_ssee","cut_flow_ssee",50,0,50,4,weight_);
      if (isGenSSmm) makeFillHisto1D<TH1F,int>("cut_flow_ssmm","cut_flow_ssmm",50,0,50,4,weight_);

      DilepHyp hyp = (hypleps.size()==2 ? DilepHyp(hypleps[0],hypleps[1]) : DilepHyp(fobs[0],fobs[1]) );
      if (hyp.p4().mass()<8) {
	if (debug) cout<< "skip, hyp mass=" << hyp.p4().mass() <<endl;
	continue;
      }
      makeFillHisto1D<TH1F,int>("cut_flow","cut_flow",50,0,50,5,weight_);
      if (isGenSS) makeFillHisto1D<TH1F,int>("cut_flow_ss","cut_flow_ss",50,0,50,5,weight_);
      if (isGenSSee) makeFillHisto1D<TH1F,int>("cut_flow_ssee","cut_flow_ssee",50,0,50,5,weight_);
      if (isGenSSmm) makeFillHisto1D<TH1F,int>("cut_flow_ssmm","cut_flow_ssmm",50,0,50,5,weight_);

      unsigned int ac_base = analysisCategory(hyp.leadLep(),hyp.traiLep());
      passesBaselineCuts(njets, nbtag, met, ht, ac_base);
      if (ac_base==0) {
	if (debug) {
	  cout << "skip, not passing baseline cuts" << endl;
	  cout << "njets=" << njets << " nbtag=" << nbtag << " ht=" << ht << " met=" << met << endl;
	}
	continue;
      }
      makeFillHisto1D<TH1F,int>("cut_flow","cut_flow",50,0,50,6,weight_);
      if (isGenSS) makeFillHisto1D<TH1F,int>("cut_flow_ss","cut_flow_ss",50,0,50,6,weight_);
      if (isGenSSee) makeFillHisto1D<TH1F,int>("cut_flow_ssee","cut_flow_ssee",50,0,50,6,weight_);
      if (isGenSSmm) makeFillHisto1D<TH1F,int>("cut_flow_ssmm","cut_flow_ssmm",50,0,50,6,weight_);

      int br = baselineRegion(nbtag);
      unsigned int ac_sig = ac_base;
      passesSignalRegionCuts(ht, met, ac_sig);
      if (debug) cout << "ac_base=" << ac_base << " ac_sig=" << ac_sig << endl;
      int sr = ac_sig!=0 ? signalRegion(njets, nbtag, met, ht) : -1;

      //write skim here (only ss)
      if (makeSSskim) {
	if (debug) cout << "ss skim" << endl;
	if (isGenSS || hyp.charge()!=0) {
	  cms2.LoadAllBranches();
	  skim_file->cd(); 
	  skim_tree->Fill();
	}
	continue;
      }

      //3rd lepton veto
      bool leptonVeto = false;
      for (unsigned int gl=0;gl<hypleps.size();++gl) {
	for (unsigned int vl=0;vl<vetoleps.size();++vl) {
	  if ( fabs(ROOT::Math::VectorUtil::DeltaR(hypleps[gl].p4(),vetoleps[vl].p4()))<0.001 ) continue;
	  if ( hypleps[gl].pdgId()!=-vetoleps[vl].pdgId() ) continue; 
	  float mll = (hypleps[gl].p4()+vetoleps[vl].p4()).mass();
	  if ( (fabs(mll-91.2)<15&&vetoleps[vl].pt()>10. ) || (mll<12&&vetoleps[vl].pt()>5. ) ) {
	    if (debug) cout << "mll=" << mll << " l1id=" << hypleps[gl].pdgId() << " pt=" << hypleps[gl].pt() << " vlid=" << vetoleps[vl].pdgId() << " pt=" << vetoleps[vl].pt() << endl;
	    leptonVeto = true;
	    break;
	  }
	}
      }
      if (leptonVeto) {
	if (debug) cout<< "skip, 3rd lepton veto" << endl;	
	if (debug) cout<< "hyp leps id=" << hypleps[0].pdgId() << " pt=" << hypleps[0].pt() << " id=" << hypleps[1].pdgId() << " pt=" << hypleps[1].pt() << endl;
	continue;
      }
      makeFillHisto1D<TH1F,int>("cut_flow","cut_flow",50,0,50,7,weight_);
      if (isGenSS) makeFillHisto1D<TH1F,int>("cut_flow_ss","cut_flow_ss",50,0,50,7,weight_);
      if (isGenSSee) makeFillHisto1D<TH1F,int>("cut_flow_ssee","cut_flow_ssee",50,0,50,7,weight_);
      if (isGenSSmm) makeFillHisto1D<TH1F,int>("cut_flow_ssmm","cut_flow_ssmm",50,0,50,7,weight_);

      //compute fake rate in ss ttbar
      if (hyp.charge()!=0 && prefix=="ttbar") {
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
	  if (!isNumerator && fr_file) {
	    TH2F* fr_h = (TH2F*) fr_file->Get((pdgid==13?"fr_mu_gen":"fr_el_gen"));
	    float maxPt=fr_h->GetXaxis()->GetBinUpEdge(fr_h->GetXaxis()->GetNbins())-0.01; 
	    float maxEta=fr_h->GetYaxis()->GetBinUpEdge(fr_h->GetYaxis()->GetNbins())-0.01;
	    float pt = fobs[fo].pt();
	    float eta = fabs(fobs[fo].eta());
	    if (pt>maxPt) pt=maxPt;
	    if (eta>maxEta) eta=maxEta;
	    float fr = fr_h->GetBinContent(fr_h->FindBin(pt,eta));
	    //float fre = fr_h->GetBinError(fr_h->FindBin(pt,eta));
	    float frW = fr/(1.-fr);
	    //std::cout << "fake rate=" << fr << " +/- " << fre << std::endl;
	    makeFillHisto2D<TH2F,float>((pdgid==13?"fr_mu_close":"fr_el_close"),(pdgid==13?"fr_mu_close":"fr_el_close"),10,0.,50.,fobs[fo].pt(),5,0.,2.5,fabs(fobs[fo].eta()),weight_*frW);
	  }
	}
      }

      bool ht400met120 = false;
      if (ht>400. || met>120.) ht400met120 = true;//continue;//fixme

      TString ll = abs(hyp.leadLep().pdgId())==13 ? "mu" : "el";
      TString lt = abs(hyp.traiLep().pdgId())==13 ? "mu" : "el";
      if (hyp.leadLep().pt()>ptCutHigh && hyp.traiLep().pt()>ptCutHigh) {
	//make sure all cuts pass except isolation
	if (fabs(hyp.traiLep().dzPV())<0.1 && fabs(hyp.traiLep().dxyPV())<0.01 && (abs(hyp.traiLep().pdgId())==11 || fabs(hyp.traiLep().dxyPV())<0.005) ) {
	  if (isFromW(hyp.traiLep())) {
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFromWtrail_"+lt+"_relIso03","hyp_ss_foFromWtrail_"+lt+"_relIso03",10,0.,1.,hyp.traiLep().relIso03(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFromWtrail_"+lt+"_iso03","hyp_ss_foFromWtrail_"+lt+"_iso03",10,0.,20.,hyp.traiLep().relIso03()*hyp.traiLep().pt(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFromWtrail_"+lt+"_pt","hyp_ss_foFromWtrail_"+lt+"_pt",10,0.,200.,hyp.traiLep().pt(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFromWtrail_"+lt+"_dxyPV","hyp_ss_foFromWtrail_"+lt+"_dxyPV",100,0,0.1,fabs(hyp.traiLep().dxyPV()),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFromWtrail_"+lt+"_dzPV","hyp_ss_foFromWtrail_"+lt+"_dzPV",100,0,0.5,fabs(hyp.traiLep().dzPV()),weight_);
	  } else {
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFaketrail_"+lt+"_relIso03","hyp_ss_foFaketrail_"+lt+"_relIso03",10,0.,1.,hyp.traiLep().relIso03(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFaketrail_"+lt+"_iso03","hyp_ss_foFaketrail_"+lt+"_iso03",10,0.,20.,hyp.traiLep().relIso03()*hyp.traiLep().pt(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFaketrail_"+lt+"_pt","hyp_ss_foFaketrail_"+lt+"_pt",10,0.,200.,hyp.traiLep().pt(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFaketrail_"+lt+"_dxyPV","hyp_ss_foFaketrail_"+lt+"_dxyPV",100,0,0.1,fabs(hyp.traiLep().dxyPV()),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFaketrail_"+lt+"_dzPV","hyp_ss_foFaketrail_"+lt+"_dzPV",100,0,0.5,fabs(hyp.traiLep().dzPV()),weight_);
	  }
	}

	//make sure all cuts pass except isolation
	if (fabs(hyp.leadLep().dzPV())<0.1 && fabs(hyp.leadLep().dxyPV())<0.01 && (abs(hyp.leadLep().pdgId())==11 || fabs(hyp.leadLep().dxyPV())<0.005) ) {
	  if (isFromW(hyp.leadLep())) {
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFromWlead_"+ll+"_relIso03","hyp_ss_foFromWlead_"+ll+"_relIso03",10,0.,1.,hyp.leadLep().relIso03(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFromWlead_"+ll+"_iso03","hyp_ss_foFromWlead_"+ll+"_iso03",10,0.,20.,hyp.leadLep().relIso03()*hyp.leadLep().pt(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFromWlead_"+ll+"_pt","hyp_ss_foFromWlead_"+ll+"_pt",10,0.,200.,hyp.leadLep().pt(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFromWlead_"+lt+"_dxyPV","hyp_ss_foFromWlead_"+lt+"_dxyPV",100,0,0.1,fabs(hyp.leadLep().dxyPV()),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFromWlead_"+lt+"_dzPV","hyp_ss_foFromWlead_"+lt+"_dzPV",100,0,0.5,fabs(hyp.leadLep().dzPV()),weight_);
	    
	    if (hyp.leadLep().relIso03()>0.1) {	    
	      //test min dR between leptons and jets
	      int lepjetidx = -1;
	      if (debug) cout << "prompt non iso lepton with lepjets.size()=" << lepjets.size() << endl;
	      float mindr = 0.7;
	      int nlepjets_per_lep = 0;
	      for (unsigned int j=0;j<lepjets.size();++j) {
		float dr = deltaR(lepjets[j].p4(),hyp.leadLep().p4());
		if (debug) cout << "dr=" << dr << endl;
		if (dr<0.7) nlepjets_per_lep++;
		if (dr<mindr) {
		  mindr = dr;
		  lepjetidx = j;
		}
	      } 
	      makeFillHisto1D<TH1F,int>("hyp_ss_foFromWlead_"+ll+"_nlepjets","hyp_ss_foFromWlead_"+ll+"_nlepjets",10,0,10,nlepjets_per_lep,weight_);
	      makeFillHisto2D<TH2F,float>("hyp_ss_foFromWlead_"+ll+"_mindr_vs_nlepjets","hyp_ss_foFromWlead_"+ll+"_mindr_vs_nlepjets",10,0,10,nlepjets_per_lep,10,0,1,mindr,weight_);	      
	      if (mindr<0.1) makeFillHisto1D<TH1F,float>("hyp_ss_foFromWlead_"+ll+"_dPtMinDr01","hyp_ss_foFromWlead_"+ll+"_dPtMinDr01",20,-1.,1.,(lepjets[lepjetidx].pt()-hyp.leadLep().pt())/hyp.leadLep().pt(),weight_); 
	      if (lepjetidx>=0) {
		float sinA = fabs(hyp.leadLep().p4().x()*lepjets[lepjetidx].p4().y()-hyp.leadLep().p4().y()*lepjets[lepjetidx].p4().x())/(sqrt(hyp.leadLep().p4().x()*hyp.leadLep().p4().x()+hyp.leadLep().p4().y()*hyp.leadLep().p4().y())*sqrt(lepjets[lepjetidx].p4().x()*lepjets[lepjetidx].p4().x()+lepjets[lepjetidx].p4().y()*lepjets[lepjetidx].p4().y()));//fixme fabs? 
		float ptrel = hyp.leadLep().pt()*sinA;
		if (debug) cout << "sinA=" << sinA << endl;
		makeFillHisto1D<TH1F,float>("hyp_ss_foFromWlead_"+ll+"_ptrel","hyp_ss_foFromWlead_"+ll+"_ptrel",10,0.,20.,ptrel,weight_);
		makeFillHisto2D<TH2F,float>("hyp_ss_foFromWlead_"+ll+"_ptrel_vs_mindr","hyp_ss_foFromWlead_"+ll+"_ptrel_vs_mindr",10,0.,1.,mindr,10,0.,20.,ptrel,weight_);
		if (ll=="mu") makeFillHisto2D<TH2F,float>("hyp_ss_foFromWlead_"+ll+"_ptrel_vs_lepfrac","hyp_ss_foFromWlead_"+ll+"_ptrel_vs_lepfrac",10,0.,1.,pfjets_muonE()[lepjets[lepjetidx].idx()]/(lepjets[lepjetidx].p4().E()*pfjets_undoJEC()[lepjets[lepjetidx].idx()]),10,0.,20.,ptrel,weight_);
		else makeFillHisto2D<TH2F,float>("hyp_ss_foFromWlead_"+ll+"_ptrel_vs_lepfrac","hyp_ss_foFromWlead_"+ll+"_ptrel_vs_lepfrac",10,0.,1.,pfjets_electronE()[lepjets[lepjetidx].idx()]/(lepjets[lepjetidx].p4().E()*pfjets_undoJEC()[lepjets[lepjetidx].idx()]),10,0.,20.,ptrel,weight_);
		makeFillHisto1D<TH1F,float>("hyp_ss_foFromWlead_"+ll+"_sinA","hyp_ss_foFromWlead_"+ll+"_sinA",10,0.,1.,ptrel/hyp.leadLep().pt(),weight_);
	      }
	    }
	    
	  } else {
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFakelead_"+ll+"_relIso03","hyp_ss_foFakelead_"+ll+"_relIso03",10,0.,1.,hyp.leadLep().relIso03(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFakelead_"+ll+"_iso03","hyp_ss_foFakelead_"+ll+"_iso03",10,0.,20.,hyp.leadLep().relIso03()*hyp.leadLep().pt(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFakelead_"+ll+"_pt","hyp_ss_foFakelead_"+ll+"_pt",10,0.,200.,hyp.leadLep().pt(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFakelead_"+lt+"_dxyPV","hyp_ss_foFakelead_"+lt+"_dxyPV",100,0,0.1,fabs(hyp.leadLep().dxyPV()),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_foFakelead_"+lt+"_dzPV","hyp_ss_foFakelead_"+lt+"_dzPV",100,0,0.5,fabs(hyp.leadLep().dzPV()),weight_);
	    
	    if (hyp.leadLep().relIso03()>0.1) {	    
	      //test min dR between leptons and jets
	      int lepjetidx = -1;
	      if (debug) cout << "prompt non iso lepton with lepjets.size()=" << lepjets.size() << endl;
	      float mindr = 0.7;
	      int nlepjets_per_lep = 0;
	      for (unsigned int j=0;j<lepjets.size();++j) {
		float dr = deltaR(lepjets[j].p4(),hyp.leadLep().p4());
		if (debug) cout << "dr=" << dr << endl;
		if (dr<0.7) nlepjets_per_lep++;
		if (dr<mindr) {
		  mindr = dr;
		  lepjetidx = j;
		}
	      } 
	      makeFillHisto1D<TH1F,int>("hyp_ss_foFakelead_"+ll+"_nlepjets","hyp_ss_foFakelead_"+ll+"_nlepjets",10,0,10,nlepjets_per_lep,weight_);
	      makeFillHisto2D<TH2F,float>("hyp_ss_foFakelead_"+ll+"_mindr_vs_nlepjets","hyp_ss_foFakelead_"+ll+"_mindr_vs_nlepjets",10,0,10,nlepjets_per_lep,10,0,1,mindr,weight_);	      
	      if (mindr<0.1) makeFillHisto1D<TH1F,float>("hyp_ss_foFakelead_"+ll+"_dPtMinDr01","hyp_ss_foFakelead_"+ll+"_dPtMinDr01",20,-1.,1.,(lepjets[lepjetidx].pt()-hyp.leadLep().pt())/hyp.leadLep().pt(),weight_); 
	      if (lepjetidx>=0) {
		float sinA = fabs(hyp.leadLep().p4().x()*lepjets[lepjetidx].p4().y()-hyp.leadLep().p4().y()*lepjets[lepjetidx].p4().x())/(sqrt(hyp.leadLep().p4().x()*hyp.leadLep().p4().x()+hyp.leadLep().p4().y()*hyp.leadLep().p4().y())*sqrt(lepjets[lepjetidx].p4().x()*lepjets[lepjetidx].p4().x()+lepjets[lepjetidx].p4().y()*lepjets[lepjetidx].p4().y()));//fixme fabs? 
		float ptrel = hyp.leadLep().pt()*sinA;
		if (debug) cout << "sinA=" << sinA << endl;
		makeFillHisto1D<TH1F,float>("hyp_ss_foFakelead_"+ll+"_ptrel","hyp_ss_foFakelead_"+ll+"_ptrel",10,0.,20.,ptrel,weight_);
		makeFillHisto2D<TH2F,float>("hyp_ss_foFakelead_"+ll+"_ptrel_vs_mindr","hyp_ss_foFakelead_"+ll+"_ptrel_vs_mindr",10,0.,1.,mindr,10,0.,20.,ptrel,weight_);
		if (ll=="mu") makeFillHisto2D<TH2F,float>("hyp_ss_foFakelead_"+ll+"_ptrel_vs_lepfrac","hyp_ss_foFakelead_"+ll+"_ptrel_vs_lepfrac",10,0.,1.,pfjets_muonE()[lepjets[lepjetidx].idx()]/(lepjets[lepjetidx].p4().E()*pfjets_undoJEC()[lepjets[lepjetidx].idx()]),10,0.,20.,ptrel,weight_);
		else makeFillHisto2D<TH2F,float>("hyp_ss_foFakelead_"+ll+"_ptrel_vs_lepfrac","hyp_ss_foFakelead_"+ll+"_ptrel_vs_lepfrac",10,0.,1.,pfjets_electronE()[lepjets[lepjetidx].idx()]/(lepjets[lepjetidx].p4().E()*pfjets_undoJEC()[lepjets[lepjetidx].idx()]),10,0.,20.,ptrel,weight_);
		makeFillHisto1D<TH1F,float>("hyp_ss_foFakelead_"+ll+"_sinA","hyp_ss_foFakelead_"+ll+"_sinA",10,0.,1.,ptrel/hyp.leadLep().pt(),weight_);
	      }
	    }
	  } 
	}
      }

      if (hypleps.size()==2) {
	makeFillHisto1D<TH1F,int>("cut_flow","cut_flow",50,0,50,8,weight_);//check that passes ID
	if (isGenSS) makeFillHisto1D<TH1F,int>("cut_flow_ss","cut_flow_ss",50,0,50,8,weight_);
	if (isGenSSee) makeFillHisto1D<TH1F,int>("cut_flow_ssee","cut_flow_ssee",50,0,50,8,weight_);
	if (isGenSSmm) makeFillHisto1D<TH1F,int>("cut_flow_ssmm","cut_flow_ssmm",50,0,50,8,weight_);
      }
      if (hypleps.size()!=2 || hypleps[0].pt()<ptCutLow || hypleps[1].pt()<ptCutLow) {
	if (debug) {
	  cout << "skip, not passing lepton cuts" << endl; 
	  cout << "fobs size = " << fobs.size() << " pdgids=" << fobs[0].pdgId() << ", " << fobs[1].pdgId() << endl;
	  if (hypleps.size()>0) cout << "Good leptons = " << hypleps.size() << " pdgids=" << hypleps[0].pdgId() << " pt=" << hypleps[0].pt() << " eta=" << hypleps[0].eta() << endl;
	  for (unsigned int fo=0;fo<fobs.size();++fo){
	    if (abs(fobs[fo].pdgId())!=13) continue;
	    cout << "fob pt=" << fobs[fo].pt() << " eta=" << fobs[fo].eta() << endl;
	    if (mus_p4().at(fobs[fo].idx()).pt()<ptCutLow) cout << "fail pt" << endl;
	    if (isMuonFO(fobs[fo].idx())==0) cout << "fail FO" << endl;
	    if (muRelIso03(fobs[fo].idx())>1.0 ) cout << "fail loose iso" << endl;
	    if (isTightMuon(fobs[fo].idx())==0) cout << "fail tight id" << endl;
	    if (muRelIso03(fobs[fo].idx())>0.1 ) cout << "fail tight iso, iso=" << muRelIso03(fobs[fo].idx()) << endl;
	  }
	}
	continue;
      }
      makeFillHisto1D<TH1F,int>("cut_flow","cut_flow",50,0,50,9,weight_);//so this is pT
      if (isGenSS) makeFillHisto1D<TH1F,int>("cut_flow_ss","cut_flow_ss",50,0,50,9,weight_);
      if (isGenSSee) makeFillHisto1D<TH1F,int>("cut_flow_ssee","cut_flow_ssee",50,0,50,9,weight_);
      if (isGenSSmm) makeFillHisto1D<TH1F,int>("cut_flow_ssmm","cut_flow_ssmm",50,0,50,9,weight_);

      if (debug) {
	cout << "Good leptons = " << hypleps.size() << " pdgids=" << hypleps[0].pdgId() << ", " << hypleps[1].pdgId() << endl;
	cout << "Total charge = " << hyp.charge() << " m=" << hyp.p4().mass() << endl;
      }

      //fill histos
      if (makehist) {
	makeFillHisto1D<TH1F,int>("hyp_charge","hyp_charge",7,-3.5,3.5,hyp.charge(),weight_);
	makeFillHisto1D<TH1F,float>("hyp_mll","hyp_mll",100,0,1000,hyp.p4().mass(),weight_);
	makeFillHisto1D<TH1F,float>("hyp_ptll","hyp_ptll",100,0,1000,hyp.p4().pt(),weight_);
	makeFillHisto1D<TH1F,int>("hyp_njets","hyp_njets",20,0,20,njets,weight_);
	makeFillHisto1D<TH1F,int>("hyp_nbtag","hyp_nbtag",20,0,20,nbtag,weight_);
	makeFillHisto1D<TH1F,float>("hyp_ht","hyp_ht",50,0,2000,ht,weight_);
	makeFillHisto1D<TH1F,float>("hyp_met","hyp_met",50,0,500,met,weight_);

	float dPhi1 = deltaPhi(hyp.leadLep().p4().phi(),evt_pfmetPhi());
	float mt1  = mt(hyp.leadLep().pt(),met,dPhi1);
	float dPhi2 = deltaPhi(hyp.traiLep().p4().phi(),evt_pfmetPhi());
	float mt2  = mt(hyp.traiLep().pt(),met,dPhi2);
	float dPhill = deltaPhi(hyp.p4().phi(),evt_pfmetPhi());
	float mtll = mt(hyp.p4().pt(),met,dPhill);

	float mtmin = min(mt1,mt2);
	float mtmin12ll = min(mtmin,mtll);

	makeFillHisto1D<TH1F,float>("hyp_mt1","hyp_mt1",15,0,300,mt1,weight_);
	makeFillHisto1D<TH1F,float>("hyp_mt2","hyp_mt2",15,0,300,mt2,weight_);
	makeFillHisto1D<TH1F,float>("hyp_mtll","hyp_mtll",15,0,300,mtll,weight_);
	makeFillHisto1D<TH1F,float>("hyp_mtmin","hyp_mtmin",15,0,300,mtmin,weight_);
	makeFillHisto1D<TH1F,float>("hyp_mtmin12ll","hyp_mtmin12ll",15,0,300,mtmin12ll,weight_);

	makeFillHisto1D<TH1F,float>("hyp_cosdPhi1","hyp_cosdPhi1",20,-1,1,cos(dPhi1),weight_);
	makeFillHisto1D<TH1F,float>("hyp_cosdPhi2","hyp_cosdPhi2",20,-1,1,cos(dPhi2),weight_);
	makeFillHisto1D<TH1F,float>("hyp_cosdPhill","hyp_cosdPhill",20,-1,1,cos(dPhill),weight_);

	int type = -1;
	if (abs(hyp.traiLep().pdgId())==13 && abs(hyp.leadLep().pdgId())==13) type=0;
	if (abs(hyp.traiLep().pdgId())==13 && abs(hyp.leadLep().pdgId())==11) type=1;
	if (abs(hyp.traiLep().pdgId())==11 && abs(hyp.leadLep().pdgId())==13) type=2;
	if (abs(hyp.traiLep().pdgId())==11 && abs(hyp.leadLep().pdgId())==11) type=3;
	if (debug) cout << "type=" << type << endl;
	makeFillHisto1D<TH1F,int>("hyp_type","hyp_type",5,0,5,type,weight_);

	if (hyp.charge()!=0) {
	  //same sign
	  //if (hyp.leadLep().mc_id()*hyp.traiLep().mc_id()<0) continue;//FIXME
	  makeFillHisto1D<TH1F,int>("cut_flow","cut_flow",50,0,50,10,weight_);
	  if (isGenSS) makeFillHisto1D<TH1F,int>("cut_flow_ss","cut_flow_ss",50,0,50,10,weight_);
	  if (isGenSSee) makeFillHisto1D<TH1F,int>("cut_flow_ssee","cut_flow_ssee",50,0,50,10,weight_);
	  if (isGenSSmm) makeFillHisto1D<TH1F,int>("cut_flow_ssmm","cut_flow_ssmm",50,0,50,10,weight_);

	  if (debug) {
	    cout << endl << "NEW SS EVENT" << endl << endl;	  
	    cout << "lead lep id=" << hyp.leadLep().pdgId() << " pt=" << hyp.leadLep().pt() << " eta=" << hyp.leadLep().eta() << " p4=" << hyp.leadLep().p4() 
		 << " mcid=" << hyp.leadLep().mc_id() << " mcp4=" << hyp.leadLep().mc_p4() << " mother_id=" << hyp.leadLep().mc_motherid()
		 << " mc3idx=" << hyp.leadLep().mc3idx() << " mc3_id=" << hyp.leadLep().mc3_id() 
		 << " mc3_motheridx=" << hyp.leadLep().mc3_motheridx() << " mc3_mother_id=" << hyp.leadLep().mc3_motherid()
		 << " genps_id_mother()[hyp.leadLep().mc3_motheridx()]=" << genps_id_mother()[hyp.leadLep().mc3_motheridx()]
		 << endl;
	    cout << "trai lep id=" << hyp.traiLep().pdgId() << " pt=" << hyp.traiLep().pt() << " eta=" << hyp.traiLep().eta() << " p4=" << hyp.traiLep().p4() 
		 << " mcid=" << hyp.traiLep().mc_id() << " mcp4=" << hyp.traiLep().mc_p4()  << " mother_id=" << hyp.traiLep().mc_motherid()
		 << " mc3idx=" << hyp.traiLep().mc3idx() << " mc3_id=" << hyp.traiLep().mc3_id() 
		 << " mc3_motheridx=" << hyp.traiLep().mc3_motheridx() << " mc3_mother_id=" << hyp.traiLep().mc3_motherid()
		 << " genps_id_mother()[hyp.traiLep().mc3_motheridx()]=" << genps_id_mother()[hyp.traiLep().mc3_motheridx()]
		 << endl;
	    cout << "njets=" << njets << " nbtag=" << nbtag << " ht=" << ht << " met=" << met << endl;
	    cout << "BR" << br << " SR" << sr << endl;
	  }

	  makeFillHisto1D<TH1F,int>("hyp_ss_njets","hyp_ss_njets",20,0,20,njets,weight_);
	  makeFillHisto1D<TH1F,int>("hyp_ss_nbtag","hyp_ss_nbtag",20,0,20,nbtag,weight_);
	  makeFillHisto1D<TH1F,float>("hyp_ss_ht","hyp_ss_ht",50,0,2000,ht,weight_);
	  makeFillHisto1D<TH1F,float>("hyp_ss_met","hyp_ss_met",50,0,500,met,weight_);

	  makeFillHisto2D<TH2F,float>("hyp_ss_met_vs_ht","hyp_ss_met_vs_ht",40,0,2000,ht,20,0,500,met,weight_);

	  //test min dR between leptons and jets
	  float mindr = 99999.;
	  for (unsigned int j=0;j<jets.size();++j) {
	    for (unsigned int gl=0;gl<goodleps.size();++gl) {
	      float dr = deltaR(jets[j].p4(),goodleps[gl].p4());
	      if (dr<mindr) {
		mindr = dr;
	      }
	    }
	  } 
	  //test min dR between leptons and btagged jets
	  makeFillHisto1D<TH1F,float>("hyp_ss_mindRlj","hyp_ss_mindRlj",50,0,5,mindr,weight_);
	  float mindrb = 99999.;
	  for (unsigned int j=0;j<btags.size();++j) {
	    for (unsigned int gl=0;gl<goodleps.size();++gl) {
	      float dr = deltaR(btags[j].p4(),goodleps[gl].p4());
	      if (dr<mindrb) {
		mindrb = dr;
	      }
	    }
	  } 
	  makeFillHisto1D<TH1F,float>("hyp_ss_mindRlb","hyp_ss_mindRlb",50,0,5,mindrb,weight_);
	  //test min dR between jets
	  float mindrj = 99999.;
	  float mass_mindrj = 99999.;
	  for (unsigned int j=0;j<jets.size();++j) {
	    for (unsigned int i=j+1;i<jets.size();++i) {
	      float dr = deltaR(jets[j].p4(),jets[i].p4());
	      if (dr<mindrj) {
		mindrj = dr;
		mass_mindrj = (jets[j].p4()+jets[i].p4()).mass();
	      }
	    }
	  } 
	  makeFillHisto1D<TH1F,float>("hyp_ss_mindRjj","hyp_ss_mindRjj",50,0,5,mindrj,weight_);
	  makeFillHisto1D<TH1F,float>("hyp_ss_mass_mindRjj","hyp_ss_mass_mindRjj",50,0,500,mass_mindrj,weight_);
	  //test dR between leptons
	  makeFillHisto1D<TH1F,float>("hyp_ss_dRll","hyp_ss_dRll",50,0,5,deltaR(goodleps[0].p4(),goodleps[1].p4()),weight_);

	  makeFillHisto1D<TH1F,int>("hyp_ss_leadjetpt","hyp_ss_leadjetpt",50,0,1000,jets[0].pt(),weight_);
	  if (btags.size()>0) makeFillHisto1D<TH1F,int>("hyp_ss_leadbtagpt","hyp_ss_leadbtagpt",50,0,1000,btags[0].pt(),weight_);

	  makeFillHisto1D<TH1F,float>("hyp_ss_mll","hyp_ss_mll",100,0,1000,hyp.p4().mass(),weight_);
	  makeFillHisto1D<TH1F,float>("hyp_ss_ptll","hyp_ss_ptll",100,0,1000,hyp.p4().pt(),weight_);
	  makeFillHisto1D<TH1F,int>("hyp_ss_type","hyp_ss_type",5,0,5,type,weight_);

	  makeFillHisto1D<TH1F,float>("hyp_ss_ngleps","hyp_ss_ngleps",10,0,10,goodleps.size(),weight_);
	  makeFillHisto1D<TH1F,float>("hyp_ss_nfobs","hyp_ss_nfobs",10,0,10,fobs.size(),weight_);
	  makeFillHisto1D<TH1F,float>("hyp_ss_nvetos","hyp_ss_nvetos",10,0,10,vetoleps.size(),weight_);

	  makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_relIso03","hyp_ss_alltrail_"+lt+"_relIso03",100,0.,1.,hyp.traiLep().relIso03(),weight_);
	  makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_relIso03","hyp_ss_alllead_"+ll+"_relIso03",100,0.,1.,hyp.leadLep().relIso03(),weight_);

	  if (abs(hyp.traiLep().pdgId())==11) {
	    unsigned int elIdx = hyp.traiLep().idx();
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_dEtaIn","hyp_ss_alltrail_"+lt+"_dEtaIn",100,0,0.05,fabs(els_dEtaIn().at(elIdx)),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_dPhiIn","hyp_ss_alltrail_"+lt+"_dPhiIn",100,0,0.05,fabs(els_dPhiIn().at(elIdx)),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_sIEtaIEta","hyp_ss_alltrail_"+lt+"_sIEtaIEta",100,0,0.05,els_sigmaIEtaIEta_full5x5().at(elIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_hOverE","hyp_ss_alltrail_"+lt+"_hOverE",100,0,0.5,els_hOverE().at(elIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_dxyPV","hyp_ss_alltrail_"+lt+"_dxyPV",100,0,0.1,fabs(els_dxyPV().at(elIdx)),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_dzPV","hyp_ss_alltrail_"+lt+"_dzPV",100,0,0.5,fabs(els_dzPV().at(elIdx)),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_EoP","hyp_ss_alltrail_"+lt+"_EoP",100,0,0.1,fabs( (1.0/els_ecalEnergy().at(elIdx)) - (els_eOverPIn().at(elIdx)/els_ecalEnergy().at(elIdx)) ),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_convflag","hyp_ss_alltrail_"+lt+"_convflag",2,0,2,els_conv_vtx_flag().at(elIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_expinlay","hyp_ss_alltrail_"+lt+"_expinlay",5,0,5,els_exp_innerlayers().at(elIdx),weight_);
	  } else {
	    unsigned int muIdx = hyp.traiLep().idx();
	    bool isGlobal  = true;
	    bool isTracker = true;
	    if (((mus_type().at(muIdx)) & (1<<1)) == 0) isGlobal  = false;
	    if (((mus_type().at(muIdx)) & (1<<2)) == 0) isTracker = false;
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_isPf","hyp_ss_alltrail_"+lt+"_isPf",2,0,2,mus_pid_PFMuon().at(muIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_isGl","hyp_ss_alltrail_"+lt+"_isGl",2,0,2,isGlobal,weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_isTk","hyp_ss_alltrail_"+lt+"_isTk",2,0,2,isTracker,weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_chi2","hyp_ss_alltrail_"+lt+"_chi2",100,0,20,mus_gfit_chi2().at(muIdx)/mus_gfit_ndof().at(muIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_vStaH","hyp_ss_alltrail_"+lt+"_vStaH",20,0,20,mus_gfit_validSTAHits().at(muIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_nMatchSt","hyp_ss_alltrail_"+lt+"_nMatchSt",20,0,20,mus_numberOfMatchedStations().at(muIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_vPixH","hyp_ss_alltrail_"+lt+"_vPixH",10,0,10,mus_validPixelHits().at(muIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_nLay","hyp_ss_alltrail_"+lt+"_nLay",20,0,20,mus_nlayers().at(muIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_dxyPV","hyp_ss_alltrail_"+lt+"_dxyPV",100,0,0.1,fabs(mus_dxyPV().at(muIdx)),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alltrail_"+lt+"_dzPV","hyp_ss_alltrail_"+lt+"_dzPV",100,0,0.5,fabs(mus_dzPV().at(muIdx)),weight_);
	  }

	  if (abs(hyp.leadLep().pdgId())==11) {
	    unsigned int elIdx = hyp.leadLep().idx();
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_dEtaIn","hyp_ss_alllead_"+ll+"_dEtaIn",100,0,0.05,fabs(els_dEtaIn().at(elIdx)),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_dPhiIn","hyp_ss_alllead_"+ll+"_dPhiIn",100,0,0.05,fabs(els_dPhiIn().at(elIdx)),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_sIEtaIEta","hyp_ss_alllead_"+ll+"_sIEtaIEta",100,0,0.5,els_sigmaIEtaIEta_full5x5().at(elIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_hOverE","hyp_ss_alllead_"+ll+"_hOverE",100,0,0.5,els_hOverE().at(elIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_dxyPV","hyp_ss_alllead_"+ll+"_dxyPV",100,0,0.1,fabs(els_dxyPV().at(elIdx)),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_dzPV","hyp_ss_alllead_"+ll+"_dzPV",100,0,0.5,fabs(els_dzPV().at(elIdx)),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_EoP","hyp_ss_alllead_"+ll+"_EoP",100,0,0.1,fabs( (1.0/els_ecalEnergy().at(elIdx)) - (els_eOverPIn().at(elIdx)/els_ecalEnergy().at(elIdx)) ),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_convflag","hyp_ss_alllead_"+ll+"_convflag",2,0,2,els_conv_vtx_flag().at(elIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_expinlay","hyp_ss_alllead_"+ll+"_expinlay",5,0,5,els_exp_innerlayers().at(elIdx),weight_);
	  } else {
	    unsigned int muIdx = hyp.leadLep().idx();
	    bool isGlobal  = true;
	    bool isTracker = true;
	    if (((mus_type().at(muIdx)) & (1<<1)) == 0) isGlobal  = false;
	    if (((mus_type().at(muIdx)) & (1<<2)) == 0) isTracker = false;
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_isPf","hyp_ss_alllead_"+ll+"_isPf",2,0,2,mus_pid_PFMuon().at(muIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_isGl","hyp_ss_alllead_"+ll+"_isGl",2,0,2,isGlobal,weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_isTk","hyp_ss_alllead_"+ll+"_isTk",2,0,2,isTracker,weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_chi2","hyp_ss_alllead_"+ll+"_chi2",100,0,20,mus_gfit_chi2().at(muIdx)/mus_gfit_ndof().at(muIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_vStaH","hyp_ss_alllead_"+ll+"_vStaH",20,0,20,mus_gfit_validSTAHits().at(muIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_nMatchSt","hyp_ss_alllead_"+ll+"_nMatchSt",20,0,20,mus_numberOfMatchedStations().at(muIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_vPixH","hyp_ss_alllead_"+ll+"_vPixH",10,0,10,mus_validPixelHits().at(muIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_nLay","hyp_ss_alllead_"+ll+"_nLay",20,0,20,mus_nlayers().at(muIdx),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_dxyPV","hyp_ss_alllead_"+ll+"_dxyPV",100,0,0.1,fabs(mus_dxyPV().at(muIdx)),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_ss_alllead_"+ll+"_dzPV","hyp_ss_alllead_"+ll+"_dzPV",100,0,0.5,fabs(mus_dzPV().at(muIdx)),weight_);
	  }

	  if (ac_base & 1<<HighPt) {
	    if (prefix=="TTWJets" && debug) cout << left << setw(10) << ls_  << " "  << setw(10) << evt_ << " "  << setw(10) << njets << " "  << setw(10) << nbtag << " "  << setw(10) << std::setprecision(8) << ht << " "  << setw(10) << std::setprecision(8) << met << " " << sr << endl;
	    makeFillHisto1D<TH1F,int>("hyp_highpt_br","hyp_highpt_br",40,0,40,br,weight_);
	    if ( isFromWZ(hyp.traiLep()) && isFromWZ(hyp.leadLep()) ) {
	      makeFillHisto1D<TH1F,int>("hyp_highpt_br_fromWZ","hyp_highpt_br_fromWZ",40,0,40,br,weight_);
	    }
	    if (makeSyncTest) cout << Form("%1d %9d %12d\t%2d\t%+2d %5.1f\t%+2d %5.1f\t%d\t%2d\t%5.1f\t%6.1f", run_, ls_, evt_, int(vetoleps.size()), hyp.leadLep().pdgId(), hyp.leadLep().pt(), 
					   hyp.traiLep().pdgId(), hyp.traiLep().pt(),  njets, nbtag, met, ht) << endl;
	    if (sr>0) {
	      makeFillHisto1D<TH1F,int>("hyp_highpt_sr","hyp_highpt_sr",40,0,40,br,weight_);
	      makeFillHisto1D<TH1F,int>("hyp_highpt_sr","hyp_highpt_sr",40,0,40,sr,weight_);
	      if ( isFromWZ(hyp.traiLep()) && isFromWZ(hyp.leadLep()) ) {
		makeFillHisto1D<TH1F,int>("hyp_highpt_sr_fromWZ","hyp_highpt_sr_fromWZ",40,0,40,br,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highpt_sr_fromWZ","hyp_highpt_sr_fromWZ",40,0,40,sr,weight_);
	      }
	      makeFillHisto1D<TH1F,int>("hyp_highpt_njets","hyp_highpt_njets",8,0,8,njets,weight_);
	      makeFillHisto1D<TH1F,int>("hyp_highpt_nbtag","hyp_highpt_nbtag",8,0,8,nbtag,weight_);
	      makeFillHisto1D<TH1F,int>("hyp_highpt_nbtag_pt20","hyp_highpt_nbtag_pt20",8,0,8,nbtag_pt20,weight_);
	      makeFillHisto1D<TH1F,int>("hyp_highpt_nbtag_pt25","hyp_highpt_nbtag_pt25",8,0,8,nbtag_pt25,weight_);
	      makeFillHisto1D<TH1F,int>("hyp_highpt_nbtag_pt30","hyp_highpt_nbtag_pt30",8,0,8,nbtag_pt30,weight_);
	      makeFillHisto1D<TH1F,int>("hyp_highpt_nbtag_pt35","hyp_highpt_nbtag_pt35",8,0,8,nbtag_pt35,weight_);
	      makeFillHisto1D<TH1F,float>("hyp_highpt_ht","hyp_highpt_ht",13,80,600,ht,weight_);
	      makeFillHisto1D<TH1F,float>("hyp_highpt_met","hyp_highpt_met",10,0,250,met,weight_);
	      makeFillHisto1D<TH1F,float>("hyp_highpt_mll","hyp_highpt_mll",100,0,1000,hyp.p4().mass(),weight_);
	      makeFillHisto1D<TH1F,float>("hyp_highpt_ptll","hyp_highpt_ptll",100,0,1000,hyp.p4().pt(),weight_);
	      makeFillHisto1D<TH1F,int>("hyp_highpt_type","hyp_highpt_type",5,0,5,type,weight_);
	      makeFillHisto1D<TH1F,float>("hyp_highpt_ptlead","hyp_highpt_ptlead",40,0,200,hyp.leadLep().pt(),weight_);
	      makeFillHisto1D<TH1F,float>("hyp_highpt_pttrai","hyp_highpt_pttrai",40,0,200,hyp.traiLep().pt(),weight_);
	      makeFillHisto1D<TH1F,float>("hyp_highpt_etalead","hyp_highpt_etalead",10,0,2.5,fabs(hyp.leadLep().eta()),weight_);
	      makeFillHisto1D<TH1F,float>("hyp_highpt_etatrai","hyp_highpt_etatrai",10,0,2.5,fabs(hyp.traiLep().eta()),weight_);
	      makeFillHisto1D<TH1F,float>("hyp_highpt_mtmin","hyp_highpt_mtmin",15,0,300,mtmin,weight_);
	    }

	    if (mtmin>100.) {
	      makeFillHisto1D<TH1F,int>("hyp_highptmt_br","hyp_highptmt_br",40,0,40,br,weight_);
	      if (sr>0) {
		makeFillHisto1D<TH1F,int>("hyp_highptmt_sr","hyp_highptmt_sr",40,0,40,br,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highptmt_sr","hyp_highptmt_sr",40,0,40,sr,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highptmt_njets","hyp_highptmt_njets",8,0,8,njets,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highptmt_nbtag","hyp_highptmt_nbtag",8,0,8,nbtag,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highptmt_nbtag_pt20","hyp_highptmt_nbtag_pt20",8,0,8,nbtag_pt20,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highptmt_nbtag_pt25","hyp_highptmt_nbtag_pt25",8,0,8,nbtag_pt25,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highptmt_nbtag_pt30","hyp_highptmt_nbtag_pt30",8,0,8,nbtag_pt30,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highptmt_nbtag_pt35","hyp_highptmt_nbtag_pt35",8,0,8,nbtag_pt35,weight_);
		makeFillHisto1D<TH1F,float>("hyp_highptmt_ht","hyp_highptmt_ht",13,80,600,ht,weight_);
		makeFillHisto1D<TH1F,float>("hyp_highptmt_met","hyp_highptmt_met",10,0,250,met,weight_);
		makeFillHisto1D<TH1F,float>("hyp_highptmt_mll","hyp_highptmt_mll",100,0,1000,hyp.p4().mass(),weight_);
		makeFillHisto1D<TH1F,float>("hyp_highptmt_ptll","hyp_highptmt_ptll",100,0,1000,hyp.p4().pt(),weight_);
		makeFillHisto1D<TH1F,int>("hyp_highptmt_type","hyp_highptmt_type",5,0,5,type,weight_);
		makeFillHisto1D<TH1F,float>("hyp_highptmt_ptlead","hyp_highptmt_ptlead",40,0,200,hyp.leadLep().pt(),weight_);
		makeFillHisto1D<TH1F,float>("hyp_highptmt_pttrai","hyp_highptmt_pttrai",40,0,200,hyp.traiLep().pt(),weight_);
		makeFillHisto1D<TH1F,float>("hyp_highptmt_etalead","hyp_highptmt_etalead",10,0,2.5,fabs(hyp.leadLep().eta()),weight_);
		makeFillHisto1D<TH1F,float>("hyp_highptmt_etatrai","hyp_highptmt_etatrai",10,0,2.5,fabs(hyp.traiLep().eta()),weight_);
		makeFillHisto1D<TH1F,float>("hyp_highptmt_mtmin","hyp_highptmt_mtmin",15,0,300,mtmin,weight_);
	      }
	    }

	    if (ht400met120) {
	      makeFillHisto1D<TH1F,int>("hyp_highpthtmet_br","hyp_highpthtmet_br",40,0,40,br,weight_);
	      if (sr>0) {
		makeFillHisto1D<TH1F,int>("hyp_highpthtmet_sr","hyp_highpthtmet_sr",40,0,40,br,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highpthtmet_sr","hyp_highpthtmet_sr",40,0,40,sr,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highpthtmet_njets","hyp_highpthtmet_njets",8,0,8,njets,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highpthtmet_nbtag","hyp_highpthtmet_nbtag",8,0,8,nbtag,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highpthtmet_nbtag_pt20","hyp_highpthtmet_nbtag_pt20",8,0,8,nbtag_pt20,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highpthtmet_nbtag_pt25","hyp_highpthtmet_nbtag_pt25",8,0,8,nbtag_pt25,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highpthtmet_nbtag_pt30","hyp_highpthtmet_nbtag_pt30",8,0,8,nbtag_pt30,weight_);
		makeFillHisto1D<TH1F,int>("hyp_highpthtmet_nbtag_pt35","hyp_highpthtmet_nbtag_pt35",8,0,8,nbtag_pt35,weight_);
		makeFillHisto1D<TH1F,float>("hyp_highpthtmet_ht","hyp_highpthtmet_ht",13,80,600,ht,weight_);
		makeFillHisto1D<TH1F,float>("hyp_highpthtmet_met","hyp_highpthtmet_met",10,0,250,met,weight_);
		makeFillHisto1D<TH1F,float>("hyp_highpthtmet_mll","hyp_highpthtmet_mll",100,0,1000,hyp.p4().mass(),weight_);
		makeFillHisto1D<TH1F,float>("hyp_highpthtmet_ptll","hyp_highpthtmet_ptll",100,0,1000,hyp.p4().pt(),weight_);
		makeFillHisto1D<TH1F,int>("hyp_highpthtmet_type","hyp_highpthtmet_type",5,0,5,type,weight_);
		makeFillHisto1D<TH1F,float>("hyp_highpthtmet_ptlead","hyp_highpthtmet_ptlead",40,0,200,hyp.leadLep().pt(),weight_);
		makeFillHisto1D<TH1F,float>("hyp_highpthtmet_pttrai","hyp_highpthtmet_pttrai",40,0,200,hyp.traiLep().pt(),weight_);
		makeFillHisto1D<TH1F,float>("hyp_highpthtmet_etalead","hyp_highpthtmet_etalead",10,0,2.5,fabs(hyp.leadLep().eta()),weight_);
		makeFillHisto1D<TH1F,float>("hyp_highpthtmet_etatrai","hyp_highpthtmet_etatrai",10,0,2.5,fabs(hyp.traiLep().eta()),weight_);
		makeFillHisto1D<TH1F,float>("hyp_highpthtmet_mtmin","hyp_highpthtmet_mtmin",15,0,300,mtmin,weight_);
	      }

	      if (mtmin>100.) {
		makeFillHisto1D<TH1F,int>("hyp_highpthtmetmt_br","hyp_highpthtmetmt_br",40,0,40,br,weight_);
		if (sr>0) {
		  makeFillHisto1D<TH1F,int>("hyp_highpthtmetmt_sr","hyp_highpthtmetmt_sr",40,0,40,br,weight_);
		  makeFillHisto1D<TH1F,int>("hyp_highpthtmetmt_sr","hyp_highpthtmetmt_sr",40,0,40,sr,weight_);
		  makeFillHisto1D<TH1F,int>("hyp_highpthtmetmt_njets","hyp_highpthtmetmt_njets",8,0,8,njets,weight_);
		  makeFillHisto1D<TH1F,int>("hyp_highpthtmetmt_nbtag","hyp_highpthtmetmt_nbtag",8,0,8,nbtag,weight_);
		  makeFillHisto1D<TH1F,int>("hyp_highpthtmetmt_nbtag_pt20","hyp_highpthtmetmt_nbtag_pt20",8,0,8,nbtag_pt20,weight_);
		  makeFillHisto1D<TH1F,int>("hyp_highpthtmetmt_nbtag_pt25","hyp_highpthtmetmt_nbtag_pt25",8,0,8,nbtag_pt25,weight_);
		  makeFillHisto1D<TH1F,int>("hyp_highpthtmetmt_nbtag_pt30","hyp_highpthtmetmt_nbtag_pt30",8,0,8,nbtag_pt30,weight_);
		  makeFillHisto1D<TH1F,int>("hyp_highpthtmetmt_nbtag_pt35","hyp_highpthtmetmt_nbtag_pt35",8,0,8,nbtag_pt35,weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpthtmetmt_ht","hyp_highpthtmetmt_ht",13,80,600,ht,weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpthtmetmt_met","hyp_highpthtmetmt_met",10,0,250,met,weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpthtmetmt_mll","hyp_highpthtmetmt_mll",100,0,1000,hyp.p4().mass(),weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpthtmetmt_ptll","hyp_highpthtmetmt_ptll",100,0,1000,hyp.p4().pt(),weight_);
		  makeFillHisto1D<TH1F,int>("hyp_highpthtmetmt_type","hyp_highpthtmetmt_type",5,0,5,type,weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpthtmetmt_ptlead","hyp_highpthtmetmt_ptlead",40,0,200,hyp.leadLep().pt(),weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpthtmetmt_pttrai","hyp_highpthtmetmt_pttrai",40,0,200,hyp.traiLep().pt(),weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpthtmetmt_etalead","hyp_highpthtmetmt_etalead",10,0,2.5,fabs(hyp.leadLep().eta()),weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpthtmetmt_etatrai","hyp_highpthtmetmt_etatrai",10,0,2.5,fabs(hyp.traiLep().eta()),weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpthtmetmt_mtmin","hyp_highpthtmetmt_mtmin",15,0,300,mtmin,weight_);
		}
	      }
	    }

	    //test gen level btag (but starting from jets with |eta|<2.4)
	    int genbjets = 0;
	    int genbjets_pt20 = 0;
	    int genbjets_pt25 = 0;
	    int genbjets_pt30 = 0;
	    int genbjets_pt35 = 0;
	    int genbjets_pt40 = 0;
	    vector<Jet> jetsbmatch;
	    for (unsigned int j=0;j<alljets.size();++j) {
	      Jet jet = alljets[j];
	      if (abs(jet.mc3_id())==5) {
		genbjets++;
		jetsbmatch.push_back(jet);
		if (jet.pt()>20) {
		  genbjets_pt20++;
		  if (jet.pt()>25) genbjets_pt25++;
		  if (jet.pt()>30) genbjets_pt30++;
		  if (jet.pt()>35) genbjets_pt35++;
		  if (jet.pt()>40) genbjets_pt40++;
		  makeFillHisto1D<TH1F,float>("hyp_highpt_genbjets_csv","hyp_highpt_genbjets_csv",50,0,1.,jet.csv(),weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpt_genbjets_tche","hyp_highpt_genbjets_tche",50,-5.,5.,pfjets_trackCountingHighEffBJetTag()[jet.idx()],weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpt_genbjets_tchp","hyp_highpt_genbjets_tchp",50,-5.,5.,pfjets_trackCountingHighPurBJetTag()[jet.idx()],weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpt_genbjets_jbp","hyp_highpt_genbjets_jbp",50,0,10.,pfjets_jetBProbabilityBJetTag()[jet.idx()],weight_);
		}
	      } else {
		if (jet.pt()>20) {
		  makeFillHisto1D<TH1F,float>("hyp_highpt_gennobjets_csv","hyp_highpt_gennobjets_csv",50,0,1.,jet.csv(),weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpt_gennobjets_tche","hyp_highpt_gennobjets_tche",50,-5.,5.,pfjets_trackCountingHighEffBJetTag()[jet.idx()],weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpt_gennobjets_tchp","hyp_highpt_gennobjets_tchp",50,-5.,5.,pfjets_trackCountingHighPurBJetTag()[jet.idx()],weight_);
		  makeFillHisto1D<TH1F,float>("hyp_highpt_gennobjets_jbp","hyp_highpt_gennobjets_jbp",50,0,10.,pfjets_jetBProbabilityBJetTag()[jet.idx()],weight_);
		  /*
		    OBJ: TNamed	pfjets_combinedSecondaryVertexBJetTag	floats_pfJetMaker_pfjetscombinedSecondaryVertexBJetTag_CMS2.obj : 0 at: 0x7fb6f3a4bf20
		    OBJ: TNamed	pfjets_jetBProbabilityBJetTag	floats_pfJetMaker_pfjetsjetBProbabilityBJetTag_CMS2.obj : 0 at: 0x7fb6f3a4c540
		    OBJ: TNamed	pfjets_jetProbabilityBJetTag	floats_pfJetMaker_pfjetsjetProbabilityBJetTag_CMS2.obj : 0 at: 0x7fb6f3a4c600
		    OBJ: TNamed	pfjets_simpleSecondaryVertexHighEffBJetTag	floats_pfJetMaker_pfjetssimpleSecondaryVertexHighEffBJetTag_CMS2.obj : 0 at: 0x7fb6f3a4ca80
		    OBJ: TNamed	pfjets_trackCountingHighEffBJetTag	floats_pfJetMaker_pfjetstrackCountingHighEffBJetTag_CMS2.obj : 0 at: 0x7fb6f3a4cc40
		    OBJ: TNamed	pfjets_trackCountingHighPurBJetTag	floats_pfJetMaker_pfjetstrackCountingHighPurBJetTag_CMS2.obj : 0 at: 0x7fb6f3a4cd10
		  */
		}
	      }
	    }
	    makeFillHisto1D<TH1F,int>("hyp_highpt_genbjets","hyp_highpt_genbjets",8,0,8,genbjets,weight_);
	    makeFillHisto1D<TH1F,int>("hyp_highpt_genbjets_pt20","hyp_highpt_genbjets_pt20",8,0,8,genbjets_pt20,weight_);
	    makeFillHisto1D<TH1F,int>("hyp_highpt_genbjets_pt25","hyp_highpt_genbjets_pt25",8,0,8,genbjets_pt25,weight_);
	    makeFillHisto1D<TH1F,int>("hyp_highpt_genbjets_pt30","hyp_highpt_genbjets_pt30",8,0,8,genbjets_pt30,weight_);
	    makeFillHisto1D<TH1F,int>("hyp_highpt_genbjets_pt35","hyp_highpt_genbjets_pt35",8,0,8,genbjets_pt35,weight_);
	    makeFillHisto1D<TH1F,int>("hyp_highpt_genbjets_pt40","hyp_highpt_genbjets_pt40",8,0,8,genbjets_pt40,weight_);
	    std::sort(jetsbmatch.begin(),jetsbmatch.end(),jetptsort);
	    if (jetsbmatch.size()>0) makeFillHisto1D<TH1F,float>("hyp_highpt_genbjets_pt0","hyp_highpt_genbjets_pt0",40,0,200,jetsbmatch[0].pt(),weight_);
	    if (jetsbmatch.size()>1) makeFillHisto1D<TH1F,float>("hyp_highpt_genbjets_pt1","hyp_highpt_genbjets_pt1",40,0,200,jetsbmatch[1].pt(),weight_);
	    if (jetsbmatch.size()>2) makeFillHisto1D<TH1F,float>("hyp_highpt_genbjets_pt2","hyp_highpt_genbjets_pt2",40,0,200,jetsbmatch[2].pt(),weight_);
	    //now check the performance as a function of the btagging algo & discriminator

	  }
	  if (ac_base & 1<<LowPt) {
	    makeFillHisto1D<TH1F,int>("hyp_lowpt_sr","hyp_lowpt_sr",40,0,40,br,weight_);
	    if (sr>0) makeFillHisto1D<TH1F,int>("hyp_lowpt_sr","hyp_lowpt_sr",40,0,40,sr,weight_);
	    makeFillHisto1D<TH1F,int>("hyp_lowpt_njets","hyp_lowpt_njets",20,0,20,njets,weight_);
	    makeFillHisto1D<TH1F,int>("hyp_lowpt_nbtag","hyp_lowpt_nbtag",20,0,20,nbtag,weight_);
	    makeFillHisto1D<TH1F,float>("hyp_lowpt_ht","hyp_lowpt_ht",50,0,2000,ht,weight_);
	    makeFillHisto1D<TH1F,float>("hyp_lowpt_met","hyp_lowpt_met",50,0,500,met,weight_);
	    makeFillHisto1D<TH1F,float>("hyp_lowpt_mll","hyp_lowpt_mll",100,0,1000,hyp.p4().mass(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_lowpt_ptll","hyp_lowpt_ptll",100,0,1000,hyp.p4().pt(),weight_);
	    makeFillHisto1D<TH1F,int>("hyp_lowpt_type","hyp_lowpt_type",5,0,5,type,weight_);
	    makeFillHisto1D<TH1F,float>("hyp_lowpt_ptlead","hyp_lowpt_ptlead",20,0,200,hyp.leadLep().pt(),weight_);
	    makeFillHisto1D<TH1F,float>("hyp_lowpt_pttrai","hyp_lowpt_pttrai",20,0,200,hyp.traiLep().pt(),weight_);
	    if (makeSyncTest) cout << Form("%1d %9d %12d\t%2d\t%+2d %5.1f\t%+2d %5.1f\t%d\t%2d\t%5.1f\t%6.1f", run_, ls_, evt_, int(vetoleps.size()), hyp.leadLep().pdgId(), hyp.leadLep().pt(), 
					   hyp.traiLep().pdgId(), hyp.traiLep().pt(),  njets, nbtag, met, ht) << endl;
	  }
	  if (ac_base & 1<<VeryLowPt) {
	    makeFillHisto1D<TH1F,int>("hyp_verylowpt_sr","hyp_verylowpt_sr",40,0,40,br,weight_);
	    if (sr>0) makeFillHisto1D<TH1F,int>("hyp_verylowpt_sr","hyp_verylowpt_sr",40,0,40,sr,weight_);
	    if (makeSyncTest) cout << Form("%1d %9d %12d\t%2d\t%+2d %5.1f\t%+2d %5.1f\t%d\t%2d\t%5.1f\t%6.1f", run_, ls_, evt_, int(vetoleps.size()), hyp.leadLep().pdgId(), hyp.leadLep().pt(), 
					   hyp.traiLep().pdgId(), hyp.traiLep().pt(),  njets, nbtag, met, ht) << endl;
	  }

	  if (prefix=="ttbar") {
	    fakeStudy(hyp,ll,lt);
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

  if (fr_file) fr_file->Close();

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

void looper::runQCDtest(vector<Lep>& fobs, vector<Lep>& goodleps, int& njets, float& met) {

  if (fobs.size()>1) return; 
  if (fobs[0].pt()<ptCutHigh) return; 
  makeFillHisto1D<TH1F,int>("njets","njets",20,0,20,njets,weight_);
  for (unsigned int ipu=0;ipu<puInfo_bunchCrossing().size();ipu++)
    if (puInfo_bunchCrossing()[ipu]==0) 
      makeFillHisto1D<TH1F,float>("npu_true","npu_true",20,0,100,puInfo_trueNumInteractions()[ipu],weight_);
  if (njets==0) return;
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
  
}


void looper::runDYtest(vector<Lep>& fobs, vector<Lep>& goodleps, int& njets, float& met) {

  //check lepton selection efficiency
  for (unsigned int gp=0;gp<genps_id().size();++gp) {
    int pdgid = abs(genps_id()[gp]);
    if (pdgid!=13 && pdgid!=11) continue;
    if (genps_id_mother()[gp]!=23 && abs(genps_id_mother()[gp])!=24) continue;
    if (genps_status()[gp]!=1) continue;//is this needed?
    if (fabs(genps_p4()[gp].eta())>2.4) continue;
    if (genps_p4()[gp].pt()<5) continue;
    if (pdgid==11 && genps_p4()[gp].pt()<10) continue;
    //denominator
    makeFillHisto2D<TH2F,float>((pdgid==13?"ef_mu_den":"ef_el_den"),(pdgid==13?"ef_mu_den":"ef_el_den"),5,0.,100.,genps_p4()[gp].pt(),3,0.,3.0,fabs(genps_p4()[gp].eta()),weight_);
    //numerator
    for (unsigned int gl=0;gl<goodleps.size();++gl) {
      if (abs(goodleps[gl].pdgId())==pdgid && deltaR(goodleps[gl].p4(),genps_p4()[gp])<0.1) {
	//cout << "goodleps.size()=" << goodleps.size() << endl;
	makeFillHisto2D<TH2F,float>((pdgid==13?"ef_mu_num":"ef_el_num"),(pdgid==13?"ef_mu_num":"ef_el_num"),5,0.,100.,genps_p4()[gp].pt(),3,0.,3.0,fabs(genps_p4()[gp].eta()),weight_);
	//charge flip
	makeFillHisto2D<TH2F,float>((pdgid==13?"flip_mu_den":"flip_el_den"),(pdgid==13?"flip_mu_den":"flip_el_den"),5,0.,100.,genps_p4()[gp].pt(),3,0.,3.0,fabs(genps_p4()[gp].eta()),weight_);
	if (goodleps[gl].pdgId()==-genps_id()[gp]) makeFillHisto2D<TH2F,float>((pdgid==13?"flip_mu_num":"flip_el_num"),(pdgid==13?"flip_mu_num":"flip_el_num"),5,0.,100.,genps_p4()[gp].pt(),3,0.,3.0,fabs(genps_p4()[gp].eta()),weight_);
	break;
      }
    }
    //numerator for fobs
    for (unsigned int fo=0;fo<fobs.size();++fo) {
      if (abs(fobs[fo].pdgId())==pdgid && deltaR(fobs[fo].p4(),genps_p4()[gp])<0.1) {
	//cout << "fobs.size()=" << fobs.size() << endl;
	makeFillHisto2D<TH2F,float>((pdgid==13?"ef_mu_num_fo":"ef_el_num_fo"),(pdgid==13?"ef_mu_num_fo":"ef_el_num_fo"),5,0.,100.,genps_p4()[gp].pt(),3,0.,3.0,fabs(genps_p4()[gp].eta()),weight_);
	break;
      }
    }
    //numerator for tight id and fobs isolation
    for (unsigned int fo=0;fo<fobs.size();++fo) {
      if (abs(fobs[fo].pdgId())==pdgid && deltaR(fobs[fo].p4(),genps_p4()[gp])<0.1) {
	//cout << "fobs.size()=" << fobs.size() << endl;
	if (fabs(fobs[fo].dzPV())>0.1) continue;
	if (fabs(fobs[fo].dxyPV())>0.01) continue;
	if (pdgid==13 && fabs(fobs[fo].dxyPV())>0.005) continue;
	makeFillHisto2D<TH2F,float>((pdgid==13?"ef_mu_num_foiso":"ef_el_num_foiso"),(pdgid==13?"ef_mu_num_foiso":"ef_el_num_foiso"),5,0.,100.,genps_p4()[gp].pt(),3,0.,3.0,fabs(genps_p4()[gp].eta()),weight_);
	break;
      }
    }
    if (pdgid==11){
      for (unsigned int elidx=0;elidx<els_p4().size();++elidx) {
	if (deltaR(els_p4()[elidx],genps_p4()[gp])>0.1) continue;
	if (fabs(els_etaSC().at(elidx)) <= 1.479) {
	  makeFillHisto1D<TH1F,float>("el_fo_sietaieta_barrel","el_fo_sietaieta_barrel",50,0.,0.025,els_sigmaIEtaIEta_full5x5().at(elidx),weight_);
	}
      }
    }
  }

}


void looper::runWZtest(vector<Lep>& vetoleps, vector<Lep>& fobs, vector<Lep>& goodleps, int& njets, float& met) {

  weight_=1.;//fixme

  makeFillHisto1D<TH1F,float>("vetoleps_size","vetoleps_size",10,-0.5,9.5,vetoleps.size(),weight_);
  makeFillHisto1D<TH1F,float>("goodleps_size","goodleps_size",10,-0.5,9.5,goodleps.size(),weight_);

  if (goodleps.size()==2) {
    makeFillHisto1D<TH1F,float>("vetoleps_size_g2","vetoleps_size_g2",10,-0.5,9.5,vetoleps.size(),weight_);
    if (goodleps[0].charge()!=goodleps[1].charge()) return;
    makeFillHisto1D<TH1F,float>("vetoleps_size_g2ss","vetoleps_size_g2ss",10,-0.5,9.5,vetoleps.size(),weight_);

    bool noZmass = true;
    for (unsigned int gl=0;gl<goodleps.size();++gl) {
      for (unsigned int vl=0;vl<vetoleps.size();++vl) {
	if ( fabs(ROOT::Math::VectorUtil::DeltaR(goodleps[gl].p4(),vetoleps[vl].p4()))<0.001 ) continue;
	if ( goodleps[gl].pdgId()!=-vetoleps[vl].pdgId() ) continue; 
	float mll = (goodleps[gl].p4()+vetoleps[vl].p4()).mass();
	if (fabs(mll-91.2)<15) {
	  noZmass = false;
	  break;
	}
      }
    }
    if (vetoleps.size()>2) makeFillHisto1D<TH1F,float>("noZmass_v3g2ss","noZmass_v3g2ss",2,0,2,noZmass,weight_);

    if (vetoleps.size()==2 || noZmass) {

      makeFillHisto1D<TH1F,float>("fobs_size_noZv3g2ss","fobs_size_noZv3g2ss",10,-0.5,9.5,fobs.size(),weight_);
      if (fobs.size()>2) return;

      cout << "new gps" << endl;

      float mz = -999.;
      float mw = -999.;
      float maxAbsEta = -1;
      float minPt = 999999.;
      vector<unsigned int> gpidv;
      for (unsigned int gp=0;gp<genps_id().size();++gp) {
	int pdgid = abs(genps_id()[gp]);
	if (pdgid==23) {
	  //cout << "Z status=" << genps_status()[gp] << " mother=" << genps_id_mother()[gp]  << endl;
	  if (genps_status()[gp]==22) mz = genps_p4()[gp].mass();
	}
	if (pdgid==24) {
	  //cout << "W status=" << genps_status()[gp] << " mother=" << genps_id_mother()[gp]  << endl;
	  if (genps_status()[gp]==22) mw = genps_p4()[gp].mass();
	}
	if (pdgid!=13 && pdgid!=11 && pdgid!=15) continue;
	if (pdgid==11) cout << "electron status=" << genps_status()[gp] << " mother=" << genps_id_mother()[gp] << " pt=" << genps_p4()[gp].pt() << endl;
	if (pdgid==13) cout << "muon status=" << genps_status()[gp] << " mother=" << genps_id_mother()[gp] << " pt=" << genps_p4()[gp].pt() << endl;
	if (pdgid==15) cout << "tau status=" << genps_status()[gp] << " mother=" << genps_id_mother()[gp] << " pt=" << genps_p4()[gp].pt() << endl;
	if (genps_id_mother()[gp]!=23 && abs(genps_id_mother()[gp])!=24) continue;
	if (pdgid==11 && genps_status()[gp]!=1) continue;//is this needed?
	if (pdgid==13 && genps_status()[gp]!=1) continue;//is this needed?
	if (pdgid==15 && genps_status()[gp]!=2) continue;//is this needed?
	gpidv.push_back(gp);
	if (fabs(genps_p4()[gp].eta())>maxAbsEta) maxAbsEta=fabs(genps_p4()[gp].eta());
	if (genps_p4()[gp].pt()<minPt) minPt=genps_p4()[gp].pt();
      }

      if (gpidv.size()==1) {
	//cout << "size1" << endl;
	for (unsigned int gp=0;gp<genps_id().size();++gp) {
	  int pdgid = abs(genps_id()[gp]);
	  //cout << "id=" << genps_id()[gp] << " status=" << genps_status()[gp] << " mother=" << genps_id_mother()[gp] << " pt=" << genps_p4()[gp].pt() << endl;
	  if (genps_status()[gp]!=1 && genps_status()[gp]!=2) continue;
	  if (abs(genps_id_mother()[gp])>=24) continue;//Ws are ok (size=1), don't want other resonances, but somehow Z leptons don't have mother id=23
	  if ((pdgid==11 && genps_status()[gp]==1) || (pdgid==13 && genps_status()[gp]==1) || (pdgid==15 && genps_status()[gp]==2)) gpidv.push_back(gp);
	}
      }
      if (gpidv.size()>3) {
	cout << "size>3 is " << gpidv.size() << endl;
	for (unsigned int gp=0;gp<gpidv.size() && gpidv.size()>3;++gp) {
	  cout << "id=" << genps_id()[gpidv[gp]] << " status=" << genps_status()[gpidv[gp]] << " mother=" << genps_id_mother()[gpidv[gp]] << " pt=" << genps_p4()[gpidv[gp]].pt() << endl;
	  if (abs(genps_id()[gpidv[gp]])==15) {//extra leps can be duplicates from tau decays
	    gpidv.erase(gpidv.begin()+gp);
	    cout << "erased position=" << gp << " new size=" << gpidv.size() << endl;
	    gp--;
	  }
	}
      }
      if (gpidv.size()>3) {
	cout << "sorting size>3 is " << gpidv.size() << endl;
	//sort by pt and erase exceeding items
	std::sort(gpidv.begin(),gpidv.end(),ptsort);
	gpidv.erase(gpidv.begin()+3,gpidv.end());
	for (unsigned int gp=0;gp<gpidv.size();++gp) {
	  cout << "id=" << genps_id()[gpidv[gp]] << " status=" << genps_status()[gpidv[gp]] 
	       << " mother=" << genps_id_mother()[gpidv[gp]] << " pt=" << genps_p4()[gpidv[gp]].pt() << endl;
	}     
      }
      if (gpidv.size()!=3) {
	cout << "anomalous size=" << gpidv.size() << endl;
	for (unsigned int gp=0;gp<gpidv.size();++gp) {
	  cout << "sel id=" << genps_id()[gpidv[gp]] << " status=" << genps_status()[gpidv[gp]] << " mother=" << genps_id_mother()[gpidv[gp]] << " pt=" << genps_p4()[gpidv[gp]].pt() << endl;
	}
	for (unsigned int gp=0;gp<genps_id().size();++gp) {
	  cout << "gps id=" << genps_id()[gp] << " status=" << genps_status()[gp] << " mother=" << genps_id_mother()[gp] << " pt=" << genps_p4()[gp].pt() << endl;
	}
      }

      int glepfromZ = -1;
      for (unsigned int gl=0;gl<goodleps.size();++gl) {
	//easy to check that it is not from W
	if (isFromW(goodleps[gl])==0 ||goodleps[gl].mc_motherid()==23) {
	  glepfromZ=gl; 
	  break;
	}
      }

      makeFillHisto1D<TH1F,float>("gps_size_g2f2","gps_size_g2f2",10,-0.5,9.5,gpidv.size(),weight_);
      makeFillHisto1D<TH1F,float>("mz_g2f2","mz_g2f2",100,0.,200,mz,weight_);
      makeFillHisto1D<TH1F,float>("mw_g2f2","mw_g2f2",100,0.,200,mw,weight_);
      if (gpidv.size()==3) {
	makeFillHisto1D<TH1F,float>("maxAbsEta_g2f2s3","maxAbsEta_g2f2s3",100,0.,5.0,maxAbsEta,weight_);
	if (maxAbsEta<2.4) makeFillHisto1D<TH1F,float>("minPt_g2f2s3eta24","minPt_g2f2s3eta24",100,0.,200.0,minPt,weight_);
	else return;

	//match goodleps and vetoleps with genps
	cout << "check duplicates" << endl;
	int unmatches = 0;
	for (unsigned int gp=0;gp<gpidv.size();++gp) {
	  unsigned int gpid = gpidv[gp];
	  int pdgid = abs(genps_id()[gpid]);
	  /*
	    cout << "sel id=" << genps_id()[gpid] << " status=" << genps_status()[gpid] 
	    << " mother=" << genps_id_mother()[gpid] << " pt=" << genps_p4()[gpid].pt() << endl;
	  */
	  bool glmatch = false;
	  for (unsigned int gl=0;gl<goodleps.size();++gl) {
	    if (goodleps[gl].mc3idx()==int(gpid) || 
		fabs(ROOT::Math::VectorUtil::DeltaR(goodleps[gl].mc_p4(),genps_p4()[gpid]))<0.1 || 
		fabs(ROOT::Math::VectorUtil::DeltaR(goodleps[gl].p4(),genps_p4()[gpid]))<0.1) {
	      cout << "match lep and gp, pt=" << goodleps[gl].mc_p4().pt() << " " << genps_p4()[gpid].pt() << endl;
	      glmatch = true;
	    }
	  }
	  if (glmatch==0) {
	    unmatches++;
	    int pfidx = -1;
	    float mindpt = 0.2;
	    for (unsigned int pfi=0; pfi<pfcands_p4().size(); ++pfi){
	      if ( pfcands_charge()[pfi]==0 ) continue;
	      float dR = fabs(ROOT::Math::VectorUtil::DeltaR(pfcands_p4()[pfi],genps_p4()[gpid]));
	      float dpt = fabs((pfcands_p4()[pfi].pt()-genps_p4()[gpid].pt())/genps_p4()[gpid].pt());
	      if ( dR<0.1 && dpt<mindpt ) {
		pfidx=pfi;
		mindpt=dpt;
		//break;
	      }
	    }
	    makeFillHisto1D<TH1F,float>("pdgId_noGL","pdgId_noGL",20,0.,20,pdgid,weight_);
	    makeFillHisto1D<TH1F,float>("momId_noGL","momId_noGL",30,0.,30,abs(genps_id_mother()[gpid]),weight_);
	    if (pdgid==15) continue;
	    if (pfidx>=0) {
	      makeFillHisto1D<TH1F,float>("deltapt_pfmatch","deltapt_pfmatch",100,-1.,1,(pfcands_p4()[pfidx].pt()-genps_p4()[gpid].pt())/genps_p4()[gpid].pt(),weight_);
	      makeFillHisto1D<TH1F,float>("pfid_pfmatch","pfid_pfmatch",250,0,250,abs(pfcands_particleId()[pfidx]),weight_);
	    } else {
	      if (pdgid==11) makeFillHisto2D<TH2F,float>("pteta_el_noPFnoGLeta24","pteta_el_noPFnoGLeta24",5,0.,50.0,genps_p4()[gpid].pt(),5,0.,2.5,fabs(genps_p4()[gpid].eta()),weight_);
	      if (pdgid==13) makeFillHisto2D<TH2F,float>("pteta_mu_noPFnoGLeta24","pteta_mu_noPFnoGLeta24",5,0.,50.0,genps_p4()[gpid].pt(),5,0.,2.5,fabs(genps_p4()[gpid].eta()),weight_);
	    }
	    if (pdgid==11) {
	      makeFillHisto1D<TH1F,float>("eta_el_noGL","eta_el_noGL",100,0.,5.0,fabs(genps_p4()[gpid].eta()),weight_);
	      if (fabs(genps_p4()[gpid].eta())<2.4) {
		makeFillHisto1D<TH1F,float>("pt_el_noGLeta24","pt_el_noGLeta24",100,0.,200.0,genps_p4()[gpid].pt(),weight_);
		if (pfidx>=0&&glepfromZ>=0) {
		  //makeFillHisto1D<TH1F,float>("mll_elpf_noGLeta24","mll_elpf_noGLeta24",40,0.,200.0,(goodleps[glepfromZ].p4()+pfcands_p4()[pfidx]).mass(),weight_);
		  float mass0 = (goodleps[0].p4()+pfcands_p4()[pfidx]).mass();
		  float mass1 = (goodleps[1].p4()+pfcands_p4()[pfidx]).mass();
		  makeFillHisto1D<TH1F,float>("mll_elpf_noGLeta24","mll_elpf_noGLeta24",40,0.,200.0,(fabs(mass0-91.2)<fabs(mass1-91.2)?mass0:mass1),weight_);
		  makeFillHisto1D<TH1F,float>("mz_elpf_noGLeta24","mz_elpf_noGLeta24",40,0.,200.0,mz,weight_);
		}
	      }
	    }		  
	    if (pdgid==13) {
	      makeFillHisto1D<TH1F,float>("eta_mu_noGL","eta_mu_noGL",100,0.,5.0,fabs(genps_p4()[gpid].eta()),weight_);
	      if (fabs(genps_p4()[gpid].eta())<2.4) {
		makeFillHisto1D<TH1F,float>("pt_mu_noGLeta24","pt_mu_noGLeta24",100,0.,200.0,genps_p4()[gpid].pt(),weight_);
		if (pfidx>=0&&glepfromZ>=0) {
		  //makeFillHisto1D<TH1F,float>("mll_mupf_noGLeta24","mll_mupf_noGLeta24",40,0.,200.0,(goodleps[glepfromZ].p4()+pfcands_p4()[pfidx]).mass(),weight_);
		  float mass0 = (goodleps[0].p4()+pfcands_p4()[pfidx]).mass();
		  float mass1 = (goodleps[1].p4()+pfcands_p4()[pfidx]).mass();
		  makeFillHisto1D<TH1F,float>("mll_mupf_noGLeta24","mll_mupf_noGLeta24",40,0.,200.0,(fabs(mass0-91.2)<fabs(mass1-91.2)?mass0:mass1),weight_);
		  makeFillHisto1D<TH1F,float>("mz_mupf_noGLeta24","mz_mupf_noGLeta24",40,0.,200.0,mz,weight_);
		}
	      }
	    }		  
	  }
		
	}
	if (unmatches>1) {
	  cout << "duplcates, unmatches=" << unmatches << endl;

	  for (unsigned int gl=0;gl<goodleps.size();++gl) {
	    cout << "good lep pt=" << goodleps[gl].pt() << " eta=" << goodleps[gl].eta() << " phi=" << goodleps[gl].p4().phi() << " mcpt=" << goodleps[gl].mc_p4().pt() << endl;
	  }

	  for (unsigned int gp=0;gp<gpidv.size();++gp) {
	    cout << "id=" << genps_id()[gpidv[gp]] << " status=" << genps_status()[gpidv[gp]] << " mother=" << genps_id_mother()[gpidv[gp]] 
		 << " pt=" << genps_p4()[gpidv[gp]].pt() << " eta=" << genps_p4()[gpidv[gp]].eta() << " phi=" << genps_p4()[gpidv[gp]].phi() << endl;
	  }		
	  //return -1;

	}

      }

    }
  }

  // makeFillHisto1D<TH1F,int>("hyp_charge","hyp_charge",7,-3.5,3.5,hyp.charge(),weight_);
  // makeFillHisto1D<TH1F,float>("hyp_mll","hyp_mll",100,0,1000,hyp.p4().mass(),weight_);
  // makeFillHisto1D<TH1F,float>("hyp_ptll","hyp_ptll",100,0,1000,hyp.p4().pt(),weight_);
  // makeFillHisto1D<TH1F,int>("hyp_njets","hyp_njets",20,0,20,njets,weight_);
  // makeFillHisto1D<TH1F,int>("hyp_nbtag","hyp_nbtag",20,0,20,nbtag,weight_);
  // makeFillHisto1D<TH1F,float>("hyp_ht","hyp_ht",50,0,2000,ht,weight_);
  // makeFillHisto1D<TH1F,float>("hyp_met","hyp_met",50,0,500,met,weight_);

  // int type = -1;
  // if (abs(hyp.traiLep().pdgId())==13 && abs(hyp.leadLep().pdgId())==13) type=0;
  // if (abs(hyp.traiLep().pdgId())==13 && abs(hyp.leadLep().pdgId())==11) type=1;
  // if (abs(hyp.traiLep().pdgId())==11 && abs(hyp.leadLep().pdgId())==13) type=2;
  // if (abs(hyp.traiLep().pdgId())==11 && abs(hyp.leadLep().pdgId())==11) type=3;
  // makeFillHisto1D<TH1F,int>("hyp_type","hyp_type",5,0,5,type,weight_);

}

void looper::fakeStudy(DilepHyp& hyp, TString& ll, TString& lt) {

  bool isLeadPrompt = 0;
  bool isTrailPrompt = 0;
  int leadType = -1;
  int trailType = -1;
	    
  if ( !isFromW(hyp.traiLep()) ) {
    //check fakes (mother not W)
	      
    makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_mc","hyp_ss_faketrail_"+lt+"_mc",11001,-5500.5,5500.5,hyp.traiLep().mc_id(),weight_);
    makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_mc_mother","hyp_ss_faketrail_"+lt+"_mc_mother",11001,-5500.5,5500.5,hyp.traiLep().mc_motherid(),weight_);
    makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_mc3","hyp_ss_faketrail_"+lt+"_mc3",11001,-5500.5,5500.5,hyp.traiLep().mc3_id(),weight_);
    makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_mc3_mother","hyp_ss_faketrail_"+lt+"_mc3_mother",11001,-5500.5,5500.5,hyp.traiLep().mc3_motherid(),weight_);
	      
    makeFillHisto1D<TH1F,float>("hyp_ss_faketrail_"+lt+"_relIso03","hyp_ss_faketrail_"+lt+"_relIso03",100,0.,1.,hyp.traiLep().relIso03(),weight_);
	      
    if (isFromB(hyp.traiLep()) ) trailType=FakeB;
    else if (isFromC(hyp.traiLep()) ) trailType=FakeC;
    else if (isFromLight(hyp.traiLep())) trailType=FakeLightTrue;
    else if (isFromLightFake(hyp.traiLep())) trailType=FakeLightFake;
    else if (hyp.traiLep().mc_id()==22 && hyp.traiLep().mc_p4().pt()>hyp.traiLep().pt()/2.) {
      trailType=FakeHiPtGamma;
      makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_hiptgamma_mc","hyp_ss_faketrail_"+lt+"_hiptgamma_mc",11001,-5500.5,5500.5,hyp.traiLep().mc_id(),weight_);
      makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_hiptgamma_mc_mother","hyp_ss_faketrail_"+lt+"_hiptgamma_mc_mother",11001,-5500.5,5500.5,hyp.traiLep().mc_motherid(),weight_);
      makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_hiptgamma_mc3","hyp_ss_faketrail_"+lt+"_hiptgamma_mc3",11001,-5500.5,5500.5,hyp.traiLep().mc3_id(),weight_);
      makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_hiptgamma_mc3_mother","hyp_ss_faketrail_"+lt+"_hiptgamma_mc3_mother",11001,-5500.5,5500.5,hyp.traiLep().mc3_motherid(),weight_);
      makeFillHisto1D<TH1F,float>("hyp_ss_faketrail_"+lt+"_hiptgamma_pt","hyp_ss_faketrail_"+lt+"_hiptgamma_pt",50,0,100,hyp.traiLep().pt(),weight_);
      makeFillHisto1D<TH1F,float>("hyp_ss_faketrail_"+lt+"_hiptgamma_ptmc","hyp_ss_faketrail_"+lt+"_hiptgamma_ptmc",50,0,100,hyp.traiLep().mc_p4().pt(),weight_);
      makeFillHisto1D<TH1F,float>("hyp_ss_faketrail_"+lt+"_hiptgamma_deltaptoverpt","hyp_ss_faketrail_"+lt+"_hiptgamma_deltaptoverpt",100,-0.5,0.5,(hyp.traiLep().pt()-hyp.traiLep().mc_p4().pt())/hyp.traiLep().mc_p4().pt(),weight_);
    } else  {
      trailType=FakeUnknown;
      if (hyp.traiLep().mc_id()==22 && hyp.traiLep().mc_p4().pt()<1.) {
	makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_category","hyp_ss_faketrail_"+lt+"_category",End,0,End,FakeLowPtGamma,weight_);
	makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_lowptgamma_mc","hyp_ss_faketrail_"+lt+"_lowptgamma_mc",11001,-5500.5,5500.5,hyp.traiLep().mc_id(),weight_);
	makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_lowptgamma_mc_mother","hyp_ss_faketrail_"+lt+"_lowptgamma_mc_mother",11001,-5500.5,5500.5,hyp.traiLep().mc_motherid(),weight_);
	makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_lowptgamma_mc3","hyp_ss_faketrail_"+lt+"_lowptgamma_mc3",11001,-5500.5,5500.5,hyp.traiLep().mc3_id(),weight_);
	makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_lowptgamma_mc3_mother","hyp_ss_faketrail_"+lt+"_lowptgamma_mc3_mother",11001,-5500.5,5500.5,hyp.traiLep().mc3_motherid(),weight_);
	makeFillHisto1D<TH1F,float>("hyp_ss_faketrail_"+lt+"_category_lowptgamma_pt","hyp_ss_faketrail_"+lt+"_category_lowptgamma_pt",50,0,100,hyp.traiLep().pt(),weight_);
      } else if (hyp.traiLep().mc_id()==-9999 && hyp.traiLep().mc_motherid()==-9999 && hyp.traiLep().mc3_id()==-9999 && hyp.traiLep().mc3_motherid()==-9999) {
	makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_category","hyp_ss_faketrail_"+lt+"_category",End,0,End,All9999,weight_);
	makeFillHisto1D<TH1F,float>("hyp_ss_faketrail_"+lt+"_category_all9999_pt","hyp_ss_faketrail_"+lt+"_category_all9999_pt",50,0,100,hyp.traiLep().pt(),weight_);
      }
      else {
	makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_category","hyp_ss_faketrail_"+lt+"_category",End,0,End,Other,weight_);
	makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_other_mc","hyp_ss_faketrail_"+lt+"_other_mc",11001,-5500.5,5500.5,hyp.traiLep().mc_id(),weight_);
	makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_other_mc_mother","hyp_ss_faketrail_"+lt+"_other_mc_mother",11001,-5500.5,5500.5,hyp.traiLep().mc_motherid(),weight_);
	makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_other_mc3","hyp_ss_faketrail_"+lt+"_other_mc3",11001,-5500.5,5500.5,hyp.traiLep().mc3_id(),weight_);
	makeFillHisto1D<TH1F,int>("hyp_ss_faketrail_"+lt+"_other_mc3_mother","hyp_ss_faketrail_"+lt+"_other_mc3_mother",11001,-5500.5,5500.5,hyp.traiLep().mc3_motherid(),weight_);
	makeFillHisto1D<TH1F,float>("hyp_ss_faketrail_"+lt+"_category_other_pt","hyp_ss_faketrail_"+lt+"_category_other_pt",50,0,100,hyp.traiLep().pt(),weight_);
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
	      
    makeFillHisto1D<TH1F,int>("hyp_ss_fakelead_"+ll+"_mc","hyp_ss_fakelead_"+ll+"_mc",11001,-5500.5,5500.5,hyp.leadLep().mc_id(),weight_);
    makeFillHisto1D<TH1F,int>("hyp_ss_fakelead_"+ll+"_mc_mother","hyp_ss_fakelead_"+ll+"_mc_mother",11001,-5500.5,5500.5,hyp.leadLep().mc_motherid(),weight_);
    makeFillHisto1D<TH1F,int>("hyp_ss_fakelead_"+ll+"_mc3","hyp_ss_fakelead_"+ll+"_mc3",11001,-5500.5,5500.5,hyp.leadLep().mc3_id(),weight_);
    makeFillHisto1D<TH1F,int>("hyp_ss_fakelead_"+ll+"_mc3_mother","hyp_ss_fakelead_"+ll+"_mc3_mother",11001,-5500.5,5500.5,hyp.leadLep().mc3_motherid(),weight_);

    makeFillHisto1D<TH1F,float>("hyp_ss_fakelead_"+lt+"_relIso03","hyp_ss_fakelead_"+lt+"_relIso03",100,0.,1.,hyp.leadLep().relIso03(),weight_);
	      
    if (isFromB(hyp.leadLep()) ) leadType=FakeB;
    else if (isFromC(hyp.leadLep()) ) leadType=FakeC;
    else if (isFromLight(hyp.leadLep())) leadType=FakeLightTrue;
    else if (isFromLightFake(hyp.leadLep())) leadType=FakeLightFake;
    else if (hyp.leadLep().mc_id()==22 && hyp.leadLep().mc_p4().pt()>hyp.leadLep().pt()/2.) leadType=FakeHiPtGamma;
    else {
      leadType=FakeUnknown;
      if (hyp.leadLep().mc_id()==22 && hyp.leadLep().mc_p4().pt()<1.) 
	makeFillHisto1D<TH1F,int>("hyp_ss_fakelead_"+ll+"_category","hyp_ss_fakelead_"+ll+"_category",End,0,End,FakeLowPtGamma,weight_);
      else if (hyp.leadLep().mc_id()==-9999 && hyp.leadLep().mc_motherid()==-9999 && hyp.leadLep().mc3_id()==-9999 && hyp.leadLep().mc3_motherid()==-9999)
	makeFillHisto1D<TH1F,int>("hyp_ss_fakelead_"+ll+"_category","hyp_ss_fakelead_"+ll+"_category",End,0,End,All9999,weight_);
      else {
	makeFillHisto1D<TH1F,int>("hyp_ss_fakelead_"+ll+"_category","hyp_ss_fakelead_"+ll+"_category",End,0,End,Other,weight_);
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
  if (trailType==PromptWS) makeFillHisto1D<TH1F,int>("hyp_fliptype","hyp_fliptype",3,0,3,1,weight_);
  else if (leadType==PromptWS) makeFillHisto1D<TH1F,int>("hyp_fliptype","hyp_fliptype",3,0,3,2,weight_);
  else makeFillHisto1D<TH1F,int>("hyp_fliptype","hyp_fliptype",3,0,3,0,weight_);
	    
  if (hyp.leadLep().mc_id()*hyp.traiLep().mc_id()<0) makeFillHisto1D<TH1F,int>("hyp_ismcss","hyp_ismcss",3,0,3,0,weight_);
  else makeFillHisto1D<TH1F,int>("hyp_ismcss","hyp_ismcss",3,0,3,1,weight_);
	    
  if (isLeadPrompt && isTrailPrompt) {
    int leadprompttype = -1;
    if (leadType==0 && trailType==0) {
      leadprompttype=0;
      //cout << "UNEXPECTED DOUBLE PROMPT" << endl;
    } else if (leadType==1 || trailType==1) leadprompttype=1;
    else if (leadType==2 || trailType==2) leadprompttype=1;
    else  leadprompttype=3;
    makeFillHisto1D<TH1F,int>("hyp_leadprompttype","hyp_leadprompttype",5,0,5,leadprompttype,weight_);
  }


}
