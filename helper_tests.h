#include <vector>

class TFile;
class TString;
class looper;
struct Lep;
struct Jet;
struct DilepHyp;

namespace tests{

  void runQCDtest( looper* loop, float& weight_, std::vector<Lep>& fobs, std::vector<Lep>& goodleps, int& njets, float& met);
  void runDYtest( looper* loop, float& weight_, std::vector<Lep>& fobs, std::vector<Lep>& goodleps, int& njets, float& met);
  void runWZtest( looper* loop, float& weight_, std::vector<Lep>& vetoleps, std::vector<Lep>& fobs, std::vector<Lep>& goodleps, int& njets, float& met);
  void fakeStudy( looper* loop, float& weight_, DilepHyp& hyp, TString& ll, TString& lt );

  void testLepIdFailMode( looper* loop, float& weight_, std::vector<Lep>& fobs );

  void computeFakeRateAndClosure( looper* loop, float& weight_, std::vector<Lep>& fobs, std::vector<Lep>& goodleps, TFile* fr_file );

  void testPtRel( looper* loop, float& weight_, DilepHyp& hyp, std::vector<Jet>& lepjets, TString& ll, TString& lt );

  void testBtag( looper* loop, float& weight_, std::vector<Jet>& alljets );

  void testLeptonIdVariables( looper* loop, float& weight_, DilepHyp& hyp, TString& ll, TString& lt );

  void makeSRplots( looper* loop, float& weight_, TString label, int& br, int& sr, DilepHyp& hyp, 
		    float& ht, float& met, float& mtmin, int& type, std::vector<Lep>& goodleps, std::vector<Lep>& fobs, 
		    std::vector<Lep>& vetoleps, std::vector<Jet>& jets, std::vector<Jet>& alljets, std::vector<Jet>& btags, 
		    TString& ll, TString& lt );
}
