#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"

using namespace ROOT;
using namespace ROOT::VecOps;


float deltaPhi(float phi1, float phi2) {                                                        
  float result = phi1 - phi2;
  while (result > float(M_PI)) result -= float(2*M_PI);
  while (result <= -float(M_PI)) result += float(2*M_PI);
  return result;
}

float deltaR2(float eta1, float phi1, float eta2, float phi2) {
  float deta = std::abs(eta1-eta2);
  float dphi = deltaPhi(phi1,phi2);
  return deta*deta + dphi*dphi;
}

float deltaR(float eta1, float phi1, float eta2, float phi2) {
  return std::sqrt(deltaR2(eta1,phi1,eta2,phi2));
}

ROOT::VecOps::RVec<std::pair<int,int>> CreateTPPair(ROOT::VecOps::RVec<Int_t> &Muon_charge,ROOT::VecOps::RVec<Int_t> &isProbe, ROOT::VecOps::RVec<Int_t> &isTag, ROOT::VecOps::RVec<Bool_t> &isTriggeredMuon, ROOT::VecOps::RVec<Bool_t> &isGenMatchedMuon ){
  ROOT::VecOps::RVec<std::pair<int,int>> TP_pairs;
  for(int iLep1=0; iLep1<Muon_charge.size();iLep1++){
    //if(!isInAcceptance[iLep1]) continue;
    if(!isTag[iLep1]) continue;
    if(!isTriggeredMuon[iLep1]) continue;
    if(!isGenMatchedMuon[iLep1]) continue;




    for(int iLep2=0; iLep2<Muon_charge.size(); iLep2++){
      if (iLep2==iLep1) continue;
      if(!isProbe[iLep2]) continue;
      if(!isTriggeredMuon[iLep2]) continue;
      if(!isGenMatchedMuon[iLep2]) continue;
      //#if(Muon_charge[iLep1] == Muon_charge[iLep2]) continue;

      std::pair<int,int> TP_pair = std::make_pair(iLep1,iLep2); 
      TP_pairs.push_back(TP_pair);
    }

  }
  return TP_pairs;
}

ROOT::VecOps::RVec<Bool_t> hasTriggerMatch(ROOT::VecOps::RVec<Float_t> &Muon_eta, ROOT::VecOps::RVec<Float_t> &Muon_phi, ROOT::VecOps::RVec<Int_t> &TrigObj_id,ROOT::VecOps::RVec<Float_t> &TrigObj_pt,ROOT::VecOps::RVec<Float_t> & TrigObj_l1pt, ROOT::VecOps::RVec<Float_t> & TrigObj_l2pt, ROOT::VecOps::RVec<Int_t> &TrigObj_filterBits, ROOT::VecOps::RVec<Float_t> &TrigObj_eta, ROOT::VecOps::RVec<Float_t> &TrigObj_phi){
  ROOT::VecOps::RVec<Bool_t> TriggerMatch;
  for (int iMuon = 0; iMuon<Muon_eta.size(); iMuon++ ){
    bool hasTrigMatch = false;
    for (unsigned int iTrig=0; iTrig<TrigObj_id.size(); ++iTrig){
      if (TrigObj_id[iTrig]  != 13 ) continue;
      if (TrigObj_pt[iTrig]   < 24.) continue;
      if (TrigObj_l1pt[iTrig] < 22.) continue;
      if (! (( TrigObj_filterBits[iTrig] & 8) || (TrigObj_l2pt[iTrig] > 10. && (TrigObj_filterBits[iTrig] & 2) )) ) continue;
      if (deltaR(Muon_eta[iMuon], Muon_phi[iMuon], TrigObj_eta[iTrig], TrigObj_phi[iTrig]) < 0.3) {
	hasTrigMatch = true;
	break;
      }
    }
    TriggerMatch.push_back(hasTrigMatch);
  }
  return TriggerMatch;
}

ROOT::VecOps::RVec<Bool_t> hasGenMatch(ROOT::VecOps::RVec<Int_t> &GenPart_pdgId, ROOT::VecOps::RVec<Int_t> &GenPart_status,ROOT::VecOps::RVec<Int_t> &GenPart_statusFlags,ROOT::VecOps::RVec<Float_t> &GenPart_eta,ROOT::VecOps::RVec<Float_t> &GenPart_phi,ROOT::VecOps::RVec<Float_t> &Muon_eta, ROOT::VecOps::RVec<Float_t> &Muon_phi){
  ROOT::VecOps::RVec<Bool_t> isGenMatched;  
  for(int iMuon=0;iMuon<Muon_eta.size();iMuon++){
    bool matchMC = false;   
    float mcmatch_tmp_dr = 999.;
    for(unsigned int iGen=0; iGen<GenPart_pdgId.size(); iGen++){
      if ( abs(GenPart_pdgId[iGen])==13 &&  GenPart_status[iGen] == 1 && (GenPart_statusFlags[iGen] & 1)  ) {
	if (deltaR(Muon_eta[iMuon], Muon_phi[iMuon], GenPart_eta[iGen], GenPart_phi[iGen]) < mcmatch_tmp_dr){
	  mcmatch_tmp_dr  = deltaR(Muon_eta[iMuon], Muon_phi[iMuon], GenPart_eta[iGen], GenPart_phi[iGen]);
	  //truePt     = GenPart_pt[ii];
	  //trueEta    = GenPart_eta[ii];
	  //trueCharge = GenPart_charge[ii];
	} 
      }
    }
    if (mcmatch_tmp_dr < 0.1) matchMC = true;
    isGenMatched.push_back(matchMC);
  }
  return isGenMatched;
}



ROOT::VecOps::RVec<Float_t> getTPmass(ROOT::VecOps::RVec<std::pair<int,int>> TPPairs, ROOT::VecOps::RVec<Float_t> &Muon_pt,ROOT::VecOps::RVec<Float_t> &Muon_eta,ROOT::VecOps::RVec<Float_t> &Muon_phi){
  ROOT::VecOps::RVec<Float_t> TPMass;
  for (int i=0;i<TPPairs.size();i++){
    std::pair<int,int> TPPair = TPPairs.at(i);
    int tag_index = TPPair.first;
    int probe_index = TPPair.second;
    TLorentzVector tagLep(0,0,0,0);
    tagLep.SetPtEtaPhiM(Muon_pt[tag_index],Muon_eta[tag_index],Muon_phi[tag_index],0);
    TLorentzVector probeLep(0,0,0,0);
    probeLep.SetPtEtaPhiM(Muon_pt[probe_index],Muon_eta[probe_index],Muon_phi[probe_index],0);
    TPMass.push_back((tagLep+probeLep).M());        
  }
  return TPMass;
}


ROOT::VecOps::RVec<Float_t> getVariables(ROOT::VecOps::RVec<std::pair<int,int>> TPPairs, ROOT::VecOps::RVec<Float_t> &Muon_variable, float option /*1 for tag and 2 for probe*/){
  ROOT::VecOps::RVec<Float_t> Variables;
  for (int i=0;i<TPPairs.size();i++){
    std::pair<int,int> TPPair = TPPairs.at(i);
    int index; 
    if (option==1) index = TPPair.first;
    else if (option==2) index = TPPair.second;
    float variable = Muon_variable[index];
    Variables.push_back(variable);
  }
  return Variables;
}


ROOT::VecOps::RVec<Int_t> getVariables(ROOT::VecOps::RVec<std::pair<int,int>> TPPairs, ROOT::VecOps::RVec<Int_t> &Muon_variable,float option /*1 for tag and 2 for probe*/){
  ROOT::VecOps::RVec<Int_t> Variables;
  for (int i=0;i<TPPairs.size();i++){
    std::pair<int,int> TPPair = TPPairs.at(i);
    int index;
    if(option==1) index = TPPair.first;
    else if (option==2) index = TPPair.second;
    float variable = Muon_variable[index];
    Variables.push_back(variable);
  }
  return Variables;
}

ROOT::VecOps::RVec<Bool_t> getVariables(ROOT::VecOps::RVec<std::pair<int,int>> TPPairs, ROOT::VecOps::RVec<Bool_t> &Muon_variable,float option /*1 for tag and 2 for probe*/){
  ROOT::VecOps::RVec<Bool_t> Variables;
  for (int i=0;i<TPPairs.size();i++){
    std::pair<int,int> TPPair = TPPairs.at(i);
    int index;
    if(option==1) index = TPPair.first;
    else if (option==2) index = TPPair.second;
    float variable = Muon_variable[index];
    Variables.push_back(variable);
  }
  return Variables;
}

ROOT::VecOps::RVec<Bool_t> isOS(ROOT::VecOps::RVec<std::pair<int,int>> TPPairs,ROOT::VecOps::RVec<Bool_t> Muon_charge){
  ROOT::VecOps::RVec<Bool_t> isOS;
  for (int i=0;i<TPPairs.size();i++){
    std::pair<int,int> TPPair = TPPairs.at(i);
    int  tag_index = TPPair.first;
    int probe_index = TPPair.second;
    bool os = false;
    if (Muon_charge[tag_index] == Muon_charge[probe_index]) os = true;
    isOS.push_back(os);
  }
  return isOS;
}



float clipGenWeight(float genWeight){
  float sign =  genWeight/ abs(genWeight);
  return sign;
}
