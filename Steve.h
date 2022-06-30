#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"

using namespace ROOT;
using namespace ROOT::VecOps;


/*class Probe{
  public:
  double pt;
  double eta;
  double phi;

  int charge;  
  int original_index;
  };*/


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

RVec<Int_t> CreateProbes_Muon(RVec<Float_t> &Muon_pt, RVec<Float_t> &Muon_standalonePt,RVec<Float_t> &Muon_eta,RVec<Float_t> &Muon_phi,RVec<Int_t> &Muon_charge, RVec<Bool_t> &Muon_mediumId, RVec<Float_t> &Muon_dxybs,RVec<Bool_t> &Muon_isGlobal){
  RVec<Int_t> Probe_Muons;
  for(int i=0;i<Muon_pt.size();i++){
    if(Muon_pt.at(i) < 15 || Muon_standalonePt.at(i) < 15 || abs(Muon_eta.at(i)) > 2.4 || !Muon_mediumId.at(i) 
	|| abs(Muon_dxybs.at(i)) > 0.05 || !Muon_isGlobal.at(i)) continue;
    //Probe probe_muon;
    //probe_muon.pt = Muon_pt.at(i);
    //probe_muon.eta = Muon_eta.at(i);
    //probe_muon.phi = Muon_phi.at(i);
    //probe_muon.charge = Muon_charge.at(i);
    //probe_muon.original_index = i;
    Probe_Muons.push_back(i);
  }
  return Probe_Muons;
}

RVec<Int_t> CreateProbes_Track(RVec<Float_t> &Track_pt, RVec<Float_t> &Track_eta,RVec<Float_t> &Track_phi,RVec<Int_t> &Track_charge, RVec<Int_t> &Track_trackOriginalAlgo){

  RVec<Int_t> Probe_Tracks;
  for(int i=0;i<Track_pt.size();i++){
    if(Track_pt.at(i) < 15. || Track_trackOriginalAlgo.at(i) == 13 || Track_trackOriginalAlgo.at(i) == 14) continue;
    //Probe probe_track;
    //probe_track.pt = Track_pt.at(i);
    //probe_track.eta = Track_eta.at(i);
    //probe_track.phi = Track_phi.at(i);
    //probe_track.charge = Track_charge.at(i);
    //probe_track.original_index = i;
    Probe_Tracks.push_back(i);
  }
  return Probe_Tracks;

}

RVec<std::pair<int,int>> CreateTPPair(RVec<Int_t> &Muon_charge, RVec<Int_t> &isTag, RVec<Bool_t> &isTriggeredMuon, RVec<Bool_t> &isGenMatchedMuon, RVec<Int_t> &Probe_Candidates, RVec<Bool_t> &isGenMatchedProbe ){
  RVec<std::pair<int,int>> TP_pairs;
  for(int iLep1=0; iLep1<Muon_charge.size();iLep1++){
    //if(!isInAcceptance[iLep1]) continue;
    if(!isTag[iLep1]) continue;
    if(!isTriggeredMuon[iLep1]) continue;
    if(!isGenMatchedMuon[iLep1]) continue;




    for(int iLep2=0; iLep2<Probe_Candidates.size(); iLep2++){
      int probe_candidate = Probe_Candidates.at(iLep2);
      //if (iLep2==iLep1) continue;
      //if(!isProbe[iLep2]) continue;
      //if(!isTriggeredMuon[iLep2]) continue;
      if(!isGenMatchedProbe[probe_candidate]) continue;
      //#if(Muon_charge[iLep1] == Muon_charge[iLep2]) continue;

      std::pair<int,int> TP_pair = std::make_pair(iLep1,probe_candidate); 
      TP_pairs.push_back(TP_pair);
    }

  }
  return TP_pairs;
}

RVec<Bool_t> hasStandAloneOrGlobalMatch(RVec<Float_t> &Track_eta, RVec<Float_t> &Track_phi,RVec<Float_t> &Muon_eta, RVec<Float_t> &Muon_phi,RVec<Bool_t> &Muon_isStandalone,RVec<Bool_t> &Muon_isGlobal){
  RVec<Int_t> hasStandAloneOrGlobalMatch;
  for(int iTrack=0;iTrack<Track_eta.size();iTrack++){
    bool has_match = false;
    for (int iMuon=0; iMuon<Muon_eta.size(); ++iMuon){
      if (!(Muon_isStandalone[iMuon] || Muon_isGlobal[iMuon])) continue;
      if (deltaR(Muon_eta[iMuon], Muon_phi[iMuon], Track_eta[iTrack], Track_phi[iTrack]) < 0.01) {has_match = 1; break;}
    } 
    hasStandAloneOrGlobalMatch.push_back(has_match);
  }
  return hasStandAloneOrGlobalMatch;
}


RVec<Float_t> trackMuonDR(RVec<Float_t> &Track_eta, RVec<Float_t> &Track_phi,RVec<Float_t> &Muon_eta,RVec<Float_t> &Muon_phi){
  RVec<Float_t> trackMuonDR; 
  for(int iTrack=0;iTrack<Track_eta.size();iTrack++){
    float dr = 999.;
    float tmp_dr  = 999.;

    for (unsigned int iMuon=0; iMuon<Muon_eta.size(); ++iMuon){
      tmp_dr  = deltaR(Muon_eta.at(iMuon), Muon_phi.at(iMuon), Track_eta.at(iTrack), Track_phi.at(iTrack));
      if (tmp_dr < dr) dr = tmp_dr;
    }
    trackMuonDR.push_back(dr);
  }  
  return trackMuonDR;
}


RVec<Float_t> trackStandaloneDR(RVec<Float_t> &Track_eta, RVec<Float_t> &Track_phi,RVec<Float_t> &Muon_standaloneEta,RVec<Float_t> &Muon_standalonePhi){
   RVec<Float_t> trackStandaloneDR;
   for(int iTrack=0;iTrack<Track_eta.size();iTrack++){
     float dr = 999.;
     float tmp_dr  = 999.;
     
     for (unsigned int iMuon=0; iMuon<Muon_standaloneEta.size(); ++iMuon){
       tmp_dr  = deltaR(Muon_standaloneEta.at(iMuon), Muon_standalonePhi.at(iMuon), Track_eta.at(iTrack), Track_phi.at(iTrack));
       if (tmp_dr < dr) dr = tmp_dr;
     } 
     trackStandaloneDR.push_back(dr);
   }
   return trackStandaloneDR;
 }

RVec<Bool_t> hasTriggerMatch(RVec<Float_t> &Muon_eta, RVec<Float_t> &Muon_phi, RVec<Int_t> &TrigObj_id,RVec<Float_t> &TrigObj_pt,RVec<Float_t> & TrigObj_l1pt, RVec<Float_t> & TrigObj_l2pt, RVec<Int_t> &TrigObj_filterBits, RVec<Float_t> &TrigObj_eta, RVec<Float_t> &TrigObj_phi){
  RVec<Bool_t> TriggerMatch;
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

RVec<Bool_t> hasGenMatch(RVec<Int_t> &GenPart_pdgId, RVec<Int_t> &GenPart_status,RVec<Int_t> &GenPart_statusFlags,RVec<Float_t> &GenPart_eta,RVec<Float_t> &GenPart_phi,RVec<Float_t> &Cand_eta, RVec<Float_t> &Cand_phi){
  RVec<Bool_t> isGenMatched;  
  for(int iCand=0;iCand<Cand_eta.size();iCand++){
    bool matchMC = false;   
    float mcmatch_tmp_dr = 999.;
    for(unsigned int iGen=0; iGen<GenPart_pdgId.size(); iGen++){
      if ( abs(GenPart_pdgId[iGen])==13 &&  GenPart_status[iGen] == 1 && (GenPart_statusFlags[iGen] & 1)  ) {
	if (deltaR(Cand_eta[iCand], Cand_phi[iCand], GenPart_eta[iGen], GenPart_phi[iGen]) < mcmatch_tmp_dr){
	  mcmatch_tmp_dr  = deltaR(Cand_eta[iCand], Cand_phi[iCand], GenPart_eta[iGen], GenPart_phi[iGen]);
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



RVec<Float_t> getTPmass(RVec<std::pair<int,int>> TPPairs, RVec<Float_t> &Muon_pt,RVec<Float_t> &Muon_eta,RVec<Float_t> &Muon_phi, RVec<Float_t> &Cand_pt, RVec<Float_t> &Cand_eta, RVec<Float_t> &Cand_phi){
  RVec<Float_t> TPMass;
  for (int i=0;i<TPPairs.size();i++){
    std::pair<int,int> TPPair = TPPairs.at(i);
    int tag_index = TPPair.first;
    int probe_index = TPPair.second;
    TLorentzVector tagLep(0,0,0,0);
    tagLep.SetPtEtaPhiM(Muon_pt[tag_index],Muon_eta[tag_index],Muon_phi[tag_index],0);
    TLorentzVector probeLep(0,0,0,0);
    probeLep.SetPtEtaPhiM(Cand_pt[probe_index],Cand_eta[probe_index],Cand_phi[probe_index],0);
    TPMass.push_back((tagLep+probeLep).M());        
  }
  return TPMass;
}


RVec<Float_t> getVariables(RVec<std::pair<int,int>> TPPairs, RVec<Float_t> &Cand_variable, float option /*1 for tag and 2 for probe*/){
  RVec<Float_t> Variables;
  for (int i=0;i<TPPairs.size();i++){
    std::pair<int,int> TPPair = TPPairs.at(i);
    float variable; 
    if (option==1) variable = Cand_variable.at(TPPair.first);
    else if (option==2) variable = Cand_variable.at(TPPair.second);
    Variables.push_back(variable);
  }
  return Variables;
}

RVec<Int_t> getVariables(RVec<std::pair<int,int>> TPPairs, RVec<Int_t> &Cand_variable,float option /*1 for tag and 2 for probe*/){
  RVec<Int_t> Variables;
  for (int i=0;i<TPPairs.size();i++){
    std::pair<int,int> TPPair = TPPairs.at(i);
    int variable;
    if(option==1) variable = Cand_variable.at(TPPair.first);
    else if (option==2) variable = Cand_variable.at(TPPair.second);
    Variables.push_back(variable);
  }
  return Variables;
}

RVec<Bool_t> getVariables(RVec<std::pair<int,int>> TPPairs, RVec<Bool_t> &Cand_variable,float option /*1 for tag and 2 for probe*/){
  RVec<Bool_t> Variables;
  for (int i=0;i<TPPairs.size();i++){
    std::pair<int,int> TPPair = TPPairs.at(i);
    bool variable;
    if(option==1) variable = Cand_variable.at(TPPair.first);
    else if (option==2) variable = Cand_variable.at(TPPair.second);
    Variables.push_back(variable);
  }
  return Variables;
}

RVec<Bool_t> isOS(RVec<std::pair<int,int>> TPPairs,RVec<Int_t> Muon_charge, RVec<Int_t> Cand_charge){
  RVec<Bool_t> isOS;
  for (int i=0;i<TPPairs.size();i++){
    std::pair<int,int> TPPair = TPPairs.at(i);
    int  tag_index = TPPair.first;
    int probe_index = TPPair.second;
    bool os = false;
    if (Muon_charge[tag_index] == Cand_charge[probe_index]) os = true;
    isOS.push_back(os);
  }
  return isOS;
}



float clipGenWeight(float genWeight){
  float sign =  genWeight/ abs(genWeight);
  return sign;
}
