import ROOT




#ROOT.gInterpreter.ProcessLine(".O3")
ROOT.ROOT.EnableImplicitMT()
#ROOT.gInterpreter.Declare('#include "somestuff.h"')


CreateTPPair_code = '''
using namespace ROOT::VecOps;
ROOT::VecOps::RVec<std::pair<int,int>> CreateTPPair(ROOT::VecOps::RVec<Int_t> &Muon_charge,ROOT::VecOps::RVec<Int_t> &isInAcceptance, ROOT::VecOps::RVec<Int_t> &isTag ){
    ROOT::VecOps::RVec<std::pair<int,int>> TP_pairs;
    for(int iLep1=0; iLep1<Muon_charge.size();iLep1++){
        if(!isInAcceptance[iLep1]) continue;
        if(!isTag[iLep1]) continue;

        
        

        for(int iLep2=0; iLep2<Muon_charge.size(); iLep2++){
            if (iLep2==iLep1) continue;
            if(!isInAcceptance[iLep2]) continue;
            if(Muon_charge[iLep1] == Muon_charge[iLep2]) continue;

            std::pair<int,int> TP_pair = std::make_pair(iLep1,iLep2); 
            TP_pairs.push_back(TP_pair);
        }

    }
    return TP_pairs;
}

'''
ROOT.gInterpreter.Declare(CreateTPPair_code)

get_TP_variables = '''

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


ROOT::VecOps::RVec<Float_t> getTagVariables(ROOT::VecOps::RVec<std::pair<int,int>> TPPairs, ROOT::VecOps::RVec<Float_t> &Muon_variable){
    ROOT::VecOps::RVec<Float_t> TagVariables;
    for (int i=0;i<TPPairs.size();i++){
        std::pair<int,int> TPPair = TPPairs.at(i);
        int tag_index = TPPair.first;
        float tag_variable = Muon_variable[tag_index];
        TagVariables.push_back(tag_variable);
    }
    return TagVariables;
}

ROOT::VecOps::RVec<Int_t> getTagVariables(ROOT::VecOps::RVec<std::pair<int,int>> TPPairs, ROOT::VecOps::RVec<Int_t> &Muon_variable){
    ROOT::VecOps::RVec<Int_t> TagVariables;
    for (int i=0;i<TPPairs.size();i++){
        std::pair<int,int> TPPair = TPPairs.at(i);
        int tag_index = TPPair.first;
        float tag_variable = Muon_variable[tag_index];
        TagVariables.push_back(tag_variable);
    }
    return TagVariables;
}

ROOT::VecOps::RVec<Float_t> getProbeVariables(ROOT::VecOps::RVec<std::pair<int,int>> TPPairs, ROOT::VecOps::RVec<Float_t> &Muon_variable){
    ROOT::VecOps::RVec<Float_t> ProbeVariables;
    for (int i=0;i<TPPairs.size();i++){
        std::pair<int,int> TPPair = TPPairs.at(i);
        int probe_index = TPPair.second;
        float probe_variable = Muon_variable[probe_index];
        ProbeVariables.push_back(probe_variable);
    }
    return ProbeVariables;
}


ROOT::VecOps::RVec<Int_t> getProbeVariables(ROOT::VecOps::RVec<std::pair<int,int>> TPPairs, ROOT::VecOps::RVec<Int_t> &Muon_variable){
    ROOT::VecOps::RVec<Int_t> ProbeVariables;
    for (int i=0;i<TPPairs.size();i++){
        std::pair<int,int> TPPair = TPPairs.at(i);
        int probe_index = TPPair.second;
        float probe_variable = Muon_variable[probe_index];
        ProbeVariables.push_back(probe_variable);
    }
    return ProbeVariables;
}
'''

ROOT.gInterpreter.Declare(get_TP_variables)


d = ROOT.RDataFrame("Events","/scratchnvme/rbhattac/TNP_NanoAOD/SingleMuon/test_457.root" )

d = d.Define("isInAcceptance","Muon_pt > 15 && Muon_standalonePt > 15 && abs(Muon_eta) < 2.4")

d = d.Define("isTag","Muon_pt > 25 && abs(Muon_eta) < 2.4 && Muon_pfRelIso04_all > 0.15 && abs(Muon_dxybs) > 0.05 && Muon_mediumId && Muon_isGlobal")

d = d.Define("weight","1") #for now (testing)


d = d.Define("TPPairs","CreateTPPair(Muon_charge,isInAcceptance,isTag)")

d = d.Define("TPmass","getTPmass(TPPairs,Muon_pt,Muon_eta,Muon_phi)")

d = d.Define("Probe_pt","getProbeVariables(TPPairs,Muon_pt)")

d = d.Define("Probe_eta","getProbeVariables(TPPairs,Muon_eta)")

d = d.Define("Probe_phi","getProbeVariables(TPPairs,Muon_phi)")

d = d.Define("Probe_isolation","getProbeVariables(TPPairs,Muon_pfRelIso04_all)")

pass_histogram = d.Define("Probe_pt_pass","Probe_pt[Probe_isolation<0.15 && (TPmass > 40 && TPmass<140)]").Define("Probe_eta_pass","Probe_eta[Probe_isolation<0.15 && (TPmass > 40 && TPmass<140)]").Define("TPmass_pass","TPmass[Probe_isolation<0.15 && (TPmass > 40 && TPmass<140)]").Histo3D(("Isolation_pass", "Isolation_pass", 20, 15., 35., 50, -2.5, 2.5, 100, 40., 140.),"Probe_pt_pass","Probe_eta_pass","TPmass_pass")
fail_histogram = d.Define("Probe_pt_fail","Probe_pt[Probe_isolation>0.15 && (TPmass > 40 && TPmass<140)]").Define("Probe_eta_fail","Probe_eta[Probe_isolation>0.15 && (TPmass > 40 && TPmass<140)]").Define("TPmass_fail","TPmass[Probe_isolation>0.15 && (TPmass > 40 && TPmass<140)]").Histo3D(("Isolation_fail", "Isolation_fail", 20, 15., 35., 50, -2.5, 2.5, 100, 40., 140.),"Probe_pt_fail","Probe_eta_fail","TPmass_fail")

f_out = ROOT.TFile("test.root","RECREATE")
pass_histogram.Write()
fail_histogram.Write()
f_out.Close()
