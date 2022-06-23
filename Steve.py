import ROOT
from array import array

ROOT.gInterpreter.ProcessLine(".O3")
ROOT.ROOT.EnableImplicitMT()
ROOT.gInterpreter.Declare('#include "Steve.h"')
from os import listdir
from os.path import isfile,join

import time

tstart = time.time()
cpustrat = time.process_time()

files=[]

dirnames=["/scratchnvme/rbhattac/TNP_NanoAOD/DY_postVFP_NanoV8MC_TagAndProbe/220405_123536/0000","/scratchnvme/rbhattac/TNP_NanoAOD/DY_postVFP_NanoV8MC_TagAndProbe/220405_123536/0001","/scratchnvme/rbhattac/TNP_NanoAOD/DY_postVFP_NanoV8MC_TagAndProbe/220405_123536/0002","/scratchnvme/rbhattac/TNP_NanoAOD/DY_postVFP_NanoV8MC_TagAndProbe/220405_123536/0003"]

for dirname in dirnames:
    for f in listdir(dirname):
        if f.endswith(".root"):
            files.append(join(dirname, f))


print(files)
filenames = ROOT.std.vector('string')()

for name in files: filenames.push_back(name)

d = ROOT.RDataFrame("Events",filenames )

d = d.Filter("HLT_IsoMu24 || HLT_IsoTkMu24")

d = d.Filter("PV_npvsGood >= 1")

d = d.Define("isProbe","Muon_pt > 15 && Muon_standalonePt > 15 && abs(Muon_eta) < 2.4")

d = d.Define("isTriggeredMuon","hasTriggerMatch(Muon_eta,Muon_phi,TrigObj_id,TrigObj_pt,TrigObj_l1pt,TrigObj_l2pt,TrigObj_filterBits,TrigObj_eta,TrigObj_phi)")

d = d.Define("isGenMatchedMuon","hasGenMatch(GenPart_pdgId,GenPart_status,GenPart_statusFlags,GenPart_eta,GenPart_phi,Muon_eta,Muon_phi)")

d = d.Define("isTag","Muon_pt > 25 && abs(Muon_eta) < 2.4 && Muon_pfRelIso04_all < 0.15 && abs(Muon_dxybs) < 0.05 && Muon_mediumId && Muon_isGlobal")

d = d.Define("weight","clipGenWeight(genWeight)") #for now (testing)


d = d.Define("All_TPPairs","CreateTPPair(Muon_charge,isProbe,isTag,isTriggeredMuon,isGenMatchedMuon)")

d = d.Define("All_TPmass","getTPmass(All_TPPairs,Muon_pt,Muon_eta,Muon_phi)")

d = d.Define("TPPairs","All_TPPairs[All_TPmass >40 && All_TPmass < 140]").Define("TPmass","All_TPmass[All_TPmass > 40 && All_TPmass < 140]")
d = d.Define("Is_Pair_OS","isOS(TPPairs,Muon_charge)")

d = d.Define("Probe_pt","getVariables(TPPairs,Muon_pt,2)")

d = d.Define("Probe_eta","getVariables(TPPairs,Muon_eta,2)")

d = d.Define("Probe_phi","getVariables(TPPairs,Muon_phi,2)")

d = d.Define("Probe_isolation","getVariables(TPPairs,Muon_pfRelIso04_all,2)")

d = d.Define("Probe_isGlobal","getVariables(TPPairs,Muon_isGlobal,2)")

d = d.Define("Probe_isStandalone","getVariables(TPPairs,Muon_isStandalone,2)")

d = d.Define("Probe_mediumId","getVariables(TPPairs,Muon_mediumId,2)")

d = d.Define("Probe_dxbs","getVariables(TPPairs,Muon_dxybs,2)")

##For Isolation

iso_df = d.Redefine("TPmass","TPmass[Probe_mediumId && abs(Probe_dxbs) < 0.05]").Redefine("Probe_pt","Probe_pt[Probe_mediumId && abs(Probe_dxbs) < 0.05]").Redefine("Probe_eta","Probe_eta[Probe_mediumId && abs(Probe_dxbs) < 0.05]").Redefine("Probe_isolation","Probe_isolation[Probe_mediumId && abs(Probe_dxbs) < 0.05]").Redefine("Probe_isGlobal","Probe_isGlobal[Probe_mediumId && abs(Probe_dxbs) < 0.05]")
iso_df = iso_df.Redefine("TPmass","TPmass[Probe_isGlobal]").Redefine("Probe_pt","Probe_pt[Probe_isGlobal]").Redefine("Probe_eta","Probe_eta[Probe_isGlobal]").Redefine("Probe_isolation","Probe_isolation[Probe_isGlobal]")

binning_pt = array('d',[24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])
binning_eta = array('d',[-2.4 + i*0.1 for i in range(49)])
binning_mass = array('d',[50 + i for i in range(81)])

model_pass = ROOT.RDF.TH3DModel("Isolation_pass", "Isolation_pass", 15, binning_pt, 48, binning_eta, 80, binning_mass)
model_fail = ROOT.RDF.TH3DModel("Isolation_fail", "Isolation_fail", 15, binning_pt, 48, binning_eta, 80, binning_mass)

pass_histogram = iso_df.Define("Probe_pt_pass","Probe_pt[Probe_isolation<0.15]").Define("Probe_eta_pass","Probe_eta[Probe_isolation<0.15]").Define("TPmass_pass","TPmass[Probe_isolation<0.15]").Histo3D(model_pass,"Probe_pt_pass","Probe_eta_pass","TPmass_pass","weight")
fail_histogram = iso_df.Define("Probe_pt_fail","Probe_pt[Probe_isolation>0.15]").Define("Probe_eta_fail","Probe_eta[Probe_isolation>0.15]").Define("TPmass_fail","TPmass[Probe_isolation>0.15]").Histo3D(model_fail,"Probe_pt_fail","Probe_eta_fail","TPmass_fail","weight")

f_out = ROOT.TFile("Steve.root","RECREATE")
pass_histogram.Write()
fail_histogram.Write()
f_out.Close()

elapsed = time.time() - tstart
elapsed_cpu = time.process_time() - cpustrat
print('Execution time:', elapsed, 'seconds')
print('CPU Execution time:', elapsed_cpu , 'seconds')
