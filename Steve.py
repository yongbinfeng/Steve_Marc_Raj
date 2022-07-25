import ROOT
from array import array
import json

ROOT.gInterpreter.ProcessLine(".O3")
ROOT.ROOT.EnableImplicitMT()
ROOT.gInterpreter.Declare('#include "Steve.h"')
from os import listdir
from os.path import isfile,join
from os import walk

import time

import argparse


def make_jsonhelper(filename):
    with open(filename) as jsonfile:
        jsondata = json.load(jsonfile)
    
    runs = []
    firstlumis = []
    lastlumis = []
    
    for run,lumipairs in jsondata.items():
        for lumipair in lumipairs:
            runs.append(int(run))
            firstlumis.append(int(lumipair[0]))
            lastlumis.append(int(lumipair[1]))
    
    jsonhelper = ROOT.JsonHelper(runs, firstlumis, lastlumis)
    
    return jsonhelper

parser = argparse.ArgumentParser()

parser.add_argument("-e","--efficiency", help="1 for reco,2 for trigger, 3 for isolation",
                    type=int)
parser.add_argument("-i","--input_path", help="path of the input root files",
                    type=str)

parser.add_argument("-o","--output_file", help="name of the output root file",
                    type=str)

parser.add_argument("-d","--isData", help="Pass 0 for MC, 1 for Data, default is 0",
                    type=int, default=0)


args = parser.parse_args()
tstart = time.time()
cpustrat = time.process_time()


if args.isData == 1:
    histo_name= "RunGtoH"
else:
    histo_name = "DY_postVFP"

files=[]

#dirnames=["/scratchnvme/rbhattac/TNP_NanoAOD/DY_postVFP_NanoV8MC_TagAndProbe/220405_123536/0000","/scratchnvme/rbhattac/TNP_NanoAOD/DY_postVFP_NanoV8MC_TagAndProbe/220405_123536/0001","/scratchnvme/rbhattac/TNP_NanoAOD/DY_postVFP_NanoV8MC_TagAndProbe/220405_123536/0002","/scratchnvme/rbhattac/TNP_NanoAOD/DY_postVFP_NanoV8MC_TagAndProbe/220405_123536/0003"]

#for dirname in dirnames:
#    for f in listdir(dirname):
#        if f.endswith(".root"):
#            files.append(join(dirname, f))

for root, dirnames, filenames in walk(args.input_path):
     for filename in filenames:
          if '.root' in filename:
              files.append(join(root, filename))



#print(files)
filenames = ROOT.std.vector('string')()

for name in files: filenames.push_back(name)

d = ROOT.RDataFrame("Events",filenames )

f_out = ROOT.TFile(args.output_file,"RECREATE")

binning_pt = array('d',[24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])
binning_eta = array('d',[-2.4 + i*0.1 for i in range(49)])
binning_mass = array('d',[50 + i for i in range(81)])

##General Cuts
d = d.Filter("HLT_IsoMu24 || HLT_IsoTkMu24","HLT Cut")

d = d.Filter("PV_npvsGood >= 1","NVtx Cut")

if (args.isData == 1):
    jsonhelper = make_jsonhelper("Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt")
    d = d.Filter(jsonhelper,["run","luminosityBlock"],"jsonfilter")

## Weights

if(args.isData == 1):
    d = d.Define("weight","1")
else:
    d = d.Define("gen_weight","clipGenWeight(genWeight)") #for now (testing)

    d = d.Define("pu_weight","puw_2016(Pileup_nTrueInt,2)")

    d = d.Define("weight","gen_weight*pu_weight")

## For Tag Muons
d = d.Define("isTriggeredMuon","hasTriggerMatch(Muon_eta,Muon_phi,TrigObj_id,TrigObj_pt,TrigObj_l1pt,TrigObj_l2pt,TrigObj_filterBits,TrigObj_eta,TrigObj_phi)")

if(args.isData == 1):
    d = d.Define("isGenMatchedMuon","createTrues(nMuon)")
else: 
    d = d.Define("isGenMatchedMuon","hasGenMatch(GenPart_pdgId,GenPart_status,GenPart_statusFlags,GenPart_eta,GenPart_phi,Muon_eta,Muon_phi)")

d = d.Define("isTag","Muon_pt > 25 && abs(Muon_eta) < 2.4 && Muon_pfRelIso04_all < 0.15 && abs(Muon_dxybs) < 0.05 && Muon_mediumId && Muon_isGlobal")


## Tracks for reco efficiency
if(args.efficiency == 1):
    d = d.Define("trackHasStandAloneorGlobalMatch","hasStandAloneOrGlobalMatch(Track_eta,Track_phi,Muon_eta,Muon_phi,Muon_isStandalone,Muon_isGlobal)")

    if(args.isData == 1):
        d = d.Define("isGenMatchedTrack","createTrues(nTrack)")
    else:
        d = d.Define("isGenMatchedTrack","hasGenMatch(GenPart_pdgId,GenPart_status,GenPart_statusFlags,GenPart_eta,GenPart_phi,Track_eta,Track_phi)")

    d = d.Define("trackMuonDR","trackMuonDR(Track_eta,Track_phi,Muon_eta,Muon_phi)")

    d = d.Define("trackStandaloneDR","trackStandaloneDR(Track_eta,Track_phi,Muon_standaloneEta,Muon_standalonePhi)")
    d = d.Define("Probe_Tracks","CreateProbes_Track(Track_pt,Track_eta,Track_phi,Track_charge,Track_trackOriginalAlgo)")

    d = d.Define("All_TPPairs","CreateTPPair(Muon_charge,isTag,isTriggeredMuon,isGenMatchedMuon,Probe_Tracks,isGenMatchedTrack)")

    d = d.Define("All_TPmass","getTPmass(All_TPPairs,Muon_pt,Muon_eta,Muon_phi,Track_pt,Track_eta,Track_phi)")

    d = d.Define("TPPairs","All_TPPairs[All_TPmass >40 && All_TPmass < 140]").Define("TPmass","All_TPmass[All_TPmass > 40 && All_TPmass < 140]")

    d = d.Define("Is_Pair_OS","isOS(TPPairs,Muon_charge,Track_charge)")

    d = d.Define("Probe_pt","getVariables(TPPairs,Track_pt,2)")

    d = d.Define("Probe_eta","getVariables(TPPairs,Track_eta,2)")

    d = d.Define("Probe_phi","getVariables(TPPairs,Track_phi,2)")

    d = d.Define("Probe_StandaloneDR","getVariables(TPPairs,trackStandaloneDR,2)")

    
    model_pass_reco = ROOT.RDF.TH3DModel("pass_mu_"+histo_name, "Reco_pass",len(binning_mass)-1, binning_mass, len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
    model_fail_reco = ROOT.RDF.TH3DModel("fail_mu_"+histo_name, "Reco_fail",len(binning_mass)-1, binning_mass, len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)

    pass_histogram_reco = d.Define("Probe_pt_pass","Probe_pt[Probe_StandaloneDR<0.3]").Define("Probe_eta_pass","Probe_eta[Probe_StandaloneDR<0.3]").Define("TPmass_pass","TPmass[Probe_StandaloneDR<0.3]").Histo3D(model_pass_reco,"TPmass_pass","Probe_pt_pass","Probe_eta_pass","weight")

    fail_histogram_reco = d.Define("Probe_pt_fail","Probe_pt[Probe_StandaloneDR>0.3]").Define("Probe_eta_fail","Probe_eta[Probe_StandaloneDR>0.3]").Define("TPmass_fail","TPmass[Probe_StandaloneDR>0.3]").Histo3D(model_fail_reco,"TPmass_fail","Probe_pt_fail","Probe_eta_fail","weight")

    pass_histogram_reco.Write()
    fail_histogram_reco.Write()


## Muons for all other efficiency
else:
    d = d.Define("Probe_Muons","CreateProbes_Muon(Muon_pt,Muon_standalonePt,Muon_eta,Muon_phi,Muon_charge,Muon_mediumId,Muon_dxybs,Muon_isGlobal)")


    #d = d.Define("isProbe","Muon_pt > 15 && Muon_standalonePt > 15 && abs(Muon_eta) < 2.4 && Muon_mediumId && abs(Muon_dxybs) < 0.05 && Muon_isGlobal")



    d = d.Define("All_TPPairs","CreateTPPair(Muon_charge,isTag,isTriggeredMuon,isGenMatchedMuon,Probe_Muons,isGenMatchedMuon)")

    d = d.Define("All_TPmass","getTPmass(All_TPPairs,Muon_pt,Muon_eta,Muon_phi,Muon_pt,Muon_eta,Muon_phi)")

    d = d.Define("TPPairs","All_TPPairs[All_TPmass >40 && All_TPmass < 140]").Define("TPmass","All_TPmass[All_TPmass > 40 && All_TPmass < 140]")
    d = d.Define("Is_Pair_OS","isOS(TPPairs,Muon_charge,Muon_charge)")

    d = d.Define("Probe_pt","getVariables(TPPairs,Muon_pt,2)")

    d = d.Define("Probe_eta","getVariables(TPPairs,Muon_eta,2)")

    d = d.Define("Probe_phi","getVariables(TPPairs,Muon_phi,2)")

    d = d.Define("Probe_isTriggered","getVariables(TPPairs,isTriggeredMuon,2)")

    d = d.Define("Probe_isolation","getVariables(TPPairs,Muon_pfRelIso04_all,2)")

    d = d.Define("Probe_isGlobal","getVariables(TPPairs,Muon_isGlobal,2)")

    d = d.Define("Probe_isStandalone","getVariables(TPPairs,Muon_isStandalone,2)")

    d = d.Define("Probe_mediumId","getVariables(TPPairs,Muon_mediumId,2)")

    d = d.Define("Probe_dxybs","getVariables(TPPairs,Muon_dxybs,2)")


    #iso_df = d.Redefine("TPmass","TPmass[Probe_mediumId && abs(Probe_dxybs) < 0.05 && Probe_isTriggered]").Redefine("Probe_pt","Probe_pt[Probe_mediumId && abs(Probe_dxybs) < 0.05 && Probe_isTriggered]").Redefine("Probe_eta","Probe_eta[Probe_mediumId && abs(Probe_dxybs) < 0.05 && Probe_isTriggered]").Redefine("Probe_isolation","Probe_isolation[Probe_mediumId && abs(Probe_dxybs) < 0.05 && Probe_isTriggered]").Redefine("Probe_isGlobal","Probe_isGlobal[Probe_mediumId && abs(Probe_dxybs) < 0.05 && Probe_isTriggered]")
    #iso_df = iso_df.Redefine("TPmass","TPmass[Probe_isGlobal]").Redefine("Probe_pt","Probe_pt[Probe_isGlobal]").Redefine("Probe_eta","Probe_eta[Probe_isGlobal]").Redefine("Probe_isolation","Probe_isolation[Probe_isGlobal]")


    # For Trigger
    if(args.efficiency == 2):
        model_pass_trig = ROOT.RDF.TH3DModel("pass_mu_"+histo_name, "Trigger_pass",len(binning_mass)-1, binning_mass, len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        model_fail_trig = ROOT.RDF.TH3DModel("fail_mu_"+histo_name, "Trigger_fail",len(binning_mass)-1, binning_mass, len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)

        pass_histogram_trig = d.Define("Probe_pt_pass","Probe_pt[Probe_isTriggered]").Define("Probe_eta_pass","Probe_eta[Probe_isTriggered]").Define("TPmass_pass","TPmass[Probe_isTriggered]").Histo3D(model_pass_trig,"TPmass_pass","Probe_pt_pass","Probe_eta_pass","weight")

        fail_histogram_trig = d.Define("Probe_pt_fail","Probe_pt[!Probe_isTriggered]").Define("Probe_eta_fail","Probe_eta[!Probe_isTriggered]").Define("TPmass_fail","TPmass[!Probe_isTriggered]").Histo3D(model_fail_trig,"TPmass_fail","Probe_pt_fail","Probe_eta_fail","weight")

        pass_histogram_trig.Write()
        fail_histogram_trig.Write()

     ##For Isolation

    if(args.efficiency == 3):
         model_pass_iso = ROOT.RDF.TH3DModel("pass_mu_"+histo_name, "Isolation_pass",len(binning_mass)-1, binning_mass, len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
         model_fail_iso = ROOT.RDF.TH3DModel("fail_mu_"+histo_name, "Isolation_fail",len(binning_mass)-1, binning_mass, len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
    
         d = d.Redefine("TPmass","TPmass[Probe_isTriggered]").Redefine("Probe_pt","Probe_pt[Probe_isTriggered]").Redefine("Probe_eta","Probe_eta[Probe_isTriggered]").Redefine("Probe_isolation","Probe_isolation[Probe_isTriggered]")
 
         pass_histogram_iso = d.Define("Probe_pt_pass","Probe_pt[Probe_isolation<0.15]").Define("Probe_eta_pass","Probe_eta[Probe_isolation<0.15]").Define("TPmass_pass","TPmass[Probe_isolation<0.15]").Histo3D(model_pass_iso,"TPmass_pass","Probe_pt_pass","Probe_eta_pass","weight")
         fail_histogram_iso = d.Define("Probe_pt_fail","Probe_pt[Probe_isolation>0.15]").Define("Probe_eta_fail","Probe_eta[Probe_isolation>0.15]").Define("TPmass_fail","TPmass[Probe_isolation>0.15]").Histo3D(model_fail_iso,"TPmass_fail","Probe_pt_fail","Probe_eta_fail","weight")

         pass_histogram_iso.Write()
         fail_histogram_iso.Write()







f_out.Close()

print(d.Report().Print())

elapsed = time.time() - tstart
elapsed_cpu = time.process_time() - cpustrat
print('Execution time:', elapsed, 'seconds')
print('CPU Execution time:', elapsed_cpu , 'seconds')
