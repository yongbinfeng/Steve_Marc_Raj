#!/usr/bin/env python3

import os, re, copy, math, array

import argparse
import sys

def safeSystem(cmd, dryRun=False, quitOnFail=True):
    print(cmd)
    if not dryRun:
        res = os.system(cmd)
        if res:
            print('-'*30)
            print("safeSystem(): error occurred when executing the following command. Aborting")
            print(cmd)
            print('-'*30)
            if quitOnFail:
                quit()
        return res
    else:
        return 0


if __name__ == "__main__":    

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--indir',  default=None, type=str, required=True,
                        help='Input directory with the root files inside (common path for data and MC)')
    parser.add_argument('-o','--outdir', default=None, type=str, required=True,
                        help='Output directory to store all root files')
    parser.add_argument('-d',  '--dryRun', action='store_true',
                        help='Do not execute commands, just print them')
    parser.add_argument('-r',  '--run', default="all", type=str, choices=["data", "mc", "all"],
                        help='Choose what to run, either data or MC, or both')
    args = parser.parse_args()

    outdir = args.outdir
    if not outdir.endswith("/"):
        outdir += "/"
    indir = args.indir
    if not indir.endswith("/"):
        indir += "/"

    if not os.path.exists(outdir):
        print(f"Creating folder {outdir}")
        safeSystem(f"mkdir -p {outdir}", dryRun=False)

    inputdir_data = "SingleMuon/"
    inputdir_mc   = "DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/"
        
    workingPoints = { 1: "reco",
                      2: "tracking",
                      3: "idip",
                      4: "trigger",
                      5: "iso",
                      6: "isonotrig"
    }

    toRun = []
    if args.run in ["all", "data"]:
        toRun.append("data")
    if args.run in ["all", "mc"]:
        toRun.append("mc")

    for xrun in toRun:
        isdata = 0 if xrun == "mc" else 1
        inpath = indir + (inputdir_data if isdata else inputdir_mc)
        for wp in workingPoints.keys():
            charges = [-1, 1] if workingPoints[wp] == "trigger" else [0]
            for ch in charges:
                step = workingPoints[wp]
                if ch:
                    step += "plus" if ch == 1 else "minus"
                outfile = f"{outdir}tnp_{step}_vertexWeights_oscharge.root"
                cmd = f"python Steve.py -i {inpath} -o {outfile} -d {isdata} -e {wp} -c {ch} -vpw"
                print("")
                print(f"Running for {xrun} and {step} efficiency")
                safeSystem(cmd, dryRun=args.dryRun)
                print("")
