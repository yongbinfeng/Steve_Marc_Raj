#!/usr/bin/env python3

# examples
#
# all working points with default options
# python runAll.py -i /scratch/shared/NanoAOD/Tnp_NanoV9/TNP/ -o testAll
#
# only mc, and only steps 1, 4, 6
# python runAll.py -i /scratch/shared/NanoAOD/Tnp_NanoV9/TNP/ -o testAll -r mc -s 1 4 6
#
# use -m to merge all output files into a single one
# use -d to test the command, without running them automatically


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

workingPoints = { 1: "reco",
                  2: "tracking",
                  3: "idip",
                  4: "trigger",
                  5: "iso",
                  6: "isonotrig"
}

    
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
    parser.add_argument('-m',  '--merge', action='store_true',
                        help='Merge root files in a new one')
    parser.add_argument('-nw', '--noVertexPileupWeight', action='store_true',
                        help='Do not use weights for vertex z position')
    parser.add_argument('-nos', '--noOppositeCharge', action='store_true',
                        help='Do not require opposite sign charge for tag-probe pairs (note that tracking never uses it by default)')
    parser.add_argument('-s','--steps', default=None, nargs='*', type=int, choices=list(workingPoints.keys()),
                        help='Default runs all working points, but can choose only some if needed')
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
        
    toRun = []
    if args.run in ["all", "data"]:
        toRun.append("data")
    if args.run in ["all", "mc"]:
        toRun.append("mc")

    outfiles = [] # store names of output files so to merge them if needed

    postfix = "vertexWeights{v}_oscharge{c}".format(v="0" if args.noVertexPileupWeight  else "1",
                                                    c="0" if args.noOppositeCharge      else "1")

    for xrun in toRun:

        isdata = 0 if xrun == "mc" else 1
        inpath = indir + (inputdir_data if isdata else inputdir_mc)
        for wp in workingPoints.keys():
            if args.steps and wp not in args.steps:
                continue
            charges = [-1, 1] if workingPoints[wp] == "trigger" else [0]
            for ch in charges:
                step = workingPoints[wp]
                if ch:
                    step += "plus" if ch == 1 else "minus"
                outfile = f"{outdir}tnp_{step}_{xrun}_{postfix}.root"
                outfiles.append(outfile)
                cmd = f"python Steve.py -i {inpath} -o {outfile} -d {isdata} -e {wp} -c {ch}"
                if not args.noVertexPileupWeight:
                    cmd += " -vpw"
                if args.noOppositeCharge:
                    cmd += " -nos"
                print("")
                print(f"Running for {xrun} and {step} efficiency")
                safeSystem(cmd, dryRun=args.dryRun)
                print("")

                
    if args.merge:
        mergedFile = f"{outdir}tnp_all_{postfix}.root"
        sourcefiles = " ".join(outfiles)
        haddcmd = f"hadd -f {mergedFile} {sourcefiles}"
        print("")
        print(f"Merging root files with hadd")
        safeSystem(haddcmd, dryRun=args.dryRun)
        print("")
