# Steve_Marc_Raj
Tool to produce the histograms, which will fed to Marc's fitting tool 


This is a simple rdf code.

Main code - Steve.py

Steve.h contains the functions which are needed to create the TP pairs

To run - use the singularity shell of Josh - 
```
singularity run /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling\:latest
```

Enable access to the EOS - 
```
kinit USERNAME@CERN.CH
```
for accessing the eos files

then "python Steve.py -h". E.g.,

```
python Steve.py -i inputs/dy_mumu_LowPU.txt -o test.root --noVertexPileupWeight -e 1 --lowPU
python Steve.py -i inputs/data_mumu_LowPU.txt -o test_data.root --noVertexPileupWeight -e 1 --lowPU -d 1
```
