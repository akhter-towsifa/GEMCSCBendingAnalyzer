# GEMCSCBendingAnalyzer

## how to check out cmssw and this package
```cmsrel CMSSW_12_4_6```

or if you want to name your work environment differently:

```scram p -n Your_Choice_of_Name CMSSW CMSSW_12_4_6```

```
cd CMSSW_12_4_6/src/
cmsenv
git cms-init
git clone https://github.com/akhter-towsifa/GEMCSCBendingAnalyzer.git -b 12_X_alcaReco
scram b -j 8
```

## Residual Packages
```cd GEMCSCBendingAnalyzer/GEM_Alignment/test/```
### GE1/1 Analyzer
```cmsRun run_GE11ana.py```
### ME1/1 Analyzer
```cmsRun run_ME11ana.py```
### GE1/1 and ME1/1 Analyzer
```cmsRun run_both_analyzers.py```

## GEM DB Maker
![alignment_cfg_flowchart](https://github.com/akhter-towsifa/GEMCSCBendingAnalyzer/assets/51368122/1539cc93-a62a-4508-9370-bdeb98358448)

```
cd GEMCSCBendingAnalyzer/GEM_Alignment/test/
cmsRun GEMAlDBWriter_cfg.py
```

## GEM_fitter.cpp: package for creating alignment estimates

```
cd GEMCSCBendingAnalyzer/GEM_Alignment/script/standAloneGemAlignment
./run_3DOF_Fitter.sh
```






The package is inherited from Jason Lee's MuonPerformance and it is used for GEM related analysis at TAMU group. The major target is to analyze GEM-CSC bending angle in real data from CMS.

The previous working version (extensively tested with cosmic data) can be found at https://github.com/aebid/GEMCSCBendingAnalyzer
