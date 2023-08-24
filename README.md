# GEMCSCBendingAnalyzer

## how to check out cmssw and this package
```
cmsrel CMSSW_12_4_6 
```
or if you want to name your work environment differently: 

```scram p -n Your_Choice_of_Name CMSSW CMSSW_12_4_6```

```
cd CMSSW_12_4_6/src/
cmsenv
git cms-init
git clone https://github.com/akhter-towsifa/GEMCSCBendingAnalyzer.git -b 12_4_0
scram b -j 8
```

## Residual Packages
```
cd GEMCSCBendingAnalyzer/GEM_Alignment/test/
```
### GE1/1 Analyzer
```
cmsRun run_GE11ana.py
```
### ME1/1 Analyzer
```
cmsRun run_ME11ana.py
```
### GE1/1 and ME1/1 Analyzer
```
cmsRun run_both_analyzers.py
```
## GEM DB Maker
```
cd GEMCSCBendingAnalyzer/GEM_Alignment/test/
cmsRun GEMAlDBWriter_cfg.py
```
## GEM_fitter.cpp: package for creating alignment estimates
```
cd GEMCSCBendingAnalyzer/GEM_Alignment/script/standAloneGemAlignment
./run_3DOF_Fitter.sh
```
![alignment_cfg_flowchart](https://github.com/akhter-towsifa/GEMCSCBendingAnalyzer/assets/51368122/e4141fa9-64f3-4be7-b0d0-6a545466a1dd)




The package is inherited from Jason Lee's MuonPerformance and it is used for GEM related analysis at TAMU group. The major target is to analyze GEM-CSC bending angle in real data from CMS.

The previous working version (extensively tested with cosmic data) can be found at https://github.com/aebid/GEMCSCBendingAnalyzer

