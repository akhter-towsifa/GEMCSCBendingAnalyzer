# GEMCSCBendingAnalyzer

## how to check out cmssw and this package
cmsrel CMSSW_12_4_6

cd CMSSW_12_4_6/src/

cmsenv

git cms-init

git clone https://github.com/akhter-towsifa/GEMCSCBendingAnalyzer.git -b 12_4_6_alcaReco

scram b -j 8

## Residual Packages
cd GEMCSCBendingAnalyzer/GEM_Alignment/test/
### GE1/1 Analyzer
cmsRun run_GE11ana.py
### ME1/1 Analyzer
cmsRun run_ME11ana.py
### GE1/1 and ME1/1 Analyzer
cmsRun run_both_analyzers.py

## GEM DB Maker
cd GEMCSCBendingAnalyzer/GEM_Alignment/test/

cmsRun GEMAlDBWriter_cfg.py

## GEM_fitter.cpp: package for creating alignment estimates

cd GEMCSCBendingAnalyzer/GEM_Alignment/script/standAloneGemAlignment

./run_3DOF_Fitter.sh







The packaged is inherited from Jason Lee's MuonPerformance and it is used for GEM related analysis at TAMU group. The major target is to analyse
GEM-CSC bending angle in real data from CMS.
