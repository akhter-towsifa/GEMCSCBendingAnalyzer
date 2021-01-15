# GEMCSCBendingAnalyzer

## how to check out cmssw and this package
cmsrel CMSSW_11_1_0  

cd CMSSW_11_1_0/src/

cmsenv

git cms-init

git clone https://github.com/aebid/GEMCSCBendingAnalyzer.git

scram b -j 9


## analyser.cc: package for GEM residual 
### run on data
cmsRun analyser.py

The packaged is inherited from Jason Lee's MuonPerformance and it is used for GEM related analysis at TAMU group. The major target is to analyse
GEM-CSC bending angle in real data from CMS. 
