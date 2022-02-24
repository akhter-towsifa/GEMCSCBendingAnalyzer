# GEMCSCBendingAnalyzer

## how to check out cmssw and this package
cmsrel CMSSW_12_0_3 

cd CMSSW_12_0_3/src/

cmsenv

git cms-init

git clone https://github.com/aebid/GEMCSCBendingAnalyzer.git

scram b -j 8

## analyzer.cc: package for GEM residual 
### run on data
cd test
cmsRun analyzer.py

## GEM_fitter.cpp: package for creating alignment estimates
### run on output form analyzer
cd script/standAloneGemAlignment
./run_GEM_fitter.sh

Will create a .csv with the recommended alignment values on each GEM chamber, options to include 3dof (dx, dy, dphiz)



The packaged is inherited from Jason Lee's MuonPerformance and it is used for GEM related analysis at TAMU group. The major target is to analyse
GEM-CSC bending angle in real data from CMS. 
