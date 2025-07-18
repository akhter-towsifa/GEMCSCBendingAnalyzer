# Welcome to GEM CSC Bending Analyzer

The code under this repository is used for offline GEM alignment in CMS experiment at CERN

*  Devin's masters defense: [Alignment of the Compact Muon Solenoid’s new Gas Electron Multiplier detector](https://indico.cern.ch/event/1179363/)
*  Towsifa's masters defense: [First Measurement of the Muon pT Dependent Bending Angle between the GEM and CSC Subdetectors using Run 3 Data](https://indico.cern.ch/event/1266980/)

For a quick start, follow the instructions below. But for a better understanding, explore the other pages on this website.

## Set up CMS Software environment and check out this alignment analyzer package

    cmsrel CMSSW_15_0_1 

or if you want to name your work environment differently: 

    scram p -n Your_Choice_of_Name CMSSW_15_0_1

Then

    cd CMSSW_15_0_1/src/
    cmsenv
    git cms-init

    git clone https://github.com/akhter-towsifa/GEMCSCBendingAnalyzer.git -b 15_X_RawReco

    scram b -j 8

* Note: there may be some compilation warnings related to `ME11ana.cc` file. This file is a legacy one (some details below).

## Residual Packages
Located under `GEM_Alignment`. To run the relevant alignment configuration files, access the `test` folder and set up proxy for the server if necessary: 

    cd GEMCSCBendingAnalyzer/GEM_Alignment/test/
    voms-proxy-init --valid 192:00 --voms cms #setting up the proxy

### GE1/1 Analyzer configuration
The associated code `analyzer.cc` is located in the `plugins` folder. To run the GEM alignment code locally,

    cmsRun run_GE11ana.py

### ME1/1 Analyzer (Legacy)
This is a legacy code that was initially developed for standalone tracker-based CSC alignment for station 1 ring 1 chambers, but not updated or used any longer.

    cmsRun run_ME11ana.py

### GE1/1 and ME1/1 Analyzer (Legacy)

    cmsRun run_both_analyzers.py

## GEM DB Maker
The associated code `GEMAlDBWriter.cc` is located in the `plugins` folder.
This code and associated configuration file creates the necessary geometry file in .db format. The input file is in csv format which found by running the alignment estimator (Minuit package)

    cd GEMCSCBendingAnalyzer/GEM_Alignment/test/
    cmsRun GEMAlDBWriter_cfg.py

## GEM_fitter.cpp: package for creating alignment estimates
```
cd GEMCSCBendingAnalyzer/GEM_Alignment/script/standAloneGemAlignment
./run_3DOF_Fitter.sh
```
![alignment_cfg_flowchart](https://github.com/akhter-towsifa/GEMCSCBendingAnalyzer/assets/51368122/e4141fa9-64f3-4be7-b0d0-6a545466a1dd)



## Resources
The package is inherited from Jason Lee's MuonPerformance and it is used for GEM related analysis at TAMU group. The major target is to analyze GEM-CSC bending angle in real data from CMS.

The previous working version (extensively tested with cosmic data) can be found at [https://github.com/aebid/GEMCSCBendingAnalyzer](https://github.com/aebid/GEMCSCBendingAnalyzer)

GEM alignment draft TDR for reference: [DN-24-012](https://gitlab.cern.ch/tdr/notes/DN-24-012)

* #### Devin's masters defense: [Alignment of the Compact Muon Solenoid’s new Gas Electron Multiplier detector](https://indico.cern.ch/event/1179363/)
* #### Towsifa' masters defense: [First Measurement of the Muon pT Dependent Bending Angle between the GEM and CSC Subdetectors using Run 3 Data](https://indico.cern.ch/event/1266980/)