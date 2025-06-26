## Set up CMS Software environment
Usually we try to pick a stable CMSSW version (example below). One other way to decide which CMSSW version to choose is to look under [CMS Data Aggregation System](https://cmsweb.cern.ch/das/). For a desired dataset (given era and datastream type), check the release version of that dataset.

    cmsrel CMSSW_15_0_1 

or if you want to name your work environment differently: 

    scram p -n Your_Choice_of_Name CMSSW_15_0_1

##  Check out this alignment analyzer package
Once the CMSSW version has been selected, set up the environment in the `src` folder and clone the github repository:

    cd CMSSW_15_0_1/src/
    cmsenv
    git cms-init

    git clone https://github.com/akhter-towsifa/GEMCSCBendingAnalyzer.git -b 15_X_RawReco

    scram b -j 8

* Note: there may be some compilation warnings related to `ME11ana.cc` file. This file is a legacy one initially created to check standalone CSC alignment (especially in station 1 ring 1). But this code is no longer utilized, and has some outdated branches that raise compilation warnings.

