# Running the GEM Analyzer

To run the GEM analyzer code, the `run_GE11ana.py` configuration file is utilized. If using a test file from the [CMS DAS](https://cmsweb.cern.ch/das/) website, the proxy will need to be set up prior to running the alignment configuration.

    cd GEMCSCBendingAnalyzer/GEM_Alignment/test/
    voms-proxy-init --valid 192:00 --voms cms #setting up the proxy
    
    cmsRun run_GE11ana.py

The output root file from this local test run can be checked by default in the `test` folder, unless the output path is changed inside the configuration file.

## GEM analyzer configuration file

The configuration file is set up in six different parts to run the EDAnalyzer code.

* CMS Framework setup: For a particular run or phase era (e.g. `Run3` or `Phase2`) the `analyzer` process is created.

* Loading necessary modules: Several modules are loaded (e.g. `MagneticField_AutoFromDBCurrent`)

* Muon System Alignment: depending on what alignment is required for analysis, the configuration is set up in this part. The `misalign` flag directs whether geometry files (sqlite db files) will be given by the user for some specific scenario. For example, the flag `misalign` can be set to `True` and then the user can give a GEM sqlite file that has the records `GEMAlignmentRcd` and `GEMAlignmentErrorExtendedRcd`. If no geometry file is given, it follows the geometry available in global tag.

    + `process.GEMGeometryESModule.applyAlignment` and/or `process.CSCGeometryESModule.applyAlignment` if set to `True`, follows either the sqlite file alignment or the global tag alignment (the flow chart below)
    + if these options are explicitly set to `False` then the geometry for the muon subdetector is assumed to be in ideal configuration.

    ![alignment_cfg_flowchart](https://github.com/akhter-towsifa/GEMCSCBendingAnalyzer/assets/51368122/e4141fa9-64f3-4be7-b0d0-6a545466a1dd)

* Global Tag: specify the global tag to be used. List of some global tags can be found  under [autoCond configuration](https://github.com/cms-sw/cmssw/blob/master/Configuration/AlCa/python/autoCond.py) or can be searched on [CMS Conditions DB Browser](https://cms-conddb.cern.ch/cmsDbBrowser/index/Prod). The autoCond configuration has some description of the listed global tags. For example, the global tag key `run3_data_prompt` points to a prompt Global Tag condition during data-taking period. 

* Event processing, input file, output file: the number of events to be processed is specified. If `maxEvents` is set to `-1`, all events are processed. The input file can be appended as a string of `file:/path/to/file/fileName.root` or `root://node-file-location//path.root`. The output file can be produced in local directory, or set to cernbox location.

* EDAnalyzer: The `analyzer` EDAnalyzer module takes input arguments to be utilized as tokens in the main `analyzer.cc` code. For example, `muons`. The available input module/label can be found by running `edmDumpEventContent` command in lxplus.

        $ edmDumpEventContent root://cms-xrd-global.cern.ch//store/data/path/to/root/file/as/found/on/cms-das/file.root
        Or
        $ edmDumpEventContent root:///eos/path/to/file/file.root

    The available modules will then be displayed if the file is accessible. `ALCA-RECO` files have some more skimmed modules compared to `Raw-Reco` files. In the example of `muons`, for `ALCA-RECO`, it was recommended to use `muons = cms.InputTag("ALCARECOMuAlCalIsolatedMu:SelectedMuons")`.

    + The other arguments are some boolean arguments that need to be set depending on which propagation methods the user is interested in. The `debug` option highly recommended when doing local test runs as it prints out several debugging statements per event. Remember to turn off the `debug` statement when submitting crab jobs as too many print-outs can cause the size of the log files to be too large and thus cause crab job failures.

    + The last part sets the execution path for running the analyzer. Note: if want to run extra steps before `analyzer`, once can add the paths, i.e. `process.p = cms.Path(process.A + process.B + .... + process.analyzer)`. The sequence is to be determined by the user.


## GEM analyzer CRAB configuration file

CRAB is a utility tool to submit CMSSW jobs and distribute these jobs on computing resources. The `crab3_run_GE11ana.py` can be modified and used to submit jobs on the server. One can use condor to submit jobs as well (there are some examples available under `GEM_Alignment/condor` folder which are no longer maintained).

Follow the [Twiki for CRAB Configuration file](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile), to understand different arguments and configure the file to your liking. Below are some standard changes to be made for GEM alignment.

* the General `requestName` should be set so that the jobs can be identified by the user. For example, `Run2025C_muon0_ZMu_150X_dataRun3_Prompt_v1_aligned_trackerprop`.

* the JobType psetName should be set to the GEM configuration file (remember to turn off `debug` under the GEM configuration file). 

* if passing an sqlite file (e.g. GEM geometry db file), set the `misalign` flag to `True` in both GEM configuration file and this CRAB configuration file. Make sure to put the sqlite db file under both configuration files. The input db file should be accessible by CRAB job. The safer option is to copy the sqlite file in the same directory as the CRAB configuration file.

    Note: user/work public directory and cernbox public directories are not accessible to read from by CRAB. CRAB can transfer output files to cernbox but not retrieve it. Docker or SandBox might be necessary to bypass this issue.

* in the `Data` section, the favored `runRange`, `inputDataset`, `lumiMask`, etc. can be set. Refer to [CMS DAS](https://cmsweb.cern.ch/das/) website to find datasets and/or runRange. If user wants to have specific input root files, a list (easily readable .list or .txt file) needs to be passed where the full path of each file is specified by the user (see `userInputFiles` argument).

* `outLFNDirBase` should be in the format `/store/user/<username>[/<subdir>*]`. CRAB will create the `subdir` path if it does not exist already.

* the `storageSite` can be set to CERNBOX or CMSLPC as the user prefers (e.g. `T3_CH_CERNBOX`, `T3_US_FNALLPC`).
