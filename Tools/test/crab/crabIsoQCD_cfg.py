from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
config.General.requestName = 'QCD_RunIISpring18MiniAOD-100X_upgrade2018_realistic_v10-v1_OnlyTree'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../muonPogNtuples_miniAOD_cfg.py'
config.JobType.pyCfgParams = ['globalTag=100X_upgrade2018_realistic_v10',
                              'ntupleName=muonPOGNtuple_IsolationStudies18_DY.root',
                              'nEvents=-1',
                              'runOnMC=True',
                              'hltPathFilter=all',
                              'minMuPt=5.0',
                              'minNMu=1'
               ]
config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
config.Data.inputDataset = '/QCD_Pt-20toInf_MuEnrichedPt15_TuneCP5_13TeV_pythia8/RunIISpring18MiniAOD-100X_upgrade2018_realistic_v10-v1/MINIAODSIM'

config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob  = 1
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
# config.Data.ignoreLocality  = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.blacklist = ['T2_US_*']
#config.Site.whitelist = ['T2_IT_Bari']
