from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
config.General.requestName = 'DYJetsToLL_madgraphMLMRunIIFall18-102X_upgrade2018'
#config.General.requestName = 'DYJetsToLL_madgraphMLMRunIIFall18-102X_upgrade2018FlatPU0to70'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../muonPogNtuples_miniAOD_cfg.py'
config.JobType.pyCfgParams = ['globalTag=102X_upgrade2018_realistic_v10',
                              'ntupleName=muonPOGNtuple_IsolationStudies18_DY.root',
                              'nEvents=-1',
                              'runOnMC=True',
                              'hltPathFilter=all',
                              'minMuPt=5.0',
                              'minNMu=1'
               ]
config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall18MiniAOD-102X_upgrade2018_realistic_v12_ext1-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall18MiniAOD-FlatPU0to70_102X_upgrade2018_realistic_v12-v1/MINIAODSIM'


config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob  = 1
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
# config.Data.ignoreLocality  = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

