import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("HPSTRACKS")

process.load("FWCore.MessageService.MessageLogger_cfi")
# Load part of the CMSSW reconstruction sequence to make vertexing possible.
# We'll need the CMS geometry and magnetic field to follow the true, non-helical shapes of tracks through the detector.
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
# process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
# process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

studyroot = {
	'SUSYMC': 
		{
			'isData':False, 
			'fileName': 'file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_26_patch1/src/AOD_pi0_study/0E2AE912-1C0E-E611-86FB-002590743042.root',
			'output_rootfile_name': "out_simAOD_SUSYMC.root"
		},
	'JetHTdata': # JetHT AOD dataset=/JetHT*/*PromptReco-v2/AOD*  => dataset=/JetHT/Run2016B-PromptReco-v2/AOD
		{
			'isData': True, 
			'fileName': 'file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_20/src/Kappa/Skimming/higgsTauTau/16DA718F-DA19-E611-BCEE-02163E01376E.root',
			'output_rootfile_name': "out_AOD_JetHTdata.root"
		},
	'kappaminidata':
		{
			'isData': True, 
			'fileName': 'file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_20/src/Kappa/Skimming/higgsTauTau/miniAOD-prod_PAT.root',
			'output_rootfile_name': "out_miniAOD_kappaminidata.root"
		},
	'JetHTdata2': #ONLY ON NFS # JetHT AOD dataset=/JetHT*/*PromptReco-v2/AOD*  => dataset=/JetHT/Run2016B-PromptReco-v2/AOD
		{
			'isData': True, 
			'fileName': 'file:/nfs/dust/cms/user/glusheno/CMSSW_8_0_20/src/Kappa/Skimming/higgsTauTau/16DA718F-DA19-E611-BCEE-02163E01376E.root',
			'output_rootfile_name': "out_AOD_JetHTdata2.root"
		}
}
filekey  = 'JetHTdata'
isData = studyroot[filekey]['isData']
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(studyroot[filekey]['fileName']))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))# -1 for all events

# Suppress messages that are less important than ERRORs.
#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring("cout"),
    cout = cms.untracked.PSet(threshold = cms.untracked.string("INFO")))

print("HPSTracksLable")
process.HPSTracksLable = cms.EDProducer('HPStracksProducer',
    TauPiZeroCollectionTag = cms.InputTag("hpsPFTauProducer","pizeros","RECO"))

print("output prep")
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

print("Run")
process.p = cms.Path(process.HPSTracksLable)

process.e = cms.EndPath(process.out)
