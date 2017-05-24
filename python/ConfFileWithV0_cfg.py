import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

pName = "HPSTRACKSWITHv0"
process = cms.Process(pName)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

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
	'JetHTdataWithHPSTracks': # JetHT AOD dataset=/JetHT*/*PromptReco-v2/AOD*  => dataset=/JetHT/Run2016B-PromptReco-v2/AOD
		{
			'isData': True, 
			'fileName': 'file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_26_patch1/src/HPStracks/HPStracksProducer/myOutputFile.root',
			'output_rootfile_name': "myOutputFileAfterV0.root"
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
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))# -1 for all events
process.MessageLogger = cms.Service("MessageLogger", destinations = cms.untracked.vstring("cout"), cout = cms.untracked.PSet(threshold = cms.untracked.string("INFO")))

print("HPSTracksLable")
process.HPSTracksLable = cms.EDProducer('HPStracksProducer')##process.load("HPStracks.HPStracksProducer.HPSTracks_cfi") # gives hpsTracks
print("HPSTracksLable is now hpsTracks")
process.load("HPStracks.HPStracksProducer.HPSTracks_cfi") # gives hpsTracks

print("This is HPS Producer config file. And this is the source:")
print(process.source)

print("SecondaryVerticesFromNewV0")
process.load("RecoVertex.V0Producer.generalV0Candidates_cfi")
process.SecondaryVerticesFromNewV0 = process.generalV0Candidates.clone( 
	beamSpot = cms.InputTag('offlineBeamSpot'),
	useVertex = cms.bool(True), # By def False
	vertices = cms.InputTag('offlinePrimaryVertices'),
	# which TrackCollection to use for vertexing, example:  "generalV0Candidates","Lambda","RECO" || vector<reco::VertexCompositeCandidate>    "generalV0Candidates"       "Lambda"          "RECO"  || #The standard track collection (label "generalTracks") is not saved in the MiniAOD event content.
	trackRecoAlgorithm = cms.InputTag("hpsTracks", "HPSTracks", pName),# "HPSTracksLable", "HPSTracks", "HPSTRACKS"
	doKShorts = cms.bool(True),
	doLambdas = cms.bool(False),
	vtxDecaySigXYCut = cms.double(10),
	vtxDecaySigXYZCut = cms.double(-1.),
	innerHitPosCut = cms.double(-1.)
)

available_v0 = {'new': {
					"process_link": process.SecondaryVerticesFromNewV0,
					"collectionName": "SecondaryVerticesFromNewV0",
					"lable": "Kshort",
					"Process": "KSHORTS",
					"newv0": True},
				'old': {
					"process_link": 1,
					"collectionName": "generalV0Candidates",
					"lable": "Kshort",
					"Process": "RECO",
					"newv0": False}
					}
which_v0 = available_v0['new']

print("output prep")
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFileWITHv0Ks.root')
)

print("Run")
#process.p = cms.Path( process.SecondaryVerticesFromNewV0 if filekey == "JetHTdataWithHPSTracks" else process.HPSTracksLable * process.SecondaryVerticesFromNewV0)
process.p = cms.Path( process.hpsTracks * process.SecondaryVerticesFromNewV0)# OR: process.p = cms.Path( process.HPSTracksLable * process.SecondaryVerticesFromNewV0)

process.e = cms.EndPath(process.out)
