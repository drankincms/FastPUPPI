import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("RESP", eras.Phase2C9_trigger)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load("L1Trigger.TrackFindingTracklet.Tracklet_cfi") 
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/Phase2HLTTDRWinter20DIGI/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3-v2/110000/005E74D6-B50E-674E-89E6-EAA9A617B476.root',
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    skipBadFiles = cms.untracked.bool(True),
    inputCommands = cms.untracked.vstring("keep *", 
        "drop l1tHGCalTowerMapBXVector_hgcalTriggerPrimitiveDigiProducer_towerMap_HLT",
        "drop *_hgcalTriggerPrimitiveDigiProducer_*_*",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT",
        "drop l1tEMTFHit2016s_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016s_simEmtfDigis__HLT")
)

process.load("L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff")

process.TTTracksFromTrackletEmulation.Hnpar = 5
process.pfTracksFromL1Tracks.nParam = 5

process.l1pfProducerRegionsBarrel = process.l1pfProducerBarrel.clone(regions = cms.VPSet(
        cms.PSet(
            etaBoundaries = cms.vdouble(-1.5,-1.0,-0.5,0.,0.5,1.,1.5),
            etaExtra = cms.double(0.25),
            phiExtra = cms.double(0.25),
            phiSlices = cms.uint32(9),
            caloNMax = cms.uint32(15),
            emcaloNMax = cms.uint32(13),
            trackNMax = cms.uint32(22),
            muonNMax = cms.uint32(2),
            pfNMax = cms.uint32(20),
            puppiNMax = cms.uint32(20),
        )
    )
)
process.l1pfProducerRegionsHGCal = process.l1pfProducerHGCal.clone(regions = cms.VPSet(
        cms.PSet(
            etaBoundaries = cms.vdouble(-2.5, -1.5),
            etaExtra = cms.double(0.25),
            phiExtra = cms.double(0.25),
            phiSlices = cms.uint32(9),
            pfNMax = cms.uint32(30),
            puppiNMax = cms.uint32(30),
        ),
        cms.PSet(
            etaBoundaries = cms.vdouble(1.5, 2.5),
            etaExtra = cms.double(0.25),
            phiExtra = cms.double(0.25),
            phiSlices = cms.uint32(9),
            pfNMax = cms.uint32(30),
            puppiNMax = cms.uint32(30),
        )
    )
)

process.l1pfCandidatesRegions = cms.EDProducer("L1TPFCandMultiMerger",
    pfProducers = cms.VInputTag(
        cms.InputTag("l1pfProducerRegionsBarrel"), 
        cms.InputTag("l1pfProducerRegionsHGCal"),
        #cms.InputTag("l1pfProducerRegionsHGCalNoTK"),
        #cms.InputTag("l1pfProducerRegionsHF")
    ),
    labelsToMerge = cms.vstring("PF", "Puppi"),
)

process.runPF = cms.Sequence( 
    process.l1ParticleFlow +
    process.l1pfProducerRegionsBarrel +
    process.l1pfProducerRegionsHGCal +
    process.l1pfCandidatesRegions
)
process.L1TkPrimaryVertex.nVtx = cms.int32(5)

process.ntuple0 = cms.EDAnalyzer("L1PFCompare",
    #emcalo = cms.InputTag("l1pfProducerHGCal","EmCalo"),
    #calo = cms.InputTag("l1pfProducerHGCal","Calo"),
    emcalo = cms.InputTag(""),
    egcalo = cms.InputTag(""),
    calo = cms.InputTag(""),
    #pf = cms.InputTag("l1pfProducer:PF"),
    pf = cms.InputTag("l1pfCandidatesRegions:PF"),
    #pup = cms.InputTag("l1pfProducer:Puppi"),
    pup = cms.InputTag("l1pfCandidatesRegions:Puppi"),
    vtx = cms.InputTag("L1TkPrimaryVertex"),
    generator = cms.InputTag('genParticles'),
    minPt = cms.double(2.),
    maxEta = cms.double(3.),
    maxN = cms.uint32(9999),
    genIDs = cms.vint32(5,-5,-4,4),
    addGenIDs = cms.vint32(),
    genStatuses = cms.vint32(23,23,23,23),
)

process.p = cms.Path(
    process.offlineBeamSpot +
    process.TTTracksFromTrackletEmulation +
    process.SimL1Emulator + 
    process.runPF + 
    process.ntuple0)

process.TFileService = cms.Service("TFileService", fileName = cms.string("pfTuple.root"))
