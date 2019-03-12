import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("RESP", eras.Phase2C4_trigger)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:inputs104X.root'),
    #fileNames = cms.untracked.vstring('file:/eos/cms/store/cmst3/user/gpetrucc/l1tr/105X/NewInputs104X/010319/ParticleGun_PU0/inputs104X_ParticleGun_PU0_job1.root'),
    #fileNames = cms.untracked.vstring('root://eoscms.cern.ch//store/cmst3/group/hzz/gpetrucc/tmp/prod104X/more-010319/ParticleGunPt0p5To5_PU200/inputs104X_ParticleGunPt0p5To5_PU200_job18.root'),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    skipBadFiles = cms.untracked.bool(True)
)

process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

process.load("L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff")
process.load("L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_split_cff")

process.pfClustersFromL1EGClustersRaw    = process.pfClustersFromL1EGClusters.clone(corrector = "")
process.pfClustersFromHGC3DClustersRaw   = process.pfClustersFromHGC3DClusters.clone(corrector = "")
process.pfClustersFromHGC3DClustersEMRaw = process.pfClustersFromHGC3DClustersRaw.clone(emOnly = True, etMin = 0.)

process.pfClustersFromCombinedCaloHCalUnclust = process.pfClustersFromCombinedCaloHCal.clone(
    ecalCandidates = []
)
process.runPF = cms.Sequence( 
    process.l1ParticleFlow_proper + # excludes the prerequisites (3D clusters and L1EG clusters)
    process.pfClustersFromL1EGClustersRaw +
    process.pfClustersFromHGC3DClustersRaw +
    process.pfClustersFromHGC3DClustersEMRaw
    + process.pfClustersFromCombinedCaloHCalUnclust
    + process.pfClustersFromHGC3DClustersEM 
    + process.pfClustersFromL1EGClusters 
    + process.pfClustersFromCombinedCalo 
    + process.pfTracksFromL1Tracks
    + process.l1pfProducer
    + process.l1PuppiForMET
)

#FROM MET TOOLS
class GatherAllModulesVisitor(object):
    """Visitor that travels within a cms.Sequence, and returns a list of objects of type gatheredInance(e.g. modules) that have it"""
    def __init__(self, gatheredInstance=cms._Module):
        self._modules = []
        self._gatheredInstance= gatheredInstance
    def enter(self,visitee):
        if isinstance(visitee,self._gatheredInstance):
            self._modules.append(visitee)
    def leave(self,visitee):
        pass
    def modules(self):
        return self._modules

def listModules(sequence):
    visitor = GatherAllModulesVisitor(gatheredInstance=cms._Module)
    sequence.visit(visitor)
    return visitor.modules()

def listSequences(sequence):
    visitor = GatherAllModulesVisitor(gatheredInstance=cms.Sequence)
    sequence.visit(visitor)
    return visitor.modules()

def __labelsInSequence(process, sequenceLabel, postfix="", keepPostFix=False):
    position = -len(postfix)
    if keepPostFix: 
        position = None

    result = [ m.label()[:position] for m in listModules( getattr(process,sequenceLabel+postfix))]
    result.extend([ m.label()[:position] for m in listSequences( getattr(process,sequenceLabel+postfix))]  )
    if postfix == "":
        result = [ m.label() for m in listModules( getattr(process,sequenceLabel+postfix))]
        result.extend([ m.label() for m in listSequences( getattr(process,sequenceLabel+postfix))]  )
    return result



process.ntuple0 = cms.EDAnalyzer("L1PFCompare",
    #pf = cms.InputTag("l1pfProducer:PF"),
    pf = cms.InputTag("l1pfCandidates:PF"),
    #pup = cms.InputTag("l1pfProducer:Puppi"),
    pup = cms.InputTag("l1pfCandidates:Puppi"),
    generator = cms.InputTag('genParticles'),
    minPt = cms.double(2.),
    maxEta = cms.double(6.),
    maxN = cms.uint32(9999),
)

process.l1ParticleFlow1 = cms.Sequence()
for l_ in  __labelsInSequence(process,"l1ParticleFlow_proper"):
    if (not hasattr(getattr(process,l_),"_seq")):
        newl = getattr(process,l_).clone()
        setattr(process,l_+"1",newl)
        process.l1ParticleFlow1 += getattr(process,l_+"1")

process.l1ParticleFlow1.replace("l1pfProducerBarrel","l1pfProducerBarrel1")
process.l1ParticleFlow1.replace("l1pfProducerHGCal","l1pfProducerHGCal1")
process.l1ParticleFlow1.replace("l1pfProducerHF","l1pfProducerHF1")
process.l1pfProducerBarrel1.vtxNum = cms.untracked.uint32(10)
process.l1pfProducerHGCal1.vtxNum = cms.untracked.uint32(10)
process.l1pfProducerHF1.vtxNum = cms.untracked.uint32(10)
process.l1pfCandidates1.pfProducers = cms.VInputTag(
        cms.InputTag("l1pfProducerBarrel1"),
        cms.InputTag("l1pfProducerHGCal1"),
        cms.InputTag("l1pfProducerHF1")
    )

process.ntuple1 = process.ntuple0.clone(pf = cms.InputTag("l1pfCandidates1:PF"), pup = cms.InputTag("l1pfCandidates1:Puppi"))

process.p = cms.Path(process.runPF + process.ntuple0 + process.l1ParticleFlow1 + process.ntuple1)

process.TFileService = cms.Service("TFileService", fileName = cms.string("pfTuple.root"))
