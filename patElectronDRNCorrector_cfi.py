import FWCore.ParameterSet.Config as cms

#from RecoEgamma.EgammaTools.patElectronDRNCorrectionProducer_cfi import patElectronDRNCorrectionProducer
"""
patElectronsDRN = patElectronDRNCorrectionProducer.clone(
                            particleSource = 'selectedPatElectrons',
                            #particleSource = 'gedGsfElectrons',
                            rhoName = 'fixedGridRhoFastjetAll',
                            Client = patElectronDRNCorrectionProducer.Client.clone(
                              mode = 'Async',
                              allowedTries = 1,
                              #modelName = 'electronObjectEnsemble',
                              modelName = 'photonObjectCombined',
                              #modelConfigPath = 'RecoEgamma/EgammaElectronProducers/data/models/electronObjectEnsemble/config.pbtxt',
                              modelConfigPath = 'RecoEgamma/EgammaElectronProducers/data/models/photonObjectCombined/config.pbtxt',
                              timeout = 10
                            )
    )
"""
from RecoEgamma.EgammaTools.gsfElectronDRNCorrectionProducer_cfi import gsfElectronDRNCorrectionProducer
gsfElectronsDRN = gsfElectronDRNCorrectionProducer.clone(
                            #particleSource = 'selectedPatElectrons',
                            particleSource = 'gedGsfElectrons',
                            rhoName = 'fixedGridRhoFastjetAll',
                            Client = gsfElectronDRNCorrectionProducer.Client.clone(
                              mode = 'Async',
                              allowedTries = 1,
                              #modelName = 'photonObjectCombined',
                              #modelName = 'electronObjectEnsemble',
                              modelName = 'electron_model',
                              modelConfigPath = 'RecoEgamma/EgammaElectronProducers/data/models/electron_model/config.pbtxt',
                              #modelConfigPath = 'RecoEgamma/EgammaPhotonProducers/data/models/photonObjectCombined/config.pbtxt',
                              timeout = 10
                            )
    )
print("Runnning gsfElectronDRN")
print("model : RecoEgamma/EgammaElectronProducers/data/models/electron_model/config.pbtxt")
