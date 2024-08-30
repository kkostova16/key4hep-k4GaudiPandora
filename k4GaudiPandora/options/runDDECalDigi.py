#
# Copyright (c) 2020-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import DDECalDigi

from Configurables import GeoSvc
from Configurables import UniqueIDGenSvc
from Configurables import RootHistSvc
from Configurables import Gaudi__Histograming__Sink__Root as RootHistoSink
import os

id_service = UniqueIDGenSvc("UniqueIDGenSvc")

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.environ["K4GEO"]+"/FCCee/CLD/compact/CLD_o2_v06/CLD_o2_v06.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

calodigi = [DDECalDigi("ECALBarrelDigi"),
            DDECalDigi("ECALEndcapDigi")]

inputcollections = [["ECalBarrelCollection"],
                    ["ECalEndcapCollection"]]

outputcollections = [["ECALBarrel"],
                     ["ECALEndcap"]]

relcollections = [["RelationCaloHitECALBarrel"],
                  ["RelationCaloHitECALEndcap"]]

#set properties
for calodigicol, inputcol, outputcol, relcol in zip(calodigi, inputcollections, outputcollections, relcollections):
    
    calodigicol.InputCaloHitCollection = inputcol           # "ECalBarrelCollection","ECalEndcapCollection"
    calodigicol.OutputCaloHitCollection = outputcol
    calodigicol.RelCollection = relcol
    
    # digitazing parameters for ECAL
    calodigicol.ECALThreshold = 5.0e-5
    calodigicol.ECALThresholdUnit = "GeV"
    calodigicol.ECALLayers = [41,100]
    calodigicol.CalibrECAL = [37.5227197175, 37.5227197175]
    calodigicol.IfDigitalEcal = 0
    calodigicol.MapsEcalCorrection = 0
    calodigicol.ECALGapCorrection = 1
    calodigicol.ECALGapCorrection = 1
    calodigicol.ECALEndcapCorrectionFactor = 1.03245503522
    calodigicol.ECALGapCorrectionFactor = 1.0
    calodigicol.ECALModuleGapCorrectionFactor = 0.0

    # timing parameters for ECAL
    calodigicol.UseEcalTiming = 1
    calodigicol.ECALCorrectTimesForPropagation = 1 
    calodigicol.ECALTimeWindowMin = -1.0
    calodigicol.ECALEndcapTimeWindowMax = 10.0
    calodigicol.ECALBarrelTimeWindowMax = 10.0
    calodigicol.ECALDeltaTimeHitResolution = 10.0
    calodigicol.ECALTimeResolution = 10.0
    calodigicol.ECALSimpleTimingCut = True

    # parameters for extra ECAL digitization effects
    calodigicol.CalibECALMIP = 1.0e-4
    calodigicol.ECALApplyRealisticDigi = 0
    calodigicol.ECAL_PPD_PE_per_MIP = 7.0
    calodigicol.ECAL_PPD_N_Pixels = 10000
    calodigicol.ECAL_PPD_N_Pixels_uncertainty = 0.05
    calodigicol.ECAL_miscalibration_uncorrel = 0.0
    calodigicol.ECAL_miscalibration_uncorrel_memorise = False
    calodigicol.ECAL_miscalibration_correl = 0.0
    calodigicol.ECAL_deadCellRate = 0.0
    calodigicol.ECAL_deadCell_memorise = False
    calodigicol.ECAL_strip_absorbtionLength = 1.0e6
    calodigicol.ECAL_pixel_spread = 0.05
    calodigicol.ECAL_elec_noise_mips = 0.0
    calodigicol.energyPerEHpair = 3.6
    calodigicol.ECAL_maxDynamicRange_MIP = 2500.0
    calodigicol.StripEcal_default_nVirtualCells = 9
    calodigicol.ECAL_default_layerConfig = "000000000000000"

#
iosvc = IOSvc()
iosvc.input = "../simulation/sim_partgun_1000.root"
iosvc.output = "../outputfiles/DDCaloDigi/outputECalDigi_Gaudi.root"

hps = RootHistSvc("HistogramPersistencySvc")
root_hist_svc = RootHistoSink("RootHistoSink")
root_hist_svc.FileName = "../outputfiles/DDCaloDigi/ddecaldigi_hist.root"

ApplicationMgr(TopAlg=[calodigi[0],
                       calodigi[1]],
               EvtSel="NONE",
               EvtMax=-1,
               ExtSvc=[EventDataSvc("EventDataSvc"), root_hist_svc],
               OutputLevel=INFO,
               )
