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
from Configurables import DDHCalDigi

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

calodigi = [DDHCalDigi("HCALBarrelDigi"),
            DDHCalDigi("HCALEndcapDigi"),
            DDHCalDigi("HCALRingDigi")]

inputcollections = [["HCalBarrelCollection"],
                    ["HCalEndcapCollection"],
                    ["HCalRingCollection"]]

outputcollections = [["HCALBarrel"],
                     ["HCALEndcap"],
                     ["HCALRing"]]

relcollections = [["RelationCaloHitHCALBarrel"],
                  ["RelationCaloHitHCALEndcap"],
                  ["RelationCaloHitHCALRing"]]

#set properties
for calodigicol, inputcol, outputcol, relcol in zip(calodigi, inputcollections, outputcollections, relcollections):
    
    calodigicol.InputCaloHitCollection = inputcol           # "HCalBarrelCollection","HCalEndcapCollection","HCalRingCollection"
    calodigicol.OutputCaloHitCollection = outputcol
    calodigicol.RelCollection = relcol
    
    # digitazing parameters for HCAL
    calodigicol.HCALThreshold = [0.00025]
    calodigicol.HCALThresholdUnit = "GeV"
    calodigicol.HCALLayers = [100]
    calodigicol.CalibrHCALBarrel = [45.9956826061]
    calodigicol.CalibrHCALEndcap = [46.9252540291]
    calodigicol.CalibrHCALOther = [57.4588011802]
    calodigicol.IfDigitalHcal = 0
    calodigicol.HCALEndcapCorrectionFactor = 1.000
    calodigicol.HCALModuleGapCorrectionFactor = 0.5

    # timing parameters for HCAL
    calodigicol.UseHcalTiming = 1
    calodigicol.HCALCorrectTimesForPropagation = 1
    calodigicol.HCALTimeWindowMin = -1.0
    calodigicol.HCALEndcapTimeWindowMax = 10.0
    calodigicol.HCALBarrelTimeWindowMax = 10.0
    calodigicol.HCALDeltaTimeHitResolution = 10.0
    calodigicol.HCALTimeResolution = 10.0
    calodigicol.HCALSimpleTimingCut = True

    # parameters for extra HCAL digitization effects
    calodigicol.CalibHCALMIP = 1.0e-4
    calodigicol.HCALApplyRealisticDigi = 0
    calodigicol.HCAL_PPD_PE_per_MIP = 10.0
    calodigicol.HCAL_PPD_N_Pixels = 400
    calodigicol.HCAL_PPD_N_Pixels_uncertainty = 0.05
    calodigicol.HCAL_miscalibration_uncorrel = 0.0
    calodigicol.HCAL_miscalibration_uncorrel_memorise = False
    calodigicol.HCAL_miscalibration_correl = 0.0
    calodigicol.HCAL_deadCellRate = 0.0
    calodigicol.HCAL_deadCell_memorise = False 
    calodigicol.HCAL_pixel_spread = 0.0
    calodigicol.HCAL_elec_noise_mips = 0.0
    calodigicol.HCAL_maxDynamicRange_MIP = 200.0

#
iosvc = IOSvc()
iosvc.input = "../simulation/sim_partgun_1000.root"
iosvc.output = "../outputfiles/DDCaloDigi/outputHCalDigi_Gaudi.root"

hps = RootHistSvc("HistogramPersistencySvc")
root_hist_svc = RootHistoSink("RootHistoSink")
root_hist_svc.FileName = "../outputfiles/DDCaloDigi/ddhcaldigi_hist.root"

ApplicationMgr(TopAlg=[calodigi[0],
                       calodigi[1],
                       calodigi[2]],
               EvtSel="NONE",
               EvtMax=-1,
               ExtSvc=[EventDataSvc("EventDataSvc"), root_hist_svc],
               OutputLevel=INFO,
               )
