/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef DDSCCALODIGI_H
#define DDSCCALODIGI_H 1

#include "Gaudi/Property.h"
#include "Gaudi/Accumulators/RootHistogram.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "k4FWCore/Transformer.h"
#include "k4FWCore/MetaDataHandle.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "CalorimeterHitType.h"
#include "DDScintillatorPpdDigi.h"
#include "DDSiliconDigi.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/RandGauss.h"

#include <random>
#include <string>
#include <vector>



/** === DDCaloDigi Processor === <br>
 *  Simple calorimeter digitizer Processor. <br>
 *  Ported from ILDCaloDigi to use DD4hep
 *  Takes SimCalorimeterHit Collections and <br>
 *  produces CalorimeterHit Collections. <br>
 *  Simulated energy depositions in active <br>
 *  layers of calorimeters are <br>
 *  converted into physical energy. This is done <br>
 *  taking into account sampling fractions of <br>
 *  ECAL and HCAL. <br>
 *  User has to specify ECAL and HCAL SimCalorimeterHit <br>
 *  collections with processor parameters <br>
 *  HCALCollections and ECALCollections. <br>
 *  The names of the output CalorimeterHit Collections <br>
 *  are specified with processor parameters <br>
 *  ECALOutputCollection and HCALOutputCollection. <br>
 *  Conversion factors for ECAL and HCAL <br>
 *  are specified via processor parameters  <br>
 *  CalibrECAL and CalibrHCAL. <br>
 *  It should be noted that ECAL and HCAL may consist <br>
 *  of several sections with different sampling fractions. <br>
 *  To handle this situation, calibration coefficients for <br>
 *  ECAL and HCAL are passed as arrays of floats with each element <br>
 *  in this array corresponding to certain section with <br>
 *  a given sampling fraction. <br>
 *  List of layer numbers terminating each section are given through <br>
 *  processor parameters ECALLayers and HCALLayers <br>
 *  There is an option to perform digitization of <br>
 *  both ECAL and HCAL in a digital mode. <br>
 *  Digital mode is activated by  <br>
 *  setting processor parameters <br>
 *  IfDigitalEcal / IfDigitalHcal to 1. <br>
 *  In this case CalibrECAL / CalibrHCAL will  <br>
 *  convert the number of hits into physical energy. <br>
 *  Thresholds on hit energies in ECAL and HCAL <br>
 *  are set with processor parameters <br>
 *  ECALThreshold and HCALThreshold.  <br>
 *  Relations between CalorimeterHits and SimCalorimeterHits <br>
 *  are held in the corresponding relation collection. <br>
 *  The name of this relation collection is specified <br>
 *  via processor parameter RelationOutputCollection. <br>
 *  <h4>Input collections and prerequisites</h4>
 *  SimCalorimeterHit collections <br>
 *  <h4>Output</h4>
 *  CalorimeterHit collections for ECal and HCal. <br>
 *  Collection of relations <br>
 *  between CalorimeterHits and SimCalorimeterHits. <br>
 *  For ECal Calorimeter hits the variable type is set to 0, <br>
 *  whereas for HCal Calorimeter hits the type is set to 1 <br>
 *  @author A. Raspereza (DESY) <br>
 *  @author M. Thomson (DESY) <br>
 *  @version $Id$ <br>
 */

const int MAX_LAYERS = 200;
const int MAX_STAVES =  16;

using retType = std::tuple<edm4hep::CalorimeterHitCollection,
                           edm4hep::CaloHitSimCaloHitLinkCollection>;

using DDCaloDigi_t = k4FWCore::MultiTransformer<retType(
                              const edm4hep::SimCalorimeterHitCollection&,
                              const edm4hep::EventHeaderCollection&)>;

struct DDScCaloDigi final
  : DDCaloDigi_t {
  DDScCaloDigi(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  retType operator()(
        const edm4hep::SimCalorimeterHitCollection& simCaloHits,
        const edm4hep::EventHeaderCollection& headers) const;

  private:
 
  Gaudi::Property<bool> m_inputColIsECAL{this, "InputColIsECAL", {true}, "If input collection is ECAL -- true, or HCAL -- false"};
  
  // digitazing parameters for ECAL and HCAL  
  Gaudi::Property<float> m_threshold{this, "Threshold", {5.0e-5}, "Threshold for ECAL/HCAL Hits in GeV"};
  Gaudi::Property<std::vector<float>> m_thresholdVec{this, "ThresholdVec", {0.00004}, "Threshold vector for HCAL hits in GeV"};
  Gaudi::Property<std::string> m_unitThreshold{this, "ThresholdUnit", {"GeV"}, "Unit for Threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants"};
  Gaudi::Property<std::vector<int>> m_calLayers{this, "CALLayers", {20, 100}, "Index of ECAL/HCAL Layers"};
  Gaudi::Property<std::vector<float>> m_calibrCoeff{this, "Calibration", {40.91, 81.81}, "Calibration coefficients for ECAL/HCAL"};
  Gaudi::Property<std::vector<float>> m_calibrCoeffBarrel{this, "CalibrBarrel", {0.0}, "Calibration coefficients for Barrel"};
  Gaudi::Property<std::vector<float>> m_calibrCoeffEndcap{this, "CalibrEndcap", {0.0}, "Calibration coefficients for Endcap"};
  Gaudi::Property<std::vector<float>> m_calibrCoeffOther{this, "CalibrOther", {0.0}, "Calibration coefficients for Other"};
  Gaudi::Property<int> m_digitalCalo{this, "IfDigitalCalo", {0}, "Digital ECAL/HCAL"};
  Gaudi::Property<int> m_mapsCalCorrection{this, "MapsCalCorrection", {0}, "ECAL correction for theta dependency of calibration for MAPS"};
  Gaudi::Property<int> m_ecalGapCorrection{this, "ECalGapCorrection", {1}, "Correct for Calorimeter gaps"};
  Gaudi::Property<float> m_endcapCorrectionFactor{this, "EndcapCorrectionFactor", {1.025}, "Energy correction for ECAL/HCAL endcap"};
  Gaudi::Property<float> m_ecalGapCorrectionFactor{this, "EcalGapCorrectionFactor", {1.0}, "Factor applied to gap correction"};
  Gaudi::Property<float> m_ecalModuleGapCorrectionFactor{this, "ECalModuleGapCorrectionFactor", {0.5}, "Factor applied to module gap correction"};
  
  
  // timing parameters
  Gaudi::Property<int> m_useTiming{this, "UseTiming", {0}, "Use ECAL/HCAL hit times"};
  Gaudi::Property<int> m_correctTimesForPropagation{this, "CorrectTimesForPropagation", {0}, "Correct ECAL/HCAL hit times for propagation: radial distance/c"};
  Gaudi::Property<float> m_timeWindowMin{this, "TimeWindowMin", {-10.0}, "Time Window minimum time in ns"};
  Gaudi::Property<float> m_endcapTimeWindowMax{this, "EndcapTimeWindowMax", {100.0}, "Endcap Time Window maximum time in ns"};
  Gaudi::Property<float> m_barrelTimeWindowMax{this, "BarrelTimeWindowMax", {100.0}, "Barrel Time Window maximum time in ns"};
  Gaudi::Property<float> m_deltaTimeHitResolution{this, "DeltaTimeHitResolution", {10.0}, "Minimum Delta Time in ns for resolving two hits"};
  Gaudi::Property<float> m_timeResolution{this, "TimeResolution", {10.0}, "Time Resolution used to smear hit times"};
  Gaudi::Property<bool> m_simpleTimingCut{this, "SimpleTimingCut", {true}, 
    "Use simple time window cut on hit times? If false: use original hit-time clustering algorithm. If true: use time window defined by BarrelTimeWindowMin and BarrelTimeWindowMax"};
  
  // parameters for extra digitization effects
  //Gaudi::Property<float> m_calibMip{this, "CalibMIP", {1.0e-4}, "Calibration to convert deposited energy to MIPs"};
  
  // !!!
  Gaudi::Property<int> m_applyEcalDigi{this, "ECALApplyRealisticDigi", {0}, "Apply realistic digitisation to ECAL hits? (0=none, 1=silicon, 2=scintillator)"};

  //Gaudi::Property<float> m_PPD_pe_per_mip{this, "PPD_PE_per_MIP", {7.0}, "# photoelectrons per MIP (in scintillator): used to poisson smear #PEs if > 0"};
  //Gaudi::Property<int> m_PPD_n_pixels{this, "PPD_N_Pixels", {10000}, "Total number of MPPC/SiPM pixels for implementation of saturation effect"};
  //Gaudi::Property<float> m_misCalibNpix{this, "PPD_N_Pixels_uncertainty", {0.05}, "Fractional uncertainty of effective total number of MPPC/SiPM pixels"};
  Gaudi::Property<float> m_misCalib_uncorrel{this, "Miscalibration_uncorrel", {0.0}, "Uncorrelated random gaussian miscalibration (as a fraction: 1.0 = 100%)"};
  Gaudi::Property<bool> m_misCalib_uncorrel_keep{this, "Miscalibration_uncorrel_memorise", {false}, "Store uncorrelated miscalbrations in memory? (WARNING: Can take a lot of memory if used...)"};
  Gaudi::Property<float> m_misCalib_correl{this, "Miscalibration_correl", {0.0}, "Correlated random gaussian miscalibration (as a fraction: 1.0 = 100%)"};
  Gaudi::Property<float> m_deadCellFraction{this, "DeadCellRate", {0.0}, "Random dead cell fraction (as a fraction: [0;1])"};
  Gaudi::Property<bool> m_deadCell_keep{this, "DeadCell_memorise", {false}, "Store dead cells in memory? (WARNING: Can take a lot of memory if used...)"};

  Gaudi::Property<float> m_strip_abs_length{this, "StripAbsorbtionLength", {1000000.0}, "Length scale for absorbtion along scintillator strip (mm)"};

  //Gaudi::Property<float> m_pixSpread{this, "PixelSpread", {0.05}, "Variation of MPPC/SiPM pixels capacitance in ECAL/HCAL (as a fraction: 0.01=1%)"};
  //Gaudi::Property<float> m_elecNoise{this, "ElecNoise_MIPs", {0.0}, "Typical electronics noise for ECAL/HCAL (in MIP units)"};

  //Gaudi::Property<float> m_ehEnergy{this, "energyPerEHpair", {3.6}, "Energy required to create e-h pair in silicon (in eV)"};
  //Gaudi::Property<float> m_maxDynRangeMIP{this, "MaxDynamicRange_MIP", {2500.0}, "Maximum of electronis dynamic range for ECAL/HCAL (in MIPs)"};
  Gaudi::Property<int> m_strip_default_nVirt{this, "Strip_default_nVirtualCells", {9}, "Default number of virtual cells (used if not found in gear file)"};
  Gaudi::Property<std::string> m_deafult_layer_config{this, "DefaultLayerConfig", {"000000000000000"}, "Default ECAL/HCAL layer configuration (used if not found in gear file"};
  
  // !!!
  Gaudi::Property<int> m_applyHcalDigi{this, "HCALApplyRealisticDigi", {0}, "Apply realistic digitisation to HCAL hits? (0=none, 1=silicon, 2=scintillator)"};

  //Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      //"The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  
  SmartIF<IGeoSvc>         m_geoSvc;
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  //const float slop = 0.25; // (mm)
  
  virtual void fillECALGaps(std::vector<edm4hep::MutableCalorimeterHit*> m_calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS],
			    std::vector<int> m_calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS]) const;

  float digitalCalibCoeff(CHT::Layout, float energy) const;
  float analogueCalibCoeff(CHT::Layout, int layer) const;

  float EnergyDigi(float energy, int id) const;

  //edm4hep::SimCalorimeterHitCollection combineVirtualStripCells(edm4hep::SimCalorimeterHitCollection const& col, bool isBarrel, int stripOrientation ) const;

  int getNumberOfVirtualCells() const;
  std::vector<std::pair <int, int>> getLayerConfig() const;
  void checkConsistency(std::string colName, int layer) const;
  std::pair<int, int> getLayerProperties(std::string const& colName, int layer) const;
  int getStripOrientationFromColName(std::string const& colName) const;


  int nRun = 0;
  int nEvt = 0;

  
  float m_zOfEcalEndcap = 0.0;
  float m_barrelPixelSizeT[MAX_LAYERS];
  float m_barrelPixelSizeZ[MAX_LAYERS];
  float m_endcapPixelSizeX[MAX_LAYERS];
  float m_endcapPixelSizeY[MAX_LAYERS];
  float m_barrelStaveDir[MAX_STAVES][2];

  std::unique_ptr<DDScintillatorPpdDigi> m_scDigi{};
  std::unique_ptr<DDSiliconDigi> m_siDigi{};

  // internal variables
  int m_strip_virt_cells = 999;
  mutable int m_countWarnings = 0;
  std::string m_ecalLayout = "";

  CLHEP::MTwistEngine *m_randomEngineDeadCell = NULL;
  float m_event_correl_miscalib = CLHEP::RandGauss::shoot(1.0, m_misCalib_correl);


  std::map <int, float> m_cell_miscalibs{};
  std::map <int, bool> m_cell_dead{};


  enum {
    SQUARE,
    STRIP_ALIGN_ALONG_SLAB,
    STRIP_ALIGN_ACROSS_SLAB,
    SIECAL=0,
    SCECAL
  };
};

DECLARE_COMPONENT(DDScCaloDigi)
#endif
