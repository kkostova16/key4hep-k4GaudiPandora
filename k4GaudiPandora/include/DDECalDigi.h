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
#ifndef DDECALDIGI_H
#define DDECALDIGI_H 1

#include "DDBaseCaloDigi.h"

#include "Gaudi/Property.h"
#include "Gaudi/Accumulators/RootHistogram.h"


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


struct DDECalDigi final : public DDBaseCaloDigi {
  DDECalDigi(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  retType operator()(
        const edm4hep::SimCalorimeterHitCollection& simCaloHits,
        const edm4hep::EventHeaderCollection& headers) const;

  private:
  
  // digitazing parameters for ECAL
  Gaudi::Property<float> m_thresholdEcal{this, "ECALThreshold", {5.0e-5}, "Threshold for ECAL Hits in GeV"};
  Gaudi::Property<std::string> m_unitThresholdEcal{this, "ECALThresholdUnit", {"GeV"}, "Unit for ECAL Threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants"};
  Gaudi::Property<std::vector<int>> m_ecalLayers{this, "ECALLayers", {20, 100}, "Index of ECal Layers"};
  Gaudi::Property<std::vector<float>> m_calibrCoeffEcal{this, "CalibrECAL", {40.91, 81.81}, "Calibration coefficients for ECAL"};
  Gaudi::Property<int> m_digitalEcal{this, "IfDigitalEcal", {0}, "Digital Ecal"};
  Gaudi::Property<int> m_mapsEcalCorrection{this, "MapsEcalCorrection", {0}, "Ecal correction for theta dependency of calibration for MAPS"};
  Gaudi::Property<int> m_ecalGapCorrection{this, "ECALGapCorrection", {1}, "Correct for ECAL gaps"};
  Gaudi::Property<float> m_ecalEndcapCorrectionFactor{this, "ECALEndcapCorrectionFactor", {1.025}, "Energy correction for ECAL endcap"};
  Gaudi::Property<float> m_ecalGapCorrectionFactor{this, "ECALGapCorrectionFactor", {1.0}, "Factor applied to gap correction"};
  Gaudi::Property<float> m_ecalModuleGapCorrectionFactor{this, "ECALModuleGapCorrectionFactor", {0.5}, "Factor applied to module gap correction ECAL"};
  
  // histograms that will be filled --- <1> or <2> is the number of dimensions of the histogram (1D or 2D)
  mutable Gaudi::Accumulators::RootHistogram<2> fEcal{this, "fEcal", "Ecal time profile", {1000, 0., 1000.0}, {1000, 0., 1000.0}};
  mutable Gaudi::Accumulators::RootHistogram<2> fEcalC{this, "fEcalC", "Ecal time profile cor", {1000, 0., 1000.0}, {1000, 0., 1000.0}};
  mutable Gaudi::Accumulators::RootHistogram<2> fEcalC1{this, "fEcalC1", "Ecal time profile cor", {100, 0., 1000.0}, {1000, 0., 1000.0}};
  mutable Gaudi::Accumulators::RootHistogram<2> fEcalC2{this, "fEcalC2", "Ecal time profile cor", {10, 0., 1000.0}, {1000, 0., 1000.0}};

  mutable Gaudi::Accumulators::RootHistogram<2> fEcalLayer1 {this, "fEcalLayer1", "Ecal layer 1 map", {{1800, -4500., 4500.0}, {1800, -4500, 4500.}}};
  mutable Gaudi::Accumulators::RootHistogram<2> fEcalLayer11{this, "fEcalLayer11", "Ecal layer 11 map", {{1800, -4500., 4500.0}, {1800, -4500, 4500.}}};
  mutable Gaudi::Accumulators::RootHistogram<2> fEcalLayer21{this, "fEcalLayer21", "Ecal layer 21 map", {{1800, -4500., 4500.0}, {1800, -4500, 4500.}}};

  mutable Gaudi::Accumulators::RootHistogram<1> fEcalRLayer1   {this, "fEcalRLayer1", "Ecal R layer 1", {100, 0., 5000.0}};
  mutable Gaudi::Accumulators::RootHistogram<1> fEcalRLayer11  {this, "fEcalRLayer11", "Ecal R layer 11", {100, 0., 5000.0}};
  mutable Gaudi::Accumulators::RootHistogram<1> fEcalRLayer21  {this, "fEcalRLayer21", "Ecal R layer 21", {100, 0., 5000.0}};
  mutable Gaudi::Accumulators::RootHistogram<2> fEcalRLayerNorm{this, "fEcalRLayerNorm", "Ecal R layer Norm", {100, 0., 5000.0}, {100, 0., 5000.0}}; 
  
  // timing parameters for ECAL
  Gaudi::Property<int> m_useEcalTiming{this, "UseEcalTiming", {0}, "Use ECAL hit times"};
  Gaudi::Property<int> m_ecalCorrectTimesForPropagation{this, "ECALCorrectTimesForPropagation", {0}, "Correct ECAL hit times for propagation: radial distance/c"};
  Gaudi::Property<float> m_ecalTimeWindowMin{this, "ECALTimeWindowMin", {-10.0}, "ECAL Time Window minimum time in ns"};
  Gaudi::Property<float> m_ecalEndcapTimeWindowMax{this, "ECALEndcapTimeWindowMax", {100.0}, "ECAL Endcap Time Window maximum time in ns"};
  Gaudi::Property<float> m_ecalBarrelTimeWindowMax{this, "ECALBarrelTimeWindowMax", {100.0}, "ECAL Barrel Time Window maximum time in ns"};
  Gaudi::Property<float> m_ecalDeltaTimeHitResolution{this, "ECALDeltaTimeHitResolution", {10.0}, "ECAL Minimum Delta Time in ns for resolving two hits"};
  Gaudi::Property<float> m_ecalTimeResolution{this, "ECALTimeResolution", {10.0}, "ECAL Time Resolution used to smear hit times"};
  Gaudi::Property<bool> m_ecalSimpleTimingCut{this, "ECALSimpleTimingCut", {true}, "Use simple time window cut on hit times? If false: use original hit-time clustering algorithm. If true: use time window defined by ECALBarrelTimeWindowMin and ECALBarrelTimeWindowMax"};
  
  // parameters for extra ECAL digitization effects
  Gaudi::Property<float> m_calibEcalMip{this, "CalibECALMIP", {1.0e-4}, "Calibration to convert ECAL deposited energy to MIPs"};
  Gaudi::Property<int> m_applyEcalDigi{this, "ECALApplyRealisticDigi", {0}, "Apply realistic digitisation to ECAL hits? (0=none, 1=silicon, 2=scintillator)"};
  Gaudi::Property<float> m_ecal_PPD_pe_per_mip{this, "ECAL_PPD_PE_per_MIP", {7.0}, "# photoelectrons per MIP (scintillator): used to poisson smear #PEs if > 0"};
  Gaudi::Property<int> m_ecal_PPD_n_pixels{this, "ECAL_PPD_N_Pixels",{10000},"ECAL total number of MPPC/SiPM pixels for implementation of saturation effect"};
  Gaudi::Property<float> m_ecal_misCalibNpix{this, "ECAL_PPD_N_Pixels_uncertainty", {0.05}, "ECAL fractional uncertainty of effective total number of MPPC/SiPM pixels"};
  Gaudi::Property<float> m_misCalibEcal_uncorrel{this, "ECAL_miscalibration_uncorrel", {0.0}, "Uncorrelated ECAL random gaussian miscalibration (as a fraction: 1.0 = 100%)"};
  Gaudi::Property<bool> m_misCalibEcal_uncorrel_keep{this, "ECAL_miscalibration_uncorrel_memorise", {false}, "Store uncorrelated ECAL miscalbrations in memory? (WARNING: Can take a lot of memory if used...)"};
  Gaudi::Property<float> m_misCalibEcal_correl{this, "ECAL_miscalibration_correl", {0.0}, "Correlated ECAL random gaussian miscalibration (as a fraction: 1.0 = 100%)"};
  Gaudi::Property<float> m_deadCellFractionEcal{this, "ECAL_deadCellRate", {0.0}, "ECAL random dead cell fraction (as a fraction: [0;1])"};
  Gaudi::Property<bool> m_deadCellEcal_keep{this, "ECAL_deadCell_memorise", {false}, "Store dead ECAL cells in memory? (WARNING: Can take a lot of memory if used...)"};
  Gaudi::Property<float> m_strip_abs_length{this, "ECAL_strip_absorbtionLength", {1000000.0}, "Length scale for absorbtion along scintillator strip (mm)"};
  Gaudi::Property<float> m_ecal_pixSpread{this, "ECAL_pixel_spread", {0.05}, "Variation of MPPC/SiPM pixels capacitance in ECAL (as a fraction: 0.01=1%)"};
  Gaudi::Property<float> m_ecal_elec_noise{this, "ECAL_elec_noise_mips", {0.0}, "Typical electronics noise for ECAL (in MIP units)"};
  Gaudi::Property<float> m_ehEnergy{this, "energyPerEHpair", {3.6}, "Energy required to create e-h pair in silicon (in eV)"};
  Gaudi::Property<float> m_ecalMaxDynMip{this, "ECAL_maxDynamicRange_MIP", {2500.0}, "Maximum of electronis dynamic range for ECAL (in MIPs)"};
  Gaudi::Property<int> m_ecalStrip_default_nVirt{this, "StripEcal_default_nVirtualCells", {9}, "Default number of virtual cells (used if not found in gear file)"};
  Gaudi::Property<std::string> m_ecal_deafult_layer_config{this, "ECAL_default_layerConfig", {"000000000000000"}, "Default ECAL layer configuration (used if not found in gear file"};
  
  virtual void fillECALGaps(std::vector<edm4hep::MutableCalorimeterHit*> m_calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS],
			    std::vector<int> m_calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS]) const;

  float digitalEcalCalibCoeff(int layer) const;
  float analogueEcalCalibCoeff(int layer) const;

  float ecalEnergyDigi(float energy, int id) const;

  float siliconDigi(float energy) const;

  std::vector<std::pair <int, int>> getLayerConfig() const;
  std::pair<int, int> getLayerProperties(std::string const& colName, int layer) const;
  void checkConsistency(std::string colName, int layer) const;

  
  float m_zOfEcalEndcap = 0.0;
 
  std::string m_ecalLayout = "";

  CLHEP::MTwistEngine *m_randomEngineDeadCellEcal = NULL;
  float m_event_correl_miscalib_ecal = CLHEP::RandGauss::shoot(1.0, m_misCalibEcal_correl);

  std::map <int, float> m_ECAL_cell_miscalibs{};
  std::map <int, bool> m_ECAL_cell_dead{};
 
};

DECLARE_COMPONENT(DDECalDigi)
#endif