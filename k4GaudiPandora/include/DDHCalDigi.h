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
#ifndef DDHCALDIGI_H
#define DDHCALDIGI_H 1

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


struct DDHCalDigi final : public DDBaseCaloDigi {
  DDHCalDigi(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  retType operator()(
        const edm4hep::SimCalorimeterHitCollection& simCaloHits,
        const edm4hep::EventHeaderCollection& headers) const;

  private:
  
  // digitazing parameters for HCAL
  Gaudi::Property<std::vector<float>> m_thresholdHcal{this, "HCALThreshold", {0.00004}, "Unit for ECAL Threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants"};
  Gaudi::Property<std::string> m_unitThresholdHcal{this, "HCALThresholdUnit", {"GeV"}, "Unit for HCAL Threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants"};
  Gaudi::Property<std::vector<int>> m_hcalLayers{this, "HCALLayers", {100}, "Index of HCal Layers"};
  Gaudi::Property<std::vector<float>> m_calibrCoeffHcalBarrel{this, "CalibrHCALBarrel", {0.0}, "Calibration coefficients for Barrel HCAL"};
  Gaudi::Property<std::vector<float>> m_calibrCoeffHcalEndcap{this, "CalibrHCALEndcap", {0.0}, "Calibration coefficients for Endcap HCAL"};
  Gaudi::Property<std::vector<float>> m_calibrCoeffHcalOther{this, "CalibrHCALOther", {0.0}, "Calibration coefficients for Other HCAL"};
  Gaudi::Property<int> m_digitalHcal{this, "IfDigitalHcal", {0}, "Digital Hcal"};
  Gaudi::Property<int> m_hcalGapCorrection{this, "HCALGapCorrection", {1}, "Correct for HCAL gaps"};
  Gaudi::Property<float> m_hcalEndcapCorrectionFactor{this, "HCALEndcapCorrectionFactor", {1.025}, "Energy correction for HCAL endcap"};
  Gaudi::Property<float> m_hcalModuleGapCorrectionFactor{this, "HCALModuleGapCorrectionFactor", {0.5}, "Factor applied to module gap correction HCAL"};
    
  // histograms that will be filled --- <1> or <2> is the number of dimensions of the histogram (1D or 2D)
  mutable Gaudi::Accumulators::RootHistogram<1> fHcal{this, "fHcal", "Hcal time profile", {1000, 0., 1000.0}};
  mutable Gaudi::Accumulators::RootHistogram<1> fHcalC{this, "fHcalC", "Hcal time profile cor", {1000, 0., 1000.0}};
  mutable Gaudi::Accumulators::RootHistogram<1> fHcalC1{this, "fHcalC1", "Hcal time profile cor", {100, 0., 1000.0}};
  mutable Gaudi::Accumulators::RootHistogram<1> fHcalC2{this, "fHcalC2", "Hcal time profile cor", {10, 0., 1000.0}};

  mutable Gaudi::Accumulators::RootHistogram<2> fHcalCvsE{this, "fHcalCvsE", "Hcal time profile cor", {{100, 0., 500.0}, {100,0.,10.}}};

  mutable Gaudi::Accumulators::RootHistogram<2> fHcalLayer1 {this, "fHcalLayer1", "Hcal layer 1 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::RootHistogram<2> fHcalLayer11{this, "fHcalLayer11", "Hcal layer 11 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::RootHistogram<2> fHcalLayer21{this, "fHcalLayer21", "Hcal layer 21 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::RootHistogram<2> fHcalLayer31{this, "fHcalLayer31", "Hcal layer 31 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::RootHistogram<2> fHcalLayer41{this, "fHcalLayer41", "Hcal layer 41 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::RootHistogram<2> fHcalLayer51{this, "fHcalLayer51", "Hcal layer 51 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::RootHistogram<2> fHcalLayer61{this, "fHcalLayer61", "Hcal layer 61 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::RootHistogram<2> fHcalLayer71{this, "fHcalLayer71", "Hcal layer 71 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};

  mutable Gaudi::Accumulators::RootHistogram<1> fHcalRLayer1   {this,  "fHcalRLayer1", "Hcal R layer 1", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::RootHistogram<1> fHcalRLayer11  {this, "fHcalRLayer11", "Hcal R layer 11", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::RootHistogram<1> fHcalRLayer21  {this, "fHcalRLayer21", "Hcal R layer 21", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::RootHistogram<1> fHcalRLayer31  {this, "fHcalRLayer31", "Hcal R layer 31", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::RootHistogram<1> fHcalRLayer41  {this, "fHcalRLayer41", "Hcal R layer 41", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::RootHistogram<1> fHcalRLayer51  {this, "fHcalRLayer51", "Hcal R layer 51", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::RootHistogram<1> fHcalRLayer61  {this, "fHcalRLayer61", "Hcal R layer 61", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::RootHistogram<1> fHcalRLayer71  {this, "fHcalRLayer71", "Hcal R layer 71", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::RootHistogram<2> fHcalRLayerNorm{this, "fHcalRLayerNorm", "Hcal R layer Norm", {50, 0., 5000.0}, {50, 0., 5000.0}};

  // timing parameters for HCAL
  Gaudi::Property<int> m_useHcalTiming{this, "UseHcalTiming", {1}, "Use HCAL hit times"};
  Gaudi::Property<int> m_hcalCorrectTimesForPropagation{this, "HCALCorrectTimesForPropagation", {0}, "Correct HCAL hit times for propagation: radial distance/c"};
  Gaudi::Property<float> m_hcalTimeWindowMin{this, "HCALTimeWindowMin", {-10.0}, "HCAL Time Window minimum time in ns"};
  Gaudi::Property<float> m_hcalEndcapTimeWindowMax{this, "HCALEndcapTimeWindowMax", {100.0}, "HCAL Endcap Time Window maximum time in ns"};
  Gaudi::Property<float> m_hcalBarrelTimeWindowMax{this, "HCALBarrelTimeWindowMax", {100.0}, "HCAL Barrel Time Window maximum time in ns"};
  Gaudi::Property<float> m_hcalDeltaTimeHitResolution{this, "HCALDeltaTimeHitResolution", {10.0}, "HCAL Minimum Delta Time in ns for resolving two hits"};
  Gaudi::Property<float> m_hcalTimeResolution{this, "HCALTimeResolution", {10.0}, "HCAL Time Resolution used to smear hit times"};
  Gaudi::Property<bool> m_hcalSimpleTimingCut{this, "HCALSimpleTimingCut", {true}, "Use simple time window cut on hit times? If false: use original hit-time clustering algorithm. If true: use time window defined by HCALBarrelTimeWindowMin and HCALBarrelTimeWindowMax"};
  
  // parameters for extra HCAL digitization effects
  Gaudi::Property<float> m_calibHcalMip{this, "CalibHCALMIP", {1.0e-4}, "Calibration to convert HCAL deposited energy to MIPs"};
  Gaudi::Property<int> m_applyHcalDigi{this, "HCALApplyRealisticDigi", {0}, "Apply realistic digitisation to HCAL hits? (0=none, 1=silicon, 2=scintillator)"};
  Gaudi::Property<float> m_hcal_PPD_pe_per_mip{this, "HCAL_PPD_PE_per_MIP", {7.0}, "# photo-electrons per MIP (scintillator): used to poisson smear #PEs if > 0"};
  Gaudi::Property<int> m_hcal_PPD_n_pixels{this, "HCAL_PPD_N_Pixels",{10000},"HCAL total number of MPPC/SiPM pixels for implementation of saturation effect"};
  Gaudi::Property<float> m_hcal_misCalibNpix{this, "HCAL_PPD_N_Pixels_uncertainty", {0.05}, "HCAL fractional uncertainty of effective total number of MPPC/SiPM pixels"};
  Gaudi::Property<float> m_misCalibHcal_uncorrel{this, "HCAL_miscalibration_uncorrel", {0.0}, "Uncorrelated HCAL random gaussian miscalibration (as a fraction: 1.0 = 100%)"};
  Gaudi::Property<bool> m_misCalibHcal_uncorrel_keep{this, "HCAL_miscalibration_uncorrel_memorise", {false}, "Store oncorrelated HCAL miscalbrations in memory? (WARNING: Can take a lot of memory if used...)"};
  Gaudi::Property<float> m_misCalibHcal_correl{this, "HCAL_miscalibration_correl", {0.0}, "Correlated HCAL random gaussian miscalibration (as a fraction: 1.0 = 100%)"};
  Gaudi::Property<float> m_deadCellFractionHcal{this, "HCAL_deadCellRate", {0.0}, "HCAL random dead cell fraction (as a fraction: [0;1])"};
  Gaudi::Property<bool> m_deadCellHcal_keep{this, "HCAL_deadCell_memorise", {false}, "Store dead HCAL cells in memory? (WARNING: Can take a lot of memory if used...)"};
  Gaudi::Property<float> m_hcal_pixSpread{this, "HCAL_pixel_spread", {0.05}, "Variation of MPPC/SiPM pixels capacitance in HCAL (as a fraction: 0.01=1%)"};
  Gaudi::Property<float> m_hcal_elec_noise{this, "HCAL_elec_noise_mips", {0.}, "Typical electronics noise for HCAL (in MIP units)"};
  Gaudi::Property<float> m_hcalMaxDynMip{this, "HCAL_maxDynamicRange_MIP", {2500.}, "Maximum of electronis dynamic range for HCAL (in MIPs)"};
  //Gaudi::Property<int> m_hcalStrip_default_nVirt{this, "StripHcal_default_nVirtualCells", {9}, "Default number of virtual cells (used if not found in gear file)"};
  //Gaudi::Property<std::string> m_hcal_deafult_layer_config{this, "HCAL_default_layerConfig", {"000000000000000"}, "Default HCAL layer configuration (used if not found in gear file"};
  //Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      //"The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  

  float digitalHcalCalibCoeff(CHT::Layout, float energy) const;
  float analogueHcalCalibCoeff(CHT::Layout, int layer) const;

  float hcalEnergyDigi(float energy, int id) const;

  std::string m_hcalLayout = "";
  
  CLHEP::MTwistEngine *m_randomEngineDeadCellHcal = NULL;
  float m_event_correl_miscalib_hcal = CLHEP::RandGauss::shoot(1.0, m_misCalibHcal_correl);

  std::map <int, float> m_HCAL_cell_miscalibs{};
  std::map <int, bool> m_HCAL_cell_dead{};

};

DECLARE_COMPONENT(DDHCalDigi)
#endif