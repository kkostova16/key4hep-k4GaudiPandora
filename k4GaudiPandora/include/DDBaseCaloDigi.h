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
#ifndef DDBASECALODIGI_H
#define DDBASECALODIGI_H 1


// DDBaseCaloDigi idea:
// define properties that are common for both ECalDigi and HCalDigi

#include "CalorimeterHitType.h"
#include "DDScintillatorPpdDigi.h"

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/EventHeaderCollection.h"

#include "k4FWCore/Transformer.h"
#include "k4FWCore/MetaDataHandle.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Factories.h"
#include "DDRec/MaterialManager.h"
#include "DDRec/DetectorData.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Constants.h"

#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

#include <random>
#include <string>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

using namespace std;
using namespace dd4hep;
using namespace DDSegmentation;

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

struct DDBaseCaloDigi : public DDCaloDigi_t {
  DDBaseCaloDigi(const std::string& name, ISvcLocator* svcLoc,
		    Gaudi::Functional::details::RepeatValues_<KeyValues, N_in> const& inputs,
		    Gaudi::Functional::details::RepeatValues_<KeyValues, N_out> const& outputs):
    MultiTransformer(name, svcLoc, inputs, outputs) {}

  //StatusCode initialize() override;
  //StatusCode finalize() override;

  //retType operator()(
  //      const edm4hep::SimCalorimeterHitCollection& simCaloHits,
  //      const edm4hep::EventHeaderCollection& headers) const;

  SmartIF<IGeoSvc>         m_geoSvc;
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  std::unique_ptr<DDScintillatorPpdDigi> m_scEcalDigi{};
  std::unique_ptr<DDScintillatorPpdDigi> m_scHcalDigi{};
 
  float scintillatorDigi(float energy, bool isEcal) const;

  //edm4hep::SimCalorimeterHitCollection combineVirtualStripCells(edm4hep::SimCalorimeterHitCollection const& col, bool isBarrel, int stripOrientation ) const;

  int getNumberOfVirtualCells() const;
  int getStripOrientationFromColName(std::string const& colName) const;

  int nRun = 0;
  int nEvt = 0;
  
  float m_barrelPixelSizeT[MAX_LAYERS];
  float m_barrelPixelSizeZ[MAX_LAYERS];
  float m_endcapPixelSizeX[MAX_LAYERS];
  float m_endcapPixelSizeY[MAX_LAYERS];
  float m_barrelStaveDir[MAX_STAVES][2];

  // internal variables
  int m_strip_virt_cells = 999;
  mutable int m_countWarnings = 0;

  enum {
    SQUARE,
    STRIP_ALIGN_ALONG_SLAB,
    STRIP_ALIGN_ACROSS_SLAB,
    SIECAL=0,
    SCECAL
  };
};

//DECLARE_COMPONENT(DDBaseCaloDigi)
#endif
