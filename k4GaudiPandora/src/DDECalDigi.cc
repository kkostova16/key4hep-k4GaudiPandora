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

#include "DDECalDigi.h"

// protect against rounding errors
// will not find caps smaller than this
const float slop = 0.25;  // (mm)
//const float pi = acos(-1.0); ///FIXME
//const float twopi = 2.0*pi;  ///FIXME: DD4HEP INTERFERES WITH THESE


// Forward Declaration, gets linked in from DDPandoraPFANewProcessor --- FIXME
                                                                    // now declare here, to be linked from DDPandoraPFANewProcessor
dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag = 0) {
  dd4hep::rec::LayeredCalorimeterData* theExtension = 0;
  
  dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
  const std::vector<dd4hep::DetElement>& theDetectors =
      dd4hep::DetectorSelector(mainDetector).detectors(includeFlag, excludeFlag);

  cout << " getExtension :  includeFlag: " << dd4hep::DetType(includeFlag)
                        << " excludeFlag: " << dd4hep::DetType(excludeFlag) << "  found : " << theDetectors.size()
                        << "  - first det: " << theDetectors.at(0).name() << endl;

  if (theDetectors.size() != 1) {
    stringstream es;
    es << " getExtension: selection is not unique (or empty)  includeFlag: " << dd4hep::DetType(includeFlag)
       << " excludeFlag: " << dd4hep::DetType(excludeFlag) << " --- found detectors : ";
    for (unsigned i = 0, N = theDetectors.size(); i < N; ++i) {
      es << theDetectors.at(i).name() << ", ";
    }
    throw std::runtime_error(es.str());
  }

  theExtension = theDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();

  return theExtension;
};


DDECalDigi::DDECalDigi(const std::string& aName, ISvcLocator* aSvcLoc)
    : DDBaseCaloDigi(aName, aSvcLoc,
          {
            KeyValues("InputCaloHitCollection", {"ECalBarrelCollection"}),
            KeyValues("HeaderName", {"EventHeader"})
          },
          {
            KeyValues("OutputCaloHitCollection", {"ECalorimeterHit1"}), //, "ECalorimeterHit2", "ECalorimeterHit3"}),
            KeyValues("RelCollection", {"RelationCaloHit"})
          })
{
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
    if (!m_uidSvc) {
      error() << "Unable to get UniqueIDGenSvc" << endmsg;
    }
  m_geoSvc = serviceLocator()->service("GeoSvc");  // important to initialize m_geoSvc

// // helper struct for string comparision --- Is needed???
// struct XToLower{
//   int operator() ( int ch ) {
//     return std::tolower ( ch );
//   }
// }
  
//m_geoSvc = aSvcLoc->service(m_geoSvcName); // ???
}

StatusCode DDECalDigi::initialize() {

  m_strip_virt_cells = -999;
  m_countWarnings    = 0;

  try {
    dd4hep::rec::LayeredCalorimeterData* ecalBarrelData =
        getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
                     (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));

    dd4hep::rec::LayeredCalorimeterData* ecalEndcapData =
        getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP),
                     (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));

    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& ecalBarrelLayers = ecalBarrelData->layers;
    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& ecalEndcapLayers = ecalEndcapData->layers;

    // Determine geometry of ECAL
    int symmetry   = ecalBarrelData->inner_symmetry;
    m_zOfEcalEndcap = ecalEndcapData->extent[2] / dd4hep::mm;

    // Determine ECAL polygon angles
    // Store radial vectors perpendicular to stave layers in m_ecalBarrelStaveDir
    // ASSUMES Mokka Stave numbering 0 = top, then numbering increases anti-clockwise
    if (symmetry > 1) {
      float nFoldSymmetry = static_cast<float>(symmetry);
      float phi0          = ecalBarrelData->phi0 / dd4hep::rad;
      for (int i = 0; i < symmetry; ++i) {
        float phi             = phi0 + i * dd4hep::twopi / nFoldSymmetry;
        m_barrelStaveDir[i][0] = cos(phi);
        m_barrelStaveDir[i][1] = sin(phi);
      }
    }

    for (unsigned int i = 0; i < ecalBarrelLayers.size(); ++i) {
      m_barrelPixelSizeT[i] = ecalBarrelLayers[i].cellSize0;
      m_barrelPixelSizeZ[i] = ecalBarrelLayers[i].cellSize1;
      debug() << "barrel pixel size " << i << " " << m_barrelPixelSizeT[i] << " " << m_barrelPixelSizeZ[i] << endl;
    }

    for (unsigned int i = 0; i < ecalEndcapLayers.size(); ++i) {
      m_endcapPixelSizeX[i] = ecalEndcapLayers[i].cellSize0;
      m_endcapPixelSizeY[i] = ecalEndcapLayers[i].cellSize1;
      debug() << "endcap pixel size " << i << " " << m_endcapPixelSizeX[i] << " " << m_endcapPixelSizeY[i] << endl;
    }

    m_strip_virt_cells = m_ecalStrip_default_nVirt;
       warning() << "taking number of virtual cells from steering file (FIXME!): " << m_strip_virt_cells << endl;
    m_ecalLayout = m_ecal_deafult_layer_config;
       warning() << "taking layer layout from steering file (FIXME): " << m_ecalLayout << endl;

  } catch (std::exception& e) {
      error() << "Could not get ECAL parameters from DD4hep!" << endmsg;
  }

  // Convert ECAL thresholds to GeV units
  if (m_unitThresholdEcal.value().compare("GeV") == 0) {
    // ECAL threshold unit is GeV, do nothing
  } else if (m_unitThresholdEcal.value().compare("MIP") == 0) {
    // ECAL threshold unit is MIP, convert via MIP2GeV
    m_thresholdEcal.value() *= m_calibEcalMip.value();
  } else if (m_unitThresholdEcal.value().compare("px") == 0) {
    // ECAL threshold unit is pixels, convert via MIP2GeV and lightyield
    m_thresholdEcal.value() *= m_ecal_PPD_pe_per_mip.value() * m_calibEcalMip.value();
  } else {
    error() << "Could not identify ECAL threshold unit. Please use \"GeV\", \"MIP\" or \"px\"! Aborting." << endmsg;
    assert(0);
  }

  // Set up the scintillator/MPPC digitizer
  m_scEcalDigi = std::unique_ptr<DDScintillatorPpdDigi>(new DDScintillatorPpdDigi()); 
  m_scEcalDigi->setPEperMIP(m_ecal_PPD_pe_per_mip);
  m_scEcalDigi->setCalibMIP(m_calibEcalMip);
  m_scEcalDigi->setNPix(m_ecal_PPD_n_pixels);
  m_scEcalDigi->setRandomMisCalibNPix(m_ecal_misCalibNpix);
  m_scEcalDigi->setPixSpread(m_ecal_pixSpread);
  m_scEcalDigi->setElecNoise(m_ecal_elec_noise);
  m_scEcalDigi->setElecRange(m_ecalMaxDynMip);
  cout << "ECAL sc digi:" << endl;
  m_scEcalDigi->printParameters();

  
  // Set up the random engines for ECAL dead cells: (could use a steering parameter though)
  if (m_deadCellEcal_keep) {
    m_randomEngineDeadCellEcal = new CLHEP::MTwistEngine(0, 0);
  } else {
    m_randomEngineDeadCellEcal = 0;
  }

  return StatusCode::SUCCESS;
}

retType DDECalDigi::operator()(
    const edm4hep::SimCalorimeterHitCollection& simCaloHitCollection,
    const edm4hep::EventHeaderCollection& eventHeaders) const {
    auto const& headers = eventHeaders.at(nRun); //FIXME maybe
  debug() << " process event : " << headers.getEventNumber() << " - run  " << headers.getRunNumber()
          << endmsg;  // headers[0].getRunNumber(),headers[0].getEventNumber()

  auto col = edm4hep::CalorimeterHitCollection();           // create output CalorimeterHit collection
  auto Relcol = edm4hep::CaloHitSimCaloHitLinkCollection(); // create relation collection CalorimeterHit-SimCalorimeterHit

  // copy the flags from the input collection
  //_flag.setBit(LCIO::CHBIT_LONG);
  //_flag.setBit(LCIO::RCHBIT_TIME);  //store timing on output hits.

  // decide on this event's correlated miscalibration
  //_event_correl_miscalib_ecal = CLHEP::RandGauss::shoot( 1.0, _misCalibEcal_correl );


  std::string colName = m_inputLocations[0][0].key();       // take input collection name
  cout << "looking for collection: " << colName << endl;

  if (colName.find("dummy") != string::npos) {
    debug() << "Ignoring input ECAL collection name (looks like dummy name)" << colName << endmsg;
  }
  
  //fg: need to establish the subdetetcor part here
  //    use collection name as cellID does not seem to have that information
  
  CHT::Layout caloLayout = layoutFromString(colName);

  // auto col = evt->getCollection( colName.c_str() ) ;
  //std::string initString; // = "FIXME"; // cellIDHandle.get();
  //std::string initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  
  //FIXME: take this from the metadata
  std::string initString = "system:5,side:2,module:8,stave:4,layer:9,submodule:4,x:32:-16,y:-16";
  
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);  // check!
                                                                    // check if decoder contains "layer"
    //CellIDDecoder<SimCalorimeterHit> idDecoder( col ); --- ???

    // Create new collection --- ???
      //LCCollectionVec *ecalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
      //auto ecalcol    = edm4hep::CalorimeterHit();
      //ecalcol->setFlag(_flag.getFlag());

      // if making gap corrections clear the vectors holding pointers to calhits --- ???
      std::vector<edm4hep::MutableCalorimeterHit*> m_calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS];
      std::vector<int> m_calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS];

      // deal with strips split into virtual cells
      //  if this collection is a strip which has been split into virtual cells, they need to be recombined
      int orientation = getStripOrientationFromColName(colName);
      if (orientation != SQUARE && getNumberOfVirtualCells() > 1) {
        //auto fixmeCollectionUnused = combineVirtualStripCells(inputCaloHitCollection, caloLayout == CHT::barrel, orientation);
      }
  
  
  //
  // * Reading Collections of ECAL Simulated Hits *
  //
    for (const auto& hit : simCaloHitCollection) {
      //  SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
      const int cellID = hit.getCellID();
      float energy = hit.getEnergy();
          
      // apply threshold cut
      if (energy > m_thresholdEcal) {
        unsigned int layer  = bitFieldCoder.get(cellID, "layer");
        unsigned int stave  = bitFieldCoder.get(cellID, "stave");
        unsigned int module = bitFieldCoder.get(cellID, "module");

        // check that layer and assumed layer type are compatible
        checkConsistency(colName, layer);

        // save hits by module/stave/layer if required later ???

        float x    = hit.getPosition()[0];
        float y    = hit.getPosition()[1];
        float z    = hit.getPosition()[2];
        float r    = sqrt(x * x + y * y + z * z);
        float rxy  = sqrt(x * x + y * y);
        float cost = fabs(z) / r;

        // fill ECAL Layer histograms
        if (z > 0) {
          if (layer == 1)
            ++fEcalLayer1[{x, y}];
          if (layer == 11)
            ++fEcalLayer11[{x, y}];
          if (layer == 21)
            ++fEcalLayer21[{x, y}];
          if (layer == 1)
            ++fEcalRLayer1[rxy];
          if (layer == 11)
            ++fEcalRLayer11[rxy];
          if (layer == 21)
            ++fEcalRLayer21[rxy];
        }

        // evaluate the ECAL calibration coefficient
        float calibr_coeff = 1.;
        if (m_digitalEcal) {
          calibr_coeff = this->digitalEcalCalibCoeff(layer);
          if (m_mapsEcalCorrection) {
            if (caloLayout == CHT::barrel) {
              float correction = 1.1387 - 0.068 * cost - 0.191 * cost * cost;
              calibr_coeff /= correction;
            } else {
              float correction = 0.592 + 0.590 * cost;
              calibr_coeff /= correction;
            }
          }
        } else {
        calibr_coeff = this->analogueEcalCalibCoeff(layer);
        }
        // if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)calibr_coeff *= m_ecalEndcapCorrectionFactor;
        if (caloLayout != CHT::endcap) {
          calibr_coeff *= m_ecalEndcapCorrectionFactor;  // more robust
        }

        // apply timing cut for ECAL
        if (m_useEcalTiming) {
          float ecalTimeWindowMax;
          if (caloLayout == CHT::barrel) {                // current SimHit is in barrel, use barrel timing cut
            ecalTimeWindowMax = m_ecalBarrelTimeWindowMax;
          } else {                                        // current simhit is not in barrel, use endcap timing cut
            ecalTimeWindowMax = m_ecalEndcapTimeWindowMax;
          }

          float dt = r / 300. - 0.1;
          auto  ecalSingleHit  = hit.getContributions();
          int   count       = 0; 
          float eCellInTime = 0.;
          float eCellOutput = 0.;
          const unsigned int n = ecalSingleHit.size();
          std::vector<bool> used(n, false);
          
          for (unsigned int i_t = 0; i_t < n; i_t++) {        // loop over all subhits
            float timei = ecalSingleHit[i_t].getTime();       // absolute hit timing of current subhit
            float energyi = ecalSingleHit[i_t].getEnergy();   // energy of current subhit
            float energySum = 0;

            float deltat = 0;
            if (m_ecalCorrectTimesForPropagation) {
              deltat = dt; //deltat now carries hit timing correction.
            }
            if (timei - deltat > m_ecalTimeWindowMin.value() && timei - deltat < ecalTimeWindowMax) {
              float ecor = energyi * calibr_coeff;
              eCellInTime += ecor;
            }

            if (!used[i_t]) {  //if current subhit has not been merged with previous hits already, take current hit as starting point to merge hits
              // merge with other hits?
              used[i_t] = true;
              for (unsigned int j_t = i_t + 1; j_t < n; j_t++) { //loop through all hits after current hit
                if (!used[j_t]) {
                  float timej   = ecalSingleHit[j_t].getTime();
                  float energyj = ecalSingleHit[j_t].getEnergy();
                  if (m_ecalSimpleTimingCut) {
                    float deltat_ij = m_ecalCorrectTimesForPropagation ? dt : 0;
                    if (timej - deltat_ij > m_ecalTimeWindowMin && timej - deltat_ij < ecalTimeWindowMax) {
                      energySum += energyj;
                      if (timej < timei) {
                        timei = timej;     //use earliest hit time for simpletimingcut
                      }
                    }
                  } else {
                      float deltat_ij = fabs(timei - timej); 
                      //if this subhit is close to current subhit, add this hit's energy to 
                      if (deltat_ij < m_ecalDeltaTimeHitResolution) {
                        if (energyj > energyi) {
                          timei = timej;
                        }
                          energyi += energyj;
                          used[j_t] = true;
                      }
                  }
                }
              }
              if (m_ecalSimpleTimingCut) {
                used = vector<bool>(n, true);  // mark everything as used to terminate for loop on next run
                energyi += energySum;          // fill energySum back into energyi to have rest of loop behave the same.
              }

              // variables and their behaviour at this point:
              // if SimpleTimingCut == false
              // energyi carries the sum of subhit energies within +- one ecal time resolution - the timecluster energy.
              // timei carries something vaguely similar to the central hit time of the merged subhits

              // if SimpleTimingCut == true
              // energyi carries the sum of subhit energies within timeWindowMin and timeWindowMax
              // timei carries the time of the earliest hit within this window

              if (m_digitalEcal) {
                calibr_coeff = this->digitalEcalCalibCoeff(layer);
                if (m_mapsEcalCorrection) {
                  if (caloLayout == CHT::barrel) {
                    float correction = 1.1387 - 0.068 * cost - 0.191 * cost * cost;
                    calibr_coeff /= correction;
                  } else {
                    float correction = 0.592 + 0.590 * cost;
                    calibr_coeff /= correction;
                  }
                }
              } else {
                calibr_coeff = this->analogueEcalCalibCoeff(layer);
              }
              // if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)calibr_coeff *= _ecalEndcapCorrectionFactor;
              if (caloLayout != CHT::barrel) {
                calibr_coeff *= m_ecalEndcapCorrectionFactor;  // more robust
              }

              // fill the hit time - energy histograms
              ++fEcal[{timei, energyi * calibr_coeff}];
              ++fEcalC[{timei - dt, energyi * calibr_coeff}];
              ++fEcalC1[{timei - dt, energyi * calibr_coeff}];
              ++fEcalC2[{timei - dt, energyi * calibr_coeff}];

              // apply extra energy digitisation effects
              energyi = ecalEnergyDigi(energyi, cellID); // this only uses the current subhit "timecluster"!

              if (energyi > m_thresholdEcal) { // now would be the correct time to do threshold comparison
                float timeCor = 0;
                if (m_ecalCorrectTimesForPropagation) {
                  timeCor = dt; 
                }
                timei = timei - timeCor;
                if (timei > m_ecalTimeWindowMin && timei < ecalTimeWindowMax) { // if current subhit timecluster is within specified timing window, 
                                                                                // create new CalorimeterHit and add to collections etc.
                  count++;
                  // CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
                  edm4hep::MutableCalorimeterHit calHit = col.create();
                  if (m_ecalGapCorrection != 0) {
                    m_calHitsByStaveLayer[stave][layer].push_back(&calHit);
                    m_calHitsByStaveLayerModule[stave][layer].push_back(module);
                  }
                  calHit.setCellID(cellID);

                  if (m_digitalEcal) {
                    calHit.setEnergy(calibr_coeff);
                  } else {
                    calHit.setEnergy(calibr_coeff * energyi);
                  }

                  eCellOutput += energyi * calibr_coeff;

                  calHit.setTime(timei);
                  calHit.setPosition(hit.getPosition());
                  calHit.setType(CHT(CHT::em, CHT::ecal, caloLayout, layer));
                  
                  // set relation with CaloHitSimCaloHitLinkCollection
                  auto ecaloRel = Relcol.create();
                  ecaloRel.setFrom(calHit);
                  ecaloRel.setTo(hit);
                  
                } else {
                  // if(caloLayout==CHT::barrel)std::cout << " Drop ECAL Barrel hit : " << timei << " " << calibr_coeff*energyi << std::endl;
                  }
              }
            }
          }
        } 
        else {  // if timing cut is not used
          // CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
          edm4hep::MutableCalorimeterHit calHit = col.create();
          if (m_ecalGapCorrection != 0) {
            m_calHitsByStaveLayer[stave][layer].push_back(&calHit);
            m_calHitsByStaveLayerModule[stave][layer].push_back(module);
          }
          float energyi = hit.getEnergy();

          // apply extra energy digitisation effects
          energyi = ecalEnergyDigi(energyi, cellID);

          calHit.setCellID(cellID);
          if (m_digitalEcal) {
            calHit.setEnergy(calibr_coeff);
          } else {
            calHit.setEnergy(calibr_coeff * energyi);
          }
          calHit.setTime(0);
          calHit.setPosition(hit.getPosition());
          calHit.setType(CHT(CHT::em, CHT::ecal, caloLayout, layer));
          //calHit.setRawHit(hit);

          // set relation with CaloHitSimCaloHitLinkCollection
          auto ecaloRel = Relcol.create();
          ecaloRel.setFrom(calHit);
          ecaloRel.setTo(hit);
        }  // timing if...else end

        //std::cout << hit->getTimeCont(0) << " count = " << count <<  " E ECAL = " << energyCal << " - " << eCellInTime << " - " << eCellOutput << std::endl;
      }  // energy threshold end
    } // hits loop end

    // if requested apply gap corrections in ECAL ?
    if (m_ecalGapCorrection != 0) {
      this->fillECALGaps(m_calHitsByStaveLayer, m_calHitsByStaveLayerModule);
      // add ECAL collection to event
      // ecalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
    }


    // fill normalisation of ECAL occupancy plots
    for (float x = 2.5; x < 3000; x += 5) {
      for (float y = 2.5; y < 3000; y += 5) {
        float r = sqrt(x * x + y * y);
        if (r > 235) {
          ++fEcalRLayerNorm[{r, 4.}];  
        }
      }
    }    

return std::make_tuple(std::move(col), std::move(Relcol)); 
} // end of ECAL digitization

StatusCode DDECalDigi::finalize() {

  //delete randomengines if needed
  if (m_randomEngineDeadCellEcal != 0) {
    delete m_randomEngineDeadCellEcal;
  }

return StatusCode::SUCCESS; 
}


void DDECalDigi::fillECALGaps(std::vector<edm4hep::MutableCalorimeterHit*> m_calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS],
			      std::vector<int> m_calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS]
			      ) const {
  // Loop over hits in the Barrel
  // For each layer calculated differences in hit positions
  // Look for gaps based on expected separation of adjacent hits
  // loop over staves and layers

  for (int is = 0; is < MAX_STAVES; ++is) {
    for (int il = 0; il < MAX_LAYERS; ++il) {
      if (m_calHitsByStaveLayer[is][il].size() > 1) {
        // compare all pairs of hits just once (j>i)

        for (unsigned int i = 0; i < m_calHitsByStaveLayer[is][il].size() - 1; ++i) {
          edm4hep::MutableCalorimeterHit* hiti    = m_calHitsByStaveLayer[is][il][i];
          int                             modulei = m_calHitsByStaveLayerModule[is][il][i];
          float                           xi      = hiti->getPosition()[0];
          float                           yi      = hiti->getPosition()[1];
          float                           zi      = hiti->getPosition()[2];

          for (unsigned int j = i + 1; j < m_calHitsByStaveLayer[is][il].size(); ++j) {
            edm4hep::MutableCalorimeterHit* hitj    = m_calHitsByStaveLayer[is][il][j];
            int                             modulej = m_calHitsByStaveLayerModule[is][il][j];
            float                           xj      = hitj->getPosition()[0];
            float                           yj      = hitj->getPosition()[1];
            float                           zj      = hitj->getPosition()[2];
            float                           dz      = fabs(zi - zj);
            // *** BARREL CORRECTION ***
            if (fabs(zi) < m_zOfEcalEndcap && fabs(zj) < m_zOfEcalEndcap) {
              // account for stave directions using normals
              // calculate difference in hit postions in z and along stave
              float dx = xi - xj;
              float dy = yi - yj;
              //float dt = fabs(dx*_barrelStaveDir[is][0] + dy*_barrelStaveDir[is][1]);
              float dt = sqrt(dx * dx + dy * dy);  // more generic
              // flags for evidence for gaps
              bool zgap  = false;  // in z direction
              bool tgap  = false;  // along stave
              bool ztgap = false;  // in both z and along stave
              bool mgap  = false;  // gaps between ECAL modules

              // criteria gaps in the z and t direction
              float zminm = 1.0 * m_barrelPixelSizeZ[il] - slop;
              float zmin  = 1.0 * m_barrelPixelSizeZ[il] + slop;
              float zmax  = 2.0 * m_barrelPixelSizeZ[il] - slop;
              float tminm = 1.0 * m_barrelPixelSizeT[il] - slop;
              float tmin  = 1.0 * m_barrelPixelSizeT[il] + slop;
              float tmax  = 2.0 * m_barrelPixelSizeT[il] - slop;

              // criteria for gaps
              // WOULD BE BETTER TO USE GEAR TO CHECK GAPS ARE OF EXPECTED SIZE
              if (dz > zmin && dz < zmax && dt < tminm)
                zgap = true;
              if (dz < zminm && dt > tmin && dt < tmax)
                tgap = true;
              if (dz > zmin && dz < zmax && dt > tmin && dt < tmax)
                ztgap = true;

              if (modulei != modulej) {
                if (dz > zmin && dz < 3.0 * m_barrelPixelSizeZ[il] - slop && dt < tmin)
                  mgap = true;
              }

              // found a gap now apply a correction based on area of gap/area of pixel
              if (zgap || tgap || ztgap || mgap) {
                float ecor = 1.;
                float f    = m_ecalGapCorrectionFactor;  // fudge
                if (mgap)
                  f = m_ecalModuleGapCorrectionFactor;
                if (zgap || mgap)
                  ecor = 1. + f * (dz - m_barrelPixelSizeZ[il]) / 2. / m_barrelPixelSizeZ[il];
                if (tgap)
                  ecor = 1. + f * (dt - m_barrelPixelSizeT[il]) / 2. / m_barrelPixelSizeT[il];
                if (ztgap)
                  ecor = 1. + f * (dt - m_barrelPixelSizeT[il]) * (dz - m_barrelPixelSizeZ[il]) / 4. /
                                  m_barrelPixelSizeT[il] / m_barrelPixelSizeZ[il];
                float ei = hiti->getEnergy() * ecor;
                float ej = hitj->getEnergy() * ecor;
                hiti->setEnergy(ei);
                hitj->setEnergy(ej);
              }

              // *** ENDCAP CORRECTION ***
            } else if (fabs(zi) > m_zOfEcalEndcap && fabs(zj) > m_zOfEcalEndcap && dz < 100) {
              float dx    = fabs(xi - xj);
              float dy    = fabs(yi - yj);
              bool  xgap  = false;
              bool  ygap  = false;
              bool  xygap = false;
              // criteria gaps in the z and t direction

              // x and y need to be swapped in different staves of endcap.
              float pixsizex, pixsizey;
              if (is % 2 == 1) {
                pixsizex = m_endcapPixelSizeY[il];
                pixsizey = m_endcapPixelSizeX[il];
              } else {
                pixsizex = m_endcapPixelSizeX[il];
                pixsizey = m_endcapPixelSizeY[il];
              }

              float xmin  = 1.0 * pixsizex + slop;
              float xminm = 1.0 * pixsizex - slop;
              float xmax  = 2.0 * pixsizex - slop;
              float ymin  = 1.0 * pixsizey + slop;
              float yminm = 1.0 * pixsizey - slop;
              float ymax  = 2.0 * pixsizey - slop;
              // look for gaps
              if (dx > xmin && dx < xmax && dy < yminm)
                xgap = true;
              if (dx < xminm && dy > ymin && dy < ymax)
                ygap = true;
              if (dx > xmin && dx < xmax && dy > ymin && dy < ymax)
                xygap = true;

              if (xgap || ygap || xygap) {
                // cout <<"NewLDCCaloDigi found endcap gap, adjusting energy! " << xgap << " " << ygap << " " << xygap << " , " << il << endl;
                // cout << "stave " << is <<  " layer " << il << endl;
                // cout << "  dx, dy " << dx<< " " << dy << " , sizes = " << pixsizex << " " << pixsizey << endl;
                // cout << " xmin... " << xmin << " " << xminm << " " << xmax << " ymin... " << ymin << " " << yminm << " " << ymax << endl;

                // found a gap make correction
                float ecor = 1.;
                float f    = m_ecalGapCorrectionFactor;  // fudge
                if (xgap)
                  ecor = 1. + f * (dx - pixsizex) / 2. / pixsizex;
                if (ygap)
                  ecor = 1. + f * (dy - pixsizey) / 2. / pixsizey;
                if (xygap)
                  ecor = 1. + f * (dx - pixsizex) * (dy - pixsizey) / 4. / pixsizex / pixsizey;

                // cout << "correction factor = " << ecor << endl;

                hiti->setEnergy(hiti->getEnergy() * ecor);
                hitj->setEnergy(hitj->getEnergy() * ecor);
              }
            }
          }
        }
      }
    }
  }

  return;
}

float DDECalDigi::digitalEcalCalibCoeff(int layer) const {
  float calib_coeff = 0.;
  for (unsigned int k(0); k < m_ecalLayers.size(); ++k) {
    int min, max;
    if (k == 0)
      min = 0;
    else
      min = m_ecalLayers[k-1];

    max = m_ecalLayers[k];
    if (layer >= min && layer < max) {
      calib_coeff = m_calibrCoeffEcal[k];
      break;
    }
  }
  return calib_coeff;
}

float DDECalDigi::analogueEcalCalibCoeff(int layer) const {
  float calib_coeff = 0;

  // retrieve calibration constants
  for (unsigned int k(0); k < m_ecalLayers.size(); ++k) {
    int min, max;
    if (k == 0) {
      min = 0;
    } else {
      min = m_ecalLayers[k-1];
    }
    max = m_ecalLayers[k];
    if (layer >= min && layer < max) {
      calib_coeff = m_calibrCoeffEcal[k];
      break;
    }
  }
  return calib_coeff;
}

float DDECalDigi::ecalEnergyDigi(float energy, int id) const {
  // some extra digi effects (daniel)
  // controlled by _applyEcalDigi = 0 (none), 1 (silicon), 2 (scintillator)

  // small update for time-constant uncorrelated miscalibrations. DJ, Jan 2015

  float e_out = energy;
  if (m_applyEcalDigi == 1) {
    e_out = DDECalDigi::siliconDigi(energy);  // silicon digi
  } else if (m_applyEcalDigi == 2) {
    e_out = DDECalDigi::scintillatorDigi(energy, true);  // scintillator digi
  }
  // add electronics dynamic range
  // Sept 2015: Daniel moved this to the ScintillatorDigi part, so it is applied before unfolding of sipm response
  // if (_ecalMaxDynMip>0) e_out = min (e_out, _ecalMaxDynMip*_calibEcalMip);

  // random miscalib
  if (m_misCalibEcal_uncorrel > 0) {
    float miscal(0);
    if (m_misCalibEcal_uncorrel_keep) {
      int id{0};
      if (m_ECAL_cell_miscalibs.find(id) !=
          m_ECAL_cell_miscalibs.end()) {  // this cell was previously seen, and a miscalib stored
        miscal = m_ECAL_cell_miscalibs.at(id);
      } else {  // we haven't seen this one yet, get a miscalib for it
        miscal = CLHEP::RandGauss::shoot(1.0, m_misCalibEcal_uncorrel);
	// FIXME: this is storing miscalibration globally for a run???
        // FIXME _ECAL_cell_miscalibs[id] = miscal;
      }
    } else {
      miscal = CLHEP::RandGauss::shoot(1.0, m_misCalibEcal_uncorrel);
    }
    e_out *= miscal;
  }

  if (m_misCalibEcal_correl > 0)
    e_out *= m_event_correl_miscalib_ecal;

  // random cell kill
  if (m_deadCellFractionEcal > 0) {
    if (m_deadCellEcal_keep == true) {
      int id;

      if (m_ECAL_cell_dead.find(id) != m_ECAL_cell_dead.end()) {  // this cell was previously seen
        if (m_ECAL_cell_dead.at(id) == true) {
          e_out = 0;
        }
      } else {  // we haven't seen this one yet, get a miscalib for it
        bool thisDead       = (CLHEP::RandFlat::shoot(m_randomEngineDeadCellEcal, .0, 1.0) < m_deadCellFractionEcal);
	// FIXME global map
        //_ECAL_cell_dead[id] = thisDead;
        if (thisDead == true) {
          e_out = 0;
        }
      }

    } else {
      if (CLHEP::RandFlat::shoot(0.0, 1.0) < m_deadCellFractionEcal)
        e_out = 0;
    }
  }

  return e_out;
}

float DDECalDigi::siliconDigi(float energy) const {
  // applies extra digitisation to silicon hits

  // calculate #e-h pairs
  float nehpairs = 1.0e9 * energy / m_ehEnergy;  // check units of energy! m_ehEnergy is in eV, energy in GeV

  // fluctuate it by Poisson
  float smeared_energy = energy * CLHEP::RandPoisson::shoot(nehpairs) / nehpairs;

  // limited electronics dynamic range // Daniel moved electronics dyn range to here
  if (m_ecalMaxDynMip > 0)
    smeared_energy = std::min(smeared_energy, m_ecalMaxDynMip * m_calibEcalMip);

  // add electronics noise
  if (m_ecal_elec_noise > 0)
    smeared_energy += CLHEP::RandGauss::shoot(0, m_ecal_elec_noise * m_calibEcalMip);

  return smeared_energy;
}

//FIXME should return a reference to somewhere else? or do this every event?
std::vector<std::pair<int, int>> DDECalDigi::getLayerConfig() const {
  // get the layer layout (silicon, scintillator)
  // first element of layerTypes is the preshower
  std::vector < std::pair <int, int> > m_layerTypes {};
  if (m_layerTypes.size() == 0) {
    for (std::string::size_type i = 0; i < m_ecalLayout.size(); ++i) {
      // convert each element of string to integer
      // int type = std::atoi( &ccdd ); // this is not well done (must be null-terminated)
      int etype = m_ecalLayout[i] - '0';  // this is less obvious, but works...

      switch (etype) {  // these originally defined in Mokka driver SEcalSD04
        case 0:
          m_layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          m_layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          break;
        case 1:
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          break;
        case 2:
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          break;
        case 3:
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          break;
        case 4:
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          break;
        case 5:
          m_layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          break;
        case 6:
          m_layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          break;
        case 7:
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          m_layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          break;
        case 8:
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          m_layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          break;
        default:
          error() << "ERROR, unknown layer type " << etype << endl;
      }
    }
  }
  return m_layerTypes;
}

std::pair<int, int> DDECalDigi::getLayerProperties(std::string const& colName, int layer) const {
  std::pair<int, int> thislayersetup(-99, -99);
  std::string         colNameLow(colName);
  std::transform(colNameLow.begin(), colNameLow.end(), colNameLow.begin(), ::tolower);
  if (colNameLow.find("presh") != string::npos) {  // preshower
    if (layer != 0) {
      warning() << "preshower layer with layer index = " << layer << " ??? " << endmsg;
    } else {
      thislayersetup = getLayerConfig()[layer];
    }
  } else if (colNameLow.find("ring") !=
             string::npos) {  // ecal ring (has no preshower), and is actually always all silicon
    if (layer < int(getLayerConfig().size())) {
      // thislayersetup = getLayerConfig()[layer];
      thislayersetup = std::pair<int, int>(SIECAL, SQUARE);
    } else {
      warning() << "unphysical layer number? " << layer << " " << getLayerConfig().size() << endmsg;
    }
  } else {  // endcap, barrel
    if (layer + 1 < int(getLayerConfig().size())) {
      thislayersetup = getLayerConfig()[layer + 1];
    } else {
      warning() << "unphysical layer number? " << layer << " " << getLayerConfig().size() << endmsg;
    }
  }
  return thislayersetup;
}

void DDECalDigi::checkConsistency(std::string colName, int layer) const{
  if (m_applyEcalDigi == 0 || m_countWarnings > 20)
    return;

  std::pair<int, int> thislayersetup = getLayerProperties(colName, layer);

  if (m_applyEcalDigi == 1 && thislayersetup.first != SIECAL) {
    error() << "collection: " << colName << endmsg;
    error() << "You seem to be trying to apply ECAL silicon digitisation to scintillator? Refusing!" << endmsg;
    error() << "Check setting of ECAL_apply_realistic_digi: " << m_applyEcalDigi << endmsg;
    assert(0);
    m_countWarnings++;
  }

  if (m_applyEcalDigi == 2 && thislayersetup.first != SCECAL) {
    error() << "collection: " << colName << endmsg;
    error() << "You seem to be trying to apply ECAL scintillator digitisation to silicon? Refusing!" << endmsg;
    error() << "Check setting of ECAL_apply_realistic_digi: " << m_applyEcalDigi << endmsg;
    assert(0);
    m_countWarnings++;
  }

  if (thislayersetup.second != getStripOrientationFromColName(colName)) {
    error() << "collection: " << colName << endmsg;
    error() << "Some inconsistency in strip orientation?" << endmsg;
    error() << " from collection name: " << getStripOrientationFromColName(colName) << endmsg;
    error() << " from layer config string: " << thislayersetup.second << endmsg;
    m_countWarnings++;
  }
  return;
}
