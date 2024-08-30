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

#include "DDHCalDigi.h"

// protect against rounding errors
// will not find caps smaller than this
const float slop = 0.25;  // (mm)
//const float pi = acos(-1.0); ///FIXME
//const float twopi = 2.0*pi;  ///FIXME: DD4HEP INTERFERES WITH THESE


DDHCalDigi::DDHCalDigi(const std::string& aName, ISvcLocator* aSvcLoc)
    : DDBaseCaloDigi(aName, aSvcLoc,
          {
            KeyValues("InputCaloHitCollection", {"HCalBarrelCollection"}),
            KeyValues("HeaderName", {"EventHeader"}),
          },
          {
            KeyValues("OutputCaloHitCollection", {"HCalorimeterHit1"}), //, "HCalorimeterHit2", "HCalorimeterHit3"}),
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

StatusCode DDHCalDigi::initialize() {

  m_strip_virt_cells = -999;
  m_countWarnings    = 0;

  // Convert HCAL thresholds to GeV units
  if (m_unitThresholdHcal.value().compare("GeV") == 0) {
    // HCAL threshold unit is GeV, do nothing
  } else if (m_unitThresholdHcal.value().compare("MIP") == 0) {
    // HCAL threshold unit is MIP, convert via MIP2GeV
    for (unsigned int i = 0; i < m_thresholdHcal.value().size(); i++) {
      m_thresholdHcal.value()[i] *= m_calibHcalMip.value();
    }
  } else if (m_unitThresholdHcal.value().compare("px") == 0) {
    // HCAL threshold unit is pixels, convert via MIP2GeV and lightyield
    for (unsigned int i = 0; i < m_thresholdHcal.size(); i++) {
      m_thresholdHcal[i] *= m_hcal_PPD_pe_per_mip.value() * m_calibHcalMip.value();
    }
  } else {
    error() << "Could not identify HCAL threshold unit. Please use \"GeV\", \"MIP\" or \"px\"! Aborting." << endmsg;
    assert(0);
  }

  m_scHcalDigi = std::unique_ptr<DDScintillatorPpdDigi>(new DDScintillatorPpdDigi());
  m_scHcalDigi->setPEperMIP(m_hcal_PPD_pe_per_mip);
  m_scHcalDigi->setCalibMIP(m_calibHcalMip);
  m_scHcalDigi->setNPix(m_hcal_PPD_n_pixels);
  m_scHcalDigi->setRandomMisCalibNPix(m_hcal_misCalibNpix);
  m_scHcalDigi->setPixSpread(m_hcal_pixSpread);
  m_scHcalDigi->setElecNoise(m_hcal_elec_noise);
  m_scHcalDigi->setElecRange(m_hcalMaxDynMip);
  cout << "HCAL sc digi:" << endl;
  m_scHcalDigi->printParameters();

  // Set up the random engines for HCAL dead cells: (could use a steering parameter though)
  if (m_deadCellHcal_keep) {
    m_randomEngineDeadCellHcal = new CLHEP::MTwistEngine(0, 0);
  } else {
    m_randomEngineDeadCellHcal = 0;
  }

  return StatusCode::SUCCESS;
}

retType DDHCalDigi::operator()(
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
  //_event_correl_miscalib_hcal = CLHEP::RandGauss::shoot( 1.0, _misCalibHcal_correl );


  std::string colName = m_inputLocations[0][0].key();       // take input collection name
  cout << "looking for collection: " << colName << endl;

  if (colName.find("dummy") != string::npos) {
    debug() << "Ignoring input HCAL collection name (looks like dummy name)" << colName << endmsg;
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
  // * Reading HCAL Collections of Simulated Hits *
  //
    for (const auto& hit : simCaloHitCollection) {
      float energy = hit.getEnergy();
      //preselect for SimHits with energySum>threshold. Doubtful at least, as lower energy hit might fluctuate up and still be counted
      if (energy > m_thresholdHcal[0] / 2) {
        int   cellID = hit.getCellID();
        float calibr_coeff = 1.0;
        unsigned int layer  = bitFieldCoder.get(cellID, "layer");
          
        // NOTE : for a digital HCAL this does not allow for varying layer thickness
        // with depth - would need a simple mod to make response proportional to layer thickness
        if (m_digitalHcal) {
          calibr_coeff = this->digitalHcalCalibCoeff(caloLayout, energy);
        } else {
          calibr_coeff = this->analogueHcalCalibCoeff(caloLayout, layer);
        }
        // if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)calibr_coeff*=_hcalEndcapCorrectionFactor;
        if (caloLayout != CHT::endcap)
          calibr_coeff *= m_hcalEndcapCorrectionFactor;  // more robust, is applied to ALL hits outside of barrel.

        //float energyCal = energy*calibr_coeff
        float x = hit.getPosition()[0];
        float y = hit.getPosition()[1];
        float z = hit.getPosition()[2];
        //float r = sqrt(x*x+y*y);
        if (m_useHcalTiming) {
          float hcalTimeWindowMax;
          if (caloLayout == CHT::barrel) {                //current SimHit is in barrel, use barrel timing cut
            hcalTimeWindowMax = m_hcalBarrelTimeWindowMax;
          } else {                                        //current SimHit is not in barrel, use endcap timing cut
            hcalTimeWindowMax = m_hcalEndcapTimeWindowMax;
          }

          float r = sqrt(x * x + y * y + z * z);  //this is a crude approximation. assumes initial particle originated at the very center of the detector.
          float dt         = r / 300 - 0.1;  //magic numbers! ~
          auto  hcalSingleHit = hit.getContributions();
          const unsigned int n          = hcalSingleHit.size();

          std::vector<bool> used(n, false);

          int count = 0;

          for (unsigned int i_t = 0; i_t < n; i_t++) {      // loop over all subhits
            float timei     = hcalSingleHit[i_t].getTime();    //absolute hit timing of current subhit
            float energyi   = hcalSingleHit[i_t].getEnergy();  //energy of current subhit
            float energySum = 0;
            //float deltat = 0;
            //if(_hcalCorrectTimesForPropagation)deltat=dt;
            //deltat now carries hit timing correction.
            //std::cout <<"outer:" << i << " " << n << std::endl;

            //idea of the following section:
            //if simpletimingcut == false
            //sum up hit energies which lie within one calo timing resolution to "timecluster" of current subhit
            //then treat each single timecluster as one hit over threshold and digitise separately. this means there can be more than one CalorimeterHit with the same cellIDs, but different hit times (!)
            //
            //if simpletimingcut == true
            //i'm very sorry. this is the worst code you will ever see.
            //sum up hit energies within timeWindowMin and timeWindowMax, use earliest subhit in this window as hit time for resulting calohit.
            //only one calorimeterhit will be generated from this.

            if (!used[i_t]) {  //if current subhit has not been merged with previous hits already, take current hit as starting point to merge hits
              // merge with other hits?
              used[i_t] = true;
              for (unsigned int j_t = i_t; j_t < n; j_t++) {  //loop through all hits after current hit
                //std::cout << "inner:" << i << " " << j << " " << n << std::endl;
                if (!used[j_t]) {
                  float timej   = hcalSingleHit[j_t].getTime();
                  float energyj = hcalSingleHit[j_t].getEnergy();
                  //              std::cout << " HCAL  deltat_ij : " << deltat_ij << std::endl;
                  if (m_hcalSimpleTimingCut) {
                    float deltat_ij = m_hcalCorrectTimesForPropagation ? dt : 0;
                    if (timej - deltat_ij > m_hcalTimeWindowMin.value() && timej - deltat_ij < hcalTimeWindowMax) {
                        energySum += energyj;
                        if (timej < timei) {
                          timei = timej;  //use earliest hit time for simpletimingcut
                        }
                      }
                    } else {
                      float deltat_ij = fabs(timei - timej);
                      //if this subhit is close to current subhit, add this hit's energy to timecluster
                      if (deltat_ij < m_hcalDeltaTimeHitResolution) {
                        if (energyj > energyi) {
                          //this is probably not what was intended. i guess this should find the largest hit of one timecluster and use its hittime for the cluster, but instead it compares the current hit energy to the sum of already found hit energies
                          timei = timej;
                        }
                        //std::cout << timei << " - " << timej << std::endl;
                        //std::cout << energyi << " - " << energyj << std::endl;
                        energyi += energyj;
                        used[j_t] = true;
                        //std::cout << timei << " " << energyi << std::endl;
                    }
                  }
                }
              }
              if (m_hcalSimpleTimingCut) {
                used = vector<bool>(n, true);  //mark everything as used to terminate loop. this is worse than goto. i'm sorry.
                energyi += energySum;  //fill energySum back into energyi to have rest of loop behave the same.
              }

              //variables and their behaviour at this point:
              //if SimpleTimingCut == false
              //energyi carries the sum of subhit energies within +- one hcal time resolution - the timecluster energy.
              //timei carries something vaguely similar to the central hit time of the merged subhits

              //if SimpleTimingCut == true
              //energyi carries the sum of subhit energies within timeWindowMin and timeWindowMax
              //timei carries the time of the earliest hit within this window

              // apply extra energy digitisation effects
              energyi = hcalEnergyDigi(energyi, cellID);  //this only uses the current subhit "timecluster"!

              if (energyi > m_thresholdHcal[0]) {  //now would be the correct time to do threshold comparison
                float timeCor = 0;
                if (m_hcalCorrectTimesForPropagation)
                  timeCor = dt;
                timei = timei - timeCor;
                if (timei > m_hcalTimeWindowMin &&
                    timei < hcalTimeWindowMax) {  //if current subhit timecluster is within specified timing window, create new CalorimeterHit and add to collections etc.
                  count++;
                  edm4hep::MutableCalorimeterHit calHit = col.create();
                  calHit.setCellID(cellID);
                  if (m_digitalHcal) {
                    calHit.setEnergy(calibr_coeff);
                    //eCellOutput+= calibr_coeff;
                  } else {
                    calHit.setEnergy(calibr_coeff * energyi);
                    //eCellOutput+= energyi*calibr_coeff;
                  }
                  calHit.setTime(timei);
                  calHit.setPosition(hit.getPosition());
                  calHit.setType(CHT(CHT::had, CHT::hcal, caloLayout, layer));
                  // Set relation with CaloHitSimCaloHitLinkCollection
                  auto hcaloRel = Relcol.create();
                  hcaloRel.setFrom(calHit);
                  hcaloRel.setTo(hit);
                } else {
                  //std::cout << "Drop HCAL hit : " << timei << " " << calibr_coeff*energyi << std::endl;
                }
              }
            }
          } // end loop over all subhits
        } 
        else {  // if timing cuo is not used
                // CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
          edm4hep::MutableCalorimeterHit calHit = col.create();
          calHit.setCellID(cellID);
          float energyi = hit.getEnergy();

          // apply realistic digitisation
          energyi = hcalEnergyDigi(energyi, cellID);

          if (m_digitalHcal) {
            calHit.setEnergy(calibr_coeff);
          } else {
            calHit.setEnergy(calibr_coeff * energyi);
          }
          //eCellOutput+= energyi*calibr_coeff;
          calHit.setTime(0);
          calHit.setPosition(hit.getPosition());
          calHit.setType(CHT(CHT::had, CHT::hcal, caloLayout, layer));
          //calHit.setRawHit(hit);

          auto hcaloRel = Relcol.create();
          hcaloRel.setFrom(calHit);
          hcaloRel.setTo(hit);
          }

          // std::cout << hit->getTimeCont(0) << " count = " << count <<  " EHCAL = " << energyCal << " - " << eCellInTime << " - " << eCellOutput << std::endl;
        }
      }
      // add HCAL collection to event
      // hcalcolparameters().setValue(LCIO::CellIDEncoding,initString);
      // evt->addCollection(hcalcol,_outputHcalCollections[i].c_str());

  // Create and add relation collection for ECAL/HCAL to event
  //chschcol = calohitNav.createLCCollection();
  // evt->addCollection(chschcol,_outputRelCollection.c_str());

  // _nEvt++;



return std::make_tuple(std::move(col), std::move(Relcol)); 
} // end of HCAL digitization


StatusCode DDHCalDigi::finalize() {
  
  //delete randomengines if needed
  if (m_randomEngineDeadCellHcal != 0) {
    delete m_randomEngineDeadCellHcal;
  }

return StatusCode::SUCCESS; 
}

float DDHCalDigi::digitalHcalCalibCoeff(CHT::Layout caloLayout, float energy) const {
  float        calib_coeff = 0;
  unsigned int ilevel      = 0;
  for (unsigned int ithresh = 1; ithresh < m_thresholdHcal.size(); ithresh++) {
    // Assume!!!  hit energies are stored as floats, i.e. 1, 2 or 3
    if (energy > m_thresholdHcal[ithresh])
      ilevel = ithresh;  // ilevel = 0 , 1, 2
  }

  switch (caloLayout) {
    case CHT::barrel:
      if (ilevel > m_calibrCoeffHcalBarrel.value().size() - 1) {
        error() << " Semi-digital level " << ilevel
		<< " greater than number of HCAL Calibration Constants (" << m_calibrCoeffHcalBarrel.value().size()
		<< ")" << std::endl;
      } else {
        calib_coeff = m_calibrCoeffHcalBarrel.value()[ilevel];
      }
      break;
    case CHT::endcap:
      if (ilevel > m_calibrCoeffHcalEndcap.value().size() - 1) {
        error() << " Semi-digital level " << ilevel
		<< " greater than number of HCAL Calibration Constants (" << m_calibrCoeffHcalEndcap.value().size()
		<< ")" << std::endl;
      } else {
        calib_coeff = m_calibrCoeffHcalEndcap.value()[ilevel];
      }
      break;
    case CHT::plug:
      if (ilevel > m_calibrCoeffHcalOther.value().size() - 1) {
	error() << " Semi-digital level " << ilevel
		<< " greater than number of HCAL Calibration Constants (" << m_calibrCoeffHcalOther.value().size()
		<< ")" << std::endl;
      } else {
        calib_coeff = m_calibrCoeffHcalOther.value()[ilevel];
      }
      break;
    default:
    error() << " Unknown HCAL Hit Type " << std::endl;
      break;
  }

  return calib_coeff;
}

float DDHCalDigi::analogueHcalCalibCoeff(CHT::Layout caloLayout, int layer) const {
  float calib_coeff = 0;

  for (unsigned int k(0); k < m_hcalLayers.size(); ++k) {
    int min, max;
    if (k == 0)
      min = 0;
    else
      min = m_hcalLayers[k - 1];

    max = m_hcalLayers[k];
    if (layer >= min && layer < max) {
      switch (caloLayout) {
        case CHT::barrel:
          calib_coeff = m_calibrCoeffHcalBarrel[k];
          break;
        case CHT::endcap:
          calib_coeff = m_calibrCoeffHcalEndcap[k];
          break;
        case CHT::plug:
        case CHT::ring:
          calib_coeff = m_calibrCoeffHcalOther[k];
          break;
        default:
          error() << " Unknown HCAL Hit Type " << std::endl;
          break;
      }
    }
  }

  return calib_coeff;
}

float DDHCalDigi::hcalEnergyDigi(float energy, int id) const {
  // some extra digi effects (daniel)
  // controlled by _applyHcalDigi = 0 (none), 1 (scintillator/SiPM)

  // small update for time-constant uncorrelated miscalibrations. DJ, Jan 2015

  float e_out(energy);
  if (m_applyHcalDigi == 1)
    e_out = scintillatorDigi(energy, false);  // scintillator digi

  // add electronics dynamic range
  // Sept 2015: Daniel moved this to the ScintillatorDigi part, so it is applied before unfolding of sipm response
  //  if (_hcalMaxDynMip>0) e_out = min (e_out, _hcalMaxDynMip*_calibHcalMip);

  // random miscalib
  //  if (_misCalibHcal_uncorrel>0) e_out*=CLHEP::RandGauss::shoot( 1.0, _misCalibHcal_uncorrel );
  if (m_misCalibHcal_uncorrel > 0) {
    float miscal(0);
    if (m_misCalibHcal_uncorrel_keep) {
      int id;
      if (m_HCAL_cell_miscalibs.find(id) !=
          m_HCAL_cell_miscalibs.end()) {  // this cell was previously seen, and a miscalib stored
        miscal = m_HCAL_cell_miscalibs.at(id);
      } else {  // we haven't seen this one yet, get a miscalib for it
        miscal                   = CLHEP::RandGauss::shoot(1.0, m_misCalibHcal_uncorrel);
	// FIXME: same as above
        //_HCAL_cell_miscalibs[id] = miscal;
      }
    } else {
      miscal = CLHEP::RandGauss::shoot(1.0, m_misCalibHcal_uncorrel);
    }
    e_out *= miscal;
  }

  if (m_misCalibHcal_correl > 0)
    e_out *= m_event_correl_miscalib_hcal;

  // random cell kill
  if (m_deadCellFractionHcal > 0) {
    if (m_deadCellHcal_keep == true) {
      int id;

      if (m_HCAL_cell_dead.find(id) != m_HCAL_cell_dead.end()) {  // this cell was previously seen
        if (m_HCAL_cell_dead.at(id) == true) {
          e_out = 0;
        }
      } else {  // we haven't seen this one yet, get a miscalib for it
        bool thisDead       = (CLHEP::RandFlat::shoot(m_randomEngineDeadCellHcal, .0, 1.0) < m_deadCellFractionHcal);
	// FIXME globally dead cell map?
        //FIXME _HCAL_cell_dead[id] = thisDead;
        if (thisDead == true) {
          e_out = 0;
        }
      }

    } else {
      if (CLHEP::RandFlat::shoot(0.0, 1.0) < m_deadCellFractionHcal)
        e_out = 0;
    }
  }
  return e_out;
}
