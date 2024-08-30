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

#include "DDBaseCaloDigi.h"

// protect against rounding errors
// will not find caps smaller than this
const float slop = 0.25;  // (mm)
//const float pi = acos(-1.0); ///FIXME
//const float twopi = 2.0*pi;  ///FIXME: DD4HEP INTERFERES WITH THESE



float DDBaseCaloDigi::scintillatorDigi(float energy, bool isEcal) const {
  // this applies some extra digitisation to scintillator+PPD hits (PPD=SiPM, MPPC)
  // - poisson fluctuates the number of photo-electrons according to #PEs/MIP
  // - applies PPD saturation according to #pixels
  // Daniel Jeans, Jan/Feb 2014.

  float digiEn(0);
  if (isEcal) {
    digiEn = m_scEcalDigi->getDigitisedEnergy(energy);  //CHECK!!!
  } else {
    digiEn = m_scHcalDigi->getDigitisedEnergy(energy);
  }
  return digiEn;
}

int DDBaseCaloDigi::getNumberOfVirtualCells() const {
  // if ( _strip_virt_cells < 0 ) { // not yet set, try to get from gear file
  //   // number of virtual cells in scintillator strip
  //   try {
  //     const gear::GearParameters& pMokka = Global::GEAR->getGearParameters("MokkaParameters");
  //     std::string nVirtualMokkaS = pMokka.getStringVal("Ecal_Sc_number_of_virtual_cells");
  //     _strip_virt_cells = std::atoi( nVirtualMokkaS.c_str() );
  //     streamlog_out (MESSAGE) << "INFO: got number of virtual cells from gear file: " << _strip_virt_cells << endl;
  //   } catch(gear::UnknownParameterException &e) { // not found in gear file, get from steering file parameter...
  //     streamlog_out (MESSAGE) << "INFO: taking number of virtual cells from steering file (not found in gear file): " << _ecalStrip_default_nVirt << endl;
  //     _strip_virt_cells = _ecalStrip_default_nVirt;
  //   }
  // }
  return m_strip_virt_cells;
}

int DDBaseCaloDigi::getStripOrientationFromColName(std::string const& colName) const {
  int orientation = -99;
  std::string colNameLow = colName;
  std::transform(colNameLow.begin(), colNameLow.end(), colNameLow.begin(), ::tolower);
  if (colNameLow.find("trans") != std::string::npos) {
    orientation = STRIP_ALIGN_ACROSS_SLAB;
  } else if (colNameLow.find("long") != std::string::npos) {
    orientation = STRIP_ALIGN_ALONG_SLAB;
  } else {  // assume square...
    orientation = SQUARE;
    std::cout << "WARNING, cannot guess strip orientation! for collection " << colName << std::endl;
  }
  return orientation;
}

// edm4hep::SimCalorimeterHitCollection DDBaseCaloDigi::combineVirtualStripCells(edm4hep::SimCalorimeterHitCollection const& col,
//       bool isBarrel, int stripOrientation) const {
//   // combines the virtual cells in a strip
//   // input collection is for virtual cells
//   // returns collection of strips

//   // daniel jeans, jan/feb 2014

//   // sanity check
//   if (stripOrientation == SQUARE) {
//     debug()  << "DDBaseCaloDigi::combineVirtualStripCells trying to deal with silicon strips??? I do not know "
//                             "how to do that, refusing to do anything!"
//                          << std::endl;
//     return 0;
//   }

//   //CellIDDecoder<SimCalorimeterHit> idDecoder( col );
//   // This variable is for tracking detectors only, not for calorimeters
//   // std::string initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());

//   // get the encoding string from the original collection
//   const auto initStringStrip = cellIDHandle.get();
//   dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initStringStrip);

//   // make the new output collection
//   //edm4hep::MutableSimCalorimeterHit tempstripcol = inputCaloHitCollection.create();;
//   auto stripcol = edm4hep::SimCalorimeterHitCollection();
//   //FIXME stripcol.setFlag(_flag.getFlag());
//   //FIXME: add initString/ cellIDEncoding to stripcol

//   // link between strip CellID and hits storage for strip hits
//   std::map<int, edm4hep::MutableSimCalorimeterHit*> newhits;  //SORT

//   // loop over input collection
//   int numElements = col.size();

//   // // sum energy for check
//   float tempenergysum(0);
//   for (const auto& hit : col) {
//     tempenergysum += hit.getEnergy();
//   }

//   float scTVirtLengthBar(-99);
//   float scLVirtLengthBar(-99);
//   float scTVirtLengthEnd(-99);
//   float scLVirtLengthEnd(-99);

//   for (const auto& hit : col) {  // loop over virtual cells

//     // for now, ignore preshower layer.
//     //   in "usual" (non-preshower) ecal collections, km1 = 0 is the first layer AFTER the preshower
//     int km1 = bitFieldCoder(hit)["layer"];  // FIXME

//     int gearlayer = km1 + 1;  // in "gearlayer", preshower is 0, first layer after absorber is 1

//     // after fix of MOKKA, this should be not needed
//     // int gearlayer(-99); // this is the layer index which the GEAR geometry knows about. depends on driver version...
//     // if ( _ecal_driver_version == 0 ) { // SEcal04 driver
//     // } else if ( _ecal_driver_version == 1 ) { // SEcal05 driver
//     //   gearlayer = km1+1;
//     // } else {
//     //   streamlog_out ( ERROR ) << "ERROR unknown value for ecal driver version, _ecal_driver_version=" << _ecal_driver_version << endl;
//     // }
//     // assert (gearlayer>=0);

//     int m       = bitFieldCoder(hit)["module"];
//     int s       = bitFieldCoder(hit)["stave"];
//     int i_index = bitFieldCoder(hit)["x"];
//     int j_index = bitFieldCoder(hit)["y"];

//     //cout << "ORIG: " << km1 << " " << m << " " << s << " " << i_index << " " << j_index << " : " <<
//     //  hit->getPosition()[0] << " "<< hit->getPosition()[1] << " "<< hit->getPosition()[2] << endl;

//     // should we combine virtual cells in the I or J direction?
//     //
//     // in barrel:
//     //   i is along slab (phi in barrel)
//     //   j is across slab (z in barrel)
//     //   increasing j = increasing z; increasing i = decreasing phi
//     //
//     // in endcap:
//     //     i is across slab
//     //     j is along slab
//     //     different staves are rotated around beam axis.
//     //        in +ve z endcap (m==6),
//     //          stave 0 is at -ve x and +ve y. in this stave, i is parallel to -x, j to +y
//     //          stave 3 is at +ve x and +ve y. in this stave, i is parallel to +y, j to +x
//     //        -ve z endcap (m==0) is same, just everything flipped in x.

//     bool collateAlongI =
//         isBarrel ? stripOrientation == STRIP_ALIGN_ALONG_SLAB :  // I is the index along the slab (R-phi in barrel)
//             stripOrientation == STRIP_ALIGN_ACROSS_SLAB;         // for some reason seems to be reversed in endcap...!!

//     debug() << "isBarrel = " << isBarrel << " : stripOrientation= " << stripOrientation
//                          << " gear layer = " << gearlayer << endl;

//     // let's get the length of the virtual cell in the direction of the strip
//     //    fixed rare crashes, and streamlined code to get virtualCellSizeAlongStrip. Sept 2014
//     float virtualCellSizeAlongStrip(0);

//     float* layerpixsize(0);
//     float* savedPixSize(0);
//     if (isBarrel) {
//       if (stripOrientation == STRIP_ALIGN_ACROSS_SLAB) {
//         layerpixsize = _barrelPixelSizeT;
//         savedPixSize = &scTVirtLengthBar;
//       } else {
//         layerpixsize = _barrelPixelSizeZ;
//         savedPixSize = &scLVirtLengthBar;
//       }
//     } else {
//       if (stripOrientation == STRIP_ALIGN_ACROSS_SLAB) {
//         layerpixsize = _endcapPixelSizeX;
//         savedPixSize = &scTVirtLengthEnd;
//       } else {
//         layerpixsize = _endcapPixelSizeY;
//         savedPixSize = &scLVirtLengthEnd;
//       }
//     }

//     virtualCellSizeAlongStrip = layerpixsize[gearlayer];
//     if (virtualCellSizeAlongStrip > 0) {  // looks OK
//       if (*savedPixSize < 0) {            // no saved size yet, keep this one
//         *savedPixSize = virtualCellSizeAlongStrip;
//       }
//     } else {                    // no valid info in gear file for this layer
//       if (*savedPixSize > 0) {  // use saved value from other layer
//         virtualCellSizeAlongStrip = *savedPixSize;
//       } else {  // look at previous layers for one with same orientation (extra check, sept 2014)
//         std::vector<std::pair<int, int>> layerTypes = getLayerConfig();

//         debug() << "could not get valid info from gear file..." << endl
//                              << "looking through previous layers to get a matching orientation" << endl
//                              << "this gear layer " << gearlayer << " type: " << layerTypes[gearlayer].first << " "
//                              << layerTypes[gearlayer].second << endl;

//         for (int il = gearlayer - 1; il >= 0; il--) {
//           // layer types include the preshower at posn "0"
//           // gearlayer has preshower as layer 0
//           if (layerTypes[il] == layerTypes[gearlayer]) {  // found layer with same setup
//             debug() << "found a match! " << il << " " << layerTypes[il].first << " "
//                                  << layerTypes[il].second << " : " << layerpixsize[il] << endl;
//             virtualCellSizeAlongStrip = layerpixsize[il];
//             *savedPixSize             = virtualCellSizeAlongStrip;
//             break;
//           }
//         }
//       }
//     }
//     debug() << "virtualCellSizeAlongStrip = " << virtualCellSizeAlongStrip << endl;

//     // calculate the new strip's I,J coordinates
//     int i_new = collateAlongI ? i_index / getNumberOfVirtualCells() : i_index;
//     int j_new = collateAlongI ? j_index : j_index / getNumberOfVirtualCells();

//     // apply position-dependent energy correction
//     //
//     // relative virtual cell position within strip (0 = centre)
//     //    distance between centre of virtual cell and centre of strip
//     int   relativeIndx = collateAlongI ? i_index : j_index;
//     float relativePos =
//         relativeIndx % getNumberOfVirtualCells();  // this is index within strip wrt strip end (0->_strip_virt_cells-1)
//     relativePos -= (getNumberOfVirtualCells() - 1) / 2.0;  // index wrt strip centre
//     relativePos *= virtualCellSizeAlongStrip;              // distance wrt strip centre

//     // effect of absorbtion length within scintillator
//     //     TODO: should check the polarity is consistent with mppc position, to make sure larger response nearer to mppc....
//     float energy_new = hit.getEnergy();

//     float energyNonuniformityScaling(1.);
//     if (_strip_abs_length > 0) {
//       energyNonuniformityScaling = exp(-relativePos / _strip_abs_length);
//       energy_new *= energyNonuniformityScaling;
//     }

//     // calculate ID of the new Strip // FIXME???
//     bitFieldCoder["system"] = m;  // FIXME ????
//     bitFieldCoder["module"] = m;
//     bitFieldCoder["stave"]  = s;
//     bitFieldCoder["layer"]  = km1;
//     bitFieldCoder["x"]      = i_new;
//     bitFieldCoder["y"]      = j_new;
//     int cellID              = bitFieldCoder.getCellID;  // FIXME

//     // calculate position of centre of this new strip
//     // it depends if we are in the barrel or endcap, and in which stave
//     //   (n.b. we could do this later, after checking for duplicates. however, keeping it here allows us to
//     //            check that we haven't supplied wrong nVirtualCells parameter)

//     float stripPos[3] = {0};
//     if (isBarrel) {
//       if (stripOrientation == STRIP_ALIGN_ACROSS_SLAB) {  // it's along z
//         stripPos[0] = hit.getPosition()[0];
//         stripPos[1] = hit.getPosition()[1];
//         stripPos[2] = hit.getPosition()[2] - relativePos;  // increasing j = increasing z
//       } else {                                             // it's along t
//         float phi   = atan2(_barrelStaveDir[s][1], _barrelStaveDir[s][0]);
//         stripPos[0] = hit.getPosition()[0] - relativePos * cos(phi);  // negative: increasing I = decreasing phi
//         stripPos[1] = hit.getPosition()[1] - relativePos * sin(phi);
//         stripPos[2] = hit.getPosition()[2];
//       }
//     } else {  // endcap
//       stripPos[0] = hit.getPosition()[0];
//       stripPos[1] = hit.getPosition()[1];
//       stripPos[2] = hit.getPosition()[2];  // z never changes
//       // there's almost certainly a neater way to get these different polarities in here...DJ
//       //    ....but this is, I believe, correct!
//       int xsign = m == 6 ? -1 : +1;  // endcaps are mirrored in x, not in y

//       switch (s) {
//         case 0:
//           if (collateAlongI)
//             stripPos[0] += -xsign * relativePos;
//           else
//             stripPos[1] += -relativePos;
//           break;
//         case 1:
//           if (collateAlongI)
//             stripPos[1] += +relativePos;
//           else
//             stripPos[0] += -xsign * relativePos;
//           break;
//         case 2:
//           if (collateAlongI)
//             stripPos[0] += +xsign * relativePos;
//           else
//             stripPos[1] += +relativePos;
//           break;
//         case 3:
//           if (collateAlongI)
//             stripPos[1] += -relativePos;
//           else
//             stripPos[0] += +xsign * relativePos;
//           break;
//       }

//     }  // endcap

//     //    cout << "new strip pos: " << stripPos[0] << " " << stripPos[1] << " " << stripPos[2] << endl;

//     edm4hep::MutableSimCalorimeterHit* newhit = nullptr;

//     // check if we have already seen a virtual cell from this strip
//     if (newhits.find(cellID) != newhits.end()) {  // already have one from this strip

//       edm4hep::MutableSimCalorimeterHit* htt = newhits.find(cellID)->second;  // previously made hit in this stirp
//       // check that calculated positions are same!
//       bool isOK = true;
//       for (int i = 0; i < 3; i++) {
//         if (fabs(htt->getPosition()[i] - stripPos[i]) > 1) {  // should be identical, but allow a little slop
//           isOK = false;
//           break;
//         }
//       }
//       if (!isOK) {
//         error() << "strip virtual cell merging: inconsistent position calc!" << std::endl
//                              << "please check that the parameter ECAL_strip_nVirtualCells is correctly set..."
//                              << getNumberOfVirtualCells() << std::endl
//                              << " layer = (k-1) " << km1 << " (gear) " << gearlayer << endl;

//         for (int i = 0; i < 3; i++) {
//             debug() << stripPos[i] << " ";
//         }
//             debug() << endl;

//         for (int i = 0; i < 3; i++) {
//             debug() << htt->getPosition()[i] << " ";
//         }
//             debug() << endl;

//         // std::cout << "K-1 = " << km1 ;
//         // std::cout << "; GEARlay = " << gearlayer;
//         // std::cout << "; m = " << m;
//         // std::cout << "; s = " << s;
//         // std::cout << "; i = " << i_index;
//         // std::cout << "; j = " << j_index << std::endl;

//         assert(0);
//       }

//       // set it to the existing one
//       newhit = htt;

//     } else {  // it's the first hit from this strip
//       // create a new hit for the strip
//       newhit = &stripcol.create();
//       newhits[cellID] = newhit;
//       newhit->setCellID(cellID);
//       newhit->setPosition(stripPos);
//       newhit->setEnergy(0);  // this is added when we add MC contributions
//     }

//     // add the MC contriutions
//     float eadd(0);
//     auto  singleHit = hit.getContributions();

//     for (int ij = 0; ij < singleHit.size(); ij++) {
//       edm4hep::MutableSimCalorimeterHitCollection contrib = col.create();
//       contrib.setParticle(hit.getParticle(ij));
//       contrib.setEnergy(hit.getEnergy(ij) * energyNonuniformityScaling);
//       contrib.setTime(hit.getTime(ij));
//       contrib.setPDG(hit.getPDG(ij));
//       newhit.addToContributions(contrib);
//       eadd += hit.getEnergy(ij) * energyNonuniformityScaling;
//     }

//     float esum(0);
//     auto  singlenewHit = newhit.getContributions();
//     for (int ij = 0; ij < singlenewHit.size(); ij++) {
//       esum += newhit.getEnergy(ij);
//     }
//   }  // loop over hits

//   // move the hits from the temporary storage to the output collection
//   for (std::map<int, edm4hep::MutableSimCalorimeterHit*>::iterator itt = newhits.begin(); itt != newhits.end(); itt++) {
//     stripcol->addElement(itt->second);
//   }

//   // return the new collection
//   return stripcol;
// }