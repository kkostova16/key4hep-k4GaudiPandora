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
#ifndef DDSCINTILLATORPPDDIGI_H
#define DDSCINTILLATORPPDDIGI_H

#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

class DDScintillatorPpdDigi : public Gaudi::Algorithm {
public:
  explicit DDScintillatorPpdDigi(const std::string& name, ISvcLocator* svcLoc);
  ~DDScintillatorPpdDigi() {}

  StatusCode execute(const EventContext&) const final;

  // expected # photoelectrons / MIP
  void setPEperMIP(float x) { m_PEperMIP = x; }

  // calibration factor from input hit energy to MIPs
  void setCalibMIP(float x) { m_calibMIP = x; }

  // # pixels of PPD
  void setNPix(int x) { m_Npix = x; }

  // random miscalibration of total # pixels (as a fraction of pixel number: 0.05 = 5% miscalibration)
  void setRandomMisCalibNPix(float x) { m_misCalibNpix = x; }

  // spread in pixel capacitance (as a fraction: 0.05 = 5% spread)
  void setPixSpread(float x) { m_pixSpread = x; }

  // electronics noise (in MIP units)
  void setElecNoise(float x) { m_elecNoise = x; }

  // electronics dynamic range (in MIP units)
  void setElecRange(float x) { m_elecMaxDynRange = x; }

  void printParameters();

  float getDigitisedEnergy(float energy);

  
//private:

  Gaudi::Property<float> m_PEperMIP{this, "PPD_PE_per_MIP", {7.0}, "# photoelectrons per MIP (scintillator): used to poisson smear #PEs if > 0"};
  Gaudi::Property<float> m_calibMIP{this, "CalibrationMIP", {1.0e-4}, "Calibration to convert deposited energy to MIPs"};
  Gaudi::Property<float> m_Npix{this, "PPD_NPixels", {10000}, "Total number of MPPC/SiPM pixels for implementation of saturation effect"};
  Gaudi::Property<float> m_misCalibNpix{this, "PPD_NPixels_uncertainty", {0.05}, "Fractional uncertainty of effective total number of MPPC/SiPM pixels"};
  Gaudi::Property<float> m_pixSpread{this, "PixelSpread", {0.05}, "Variation of MPPC/SiPM pixels capacitance in ECAL/HCAL (as a fraction: 0.01=1%)"};
  Gaudi::Property<float> m_elecNoise{this, "ElectronicsNoise_MIP", {0.0}, "Typical electronics noise (in MIP units)"};
  Gaudi::Property<float> m_elecMaxDynRange{this, "ElectronicsMaxDynamicRange_MIP", {2500.0}, "Maximum of electronis dynamic range (in MIP units)"};
};

#endif
