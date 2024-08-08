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
#ifndef _DDScintillatorPpdDigi_h_
#define _DDScintillatorPpdDigi_h_

class DDScintillatorPpdDigi {
public:
  DDScintillatorPpdDigi();
  ~DDScintillatorPpdDigi() {}

  // expected # photoelectrons / MIP
  void setPEperMIP(float x) { _pe_per_mip = x; }

  // calibration factor from input hit energy to MIPs
  void setCalibMIP(float x) { _calib_mip = x; }

  // #pixels of PPD
  void setNPix(int x) { _npix = x; }

  // random miscalibration of total #pixels (as a fraction of pixel number: 0.05 = 5% miscalibration)
  void setRandomMisCalibNPix(float x) { _misCalibNpix = x; }

  // spread in pixel capacitance (as a fraction: 0.05 = 5% spread)
  void setPixSpread(float x) { _pixSpread = x; }

  // electronics noise (in MIP units)
  void setElecNoise(float x) { _elecNoise = x; }

  // electronics dynamic range (in MIP units)
  void setElecRange(float x) { _elecMaxDynRange_MIP = x; }

  float getDigitisedEnergy(float energy);

  void printParameters();

private:
  float _pe_per_mip          = -99;
  float _calib_mip           = -99;
  float _npix                = -99;
  float _misCalibNpix        = 0;
  float _pixSpread           = 0;
  float _elecNoise           = 0;
  float _elecMaxDynRange_MIP = 0;
};

#endif
