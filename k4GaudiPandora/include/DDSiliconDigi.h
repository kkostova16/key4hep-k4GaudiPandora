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
#ifndef DDSILICONDIGI_H
#define DDSILICONDIGI_H 1

#include "DDScCaloDigi.h"
#include "Gaudi/Property.h"

class DDSiliconDigi final : public DDScCaloDigi {
public:
  DDSiliconDigi(const std::string& name, ISvcLocator* svcLoc);
  //~DDSiliconDigi() {}

float siliconDigi(float energy) const;

private:
Gaudi::Property<float> m_ehEnergy{this, "EnergyPerEHpair", {3.6}, "Energy required to create e-h pair in silicon (in eV)"}; 
Gaudi::Property<float> m_elecMaxDynRange{this, "ElectronicsMaxDynamicRange_MIP", {2500.0}, "Maximum of electronis dynamic range (in MIP units)"};
Gaudi::Property<float> m_elecNoise{this, "ElectronicsNoise_MIP", {0.0}, "Typical electronics noise (in MIP units)"};
Gaudi::Property<float> m_calibMIP{this, "CalibrationMIP", {1.0e-4}, "Calibration to convert deposited energy to MIPs"};
};

#endif