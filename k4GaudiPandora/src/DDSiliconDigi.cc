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
#include "DDSiliconDigi.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

DECLARE_COMPONENT(DDSiliconDigi)

DDSiliconDigi::DDSiliconDigi(const std::string& name, ISvcLocator* svcLoc) 
  : DDScCaloDigi(name, svcLoc) {}

// this applies extra digitisation to silicon hits
float DDSiliconDigi::siliconDigi(float energy) const {

  // calculate #e-h pairs
  float nehpairs = 1.0e9 * energy / m_ehEnergy;  // check units of energy! m_ehEnergy is in eV, energy in GeV

  // fluctuate it by Poisson
  float smeared_energy = energy * CLHEP::RandPoisson::shoot(nehpairs) / nehpairs;

  // limited electronics dynamic range // Daniel moved electronics dyn range to here
  if (m_elecMaxDynRange > 0)
    smeared_energy = std::min(smeared_energy, m_elecMaxDynRange * m_calibMIP);

  // add electronics noise
  if (m_elecNoise > 0)
    smeared_energy += CLHEP::RandGauss::shoot(0, m_elecNoise * m_calibMIP);

  return smeared_energy;
}