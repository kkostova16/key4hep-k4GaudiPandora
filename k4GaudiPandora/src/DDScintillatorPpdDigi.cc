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
#include "DDScintillatorPpdDigi.h"
#include <assert.h>
#include <iostream>
#include "CLHEP/Random/RandBinomial.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"
using std::cout;
using std::endl;

DECLARE_COMPONENT(DDScintillatorPpdDigi)

// this applies some extra digitisation to scintillator+PPD hits (PPD=SiPM, MPPC)
// - poisson fluctuates the number of photo-electrons according to #PEs/MIP
// - applies PPD saturation according to #pixels
// Daniel Jeans, Jan/Feb 2014.
// (split off from ILDCaloDigi Aug'14.)

DDScintillatorPpdDigi::DDScintillatorPpdDigi(const std::string& name, ISvcLocator* svcLoc) 
  : DDScCaloDigi(name, svcLoc) {}

void DDScintillatorPpdDigi::printParameters() {
  cout << "--------------------------------" << endl;
  cout << " DDScintillatorPpdDigi parameters " << endl;
  cout << " PEperMIP        = " << m_PEperMIP << endl;
  cout << " calibrationMIP  = " << m_calibMIP << endl;
  cout << " Npixels         = " << m_Npix << endl;
  cout << " misCalibNpix    = " << m_misCalibNpix << endl;
  cout << " pixSpread       = " << m_pixSpread << endl;
  cout << " elecNoise       = " << m_elecNoise << endl;
  cout << " elecMaxDynRange = " << m_elecMaxDynRange << endl;
  cout << "--------------------------------" << endl;
  return;
}

// Convert thresholds to GeV units
float DDScintillatorPpdDigi::convertThresholdUnits(std::string unitThreshold, float threshold) {

  if (unitThreshold.compare("GeV") == 0) {
    // threshold unit is GeV, do nothing
  } else if (unitThreshold.compare("MIP") == 0) {
    // threshold unit is MIP, convert via MIP2GeV
    threshold *= m_calibMIP.value();
  } else if (unitThreshold.compare("px") == 0) {
    // threshold unit is pixels, convert via MIP2GeV and lightyield
    threshold *= m_PEperMIP.value() * m_calibMIP.value();
  } else {
    error() << "Could not identify threshold unit. Please use \"GeV\", \"MIP\" or \"px\"! Aborting." << endmsg;
    assert(0);
  }
  return threshold; // in GeV units
}

float DDScintillatorPpdDigi::getDigitisedEnergy(float energy) {
  float correctedEnergy = energy;

  if (m_PEperMIP <= 0 || m_calibMIP <= 0 || m_Npix <= 0) {
    cout << "ERROR, crazy parameters for DDScintillatorPpdDigi: PE/MIP = " << m_PEperMIP << ", MIP calibration = " << m_calibMIP
         << ", # pixels = " << m_Npix << endl;
    cout << "You must specify at least the #PE/MIP, MIP calibration, and #pixels for realistic scintillator digitisation!!"
         << endl;
    cout << "Refusing to proceed!" << endl;
    assert(0);
  }
  else {
    // 1. convert energy to expected # photoelectrons (via MIPs)
    float Npe = m_PEperMIP * energy / m_calibMIP;

    //oh: commented out Daniel's digitisation model. (npe -> poisson -> saturation -> stoykov smearing).
    // Stoykov smearing used with Gaussian shape for lack of better model.
    /*
    // 2. smear according to poisson (PE statistics)
    npe = CLHEP::RandPoisson::shoot( npe );

    if (_npix>0) {
      // 3. apply average sipm saturation behaviour
      npe = _npix*(1.0 - exp( -npe/_npix ) );

      // 4. apply intrinsic sipm fluctuations (see e.g arXiv:0706.0746 [physics.ins-det])
      float alpha = npe/_npix; // frac of hit pixels
      float width = sqrt( _npix*exp(-alpha)*( 1. - (1.+alpha)*exp(-alpha) ) );

      // make an integer correction
      int dpix = int( floor( CLHEP::RandGauss::shoot(0, width) + 0.5 ) );

      npe += dpix;
    }
  */

    //AHCAL TB style digitisation: npe -> saturation -> binomial smearing
    //shown to be mathematically equivalent to Daniel's model above, but slightly faster and automatically generates correct shape instead of Gaussian approximation

    if (m_Npix > 0) {
      // apply average SiPM saturation behaviour
      Npe = m_Npix * (1.0 - exp(- Npe / m_Npix));

      //apply binomial smearing
      float p = Npe / m_Npix;                           // fraction of hit pixels on SiPM
      Npe     = CLHEP::RandBinomial::shoot(m_Npix, p);  // # photoelectrons now quantised to integer pixels
    }

    if (m_pixSpread > 0) {
      // variations in pixel capacitance
      Npe *= CLHEP::RandGauss::shoot(1., m_pixSpread / sqrt(Npe));
    }

    if (m_elecMaxDynRange > 0) {
      // limited dynamic range of readout electronics
      // Daniel moved this here, before the unfolding of saturation (September 2015)
      Npe = std::min(Npe, m_elecMaxDynRange * m_PEperMIP);
    }

    if (m_elecNoise > 0) {
      // add electronics noise
      Npe += CLHEP::RandGauss::shoot(0., m_elecNoise * m_PEperMIP);
    }

    if (m_Npix > 0) {
      // 4. unfold the saturation
      // - miscalibration of Npix
      float smearedNpix = m_misCalibNpix > 0 ? static_cast<Gaudi::Property<float>> (m_Npix * CLHEP::RandGauss::shoot(1.0, m_misCalibNpix)) : m_Npix;

    //oh: commented out daniel's implmentation of dealing with hits>smearedNpix. using linearisation of saturation-reconstruction for high amplitude hits instead.
    /*
    // - this is to deal with case when #pe is larger than #pixels (would mean taking log of negative number)
    float epsilon=1; // any small number...
    if ( npe>smearedNpix-epsilon ) npe=smearedNpix-epsilon;
    // - unfold saturation
    npe = -smearedNpix * std::log ( 1. - ( npe / smearedNpix ) );
    */

    const float r =
        0.95;  //this is the fraction of SiPM pixels fired above which a linear continuation of the saturation-reconstruction function is used. 0.95 of nPixel corresponds to a energy correction of factor ~3.

    if (Npe < r * smearedNpix) {  //current hit below linearisation threshold, reconstruct energy normally:
      Npe = -smearedNpix * std::log(1. - (Npe / smearedNpix));
    } 
    else {  //current hit is aove linearisation threshold, reconstruct using linear continuation function:
      Npe = 1 / (1 - r) * (Npe - r * smearedNpix) - smearedNpix * std::log(1 - r);
    }
  }

  // convert back to energy
  correctedEnergy = m_calibMIP * Npe / m_PEperMIP;
  }
  return correctedEnergy;
}
