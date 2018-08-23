/*
 * Copyright (c) 2009--2018, the KLFitter developer team
 *
 * This file is part of KLFitter.
 *
 * KLFitter is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * KLFitter is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with KLFitter. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef KLFITTER_PARTICLES_PHOTON_H_
#define KLFITTER_PARTICLES_PHOTON_H_

#include "KLFitter/Particle/Base.h"

namespace KLFitter {
namespace Particle {
/**
 * Implementation of a photon class. This inherits a generic
 * structure from Particle::Base, reimplements all purely virtual
 * functions and adds functionality specific to photons.
 */
class Photon : public Base {
 public:
  /**
   * Photon constructor taking name and four momentum.
   * @param name The name of the particle
   * @param p4 The four momentum
   */
  Photon(const std::string& name, const TLorentzVector& p4)
    : m_name(name)
    , m_p4(p4) {
    // empty
  }

  /// The (defaulted) destructor.
  ~Photon() = default;

  /// The (defaulted) copy constructor.
  Photon(const Photon&) = default;

  /// The (defaulted) assignment operator.
  Photon& operator=(const Photon&) = default;

  /** \name Inherited, reimplemented methods */
  /** @{ */

  /// Const accessor for the particle's name.
  const std::string& GetName() const override { return m_name; }

  /// Const accessor to the assigned identifier of the particle.
  unsigned int GetIdentifier() const override { return m_identifier; }

  /// Const accessor for the particle's four momentum.
  const TLorentzVector& GetP4() const override { return m_p4; }

  /// Non-const accessor for the particle's four momentum.
  TLorentzVector& GetP4() override { return m_p4; }

  /// Set the identifier of the particle.
  void SetIdentifier(unsigned int id) override { m_identifier = id; }

  /// Set the four momentum.
  void SetP4(const TLorentzVector& p4) override { m_p4 = p4; }

  /** @} */
  /** \name Get methods */
  /** @{ */

  /// Get the detector eta.
  double GetDetEta() const { return m_det_eta; }

  /** @} */
  /** \name Set methods */
  /** @{ */

  /// Set the detector eta.
  void SetDetEta(double val) { m_det_eta = val; }

  /** @} */

 private:
  /// String with the name of the particle.
  std::string m_name;

  /// Assigned identifier of the particle.
  unsigned int m_identifier{0};

  /// Four vector of the particle.
  TLorentzVector m_p4;

  double m_det_eta{0.};
};
}  // namespace Particles
}  // namespace KLFitter

#endif  // KLFITTER_PARTICLES_PHOTON_H_
