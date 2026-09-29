#ifndef ZAKU_CORE_DETECTOR_H
#define ZAKU_CORE_DETECTOR_H
#include <string>
#include <vector>

#include "zaku/core/Types.h"

namespace zaku {

/// The target description NEUT's card wants: one bound nucleus plus a number
/// of free protons per nucleus (H2O = 16 nucleons, 8 bound protons, 2 free protons).
struct NeutTarget {
  int numNucleons = 0;
  int numBoundProtons = 0;
  int numFreeProtons = 0;
};

/// A parametric detector: a fiducial mass of a homogeneous material. There is
/// no geometry; the fiducial mass is an assumption the user asserts and every
/// normalisation inherits it.
class Detector {
public:
  Detector() = default;
  /// Composition mass fractions must sum to 1 (within 1e-6).
  Detector(std::string name, double fiducialMass_kt, std::vector<Target> composition);

  static Detector water(double fiducialMass_kt);         ///< H2O: O-16 16/18, H-1 2/18 by mass
  static Detector argon(double fiducialMass_kt);         ///< Ar-40
  static Detector scintillator(double fiducialMass_kt);  ///< CH: C-12 12/13, H-1 1/13 by mass
  static Detector iron(double fiducialMass_kt);          ///< Fe-56
  /// water | argon | scintillator | iron (case-insensitive).
  static Detector fromMaterial(const std::string& material, double fiducialMass_kt);

  const std::string& name() const { return name_; }
  double fiducialMass_kt() const { return fiducialMass_kt_; }
  const std::vector<Target>& composition() const { return composition_; }

  /// Number of nucleons in the fiducial mass: mass[g] x N_A (molar mass ~ A g/mol,
  /// so the composition cancels), unless overridden with setNucleonsOverride.
  /// 188 kt of water gives 1.132e35.
  double nNucleons() const;
  /// Number of nuclei of one target species: mass[g] x f / A x N_A.
  double nNuclei(const Target& t) const;
  double nNucleiTotal() const;
  /// Assert the nucleon count directly (0 restores the mass-based value).
  void setNucleonsOverride(double n) { nucleonsOverride_ = n; }
  double nucleonsOverride() const { return nucleonsOverride_; }

  /// GENIE target mix: "1000080160[0.888889],1000010010[0.111111]".
  std::string genieTargetString() const;

  /// NEUT target: the heaviest nucleus plus the number of hydrogen atoms per
  /// such nucleus as free protons. Throws std::runtime_error when the
  /// composition has more than one non-hydrogen nucleus and no override was set.
  NeutTarget neutTarget() const;
  void setNeutTarget(const NeutTarget& t);

  /// Human-readable summary.
  std::string describe() const;

private:
  std::string name_ = "detector";
  double fiducialMass_kt_ = 0.0;
  std::vector<Target> composition_;
  double nucleonsOverride_ = 0.0;
  NeutTarget neutTargetOverride_;
  bool hasNeutTargetOverride_ = false;
};

}  // namespace zaku

#endif  // ZAKU_CORE_DETECTOR_H
