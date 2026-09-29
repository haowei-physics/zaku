#include "zaku/core/Detector.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "zaku/core/Units.h"

namespace zaku {

Detector::Detector(std::string name, double fiducialMass_kt, std::vector<Target> composition)
    : name_(std::move(name)), fiducialMass_kt_(fiducialMass_kt), composition_(std::move(composition)) {
  if (fiducialMass_kt_ <= 0.0) throw std::invalid_argument("Detector '" + name_ + "': fiducial mass must be positive");
  if (composition_.empty()) throw std::invalid_argument("Detector '" + name_ + "': composition is empty");
  double sum = 0.0;
  for (const auto& t : composition_) {
    if (t.A <= 0 || t.Z < 0 || t.Z > t.A)
      throw std::invalid_argument("Detector '" + name_ + "': invalid target Z=" + std::to_string(t.Z) +
                                  " A=" + std::to_string(t.A));
    if (t.massFraction < 0.0)
      throw std::invalid_argument("Detector '" + name_ + "': negative mass fraction for " + t.name());
    sum += t.massFraction;
  }
  if (std::fabs(sum - 1.0) > 1.0e-6)
    throw std::invalid_argument("Detector '" + name_ + "': mass fractions sum to " + std::to_string(sum) +
                                ", expected 1");
}

Detector Detector::water(double kt) {
  return Detector("water", kt, {Target(8, 16, 16.0 / 18.0), Target(1, 1, 2.0 / 18.0)});
}
Detector Detector::argon(double kt) { return Detector("argon", kt, {Target(18, 40, 1.0)}); }
Detector Detector::scintillator(double kt) {
  return Detector("scintillator", kt, {Target(6, 12, 12.0 / 13.0), Target(1, 1, 1.0 / 13.0)});
}
Detector Detector::iron(double kt) { return Detector("iron", kt, {Target(26, 56, 1.0)}); }

Detector Detector::fromMaterial(const std::string& material, double kt) {
  std::string m = material;
  std::transform(m.begin(), m.end(), m.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (m == "water" || m == "h2o") return water(kt);
  if (m == "argon" || m == "ar" || m == "lar") return argon(kt);
  if (m == "scintillator" || m == "ch" || m == "plastic") return scintillator(kt);
  if (m == "iron" || m == "fe" || m == "steel") return iron(kt);
  throw std::invalid_argument("Detector::fromMaterial: unknown material '" + material +
                              "' (water, argon, scintillator, iron; or give an explicit composition)");
}

double Detector::nNucleons() const {
  if (nucleonsOverride_ > 0.0) return nucleonsOverride_;
  // N = sum_i (M f_i / A_i) N_A A_i = M N_A, with the molar mass taken as A g/mol.
  double fsum = 0.0;
  for (const auto& t : composition_) fsum += t.massFraction;
  return fiducialMass_kt_ * units::kt_to_g * units::Avogadro * fsum;
}

double Detector::nNuclei(const Target& t) const {
  if (t.A <= 0) return 0.0;
  return fiducialMass_kt_ * units::kt_to_g * t.massFraction / static_cast<double>(t.A) * units::Avogadro;
}

double Detector::nNucleiTotal() const {
  double n = 0.0;
  for (const auto& t : composition_) n += nNuclei(t);
  return n;
}

std::string Detector::genieTargetString() const {
  std::ostringstream os;
  os << std::fixed << std::setprecision(6);
  bool first = true;
  for (const auto& t : composition_) {
    if (!first) os << ",";
    first = false;
    os << t.nuclearPdg() << "[" << t.massFraction << "]";
  }
  return os.str();
}

NeutTarget Detector::neutTarget() const {
  if (hasNeutTargetOverride_) return neutTargetOverride_;
  const Target* heavy = nullptr;
  double hydrogenFraction = 0.0;
  int nHeavy = 0;
  for (const auto& t : composition_) {
    if (t.Z == 1 && t.A == 1) {
      hydrogenFraction += t.massFraction;
      continue;
    }
    ++nHeavy;
    if (!heavy || t.A > heavy->A) heavy = &t;
  }
  if (nHeavy > 1) {
    throw std::runtime_error("Detector '" + name_ +
                             "': NEUT needs a single bound nucleus; the composition has " + std::to_string(nHeavy) +
                             ". Call setNeutTarget() to choose one explicitly.");
  }
  NeutTarget nt;
  if (!heavy) {
    // Pure hydrogen: NEUT models it as free protons only.
    nt.numNucleons = 1;
    nt.numBoundProtons = 0;
    nt.numFreeProtons = 1;
    return nt;
  }
  nt.numNucleons = heavy->A;
  nt.numBoundProtons = heavy->Z;
  // Free protons per heavy nucleus: (moles of H) / (moles of heavy nucleus).
  const double nH = hydrogenFraction / 1.0;
  const double nX = heavy->massFraction / static_cast<double>(heavy->A);
  nt.numFreeProtons = nX > 0.0 ? static_cast<int>(std::lround(nH / nX)) : 0;
  return nt;
}

void Detector::setNeutTarget(const NeutTarget& t) {
  neutTargetOverride_ = t;
  hasNeutTargetOverride_ = true;
}

std::string Detector::describe() const {
  std::ostringstream os;
  os << "Detector '" << name_ << "': " << fiducialMass_kt_ << " kt fiducial, composition";
  for (const auto& t : composition_) os << " " << t.name() << "(" << t.massFraction << ")";
  os << std::scientific << std::setprecision(4) << ", nucleons " << nNucleons();
  if (nucleonsOverride_ > 0.0) os << " (overridden)";
  return os.str();
}

}  // namespace zaku
