#ifndef ZAKU_DETECTOR_SMEARING_H
#define ZAKU_DETECTOR_SMEARING_H
#include <memory>
#include <string>

#include "zaku/core/Units.h"
#include "zaku/detector/Formula.h"
#include "zaku/detector/Variables.h"

class TH1;
class TH1D;

namespace zaku {

/// How wide the Gaussian smearing of one quantity is, as a function of the
/// true value (and, for formulas, any truth variable).
class ResolutionModel {
public:
  enum class Kind : int {
    None = 0,        ///< no smearing
    Fractional = 1,  ///< sigma = fraction x value
    Absolute = 2,    ///< sigma = constant (same unit as the value)
    Histogram = 3,   ///< sigma from a TH1 binned in the true value, bins clamped (legacy)
    Stochastic = 4,  ///< sigma/E = a/sqrt(E) (+) b, i.e. sigma = sqrt((a sqrt(E))^2 + (b E)^2)
    Formula = 5      ///< sigma = formula over the truth variables and `value`
  };

  ResolutionModel() = default;
  static ResolutionModel none();
  static ResolutionModel fractional(double fraction);
  static ResolutionModel absolute(double sigma);
  static ResolutionModel histogram(const std::string& file, const std::string& name, EnergyUnit units = EnergyUnit::GeV);
  static ResolutionModel histogram(TH1* h, const std::string& label = "");
  static ResolutionModel stochastic(double a, double b);
  static ResolutionModel formula(const std::string& expression);

  Kind kind() const { return kind_; }
  bool active() const { return kind_ != Kind::None; }
  double parameterA() const { return a_; }
  double parameterB() const { return b_; }
  const std::string& expression() const { return expression_; }
  const std::string& histogramLabel() const { return label_; }
  const TH1* histogramPtr() const { return hist_.get(); }

  void prepare(const VariableSet& vars, const std::string& where);
  /// sigma for a true value (the variable set supplies formula inputs).
  double sigma(double value, const VariableSet& vars) const;

  std::string describe() const;

private:
  Kind kind_ = Kind::None;
  double a_ = 0.0, b_ = 0.0;
  std::shared_ptr<TH1> hist_;
  std::string label_;
  std::string expression_;
  Formula formula_;
  int valueIndex_ = -1;
};

/// Decay-electron (Michel) tagging efficiencies used for `true_ME`/`recon_ME`.
struct DecayElectronTagging {
  double muPlus = 0.96;   ///< mu+
  double muMinus = 0.83;  ///< mu-, nuclear capture included
  double piPlus = 0.45;   ///< pi+ -> mu+ -> e+
};

/// The reconstruction resolutions of the fake detector.
///
/// The legacy chain smeared the outgoing lepton momentum with a fractional
/// resolution of 2.4% and its direction with 1.8 degrees, both through
/// histograms binned in the true value (`data/hk/smearing.root`).
/// `parametric()` rebuilds exactly those histograms; `fromLegacyFile()` reads them.
class Smearing {
public:
  Smearing();

  static Smearing parametric(double leptonMomentumFraction = 0.024, double leptonAngle_deg = 1.8);
  static Smearing fromLegacyFile(const std::string& smearingFile);
  /// The legacy histogram of sigma_p vs p: 50 bins on [0, 5] GeV, sigma = fraction x centre.
  static TH1D* legacyMomentumHistogram(double fraction, const std::string& name = "h_dp_vs_p");
  /// The legacy histogram of sigma_cos vs cos: 40 bins on [-1, 1], sigma = sin(theta) x angle.
  static TH1D* legacyAngleHistogram(double angle_deg, const std::string& name = "h_dcos_vs_cos");

  ResolutionModel& leptonMomentum() { return leptonMomentum_; }
  ResolutionModel& leptonAngle() { return leptonAngle_; }       ///< sigma of cos(theta)
  ResolutionModel& hadronicEnergy() { return hadronicEnergy_; }
  ResolutionModel& electronEnergy() { return electronEnergy_; }  ///< calorimetric lepton term, e hypothesis
  ResolutionModel& muonEnergy() { return muonEnergy_; }          ///< calorimetric lepton term, mu hypothesis
  const ResolutionModel& leptonMomentum() const { return leptonMomentum_; }
  const ResolutionModel& leptonAngle() const { return leptonAngle_; }
  const ResolutionModel& hadronicEnergy() const { return hadronicEnergy_; }
  const ResolutionModel& electronEnergy() const { return electronEnergy_; }
  const ResolutionModel& muonEnergy() const { return muonEnergy_; }
  Smearing& setLeptonMomentum(const ResolutionModel& m) { leptonMomentum_ = m; return *this; }
  Smearing& setLeptonAngle(const ResolutionModel& m) { leptonAngle_ = m; return *this; }
  Smearing& setHadronicEnergy(const ResolutionModel& m) { hadronicEnergy_ = m; return *this; }
  Smearing& setElectronEnergy(const ResolutionModel& m) { electronEnergy_ = m; return *this; }
  Smearing& setMuonEnergy(const ResolutionModel& m) { muonEnergy_ = m; return *this; }

  DecayElectronTagging& decayElectronTagging() { return tagging_; }
  const DecayElectronTagging& decayElectronTagging() const { return tagging_; }
  Smearing& setDecayElectronTagging(const DecayElectronTagging& t) { tagging_ = t; return *this; }

  /// Visible fraction of neutron kinetic energy in Ehad_true (0 = invisible).
  double neutronVisibleFraction() const { return neutronFraction_; }
  Smearing& setNeutronVisibleFraction(double f) { neutronFraction_ = f; return *this; }

  void prepare(const VariableSet& vars);

  static Smearing fromYamlString(const std::string& yaml, const std::string& baseDir = "");
  std::string toYaml() const;
  std::string describe() const;

private:
  ResolutionModel leptonMomentum_;
  ResolutionModel leptonAngle_;
  ResolutionModel hadronicEnergy_;
  ResolutionModel electronEnergy_;
  ResolutionModel muonEnergy_;
  DecayElectronTagging tagging_;
  double neutronFraction_ = 0.0;
};

}  // namespace zaku

#endif  // ZAKU_DETECTOR_SMEARING_H
