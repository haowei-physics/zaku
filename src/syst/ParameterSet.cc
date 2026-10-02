#include "zaku/syst/ParameterSet.h"

#include <TFile.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TVectorD.h>

#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace zaku {

std::string Parameter::describe() const {
  std::ostringstream os;
  os << name << ": prior " << prior << ", sigma " << sigma;
  if (prior != 0.0) os << " (" << std::fixed << std::setprecision(2) << 100.0 * sigma / std::fabs(prior) << "%)";
  os << std::defaultfloat;
  if (hasBounds) os << ", bounds [" << lower << ", " << upper << "]";
  if (fixed) os << ", fixed";
  return os.str();
}

ParameterSet::ParameterSet(const ParameterSet& o)
    : prefix_(o.prefix_), title_(o.title_), params_(o.params_), binning_(o.binning_), hasBinning_(o.hasBinning_) {
  cov_.ResizeTo(o.cov_);
  cov_ = o.cov_;
}

ParameterSet& ParameterSet::operator=(const ParameterSet& o) {
  if (this == &o) return *this;
  prefix_ = o.prefix_;
  title_ = o.title_;
  params_ = o.params_;
  binning_ = o.binning_;
  hasBinning_ = o.hasBinning_;
  cov_.ResizeTo(o.cov_);
  cov_ = o.cov_;
  return *this;
}

ParameterSet& ParameterSet::add(const Parameter& p) {
  if (p.name.empty()) throw std::invalid_argument("ParameterSet::add: a parameter needs a name");
  if (index(p.name) >= 0) throw std::invalid_argument("ParameterSet::add: parameter '" + p.name + "' already exists");
  if (p.sigma < 0.0) throw std::invalid_argument("ParameterSet::add: '" + p.name + "' has a negative sigma");
  if (p.hasBounds && p.lower > p.upper) {
    throw std::invalid_argument("ParameterSet::add: '" + p.name + "' has lower bound above the upper bound");
  }
  const int n = static_cast<int>(params_.size());
  params_.push_back(p);
  cov_.ResizeTo(n + 1, n + 1);  // existing elements are preserved, new ones are zero
  cov_(n, n) = p.sigma * p.sigma;
  return *this;
}

int ParameterSet::index(const std::string& name) const {
  for (std::size_t i = 0; i < params_.size(); ++i) {
    if (params_[i].name == name) return static_cast<int>(i);
  }
  return -1;
}

std::vector<std::string> ParameterSet::names() const {
  std::vector<std::string> out;
  for (const auto& p : params_) out.push_back(p.name);
  return out;
}

bool ParameterSet::anyBounds() const {
  for (const auto& p : params_) {
    if (p.hasBounds) return true;
  }
  return false;
}

void ParameterSet::setCovariance(const TMatrixD& cov, double symTol) {
  const int n = static_cast<int>(params_.size());
  if (cov.GetNrows() != n || cov.GetNcols() != n) {
    throw std::invalid_argument("ParameterSet::setCovariance: matrix is " + std::to_string(cov.GetNrows()) + "x" +
                                std::to_string(cov.GetNcols()) + " but the set has " + std::to_string(n) + " parameters");
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      const double scale = std::max(1.0, std::max(std::fabs(cov(i, j)), std::fabs(cov(j, i))));
      if (std::fabs(cov(i, j) - cov(j, i)) > symTol * scale) {
        throw std::invalid_argument("ParameterSet::setCovariance: matrix is not symmetric at (" + std::to_string(i) +
                                    "," + std::to_string(j) + ")");
      }
    }
  }
  cov_.ResizeTo(cov);
  cov_ = cov;
  for (int i = 0; i < n; ++i) params_[static_cast<std::size_t>(i)].sigma = std::sqrt(std::max(0.0, cov_(i, i)));
}

void ParameterSet::setCorrelation(int i, int j, double rho) {
  if (i < 0 || j < 0 || i >= static_cast<int>(params_.size()) || j >= static_cast<int>(params_.size())) {
    throw std::out_of_range("ParameterSet::setCorrelation: index out of range");
  }
  if (rho < -1.0 || rho > 1.0) throw std::invalid_argument("ParameterSet::setCorrelation: rho must be in [-1, 1]");
  if (i == j) return;
  const double c = rho * sigma(i) * sigma(j);
  cov_(i, j) = c;
  cov_(j, i) = c;
}

double ParameterSet::sigma(int i) const { return std::sqrt(std::max(0.0, cov_(i, i))); }

double ParameterSet::correlation(int i, int j) const {
  const double d = sigma(i) * sigma(j);
  return d > 0.0 ? cov_(i, j) / d : 0.0;
}

bool ParameterSet::isPositiveDefinite(double* lambdaMin, double* lambdaMax) const {
  const int n = cov_.GetNrows();
  if (n == 0) return false;
  TMatrixDSym sym(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) sym(i, j) = 0.5 * (cov_(i, j) + cov_(j, i));
  }
  TMatrixDSymEigen eig(sym);
  const TVectorD& v = eig.GetEigenValues();  // descending
  if (lambdaMin) *lambdaMin = v(n - 1);
  if (lambdaMax) *lambdaMax = v(0);
  return v(n - 1) > 0.0;
}

void ParameterSet::write(const std::string& file) const {
  if (params_.empty()) throw std::runtime_error("ParameterSet::write: the set '" + prefix_ + "' is empty");
  if (prefix_.empty()) throw std::runtime_error("ParameterSet::write: the set has no prefix");
  const int n = static_cast<int>(params_.size());
  std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "RECREATE"));
  if (!f || f->IsZombie()) throw std::runtime_error("ParameterSet::write: cannot create '" + file + "'");

  TObjArray names(n);
  names.SetOwner(true);
  TVectorD priors(n), lower(n), upper(n);
  for (int i = 0; i < n; ++i) {
    const Parameter& p = params_[static_cast<std::size_t>(i)];
    names.Add(new TObjString(p.name.c_str()));
    priors(i) = p.prior;
    lower(i) = p.hasBounds ? p.lower : p.prior;
    upper(i) = p.hasBounds ? p.upper : p.prior;
  }
  f->WriteTObject(&names, (prefix_ + "_param_names").c_str(), "SingleKey Overwrite");
  f->WriteTObject(&priors, (prefix_ + "_param_priors").c_str(), "Overwrite");
  if (anyBounds()) {
    f->WriteTObject(&lower, (prefix_ + "_param_lower").c_str(), "Overwrite");
    f->WriteTObject(&upper, (prefix_ + "_param_upper").c_str(), "Overwrite");
  }
  TMatrixD cov(cov_);  // a TMatrixD, never a TMatrixDSym
  f->WriteTObject(&cov, (prefix_ + "_param_cov").c_str(), "Overwrite");
  f->Close();
}

ParameterSet ParameterSet::read(const std::string& file, const std::string& prefix) {
  std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("ParameterSet::read: cannot open '" + file + "'");
  auto* names = f->Get<TObjArray>((prefix + "_param_names").c_str());
  auto* priors = f->Get<TVectorD>((prefix + "_param_priors").c_str());
  auto* cov = f->Get<TMatrixD>((prefix + "_param_cov").c_str());
  if (!names || !priors || !cov) {
    throw std::runtime_error("ParameterSet::read: '" + file + "' has no complete '" + prefix + "' parameter set");
  }
  auto* lower = f->Get<TVectorD>((prefix + "_param_lower").c_str());
  auto* upper = f->Get<TVectorD>((prefix + "_param_upper").c_str());
  const int n = names->GetEntries();
  if (priors->GetNrows() != n || cov->GetNrows() != n || cov->GetNcols() != n) {
    throw std::runtime_error("ParameterSet::read: inconsistent sizes in '" + file + "'");
  }
  ParameterSet s(prefix);
  for (int i = 0; i < n; ++i) {
    Parameter p;
    p.name = names->At(i)->GetName();
    p.prior = (*priors)(i);
    p.sigma = std::sqrt(std::max(0.0, (*cov)(i, i)));
    if (lower && upper && lower->GetNrows() == n && upper->GetNrows() == n) {
      p.hasBounds = true;
      p.lower = (*lower)(i);
      p.upper = (*upper)(i);
    }
    s.add(p);
  }
  s.setCovariance(*cov, 1e-6);
  return s;
}

namespace {
bool agree(double a, double b, double tol) {
  const double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
  return std::fabs(a - b) <= tol * scale;
}
}  // namespace

bool ParameterSet::sameAs(const ParameterSet& other, double tol, std::string* why) const {
  auto fail = [&](const std::string& s) {
    if (why) *why = s;
    return false;
  };
  if (params_.size() != other.params_.size()) {
    return fail("size " + std::to_string(params_.size()) + " vs " + std::to_string(other.params_.size()));
  }
  const int n = static_cast<int>(params_.size());
  // Bounds are written as full vectors (the prior standing in for an unbounded
  // parameter), so compare the effective bounds whenever either set has any.
  if (anyBounds() != other.anyBounds()) return fail("bounds present in one set only");
  for (int i = 0; i < n; ++i) {
    const Parameter& a = params_[static_cast<std::size_t>(i)];
    const Parameter& b = other.params_[static_cast<std::size_t>(i)];
    if (a.name != b.name) return fail("name " + std::to_string(i) + ": " + a.name + " vs " + b.name);
    if (!agree(a.prior, b.prior, tol)) return fail(a.name + " prior " + std::to_string(a.prior) + " vs " + std::to_string(b.prior));
    if (anyBounds()) {
      const double al = a.hasBounds ? a.lower : a.prior, au = a.hasBounds ? a.upper : a.prior;
      const double bl = b.hasBounds ? b.lower : b.prior, bu = b.hasBounds ? b.upper : b.prior;
      if (!agree(al, bl, tol) || !agree(au, bu, tol)) return fail(a.name + " bounds differ");
    }
    for (int j = 0; j < n; ++j) {
      if (!agree(cov_(i, j), other.cov_(i, j), tol)) {
        return fail("covariance (" + a.name + ", " + params_[static_cast<std::size_t>(j)].name + ") " +
                    std::to_string(cov_(i, j)) + " vs " + std::to_string(other.cov_(i, j)));
      }
    }
  }
  return true;
}

std::string ParameterSet::describe() const {
  std::ostringstream os;
  os << prefix_ << " parameter set";
  if (!title_.empty()) os << " (" << title_ << ")";
  os << ": " << params_.size() << " parameters";
  double lo = 0, hi = 0;
  if (!params_.empty()) {
    const bool pd = isPositiveDefinite(&lo, &hi);
    os << ", " << (pd ? "positive definite" : "NOT positive definite") << " (eigenvalues " << lo << " .. " << hi << ")";
  }
  if (hasBinning_) os << ", " << binning_.describe();
  os << "\n";
  for (const auto& p : params_) os << "  " << p.describe() << "\n";
  return os.str();
}

}  // namespace zaku
