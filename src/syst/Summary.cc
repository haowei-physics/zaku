#include "zaku/syst/Summary.h"

#include <TBranch.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>

#include "zaku/core/Provenance.h"
#include "zaku/detector/Shifts.h"
#include "zaku/syst/Covariances.h"

namespace zaku {

// --- SummaryEntry / SummaryTable -----------------------------------------------------

double SummaryEntry::maxFrac() const { return std::max(std::fabs(plusFrac()), std::fabs(minusFrac())); }

namespace {
std::vector<std::string> unique(const std::vector<SummaryEntry>& e, const std::string SummaryEntry::*field) {
  std::vector<std::string> out;
  for (const auto& x : e) {
    if (std::find(out.begin(), out.end(), x.*field) == out.end()) out.push_back(x.*field);
  }
  return out;
}
}  // namespace

std::vector<std::string> SummaryTable::sets() const { return unique(entries_, &SummaryEntry::set); }
std::vector<std::string> SummaryTable::parameters() const { return unique(entries_, &SummaryEntry::parameter); }
std::vector<std::string> SummaryTable::samples() const { return unique(entries_, &SummaryEntry::sample); }

const SummaryEntry* SummaryTable::find(const std::string& parameter, const std::string& sample) const {
  for (const auto& e : entries_) {
    if (e.parameter == parameter && e.sample == sample) return &e;
  }
  return nullptr;
}

std::vector<SummaryEntry> SummaryTable::top(int n) const {
  std::vector<SummaryEntry> out = entries_;
  std::stable_sort(out.begin(), out.end(), [](const SummaryEntry& a, const SummaryEntry& b) { return a.maxFrac() > b.maxFrac(); });
  if (n >= 0 && static_cast<std::size_t>(n) < out.size()) out.resize(static_cast<std::size_t>(n));
  return out;
}

std::string SummaryTable::csv() const {
  std::ostringstream os;
  os << std::setprecision(10);
  os << "set,parameter,sample,nominal,plus,minus,plus_frac,minus_frac\n";
  for (const auto& e : entries_) {
    os << e.set << "," << e.parameter << "," << e.sample << "," << e.nominal << "," << e.plus << "," << e.minus << ","
       << e.plusFrac() << "," << e.minusFrac() << "\n";
  }
  return os.str();
}

SummaryTable SummaryTable::fromCsv(const std::string& text) {
  SummaryTable t;
  std::istringstream in(text);
  std::string line;
  bool header = true;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (header) {
      header = false;
      continue;
    }
    std::vector<std::string> f;
    std::istringstream ls(line);
    std::string item;
    while (std::getline(ls, item, ',')) f.push_back(item);
    if (f.size() < 6) throw std::runtime_error("SummaryTable::fromCsv: cannot parse '" + line + "'");
    SummaryEntry e;
    e.set = f[0];
    e.parameter = f[1];
    e.sample = f[2];
    e.nominal = std::stod(f[3]);
    e.plus = std::stod(f[4]);
    e.minus = std::stod(f[5]);
    t.add(e);
  }
  return t;
}

void SummaryTable::writeCsv(const std::string& path) const {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("SummaryTable::writeCsv: cannot write '" + path + "'");
  out << csv();
}

SummaryTable SummaryTable::readCsv(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("SummaryTable::readCsv: cannot read '" + path + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  return fromCsv(ss.str());
}

std::string SummaryTable::text() const {
  const std::vector<std::string> smp = samples();
  const std::vector<std::string> par = parameters();
  std::ostringstream os;
  os << std::left << std::setw(6) << "set" << std::setw(28) << "parameter";
  for (const auto& s : smp) os << std::right << std::setw(16) << s;
  os << "\n";
  os << std::left << std::setw(6) << "" << std::setw(28) << "nominal";
  for (const auto& s : smp) {
    double n = 0.0;
    for (const auto& e : entries_) {
      if (e.sample == s) {
        n = e.nominal;
        break;
      }
    }
    os << std::right << std::setw(16) << std::fixed << std::setprecision(1) << n;
  }
  os << "\n";
  for (const auto& p : par) {
    std::string set;
    for (const auto& e : entries_) {
      if (e.parameter == p) {
        set = e.set;
        break;
      }
    }
    os << std::left << std::setw(6) << set << std::setw(28) << p;
    for (const auto& s : smp) {
      const SummaryEntry* e = find(p, s);
      if (!e) {
        os << std::right << std::setw(16) << "-";
        continue;
      }
      std::ostringstream cell;
      cell << std::showpos << std::fixed << std::setprecision(2) << 100.0 * e->plusFrac() << "/" << 100.0 * e->minusFrac() << "%";
      os << std::right << std::setw(16) << cell.str();
    }
    os << "\n";
  }
  return os.str();
}

// --- Spectrum ----------------------------------------------------------------------------

void Spectrum::reset(const std::vector<double>& e) {
  edges = e;
  contents.assign(e.size() > 1 ? e.size() - 1 : 0, 0.0);
}

int Spectrum::find(double x) const {
  if (edges.size() < 2 || x < edges.front() || x >= edges.back()) return -1;
  const auto it = std::upper_bound(edges.begin(), edges.end(), x);
  return static_cast<int>(it - edges.begin()) - 1;
}

void Spectrum::fill(double x, double w) {
  const int b = find(x);
  if (b >= 0) contents[static_cast<std::size_t>(b)] += w;
}

double Spectrum::total() const {
  double s = 0.0;
  for (double c : contents) s += c;
  return s;
}

TH1D* Spectrum::histogram(const std::string& name, const std::string& title) const {
  auto* h = new TH1D(name.c_str(), (title.empty() ? sample : title).c_str(), nBins(), edges.data());
  h->SetDirectory(nullptr);
  h->GetXaxis()->SetTitle(variable.c_str());
  for (int i = 0; i < nBins(); ++i) h->SetBinContent(i + 1, contents[static_cast<std::size_t>(i)]);
  return h;
}

// --- SystematicsSummary --------------------------------------------------------------------

namespace {

/// Sorted unique edges of the first range variable of a binning (file binnings have no axes).
bool projectedEdges(const Binning& b, std::string& variable, std::vector<double>& edges) {
  for (std::size_t v = 0; v < b.variables().size(); ++v) {
    if (!b.isRange()[v]) continue;
    std::set<double> e;
    for (const auto& bin : b.bins()) {
      e.insert(bin.low[v]);
      e.insert(bin.high[v]);
    }
    if (e.size() < 2) continue;
    variable = b.variables()[v];
    edges.assign(e.begin(), e.end());
    return true;
  }
  return false;
}

/// Edges of a variable if it appears as a range in the binning; else empty.
std::vector<double> edgesOf(const Binning& b, const std::string& variable) {
  for (std::size_t v = 0; v < b.variables().size(); ++v) {
    if (b.variables()[v] != variable || !b.isRange()[v]) continue;
    std::set<double> e;
    for (const auto& bin : b.bins()) {
      e.insert(bin.low[v]);
      e.insert(bin.high[v]);
    }
    return std::vector<double>(e.begin(), e.end());
  }
  return {};
}

const std::vector<double>& defaultEnuEdges() {
  static const std::vector<double> e = {0.0, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0};
  return e;
}

struct Universe {
  std::string set;
  std::string parameter;
  std::vector<Spectrum> dplus;    ///< per sample: plus minus nominal
  std::vector<Spectrum> dminus;
  std::vector<Spectrum> dplusEnu;   ///< shifts only: the same in Enu_rec
  std::vector<Spectrum> dminusEnu;
};

}  // namespace

struct SystematicsSummary::Impl {
  Experiment exp;
  ParameterSet xsec, flux, det, osc;
  std::vector<SampleDefinition> samples;
  std::vector<std::string> variable;            ///< per sample: the plotted variable
  std::vector<std::vector<double>> edges;       ///< per sample
  std::vector<std::vector<double>> enuEdges;    ///< per sample: Enu_rec edges for the shift covariance
  std::vector<Spectrum> nominalEnu;
  std::vector<Universe> universes;
  std::vector<int> xsecU, fluxU, detU, oscU, shiftU;   ///< universe indices per set
  std::vector<double> oscPlusP, oscMinusP;             ///< scratch
  bool finished = false;

  Impl(const Experiment& e, const ParameterSet& x, const ParameterSet& f, const ParameterSet& d, const ParameterSet& o)
      : exp(e), xsec(x), flux(f), det(d), osc(o), samples(e.samples()) {}

  void addUniverse(const std::string& set, const std::string& name, std::vector<int>& list, bool withEnu) {
    Universe u;
    u.set = set;
    u.parameter = name;
    for (std::size_t s = 0; s < samples.size(); ++s) {
      Spectrum sp;
      sp.sample = samples[s].name;
      sp.variable = variable[s];
      sp.reset(edges[s]);
      u.dplus.push_back(sp);
      u.dminus.push_back(sp);
      if (withEnu) {
        Spectrum se;
        se.sample = samples[s].name;
        se.variable = "Enu_rec";
        se.reset(enuEdges[s]);
        u.dplusEnu.push_back(se);
        u.dminusEnu.push_back(se);
      }
    }
    list.push_back(static_cast<int>(universes.size()));
    universes.push_back(u);
  }
};

SystematicsSummary::SystematicsSummary(const Experiment& exp, const ParameterSet& xsec, const ParameterSet& flux,
                                       const ParameterSet& det, const ParameterSet& osc)
    : impl_(new Impl(exp, xsec, flux, det, osc)) {
  Impl& I = *impl_;
  if (I.samples.empty()) throw std::runtime_error("SystematicsSummary: the experiment has no samples");
  for (const auto& s : I.samples) {
    std::string var;
    std::vector<double> e;
    if (!projectedEdges(s.binning, var, e)) {
      var = "Enu_rec";
      e = defaultEnuEdges();
    }
    I.variable.push_back(var);
    I.edges.push_back(e);
    std::vector<double> enu = edgesOf(s.binning, "Enu_rec");
    if (enu.size() < 2) enu = defaultEnuEdges();
    I.enuEdges.push_back(enu);
    Spectrum n;
    n.sample = s.name;
    n.variable = var;
    n.reset(e);
    nominal_.push_back(n);
    Spectrum ne;
    ne.sample = s.name;
    ne.variable = "Enu_rec";
    ne.reset(enu);
    I.nominalEnu.push_back(ne);
  }
  for (const auto& p : I.xsec.parameters()) I.addUniverse("xsec", p.name, I.xsecU, false);
  for (const auto& p : I.flux.parameters()) I.addUniverse("flux", p.name, I.fluxU, false);
  for (const auto& p : I.det.parameters()) I.addUniverse("det", p.name, I.detU, false);
  for (const auto& p : I.osc.parameters()) I.addUniverse("osc", p.name, I.oscU, false);
  for (const auto& sh : exp.shifts()) I.addUniverse("shift", sh.name, I.shiftU, true);
  if (!I.flux.empty() && !I.flux.hasBinning()) throw std::runtime_error("SystematicsSummary: the flux set has no binning");
  if (!I.det.empty() && !I.det.hasBinning()) throw std::runtime_error("SystematicsSummary: the detector set has no binning");
}

SystematicsSummary::~SystematicsSummary() = default;

void SystematicsSummary::process(BeamMode mode, const std::string& file) {
  Impl& I = *impl_;
  if (I.finished) throw std::logic_error("SystematicsSummary::process after finish()");
  std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("SystematicsSummary: cannot open '" + file + "'");
  Provenance::requireCompatible(*f, "summary: " + file);
  TTree* t = f->Get<TTree>("FakeSelect");
  if (!t) throw std::runtime_error("SystematicsSummary: no FakeSelect tree in '" + file + "'");

  // --- branches -------------------------------------------------------------------
  Int_t EventType = 0, flux_pdg = 0, PDGnu = 0, PDGLep_rec = 0;
  Bool_t isRHCb = false, isAppearance = false;
  Char_t cc = 0;
  Float_t Enu_rec = 0, Q2_rec = 0, PLep_rec = 0, CosLep_rec = 0, Theta_rec = 0, Enu_true = 0, Elep_rec = 0, Ehad_rec = 0,
          EhadFrac_rec = 0, Enu_rec_ccqe_mu = 0, Enu_rec_ccqe_e = 0, Enu_rec_calo_mu = 0, Enu_rec_calo_e = 0;
  Double_t POTWeight = 0, OscWeight = 0, SelEffWeight = 0, SampleProb = 1.0;
  auto bind = [&](const char* name, void* addr) {
    if (!t->GetBranch(name)) throw std::runtime_error(std::string("SystematicsSummary: '") + file + "' has no branch " + name);
    t->SetBranchAddress(name, addr);
  };
  bind("EventType", &EventType);
  bind("flux_pdg", &flux_pdg);
  bind("PDGnu", &PDGnu);
  bind("PDGLep_rec", &PDGLep_rec);
  bind("isRHC", &isRHCb);
  bind("isAppearance", &isAppearance);
  bind("cc", &cc);
  bind("Enu_rec", &Enu_rec);
  bind("Q2_rec", &Q2_rec);
  bind("PLep_rec", &PLep_rec);
  bind("CosLep_rec", &CosLep_rec);
  bind("Theta_rec", &Theta_rec);
  bind("Enu_true", &Enu_true);
  bind("Elep_rec", &Elep_rec);
  bind("Ehad_rec", &Ehad_rec);
  bind("EhadFrac_rec", &EhadFrac_rec);
  bind("Enu_rec_ccqe_mu", &Enu_rec_ccqe_mu);
  bind("Enu_rec_ccqe_e", &Enu_rec_ccqe_e);
  bind("Enu_rec_calo_mu", &Enu_rec_calo_mu);
  bind("Enu_rec_calo_e", &Enu_rec_calo_e);
  bind("POTWeight", &POTWeight);
  bind("OscWeight", &OscWeight);
  bind("SelEffWeight", &SelEffWeight);
  bind("SampleProb", &SampleProb);

  std::vector<TGraph*> graphs(I.xsec.size(), nullptr);
  for (std::size_t i = 0; i < I.xsec.size(); ++i) {
    const std::string& name = I.xsec.parameter(static_cast<int>(i)).name;
    TBranch* br = t->GetBranch(name.c_str());
    if (!br || !TClass::GetClass(br->GetClassName()) || !TClass::GetClass(br->GetClassName())->InheritsFrom(TGraph::Class())) {
      throw std::runtime_error("SystematicsSummary: '" + file + "' has no response graph '" + name +
                               "'; run the splines stage and then select again (--force)");
    }
    t->SetBranchAddress(name.c_str(), &graphs[i]);
  }

  // per-variable value lookup for the plotted spectra
  auto valueOf = [&](const std::string& var) -> double {
    if (var == "Enu_rec") return Enu_rec;
    if (var == "PLep_rec") return PLep_rec;
    if (var == "Theta_rec") return Theta_rec;
    if (var == "Q2_rec") return Q2_rec;
    if (var == "CosLep_rec") return CosLep_rec;
    if (var == "Elep_rec") return Elep_rec;
    if (var == "Ehad_rec") return Ehad_rec;
    if (var == "EhadFrac_rec") return EhadFrac_rec;
    if (var == "Enu_true") return Enu_true;
    throw std::runtime_error("SystematicsSummary: cannot bin in '" + var + "' (not a FakeSelect branch)");
  };

  const Oscillator nominalOsc = I.exp.oscillator();
  std::vector<Oscillator> oscPlus, oscMinus;
  for (const auto& p : I.osc.parameters()) {
    oscPlus.emplace_back(OscCovariance::moved(I.exp.oscillation(), p.name, p.sigma), I.exp.baseline_km(), I.exp.rhoYe(), 0);
    oscMinus.emplace_back(OscCovariance::moved(I.exp.oscillation(), p.name, -p.sigma), I.exp.baseline_km(), I.exp.rhoYe(), 0);
  }
  const EnergyEstimator estimator = I.exp.energyEstimator();
  const double Eb = I.exp.bindingEnergy();
  const std::vector<ShiftParameter>& shifts = I.exp.shifts();
  const std::vector<int> sampleIdx = I.exp.sampleIndices(mode);

  const long long n = t->GetEntries();
  for (long long i = 0; i < n; ++i) {
    t->GetEntry(i);
    ++entries_;
    if (EventType < 0 || EventType >= static_cast<int>(I.samples.size())) continue;
    const std::size_t s = static_cast<std::size_t>(EventType);
    if (std::find(sampleIdx.begin(), sampleIdx.end(), EventType) == sampleIdx.end()) continue;
    const double w = POTWeight * OscWeight * SelEffWeight * SampleProb;
    const double x = valueOf(I.variable[s]);
    nominal_[s].fill(x, w);
    I.nominalEnu[s].fill(Enu_rec, w);
    if (w == 0.0) continue;

    // xsec: graph at prior +- sigma
    for (std::size_t k = 0; k < I.xsec.size(); ++k) {
      const TGraph* g = graphs[k];
      if (!g || g->GetN() == 0) continue;
      const Parameter& p = I.xsec.parameter(static_cast<int>(k));
      const double lo = g->GetX()[0], hi = g->GetX()[g->GetN() - 1];
      const double xp = std::min(hi, std::max(lo, p.prior + p.sigma));
      const double xm = std::min(hi, std::max(lo, p.prior - p.sigma));
      Universe& u = I.universes[static_cast<std::size_t>(I.xsecU[k])];
      u.dplus[s].fill(x, w * (g->Eval(xp) - 1.0));
      u.dminus[s].fill(x, w * (g->Eval(xm) - 1.0));
    }
    // flux: the parameter of this (isRHC, flux_pdg, Enu_true) bin
    if (!I.flux.empty()) {
      const int k = I.flux.binning().find({isRHCb ? 1.0 : 0.0, static_cast<double>(flux_pdg), static_cast<double>(Enu_true)});
      if (k >= 0 && k < static_cast<int>(I.flux.size())) {
        const double sig = I.flux.sigma(k);
        Universe& u = I.universes[static_cast<std::size_t>(I.fluxU[static_cast<std::size_t>(k)])];
        u.dplus[s].fill(x, w * sig);
        u.dminus[s].fill(x, -w * sig);
      }
    }
    // detector normalisation of this sample
    if (!I.det.empty()) {
      const int k = I.det.binning().find({static_cast<double>(EventType)});
      if (k >= 0 && k < static_cast<int>(I.det.size())) {
        const double sig = I.det.sigma(k);
        Universe& u = I.universes[static_cast<std::size_t>(I.detU[static_cast<std::size_t>(k)])];
        u.dplus[s].fill(x, w * sig);
        u.dminus[s].fill(x, -w * sig);
      }
    }
    // oscillation: CC events re-weighted by P'/P
    if (!I.osc.empty() && cc && OscWeight > 0.0) {
      const double P = nominalOsc.probability(flux_pdg, PDGnu, Enu_true);
      if (P > 0.0) {
        for (std::size_t k = 0; k < I.osc.size(); ++k) {
          Universe& u = I.universes[static_cast<std::size_t>(I.oscU[k])];
          u.dplus[s].fill(x, w * (oscPlus[k].probability(flux_pdg, PDGnu, Enu_true) / P - 1.0));
          u.dminus[s].fill(x, w * (oscMinus[k].probability(flux_pdg, PDGnu, Enu_true) / P - 1.0));
        }
      }
    }
    // shifts: re-derive the reconstruction and move the event between bins
    if (!shifts.empty()) {
      RecoObservables ro;
      ro.PLep_rec = PLep_rec;
      ro.CosLep_rec = CosLep_rec;
      ro.Theta_rec = Theta_rec;
      ro.Ehad_rec = Ehad_rec;
      ro.Elep_rec_mu = Elep_rec;
      ro.Elep_rec_e = Elep_rec;
      ro.Enu_rec_ccqe_mu = Enu_rec_ccqe_mu;
      ro.Enu_rec_ccqe_e = Enu_rec_ccqe_e;
      ro.Enu_rec_calo_mu = Enu_rec_calo_mu;
      ro.Enu_rec_calo_e = Enu_rec_calo_e;
      const Hypothesis hyp = I.samples[s].hypothesis;
      auto shiftedValue = [&](const RecoObservables& r, const std::string& var) -> double {
        if (var == "Enu_rec") return r.enuRec(hyp, estimator, PDGLep_rec);
        if (var == "Q2_rec") return r.q2Rec(hyp, PDGLep_rec);
        if (var == "PLep_rec") return r.PLep_rec;
        if (var == "Theta_rec") return r.Theta_rec;
        if (var == "CosLep_rec") return r.CosLep_rec;
        if (var == "Elep_rec") return r.elepRec(hyp, PDGLep_rec);
        if (var == "Ehad_rec") return r.Ehad_rec;
        if (var == "EhadFrac_rec") {
          const double e = r.enuRec(hyp, estimator, PDGLep_rec);
          return e > 0.0 ? r.Ehad_rec / e : 0.0;
        }
        return valueOf(var);
      };
      for (std::size_t k = 0; k < shifts.size(); ++k) {
        Universe& u = I.universes[static_cast<std::size_t>(I.shiftU[k])];
        for (int sign = -1; sign <= 1; sign += 2) {
          RecoObservables r = ro;
          ShiftInjector::apply(r, shifts[k], sign, isRHCb, Eb, cc != 0, hyp, PDGLep_rec);
          const double xs = shiftedValue(r, I.variable[s]);
          const double es = r.enuRec(hyp, estimator, PDGLep_rec);
          Spectrum& d = sign > 0 ? u.dplus[s] : u.dminus[s];
          Spectrum& de = sign > 0 ? u.dplusEnu[s] : u.dminusEnu[s];
          d.fill(xs, w);
          d.fill(x, -w);
          de.fill(es, w);
          de.fill(Enu_rec, -w);
        }
      }
    }
  }
  t->ResetBranchAddresses();
}

SummaryTable SystematicsSummary::finish() {
  Impl& I = *impl_;
  table_ = SummaryTable();
  for (const auto& u : I.universes) {
    for (std::size_t s = 0; s < I.samples.size(); ++s) {
      SummaryEntry e;
      e.set = u.set;
      e.parameter = u.parameter;
      e.sample = I.samples[s].name;
      e.nominal = nominal_[s].total();
      e.plus = e.nominal + u.dplus[s].total();
      e.minus = e.nominal + u.dminus[s].total();
      table_.add(e);
    }
  }
  I.finished = true;
  return table_;
}

Spectrum SystematicsSummary::spectrum(const std::string& parameter, int sample, bool plus) const {
  const Impl& I = *impl_;
  if (sample < 0 || sample >= static_cast<int>(I.samples.size())) throw std::out_of_range("SystematicsSummary::spectrum: sample index");
  for (const auto& u : I.universes) {
    if (u.parameter != parameter) continue;
    Spectrum out = nominal_[static_cast<std::size_t>(sample)];
    const Spectrum& d = plus ? u.dplus[static_cast<std::size_t>(sample)] : u.dminus[static_cast<std::size_t>(sample)];
    for (std::size_t b = 0; b < out.contents.size(); ++b) out.contents[b] += d.contents[b];
    return out;
  }
  throw std::out_of_range("SystematicsSummary::spectrum: no parameter '" + parameter + "'");
}

void SystematicsSummary::writeHistograms(const std::string& file) const {
  const Impl& I = *impl_;
  std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "RECREATE"));
  if (!f || f->IsZombie()) throw std::runtime_error("SystematicsSummary::writeHistograms: cannot create '" + file + "'");
  for (std::size_t s = 0; s < I.samples.size(); ++s) {
    TDirectory* d = f->mkdir(I.samples[s].name.c_str());
    std::unique_ptr<TH1D> h(nominal_[s].histogram("nominal", I.samples[s].title + " nominal"));
    d->WriteTObject(h.get(), "nominal", "Overwrite");
    for (const auto& u : I.universes) {
      for (int sign = -1; sign <= 1; sign += 2) {
        const Spectrum sp = spectrum(u.parameter, static_cast<int>(s), sign > 0);
        const std::string name = u.parameter + (sign > 0 ? "_plus" : "_minus");
        std::unique_ptr<TH1D> hh(sp.histogram(name, I.samples[s].title + " " + u.parameter + (sign > 0 ? " +1#sigma" : " -1#sigma")));
        d->WriteTObject(hh.get(), name.c_str(), "Overwrite");
      }
    }
  }
  std::ostringstream prov;
  prov << "experiment " << I.exp.name() << "\n" << I.universes.size() << " parameters, " << I.samples.size()
       << " samples, " << entries_ << " FakeSelect entries";
  Provenance::stamp(*f, "summary", prov.str());
  f->Close();
}

ParameterSet SystematicsSummary::shiftCovariance() const {
  const Impl& I = *impl_;
  ParameterSet s("det");
  s.setTitle("detector covariance from the shift systematics (fractional, per sample and Enu_rec bin)");
  if (I.shiftU.empty()) throw std::runtime_error("SystematicsSummary::shiftCovariance: the experiment has no shift parameters");
  // rows
  std::vector<std::pair<std::size_t, std::size_t>> rows;  // (sample, bin)
  std::ostringstream txt;
  txt << "variables: EventType Enu_rec Enu_rec\n";
  for (std::size_t smp = 0; smp < I.samples.size(); ++smp) {
    const Spectrum& nom = I.nominalEnu[smp];
    for (int b = 0; b < nom.nBins(); ++b) {
      Parameter p;
      std::ostringstream name;
      name << DetectorCovariance::defaultName(I.samples[smp].name) << "_" << std::setw(2) << std::setfill('0') << b;
      p.name = name.str();
      p.prior = 1.0;
      p.sigma = 0.0;
      s.add(p);
      rows.emplace_back(smp, static_cast<std::size_t>(b));
      txt << smp << " " << std::setprecision(10) << nom.edges[static_cast<std::size_t>(b)] << " "
          << nom.edges[static_cast<std::size_t>(b) + 1] << "\n";
    }
  }
  const int n = static_cast<int>(rows.size());
  TMatrixD cov(n, n);
  for (int ui : I.shiftU) {
    const Universe& u = I.universes[static_cast<std::size_t>(ui)];
    std::vector<double> delta(static_cast<std::size_t>(n), 0.0);
    for (int r = 0; r < n; ++r) {
      const std::size_t smp = rows[static_cast<std::size_t>(r)].first, b = rows[static_cast<std::size_t>(r)].second;
      const double N = I.nominalEnu[smp].contents[b];
      if (N <= 0.0) continue;
      delta[static_cast<std::size_t>(r)] = (u.dplusEnu[smp].contents[b] - u.dminusEnu[smp].contents[b]) / (2.0 * N);
    }
    for (int r = 0; r < n; ++r) {
      for (int c = 0; c < n; ++c) cov(r, c) += delta[static_cast<std::size_t>(r)] * delta[static_cast<std::size_t>(c)];
    }
  }
  // a rank-deficient matrix is not invertible: add a floor so GUNDAM can use it
  for (int r = 0; r < n; ++r) cov(r, r) += 1e-8;
  s.setCovariance(cov, 1e-9);
  s.setBinning(Binning::fromText(txt.str(), "shift detector binning"));
  return s;
}

}  // namespace zaku
