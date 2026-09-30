// zaku - command-line front end to the Zaku library.
//
// Usage: zaku <experiment.yaml> <stage> [--production <name>] [--force]
//        zaku --version | --help
//
// Stages are added as the library grows; --help lists the ones this build has.
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "zaku/Zaku.h"

namespace {

struct StageInfo {
  const char* name;
  const char* summary;
  bool available;
};

const StageInfo kStages[] = {
    {"describe", "print the experiment as this build understands it", true},
    {"validate", "check that the experiment is complete and its inputs exist", true},
    {"status", "report which production files exist", true},
    {"list-productions", "print the production names in index order", true},
    {"list-beam-modes", "print the beam modes with POT and flux integrals", true},
    {"prepare", "write generation, weight and flatten inputs + manifests", false},
    {"prepare-generation", "write generator cards/seed files + generate manifest", false},
    {"prepare-weights", "write dial configs + weight manifest", false},
    {"prepare-flatten", "write the flatten manifest", false},
    {"recon", "apply the response, smearing and energy estimators", true},
    {"splines", "build per-event systematic splines", false},
    {"select", "assign samples, write the FakeSelect tree and rate table", true},
    {"covariances", "write the parameter-set files", false},
    {"summary", "summarise the systematics", false},
    {"gundam-config", "write the GUNDAM configuration tree + fit manifest", false},
    {"list-dials", "print the dial catalogue", false},
    {"list-models", "print the NEUT model choices", false},
    {"list-tunes", "print the GENIE tune catalogue", false},
    {"show-response", "print the response tables (or write them as YAML with --yaml)", true},
    {"show-samples", "print the analysis samples", true},
};

struct Options {
  std::string config;
  std::string stage;
  std::string production;
  bool force = false;
  bool yaml = false;
};

void printUsage() {
  std::printf("zaku %s\n\n", zaku::buildInfo().c_str());
  std::printf("Usage:\n");
  std::printf("  zaku <experiment.yaml> <stage> [--production <name>] [--force] [--yaml]\n");
  std::printf("  zaku --version | --help\n\n");
  std::printf("Stages:\n");
  for (const auto& s : kStages) {
    std::printf("  %-20s %s%s\n", s.name, s.summary, s.available ? "" : "  [not in this build]");
  }
  std::printf("\nZaku never runs NEUT, GENIE, T2KReWeight, NUISANCE or GUNDAM itself: the prepare\n");
  std::printf("stages write their inputs and a manifest of command lines for you to run.\n");
}

int run(const Options& opt) {
  zaku::Experiment exp = zaku::Experiment::fromConfig(opt.config);
  const std::string& stage = opt.stage;
  if (stage == "describe") {
    std::printf("%s", exp.describe().c_str());
    return 0;
  }
  if (stage == "validate") {
    exp.validate();
    std::printf("%s: complete\n", exp.name().c_str());
    return 0;
  }
  if (stage == "status") {
    std::printf("%s", exp.status().c_str());
    return 0;
  }
  if (stage == "list-productions") {
    for (const auto& s : exp.productions()) {
      if (!opt.production.empty() && opt.production != s.name()) continue;
      std::printf("%s\n", s.name().c_str());
    }
    return 0;
  }
  if (stage == "list-beam-modes") {
    for (zaku::BeamMode mode : exp.beamModes()) {
      std::printf("%s: POT %.4e\n", zaku::toString(mode).c_str(), exp.pot(mode));
      for (int pdg : exp.flux(mode).pdgs()) {
        std::printf("  %-8s %.6e /cm^2/POT\n", zaku::NuSpecies::fromPdg(pdg).name().c_str(),
                    exp.fluxIntegralPerPOT(mode, pdg));
      }
    }
    return 0;
  }
  if (stage == "show-response") {
    if (opt.yaml) std::printf("%s", exp.response().toYaml().c_str());
    else std::printf("%s", exp.response().describe().c_str());
    return 0;
  }
  if (stage == "show-samples") {
    const auto& samples = exp.samples();
    for (std::size_t i = 0; i < samples.size(); ++i) std::printf("[%zu] %s\n", i, samples[i].describe().c_str());
    return 0;
  }
  if (stage == "recon") {
    exp.reconstruct(opt.production, opt.force);
    return 0;
  }
  if (stage == "select") {
    exp.select(opt.force);
    return 0;
  }
  for (const auto& s : kStages) {
    if (stage == s.name && !s.available) {
      std::fprintf(stderr, "zaku: stage '%s' is not available in this build (Zaku %s)\n", s.name,
                   zaku::version().c_str());
      return 2;
    }
  }
  std::fprintf(stderr, "zaku: unknown stage '%s'\n", stage.c_str());
  printUsage();
  return 1;
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> args(argv + 1, argv + argc);
  if (args.empty() || args[0] == "--help" || args[0] == "-h") {
    printUsage();
    return args.empty() ? 1 : 0;
  }
  if (args[0] == "--version" || args[0] == "-v") {
    std::printf("%s\n", zaku::version().c_str());
    return 0;
  }
  if (args.size() < 2) {
    std::fprintf(stderr, "zaku: expected <experiment.yaml> <stage>\n");
    printUsage();
    return 1;
  }
  Options opt;
  opt.config = args[0];
  opt.stage = args[1];
  for (std::size_t i = 2; i < args.size(); ++i) {
    if (args[i] == "--force") {
      opt.force = true;
    } else if (args[i] == "--yaml") {
      opt.yaml = true;
    } else if (args[i] == "--production" && i + 1 < args.size()) {
      opt.production = args[++i];
    } else {
      std::fprintf(stderr, "zaku: unknown option '%s'\n", args[i].c_str());
      return 1;
    }
  }
  try {
    return run(opt);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "zaku: %s\n", e.what());
    return 3;
  }
}
