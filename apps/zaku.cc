// zaku - command-line front end to the Zaku library.
//
// Usage: zaku <experiment.yaml> <stage> [--sample <name>] [--force]
//        zaku --version | --help
//
// Stages are added as the library grows; the ones not available yet in this
// build are listed as such by --help.
#include <cstdio>
#include <cstring>
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
    {"prepare", "write generation, weight and flatten inputs + manifests", false},
    {"prepare-generation", "write generator cards/seed files + generate manifest", false},
    {"prepare-weights", "write dial configs + weight manifest", false},
    {"prepare-flatten", "write the flatten manifest", false},
    {"recon", "apply the response, smearing and energy estimators", false},
    {"splines", "build per-event systematic splines", false},
    {"select", "assign samples, write the FakeSelect tree and rate table", false},
    {"covariances", "write the parameter-set files", false},
    {"summary", "summarise the systematics", false},
    {"gundam-config", "write the GUNDAM configuration tree + fit manifest", false},
    {"status", "report which inputs and outputs exist", false},
    {"list-samples", "print the MC sample names in index order", false},
    {"list-dials", "print the dial catalogue", false},
    {"list-models", "print the NEUT model choices", false},
    {"list-tunes", "print the GENIE tune catalogue", false},
    {"show-response", "print the response tables", false},
};

void printUsage() {
  std::printf("zaku %s\n\n", zaku::buildInfo().c_str());
  std::printf("Usage:\n");
  std::printf("  zaku <experiment.yaml> <stage> [--sample <name>] [--force]\n");
  std::printf("  zaku --version | --help\n\n");
  std::printf("Stages:\n");
  for (const auto& s : kStages) {
    std::printf("  %-20s %s%s\n", s.name, s.summary, s.available ? "" : "  [not in this build]");
  }
  std::printf("\nZaku never runs NEUT, GENIE, T2KReWeight, NUISANCE or GUNDAM itself: the prepare\n");
  std::printf("stages write their inputs and a manifest of command lines for you to run.\n");
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
  const std::string& stage = args[1];
  for (const auto& s : kStages) {
    if (stage == s.name) {
      if (!s.available) {
        std::fprintf(stderr, "zaku: stage '%s' is not available in this build (Zaku %s)\n", s.name,
                     zaku::version().c_str());
        return 2;
      }
    }
  }
  std::fprintf(stderr, "zaku: unknown stage '%s'\n", stage.c_str());
  printUsage();
  return 1;
}
