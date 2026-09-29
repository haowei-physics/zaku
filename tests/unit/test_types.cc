#include "ZakuTest.h"
#include "zaku/core/Types.h"

int main() {
  using namespace zaku;

  // NuSpecies
  NuSpecies numu;
  ZAKU_CHECK_EQ(numu.pdg(), 14);
  ZAKU_CHECK_STR_EQ(numu.name(), "numu");
  ZAKU_CHECK_STR_EQ(numu.t2kSuffix(), "numu");
  ZAKU_CHECK_STR_EQ(numu.latex(), "#nu_{#mu}");

  const NuSpecies nuebar = NuSpecies::fromPdg(-12);
  ZAKU_CHECK(nuebar.flavour == Flavour::NuE);
  ZAKU_CHECK(nuebar.isAnti);
  ZAKU_CHECK_STR_EQ(nuebar.name(), "nuebar");
  ZAKU_CHECK_STR_EQ(nuebar.t2kSuffix(), "nueb");
  ZAKU_CHECK_STR_EQ(nuebar.latex(), "#bar{#nu}_{e}");
  ZAKU_CHECK_THROWS(NuSpecies::fromPdg(13));

  ZAKU_CHECK_EQ(NuSpecies::fromName("numubar").pdg(), -14);
  ZAKU_CHECK_EQ(NuSpecies::fromName("numub").pdg(), -14);
  ZAKU_CHECK_EQ(NuSpecies::fromName("anti-numu").pdg(), -14);
  ZAKU_CHECK_EQ(NuSpecies::fromName("nu_mu_bar").pdg(), -14);
  ZAKU_CHECK_EQ(NuSpecies::fromName("NuE").pdg(), 12);
  ZAKU_CHECK_EQ(NuSpecies::fromName("nueb").pdg(), -12);
  ZAKU_CHECK_EQ(NuSpecies::fromName("nutau").pdg(), 16);
  ZAKU_CHECK_EQ(NuSpecies::fromName("-16").pdg(), -16);
  ZAKU_CHECK_EQ(NuSpecies::fromName("14").pdg(), 14);
  ZAKU_CHECK_THROWS(NuSpecies::fromName("proton"));
  ZAKU_CHECK(NuSpecies::fromName("numu") == numu);
  ZAKU_CHECK(NuSpecies::fromName("numubar") != numu);

  // BeamMode
  ZAKU_CHECK(beamModeFromString("FHC") == BeamMode::FHC);
  ZAKU_CHECK(beamModeFromString("rhc") == BeamMode::RHC);
  ZAKU_CHECK(beamModeFromString("anumode") == BeamMode::RHC);
  ZAKU_CHECK(beamModeFromString("plus250kA") == BeamMode::FHC);
  ZAKU_CHECK_THROWS(beamModeFromString("sideways"));
  ZAKU_CHECK_STR_EQ(toString(BeamMode::RHC), "RHC");
  ZAKU_CHECK(isRHC(BeamMode::RHC));
  ZAKU_CHECK(!isRHC(BeamMode::FHC));
  ZAKU_CHECK_EQ(static_cast<int>(allBeamModes().size()), 2);

  // Target
  const Target o16(8, 16, 16.0 / 18.0);
  ZAKU_CHECK_EQ(o16.nuclearPdg(), 1000080160);
  ZAKU_CHECK_STR_EQ(o16.name(), "O16");
  ZAKU_CHECK_STR_EQ(Target(1, 1, 1.0).name(), "H1");
  ZAKU_CHECK_STR_EQ(Target(18, 40, 1.0).name(), "Ar40");
  ZAKU_CHECK_EQ(Target(18, 40, 1.0).nuclearPdg(), 1000180400);
  ZAKU_CHECK_STR_EQ(elementSymbol(26), "Fe");
  ZAKU_CHECK_STR_EQ(elementSymbol(200), "Z200");

  return zakutest::finish();
}
