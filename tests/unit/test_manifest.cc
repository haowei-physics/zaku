#include <fstream>

#include "ZakuTest.h"
#include "ZakuTestPaths.h"
#include "zaku/core/Manifest.h"

using namespace zaku;

int main() {
  const std::string scratch = zakutest::scratchDir();

  Manifest m("generate");
  m.addPreamble("source ~/setup_neut.sh");
  ManifestJob& a = m.add("numu_x_numu_2K_FHC", "neutroot2 card.toml out.root");
  a.workDir = scratch + "/gen/numu";
  a.outputs = {scratch + "/gen/numu/out.root"};
  a.comment = "2000 events";
  ManifestJob& b = m.add("done_job", "true");
  b.outputs = {scratch};  // exists

  ZAKU_CHECK_EQ(static_cast<int>(m.size()), 2);
  ZAKU_CHECK_EQ(m.nDone(), 1);
  ZAKU_CHECK_EQ(static_cast<int>(m.pending().size()), 1);
  ZAKU_CHECK_STR_EQ(m.pending()[0]->name, "numu_x_numu_2K_FHC");
  ZAKU_CHECK_STR_EQ(m.summary(), "generate: 2 jobs, 1 done");

  const std::string text = m.text();
  ZAKU_CHECK(text.find("# job: numu_x_numu_2K_FHC  (2000 events)") != std::string::npos);
  ZAKU_CHECK(text.find("cd " + scratch + "/gen/numu && neutroot2 card.toml out.root") != std::string::npos);
  ZAKU_CHECK(text.find("# env: source ~/setup_neut.sh") != std::string::npos);
  // exactly two command lines (non-comment)
  int commands = 0;
  {
    std::size_t pos = 0;
    while (pos < text.size()) {
      const std::size_t end = text.find('\n', pos);
      const std::string line = text.substr(pos, end == std::string::npos ? std::string::npos : end - pos);
      if (!line.empty() && line[0] != '#') ++commands;
      if (end == std::string::npos) break;
      pos = end + 1;
    }
  }
  ZAKU_CHECK_EQ(commands, 2);

  const std::string json = m.json();
  ZAKU_CHECK(json.find("\"stage\": \"generate\"") != std::string::npos);
  ZAKU_CHECK(json.find("\"name\": \"done_job\"") != std::string::npos);

  const std::string txt = m.write(scratch + "/manifests");
  ZAKU_CHECK_STR_EQ(txt, scratch + "/manifests/generate.txt");
  ZAKU_CHECK(zakutest::exists(scratch + "/manifests/generate.json"));
  {
    std::ifstream f(txt);
    std::string all((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    ZAKU_CHECK_STR_EQ(all, text);
  }
  ZAKU_CHECK_THROWS(Manifest("").write(scratch));

  ZAKU_CHECK_STR_EQ(Manifest::quote("plain/path.root"), "plain/path.root");
  ZAKU_CHECK_STR_EQ(Manifest::quote("with space"), "'with space'");
  ZAKU_CHECK_STR_EQ(Manifest::quote("it's"), "'it'\\''s'");
  ZAKU_CHECK_STR_EQ(Manifest::quote(""), "''");

  return zakutest::finish();
}
