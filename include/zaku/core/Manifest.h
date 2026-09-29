#ifndef ZAKU_CORE_MANIFEST_H
#define ZAKU_CORE_MANIFEST_H
#include <string>
#include <vector>

namespace zaku {

/// One external job the user has to run: a command line, the directory to run
/// it in, and the files it is expected to produce.
struct ManifestJob {
  std::string name;
  std::string command;
  std::string workDir;                ///< empty = the run work directory
  std::vector<std::string> outputs;   ///< absolute paths
  std::string comment;                ///< free text shown next to the job

  std::vector<std::string> missingOutputs() const;
  bool done() const { return missingOutputs().empty(); }
};

/// The list of jobs a prepare stage hands to the user. Zaku never runs them.
///
/// `write()` produces two files in a directory: `<stage>.txt`, one job per
/// line (preceded by a `# name` comment line) so it can be fed to `xargs`,
/// a shell loop or a batch array; and `<stage>.json` with the same content
/// plus the expected outputs, for tooling.
class Manifest {
public:
  explicit Manifest(std::string stage = "");

  const std::string& stage() const { return stage_; }
  void setStage(const std::string& s) { stage_ = s; }

  void add(const ManifestJob& job);
  ManifestJob& add(const std::string& name, const std::string& command);
  const std::vector<ManifestJob>& jobs() const { return jobs_; }
  std::size_t size() const { return jobs_.size(); }
  bool empty() const { return jobs_.empty(); }

  /// Environment lines written at the top of the text manifest (e.g. "source ~/setup_neut.sh").
  void addPreamble(const std::string& line) { preamble_.push_back(line); }
  const std::vector<std::string>& preamble() const { return preamble_; }

  /// Jobs whose outputs do not all exist yet.
  std::vector<const ManifestJob*> pending() const;
  int nDone() const;

  /// The text manifest as a string (what write() puts in <stage>.txt).
  std::string text() const;
  std::string json() const;
  /// Write <dir>/<stage>.txt and <dir>/<stage>.json; returns the .txt path.
  std::string write(const std::string& dir) const;

  /// One-line summary such as "generate: 12 jobs, 3 done".
  std::string summary() const;

  /// Shell-quote an argument if needed.
  static std::string quote(const std::string& arg);

private:
  std::string stage_;
  std::vector<std::string> preamble_;
  std::vector<ManifestJob> jobs_;
};

}  // namespace zaku

#endif  // ZAKU_CORE_MANIFEST_H
