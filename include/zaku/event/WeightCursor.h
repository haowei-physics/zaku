#ifndef ZAKU_EVENT_WEIGHTCURSOR_H
#define ZAKU_EVENT_WEIGHTCURSOR_H
#include <memory>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace zaku {

/// The knots of one systematic dial: the parameter values at which weights
/// were evaluated, and the value the events were generated at.
struct DialKnots {
  std::string name;
  std::vector<double> knots;  ///< absolute parameter values, in the order of the weight vectors
  double nominal = 0.0;       ///< generated value ("from")
  /// True if the nominal coincides with one of the knots (within 1e-9 relative).
  bool nominalIsKnot() const;
  /// Knots plus the nominal (if absent), sorted ascending.
  std::vector<double> knotsWithNominal() const;
  std::string describe() const;
};

/// A source of per-event systematic weights: for every row, one weight per
/// knot per dial. Rows are matched to flat-tree entries by EventJoin.
class WeightSource {
public:
  virtual ~WeightSource() = default;
  virtual const std::string& fileName() const = 0;
  virtual const std::vector<DialKnots>& dials() const = 0;
  std::vector<std::string> dialNames() const;
  /// Index of a dial or -1.
  int dialIndex(const std::string& name) const;

  virtual long long entries() const = 0;
  virtual long long entry() const = 0;
  virtual bool next() = 0;
  virtual bool goTo(long long entry) = 0;
  /// Position before the first row, so that next() loads row 0.
  virtual void rewind() = 0;

  /// Weights of the current row for dial `i`, one per knot.
  virtual const std::vector<double>& weights(int i) const = 0;
  const std::vector<double>& weights(const std::string& dial) const;
  /// Nominal weight of the row (1 for a consistent reweighting; 0 marks events the engine could not handle).
  virtual double cvWeight() const = 0;

  /// Key used to match the row to a flat-tree event.
  virtual int Mode() const = 0;
  virtual int PDGnu() const = 0;
  virtual double Enu() const = 0;   ///< GeV
  /// Diagnostics only (may disagree with the flat tree in multi-lepton events).
  virtual int PDGLep() const = 0;
  virtual double ELep() const = 0;
  /// Explicit event number when the format has one (GENIE), else -1.
  virtual long long eventNumber() const { return -1; }
};

/// Reader for a T2KReWeight `EventWeightGen` output (`*_Weight.flat.root`):
/// tree `T2KRWFlatTree` with one `vector<double>` per dial plus the NIWG
/// kinematics, and tree `TweakMetaData` with the knot values (`<dial>`) and
/// the generated value (`<dial>_from`). The file holds several cycles of the
/// trees; the highest is used.
class T2KRWWeightCursor : public WeightSource {
public:
  explicit T2KRWWeightCursor(const std::string& file);
  ~T2KRWWeightCursor() override;
  T2KRWWeightCursor(const T2KRWWeightCursor&) = delete;
  T2KRWWeightCursor& operator=(const T2KRWWeightCursor&) = delete;

  const std::string& fileName() const override { return fileName_; }
  const std::vector<DialKnots>& dials() const override { return dials_; }
  long long entries() const override;
  long long entry() const override { return entry_; }
  bool next() override;
  bool goTo(long long entry) override;
  void rewind() override { entry_ = -1; }
  using WeightSource::weights;
  const std::vector<double>& weights(int i) const override;
  double cvWeight() const override;
  int Mode() const override;
  int PDGnu() const override;
  double Enu() const override;
  int PDGLep() const override;
  double ELep() const override;

  /// The per-event cross-section scale stored by T2KReWeight (same as fScaleFactor).
  double scaleWeight() const;

  static bool looksLikeT2KRW(const std::string& file);

private:
  struct Buffers;
  std::string fileName_;
  std::unique_ptr<TFile> file_;
  TTree* tree_ = nullptr;
  std::unique_ptr<Buffers> buf_;
  std::vector<DialKnots> dials_;
  long long entry_ = -1;
};

}  // namespace zaku

#endif  // ZAKU_EVENT_WEIGHTCURSOR_H
