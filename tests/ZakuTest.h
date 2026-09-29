// Minimal self-contained test harness (no external framework needed).
//
//   #include "ZakuTest.h"
//   int main() {
//     ZAKU_CHECK(1 + 1 == 2);
//     ZAKU_CHECK_CLOSE(3.14159, 3.1416, 1e-4);   // relative tolerance
//     ZAKU_CHECK_THROWS(throwingCall());
//     return zakutest::finish();
//   }
#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace zakutest {

inline int& failures() {
  static int n = 0;
  return n;
}
inline int& checks() {
  static int n = 0;
  return n;
}

inline void report(bool ok, const char* what, const char* file, int line, const std::string& detail = "") {
  ++checks();
  if (ok) return;
  ++failures();
  std::fprintf(stderr, "FAIL %s:%d: %s%s%s\n", file, line, what, detail.empty() ? "" : " -- ", detail.c_str());
}

inline bool close(double a, double b, double relTol) {
  const double scale = std::fmax(std::fabs(a), std::fabs(b));
  if (scale == 0.0) return true;
  return std::fabs(a - b) <= relTol * scale;
}

inline int finish() {
  std::printf("%d checks, %d failures\n", checks(), failures());
  return failures() == 0 ? 0 : 1;
}

/// Skip the whole test (exit code 77 is CTest's SKIP_RETURN_CODE convention).
inline void skip(const std::string& why) {
  std::printf("SKIPPED: %s\n", why.c_str());
  std::exit(77);
}

}  // namespace zakutest

#define ZAKU_CHECK(expr) ::zakutest::report(static_cast<bool>(expr), #expr, __FILE__, __LINE__)

#define ZAKU_CHECK_MSG(expr, msg) ::zakutest::report(static_cast<bool>(expr), #expr, __FILE__, __LINE__, (msg))

#define ZAKU_CHECK_EQ(a, b)                                                                            \
  do {                                                                                                 \
    const auto _za = (a);                                                                              \
    const auto _zb = (b);                                                                              \
    ::zakutest::report(_za == _zb, #a " == " #b, __FILE__, __LINE__,                                   \
                       "got " + std::to_string(_za) + " vs " + std::to_string(_zb));                   \
  } while (0)

#define ZAKU_CHECK_STR_EQ(a, b)                                                                        \
  do {                                                                                                 \
    const std::string _za = (a);                                                                       \
    const std::string _zb = (b);                                                                       \
    ::zakutest::report(_za == _zb, #a " == " #b, __FILE__, __LINE__, "got '" + _za + "' vs '" + _zb + "'"); \
  } while (0)

#define ZAKU_CHECK_CLOSE(a, b, tol)                                                                    \
  do {                                                                                                 \
    const double _za = (a);                                                                            \
    const double _zb = (b);                                                                            \
    ::zakutest::report(::zakutest::close(_za, _zb, (tol)), #a " ~ " #b, __FILE__, __LINE__,            \
                       "got " + std::to_string(_za) + " vs " + std::to_string(_zb));                   \
  } while (0)

#define ZAKU_CHECK_THROWS(stmt)                                                                        \
  do {                                                                                                 \
    bool _thrown = false;                                                                              \
    try {                                                                                              \
      stmt;                                                                                            \
    } catch (...) {                                                                                    \
      _thrown = true;                                                                                  \
    }                                                                                                  \
    ::zakutest::report(_thrown, "throws: " #stmt, __FILE__, __LINE__);                                 \
  } while (0)

#define ZAKU_CHECK_NOTHROW(stmt)                                                                       \
  do {                                                                                                 \
    bool _thrown = false;                                                                              \
    std::string _what;                                                                                 \
    try {                                                                                              \
      stmt;                                                                                            \
    } catch (const std::exception& e) {                                                                \
      _thrown = true;                                                                                  \
      _what = e.what();                                                                                \
    } catch (...) {                                                                                    \
      _thrown = true;                                                                                  \
    }                                                                                                  \
    ::zakutest::report(!_thrown, "no throw: " #stmt, __FILE__, __LINE__, _what);                       \
  } while (0)
