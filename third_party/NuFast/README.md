# NuFast (vendored)

`NuFast_LBL.h` / `NuFast_LBL.cc` are P. Denton's NuFast long-baseline
oscillation probability code (https://github.com/PeterDenton/NuFast, MIT
licence, see `LICENSE`), taken from the copy used by the legacy
`mock_lbl_inputs` macros. The only changes are:

- the demonstration `main()` was removed so the code can be linked into a library;
- the declarations were split into a header;
- the file-scope `using std::vector, std::array;` became local to the namespace.

The arithmetic is untouched, so `zaku::Oscillator` reproduces the legacy
`OscWeight` bit for bit.
