// ROOT dictionary selection for libZaku.
//
// Everything in namespace zaku is exposed to Cling and PyROOT. Plain value
// structs are linked with '+' so they can be stored in ROOT files and used as
// tree branches; behavioural classes get no streamer.
#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;

#pragma link C++ namespace zaku;
#pragma link C++ namespace zaku::units;

// core/Types.h
#pragma link C++ enum zaku::Flavour;
#pragma link C++ enum zaku::BeamMode;
#pragma link C++ struct zaku::NuSpecies+;
#pragma link C++ struct zaku::Target+;
#pragma link C++ class std::vector<zaku::Target>+;
#pragma link C++ class std::vector<zaku::BeamMode>+;

// core/Units.h
#pragma link C++ enum zaku::EnergyUnit;

// core/Detector.h
#pragma link C++ struct zaku::NeutTarget+;
#pragma link C++ class zaku::Detector+;

#endif
