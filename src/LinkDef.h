// ROOT dictionary selection for libZaku.
//
// Everything in namespace zaku is exposed to Cling and PyROOT. Plain value
// structs are linked with '+' so they can be stored in ROOT files and used as
// tree branches; behavioural classes (cursors, the Experiment facade) get no
// streamer ('-').
#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;

#pragma link C++ namespace zaku;
#pragma link C++ namespace zaku::units;
#pragma link C++ namespace zaku::Layout;
#pragma link C++ namespace zaku::FluxPreparer;

// Enums are deliberately not listed: an explicit '#pragma link C++ enum' puts
// them in the rootmap, and ROOT then autoparses the whole dictionary payload
// the moment a user header declares the enum, before the library is loaded,
// which fails on the not-yet-complete headers. The enums are still fully
// usable from Cling and PyROOT through the headers.

// core/Types.h
#pragma link C++ struct zaku::NuSpecies+;
#pragma link C++ struct zaku::Target+;
#pragma link C++ class std::vector<zaku::Target>+;
#pragma link C++ class std::vector<zaku::BeamMode>+;

// core/Detector.h
#pragma link C++ struct zaku::NeutTarget+;
#pragma link C++ class zaku::Detector+;

// core/FluxSpec.h
#pragma link C++ struct zaku::FluxSpec+;
#pragma link C++ struct zaku::FluxCovarianceSpec+;

// core/Production.h
#pragma link C++ struct zaku::Production+;
#pragma link C++ class std::vector<zaku::Production>+;

// core/Manifest.h
#pragma link C++ struct zaku::ManifestJob+;
#pragma link C++ class std::vector<zaku::ManifestJob>+;
#pragma link C++ class zaku::Manifest+;

// core/Provenance.h, core/RunSettings.h
#pragma link C++ struct zaku::Provenance+;
#pragma link C++ struct zaku::RunSettings+;

// core/Experiment.h
#pragma link C++ class zaku::Experiment-;

// flux/
#pragma link C++ class zaku::FluxHistograms-;
#pragma link C++ struct zaku::FluxCovarianceBlock+;
#pragma link C++ class std::vector<zaku::FluxCovarianceBlock>+;
#pragma link C++ class zaku::FluxCovarianceInput-;

// event/
#pragma link C++ struct zaku::FlatEvent+;
#pragma link C++ struct zaku::ProductionContext+;
#pragma link C++ struct zaku::TruthVariables+;
#pragma link C++ class zaku::FlatTreeCursor-;
#pragma link C++ struct zaku::DialKnots+;
#pragma link C++ class std::vector<zaku::DialKnots>+;
#pragma link C++ class zaku::WeightSource-;
#pragma link C++ class zaku::T2KRWWeightCursor-;
#pragma link C++ struct zaku::JoinOptions+;
#pragma link C++ struct zaku::JoinStats+;
#pragma link C++ class zaku::EventJoin-;

#endif
