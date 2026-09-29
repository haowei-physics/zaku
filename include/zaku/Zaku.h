/// \file Zaku.h
/// Umbrella header: include this to get the whole public Zaku API.
//
// All Zaku headers use classic include guards rather than #pragma once: ROOT's
// autoparsing may include a header a second time through a different path
// spelling, and only macro guards survive that.
#ifndef ZAKU_ZAKU_H
#define ZAKU_ZAKU_H

#include "zaku/core/Version.h"
#include "zaku/core/Types.h"
#include "zaku/core/Units.h"
#include "zaku/core/Detector.h"
#include "zaku/core/FluxSpec.h"
#include "zaku/core/Production.h"
#include "zaku/core/Manifest.h"
#include "zaku/core/Provenance.h"
#include "zaku/core/RunSettings.h"
#include "zaku/core/Experiment.h"
#include "zaku/flux/FluxHistograms.h"
#include "zaku/flux/FluxPreparer.h"
#include "zaku/flux/FluxCovarianceInput.h"
#include "zaku/event/FlatEvent.h"
#include "zaku/event/FlatTreeCursor.h"
#include "zaku/event/WeightCursor.h"
#include "zaku/event/EventJoin.h"

#endif  // ZAKU_ZAKU_H
