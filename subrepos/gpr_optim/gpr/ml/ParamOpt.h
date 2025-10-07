#pragma once

#include <variant>

#include "ADAM.h"
#include "SCG.h"
#include "Ceres.h"

namespace funcmin {
using ParamOptimizer = std::variant<ADAM, SCG, Ceres>;
}
