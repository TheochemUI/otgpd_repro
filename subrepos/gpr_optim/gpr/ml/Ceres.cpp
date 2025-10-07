#include "Ceres.h"

namespace funcmin {

Ceres::Ceres()
{
    settings.setDefault();
    failedOptim = false;
}

Ceres::~Ceres()
{
    settings.setDefault();
}
}  // namespace funcmin
