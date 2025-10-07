#include "ADAM.h"

namespace funcmin {

ADAM::ADAM()
{
    settings.setDefault();
    failedOptim = false;
}

ADAM::~ADAM()
{
    // settings.setDefault(); // Optional cleanup
}

} /* namespace funcmin */
