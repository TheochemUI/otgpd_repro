/*
 * SCG.cpp
 *
 *  Created on: 20 Nov 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "SCG.h"

namespace funcmin {

SCG::SCG()
{
    settings.setDefault();
    failedOptim=false;
}

SCG::~SCG()
{
    settings.setDefault();
}

} /* namespace funcmin */
