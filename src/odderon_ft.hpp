//
//  odderon_ft.hpp
//  rholib
//
//  Created by Heikki Mäntysaari on 14.10.2021.
//

#ifndef odderon_ft_hpp
#define odderon_ft_hpp

#include <stdio.h>


const double VEGAS_CHISQR_TOLERANCE = 0.2;
const double MC_ERROR_TOLERANCE = 0.2;

enum QMIN_CUTOFF
{
 HARD,
 GAUSSIAN
};

const QMIN_CUTOFF qmin_ir_cutoff = HARD;

#endif /* odderon_ft_hpp */
