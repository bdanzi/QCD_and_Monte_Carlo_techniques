#ifndef courseLIB__H
#define courseLIB__H

#include "ranlxd.h"

inline double Rand()
{
    double r; 
    ranlxd(&r,1);
    return r;
}
 
#endif
