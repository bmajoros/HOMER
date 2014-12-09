/**************************************************************
Constants.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR),
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_Constants_H
#define INCL_Constants_H
using namespace std;
#include <math.h>
#include <limits.h>
#include <float.h>

//==================================================================
// see /usr/include/limits.h and /usr/include/values.h for others!!!
//==================================================================

const double POSITIVE_INFINITY    =-log(0.0);
const double NEGATIVE_INFINITY    =log(0.0);
const int    LARGEST_INTEGER      =INT_MAX;
const float  LARGEST_FLOAT        =FLT_MAX;
const float  SMALLEST_FLOAT       =FLT_MIN;
const double LARGEST_DOUBLE       =DBL_MAX;
const double SMALLEST_DOUBLE      =DBL_MIN;
const double PI                   =acos(-1.0);

inline bool isInfinity(double d)
{
  return d==NEGATIVE_INFINITY || d==POSITIVE_INFINITY;
}

inline bool isNaN(double d)
{
  return isnan(d);
}

inline bool isFinite(double d)
{
  return !(isinf(d) || isnan(d));
}

#endif
