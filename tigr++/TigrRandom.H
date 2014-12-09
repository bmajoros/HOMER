/**************************************************************
TigrRandom.H
bmajoros@tigr.org 1/1/2003

TIGR++ : C++ class template library for bioinformatics

Copyright (c) 2003, The Institute for Genomic Research (TIGR), 
Rockville, Maryland, U.S.A.  All rights reserved.
***************************************************************/
#ifndef INCL_TigrRandom_H
#define INCL_TigrRandom_H
#include <limits.h>
using namespace std;

#define INCL_RANDOM_H

const int LARGEST_RANDOM_NUMBER=SHRT_MAX;  

unsigned GetRandomSeed();

void SeedRandomizer(unsigned);

void randomize();

int RandomNumber(int n);

float Random0to1();

float RandomFloat(float range);

float RandomFloat(float from,float to);

float RandomGaussian(float min,float max,int n);

#endif
