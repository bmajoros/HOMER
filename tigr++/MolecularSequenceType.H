/****************************************************************
 MolecularSequenceType.H
 bmajoros@tigr.org

 Copyright (c)2003, The Institute for Genomic Research (TIGR),
 Rockville, Maryland, U.S.A.  All rights reserved.  This software
 is governed by the ARTISTIC LICENSE (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_MolecularSequenceType_H
#define INCL_MolecularSequenceType_H
#include <iostream>

enum MolecularSequenceType
  {
    DNA,
    PROTEIN
  };



inline ostream &operator<<(ostream &os,MolecularSequenceType t)
{
  switch(t)
    {
    case DNA: os<<"DNA"; break;
    case PROTEIN: os<<"PROTEIN"; break;
    }
  return os;
}

#endif
