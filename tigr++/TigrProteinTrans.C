#include <iostream>
#include "TigrProteinTrans.H"
using namespace std;



TigrMap<const char*,char,TigrProteinTrans::mapstrcmp> 
  TigrProteinTrans::codonMap;
TigrProteinTrans TigrProteinTrans::initializer;



TigrProteinTrans::TigrProteinTrans()
{
  setupCodonMap();
}



TigrString TigrProteinTrans::reverseComplement(const TigrString &s)
{
  int length=s.length();
  char *buf=new char[length+1];
  buf[length]='\0';
  for(int i=0,j=length-1 ; i<length ; ++i,--j)
    buf[i]=complement(s[j]);
  TigrString r=buf;
  delete [] buf;
  return r;
}



TigrString TigrProteinTrans::translate(const TigrString &transcript)
{
  TigrString translation;
  const char *pTranscript=transcript.c_str();

  int length=transcript.length();
  char buf[4];
  buf[3]='\0';
  for(int i=0 ; i+2<length ; i+=3)
    {
      for(int j=0 ; j<3 ; ++j)
	buf[j]=pTranscript[i+j];
      translation+=codonMap[buf];
    }

  return translation;
}



char TigrProteinTrans::complement(char c)
{
  switch(c)
    {
    case 'N': return 'N';
    case 'R': return 'Y';
    case 'Y': return 'R';
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    case ' ': return ' ';
    default:
      throw TigrString("Can't complement base: ")+c;
    }
}



char TigrProteinTrans::mapCodon(const char *codon)
{
  return codonMap.find(codon)==codonMap.end() ? 'X' : codonMap[codon];
}



void TigrProteinTrans::setupCodonMap()
{
  codonMap["ACA"]='T';
  codonMap["GAA"]='E';
  codonMap["CTG"]='L';
  codonMap["TGY"]='C';
  codonMap["GGY"]='G';
  codonMap["GGG"]='G';
  codonMap["CTA"]='L';
  codonMap["AGY"]='S';
  codonMap["TTC"]='F';
  codonMap["GGN"]='G';
  codonMap["GCG"]='A';
  codonMap["GAT"]='D';
  codonMap["TAT"]='Y';
  codonMap["CGT"]='R';
  codonMap["GGT"]='G';
  codonMap["AGC"]='S';
  codonMap["ACG"]='T';
  codonMap["TCT"]='S';
  codonMap["TAA"]='*';
  codonMap["TAC"]='Y';
  codonMap["AGA"]='R';
  codonMap["CGC"]='R';
  codonMap["CAC"]='H';
  codonMap["ATG"]='M';
  codonMap["GTR"]='V';
  codonMap["TCG"]='S';
  codonMap["TAR"]='*';
  codonMap["TAG"]='*';
  codonMap["CAR"]='Q';
  codonMap["CTT"]='L';
  codonMap["ATC"]='I';
  codonMap["GAC"]='D';
  codonMap["TCC"]='S';
  codonMap["CCR"]='P';
  codonMap["AAR"]='K';
  codonMap["GCC"]='A';
  codonMap["GGR"]='G';
  codonMap["CGA"]='R';
  codonMap["TGG"]='W';
  codonMap["CTN"]='L';
  codonMap["GTY"]='V';
  codonMap["CAG"]='Q';
  codonMap["GCR"]='A';
  codonMap["TTG"]='L';
  codonMap["CGG"]='R';
  codonMap["GAY"]='D';
  codonMap["TTY"]='F';
  codonMap["AGT"]='S';
  codonMap["TAY"]='Y';
  codonMap["ATT"]='I';
  codonMap["TGC"]='C';
  codonMap["CCY"]='P';
  codonMap["TAG"]='*';
  codonMap["GTA"]='V';
  codonMap["AAA"]='K';
  codonMap["CCT"]='P';
  codonMap["CCG"]='P';
  codonMap["GCN"]='A';
  codonMap["AGG"]='R';
  codonMap["GCT"]='A';
  codonMap["GGA"]='G';
  codonMap["AAC"]='N';
  codonMap["CGN"]='R';
  codonMap["GAG"]='E';
  codonMap["CAY"]='H';
  codonMap["AAY"]='N';
  codonMap["CAT"]='H';
  codonMap["CCC"]='P';
  codonMap["GTN"]='V';
  codonMap["GTG"]='V';
  codonMap["AAG"]='K';
  codonMap["GGC"]='G';
  codonMap["ATA"]='I';
  codonMap["GAR"]='E';
  codonMap["TTA"]='L';
  codonMap["CCN"]='P';
  codonMap["GCY"]='A';
  codonMap["CAA"]='Q';
  codonMap["TGT"]='C';
  codonMap["AAT"]='N';
  codonMap["YTR"]='L';
  codonMap["ACT"]='T';
  codonMap["TTR"]='L';
  codonMap["GCA"]='A';
  codonMap["TCA"]='S';
  codonMap["CCA"]='P';
  codonMap["AGY"]='R';
  codonMap["TGA"]='*';
  codonMap["TTT"]='F';
  codonMap["ACC"]='T';
  codonMap["GTT"]='V';
  codonMap["GTC"]='V';
  codonMap["TCN"]='S';
  codonMap["ATH"]='I';
  codonMap["CTC"]='L';
  codonMap["ACN"]='T';
}



Sequence *TigrProteinTrans::translate(const Sequence &source,
				      Alphabet &dnaAlphabet,
				      Alphabet &proteinAlphabet)
{
  TigrString *sourceString=source.toString(dnaAlphabet);
  TigrString destString=translate(*sourceString);
  delete sourceString;
  return new Sequence(destString,proteinAlphabet);
}



void TigrProteinTrans::translate(const Sequence &source,
				 Alphabet &dnaAlphabet,
				 Alphabet &proteinAlphabet,
				 Sequence &dest)
{
  TigrString *sourceString=source.toString(dnaAlphabet);
  TigrString destString=translate(*sourceString);
  delete sourceString;
  dest.copyFrom(destString,proteinAlphabet);
}


