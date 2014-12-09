#include <iostream>
#include "TigrStrTokenizer.H"
#include "TigrGffFeature.H"
using namespace std;



TigrGffFeature::TigrGffFeature(const TigrString &rawLine,const TigrString &substrate,
		       const TigrString &source,const TigrString &featureType,
		       int begin,int end,double score,bool hasScore,
		       char strand,int frame,bool hasFrame)
  : rawLine(rawLine), substrate(substrate), source(source), 
    featureType(featureType), begin(begin), end(end), score(score),
    hasScore(hasScore), strand(strand), frame(frame), hasFrame(hasFrame)
{
}



TigrGffFeature::TigrGffFeature(const TigrString &rawLine)
  : rawLine(rawLine), score(0)
{
  parseLine(rawLine);
}



ostream &operator<<(ostream &os,TigrGffFeature &feature)
{
  os<<feature.toGff();
  return os;
}



TigrString TigrGffFeature::toGff() const
{
  TigrString gff=
    substrate+"\t"+source+"\t"+featureType+"\t";
  gff+=TigrString(begin+1)+"\t"+TigrString(end)+"\t";
  if(hasScore) gff+=TigrString("")+score; else gff+=".";
  gff+=TigrString("\t")+strand+"\t";
  if(hasFrame) gff+=TigrString("")+frame; else gff+=".";
  gff+="\n";
  return gff;
}



TigrVector<TigrString> &TigrGffFeature::getExtraFields()
{
  return extraFields;
}



bool TigrGffFeature::hasExtraFields() const
{
  return extraFields.size()>0;
}



bool TigrGffFeature::isFramed() const
  {
    return hasFrame;
  }



bool TigrGffFeature::isScored() const
  {
    return hasScore;
  }



bool TigrGffFeature::isStranded() const
  {
    return strand!='.';
  }



char TigrGffFeature::getStrand() const
  {
    return strand;
  }



const TigrString &TigrGffFeature::getFeatureType() const
{
  return featureType;
}



const TigrString &TigrGffFeature::getRawField(int index)
{
  return allFields[index];
}



const TigrString &TigrGffFeature::getRawLine() const
{
  return rawLine;
}



const TigrString &TigrGffFeature::getSource() const
{
  return source;
}



const TigrString &TigrGffFeature::getSubstrate() const
{
  return substrate;
}



double TigrGffFeature::getScore() const
{
  return score;
}



int TigrGffFeature::getBegin() const
{
   
  return begin;
}



int TigrGffFeature::getEnd() const
{
   
  return end;
}



int TigrGffFeature::getFrame() const
  {
    return frame;
  }



void TigrGffFeature::parseLine(const TigrString &line)
{
  TigrStrTokenizer tokenizer(line);
  int fieldNum=0;
  TigrString dot=".";
  for(; tokenizer.hasMoreTokens() ; ++fieldNum)
    {
      const char *token=tokenizer.nextToken();
      allFields.push_back(token);
      switch(fieldNum)
	{
	case 0:
	  substrate=token;
	  break;
	case 1:
	  source=token;
	  break;
	case 2:
	  featureType=token;
	  break;
	case 3:
	  begin=atoi(token)-1;
	  break;
	case 4:
	  end=atoi(token);
	  break;
	case 5:
	  if(dot==token)
	    hasScore=false;
	  else
	    {
	      hasScore=true;
	      score=atof(token);
	    }
	  break;
	case 6:
	  strand=*token;
	  break;
	case 7:
	  if(dot==token)
	    hasFrame=false;
	  else
	    {
	      hasFrame=true;
	      frame=atoi(token);
	    }
	  break;
	default:
	  extraFields.push_back(token);
	}
    }
}
