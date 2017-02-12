#ifndef SELCALI_DEF_hh
#define SELCALI_DEF_hh 1

#include <TString.h>
#include <TMath.h>

namespace selCaliDef
{
  //static const TString SCDinfile = "data/data0190.root";//input file
  static const TString SCDpath = "/data/d2/CIAE_2016_BCO_DATA/datafile/";//input path
  static const TString SCDtreeName = "tree";//tree name
  static const TString SCDadcName = "madc";//branch name of adc

  static const int SCDadcNum = 11;//num of adc
  static const int SCDfaceNum = 1;//adc No.
  static const int SCDbackNum = 2;
  
  static const int SCDiZero = 0;//i strip start
  static const int SCDjZero = 0;  
  static const int SCDiNum = 32;//num of i strip
  static const int SCDjNum = 32;

  static const int SCDmaxCh = 8000;//madc : 8000; adc : 4000;

  static const Double_t deg = TMath::DegToRad();
  inline Double_t htCurve(Double_t *x,Double_t *p)
  {
    return p[0]*TMath::Cos(x[0]*deg)+p[1]*TMath::Sin(x[0]*deg);
  }
}
#endif
