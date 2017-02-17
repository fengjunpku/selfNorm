#ifndef SNDEF_hh
#define SNDEF_hh 1

#include <TString.h>
#include <TMath.h>

namespace snDefine
{
  static const TString SNpath = "/data/d2/CIAE_2016_BCO_data/datafile/";//input path
  static const TString SNtreeName = "tree";//tree name
  static const TString SNadcName = "madc";//branch name of adc

  static const int SNadcNum = 5;//num of adc
  static const int SNmadcNum = 11;//num of adc
  static const int SNfaceNum = 1;//adc No.
  static const int SNbackNum = 2;

  static const int SNmaxCh = 8000;//madc : 8000; adc : 4000;

  static const Double_t deg = TMath::DegToRad();
  inline Double_t htCurve(Double_t *x,Double_t *p)
  {
    return p[0]*TMath::Cos(x[0]*deg)+p[1]*TMath::Sin(x[0]*deg);
  }
}
#endif