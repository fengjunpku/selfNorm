#ifndef SNPIXEL_HH
#define SNPIXEL_HH 1

#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include <TGraph.h>
#include <TF1.h>
#include <TList.h>
#include <TH2F.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include "snDef.hh"
#include "JunErrors.hh"

using namespace std;
using namespace snDefine;

class snPixel
{
public:
  snPixel();
  virtual ~snPixel();
  void Init(int frontCh,int backCh);
  void Fill(double x,double y);
  void Write();
  void Clear();
  int GetN();
  int fch,bch;
  TGraph *Pixel,*fg,*fg2;

  void ReFill();
  void ReFill2();
  double k,b;
  double w;
  double p0,p1,e0,e1;
  bool Drop;
  int n2k;
};
#endif