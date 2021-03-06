#ifndef SCRAWPIXEL_HH
#define SCRAWPIXEL_HH 1

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

#include "selCali_def.hh"
#include "JunErrors.hh"

using namespace std;
using namespace selCaliDef;

class scRawPixel
{
public:
  scRawPixel();
  virtual ~scRawPixel();
  void Init(int frontCh,int backCh);
  void Fill(double x,double y);
  void Write();
  void Clear();
  int GetN();
  int fch,bch;
  TGraph *Pixel;
  TGraph *fg;
  TGraph *fg2;

  void FindLine();
  void ReFill();
  void ReFill2();
  double k,b;
  double w;
  double p0,p1,e0,e1;
  bool Drop;
  int n2k;
};
#endif
