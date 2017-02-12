#ifndef SCEXTRACT_HH
#define SCEXTRACT_HH_HH 1

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStatistic.h>

#include "JunErrors.hh"
#include "selCali_def.hh"
#include "scRawPixel.hh"

using namespace selCaliDef;

class scExtract
{
public:
  scExtract(int runNo);
  virtual ~scExtract();
  void Init();
  void Save(const char* tmpfile  = "pixels.root");
  void Write();
  void Load(int face,int back);
  bool Select(double x,double y);
  void Fit();
  void Weight();
  void Clear();
  void Output(const char* outfile = "pars.txt");

  TTree *dtree;
  scRawPixel rawPix[SCDiNum][SCDjNum];
  Double_t Kf[SCDiNum],Bf[SCDiNum],Wf[SCDiNum];
  Double_t Kb[SCDjNum],Bb[SCDjNum],Wb[SCDjNum];
  float fNode[SCDiNum],bNode[SCDjNum];

private:
  TFile *inputFile;
  int adc[SCDadcNum][32];
  int runNum;
  TStatistic *Ks,*Bs;
};
#endif
