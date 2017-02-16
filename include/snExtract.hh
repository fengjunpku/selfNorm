#ifndef SCEXTRACT_HH
#define SNEXTRACT_HH_HH 1

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStatistic.h>

#include "JunErrors.hh"
#include "selCali_def.hh"
#include "scRawPixel.hh"

using namespace selCaliDef;

class snExtract
{
public:
  snExtract(int runNo, TString tele, TString dssd);
  virtual ~snExtract();
  void Init();
  void Save(const char* tmpfile  = "pixels.root");
  void Write();
  void Load();
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
  TFile *_inputFile;
  int _adc[SCDadcNum][32];
  int _madc[SCDmadcNum][32];
  int _runNum;
  int _iNum,_jNum;// 16 or 32
  int _iZero,_jZero;
  int _iGeo, _jGeo;
  TStatistic *_Ks,*_Bs;
};
#endif