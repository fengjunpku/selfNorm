#ifndef SCEXTRACT_HH
#define SNEXTRACT_HH_HH 1

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStatistic.h>

#include "JunErrors.hh"
#include "snDef.hh"
#include "snPixel.hh"

using namespace snDefine;

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
  snPixel rawPix[32][32];
  Double_t Kf[32],Bf[32],Wf0[32],Wf1[32];
  Double_t Kb[32],Bb[32],Wb[32];
  float fNode[32],bNode[32];

private:
  TFile *_inputFile;
  int _adc[SNadcNum][32];
  int _madc[SNmadcNum][32];
  int _runNum;
  int _iNum,_jNum;// 16 or 32
  int _iZero,_jZero;
  int _iGeo, _jGeo;
  TStatistic *_Ks,*_Bs;
};
#endif