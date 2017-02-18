#ifndef SCEXTRACT_HH
#define SNEXTRACT_HH_HH 1

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

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
  void Clear();
  void Weight();
  void Output();

  TTree *dtree;
  snPixel rawPix[32][32];
  Double_t P_f[32][2],W_f[32][2],e_f[32][2];
  Double_t P_b[32][2],W_b[32][2],e_b[32][2];
  float fNode[32],bNode[32];

private:
  TFile *_inputFile;
  int _adc[SNadcNum][32];
  int _madc[SNmadcNum][32];
  int _runNum;
  int _iNum,_jNum;// 16 or 32
  int _iZero,_jZero;
  int _iGeo, _jGeo;
  TString _dname;
};
#endif