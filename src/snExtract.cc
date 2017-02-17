#include "snExtract.hh"

snExtract::snExtract(int runNo, TString tele, TString dssd)
{
  //check input
  if(!tele.EqualTo("l0")&&!tele.EqualTo("l1")&&!tele.EqualTo("l2")
    &&!tele.EqualTo("r0")&&!tele.EqualTo("r1")&&!tele.EqualTo("r2"))
    MiaoError("Tele name should be l0/l1/l2/r0/r1/r2!");
  if(!dssd.EqualTo("w1")&&!dssd.EqualTo("bb7")&&!dssd.EqualTo("ssd"))
    MiaoError("Detector name should be w1/bb7/ssd!");
  if((tele.EqualTo("l1")||tele.EqualTo("r1"))&&dssd.EqualTo("bb7"))
    MiaoError("There is not l1/r1 - bb7");
  if((tele.EqualTo("l2")||tele.EqualTo("r2"))&&dssd.EqualTo("w1"))
    MiaoError("There is not l2/r2 - w1");
  //open tele root file
  char filename[20];
  _runNum = runNo;
  sprintf(filename,"data%04d.root",_runNum);
  TString datafile = SNpath+ filename;
  _inputFile = TFile::Open(datafile.Data());
  if(!_inputFile || !_inputFile->IsOpen())
  {
    char errStr[36];
    sprintf(errStr,"Can not open input file %s, Exit !",datafile.Data());
    MiaoError(errStr);
  }
  _inputFile->GetObject("tree",dtree);
  //
  if(dssd.EqualTo("w1"))
  {
    _iNum = _jNum = 16;
    _iZero = 0;
    _jZero = 16;
    if(tele.EqualTo("l0")) {_iGeo = 10; _jGeo = 10;}
    if(tele.EqualTo("r0")) {_iGeo = 13; _jGeo = 13;}
    if(tele.EqualTo("l1")) {_iGeo = 17; _jGeo = 18;_jZero = 0;}
    if(tele.EqualTo("r1")) {_iGeo = 17; _jGeo = 18;_iZero = 16;}
  }
  if(dssd.EqualTo("bb7"))
  {
    _iNum = _jNum = 32;
    _iZero = _jZero = 0;
    if(tele.EqualTo("l0")) {_iGeo = 11; _jGeo = 12;}
    if(tele.EqualTo("r0")) {_iGeo = 14; _jGeo = 15;}
    if(tele.EqualTo("l2")) {_iGeo = 16; _jGeo = 1;}
    if(tele.EqualTo("r2")) {_iGeo = 2; _jGeo = 3;}
  }
  Init();
}

snExtract::~snExtract()
{
  _inputFile->Close();
  delete _inputFile;
  delete _Ks;
  delete _Bs;
}

void snExtract::Init()
{
  _Ks = new TStatistic("Sta_K");
  _Bs = new TStatistic("Sta_B");
  _inputFile->GetObject(SNtreeName.Data(),dtree);
  dtree->SetBranchStatus("*",0);
  dtree->SetBranchStatus("madc",1);
  dtree->SetBranchStatus("adc",1);
  dtree->SetBranchAddress("madc",_madc);
  dtree->SetBranchAddress("adc",_adc);
  for(int i=0;i<32;i++)
    for(int j=0;j<32;j++)
      rawPix[i][j].Init(i,j);
}

void snExtract::Load()
{
  cout<<"   =======Loading======="<<endl;
  if(!dtree) MiaoError("data::loop the tree is null");
  long nentries = dtree->GetEntriesFast();
  for(long ientry=0;ientry<nentries;ientry++)
  {
    if(!(ientry%100000))
      printf(" %10ld : %5.2f %%\n",ientry,(double)ientry/(double)nentries*100.);
      //cout<<ientry<<"    "<<(double)ientry/(double)nentries*100.<<" %"<<endl;
    dtree->GetEntry(ientry);
    //count the hit num
    int nhit = 0;
    int *fadc,*badc;
    if(_iGeo<10) fadc = _adc[_iGeo];
    else     fadc = _madc[_iGeo-10];
    if(_jGeo<10) badc = _adc[_jGeo];
    else     badc = _madc[_jGeo-10];
    for(int i=_iZero;i<_iNum+_iZero;i++)
    {
      for(int j=_jZero;j<_jNum+_jZero;j++)
      {
        if(fadc[i]>250&&badc[j]>250) nhit++;
      }
    }
    if(nhit != 1) continue;
    //load data
    for(int i=_iZero;i<_iNum+_iZero;i++)
    {
      for(int j=_jZero;j<_jNum+_jZero;j++)
      {
        double x,y;
        x = fadc[i];
        y = badc[j];
        if(_jGeo == 1) y*=2;
        if(Select(x,y))
          rawPix[i][j].Fill(x,y);
      }
    }
  }
  Fit();
  TFile *tmpf = new TFile("x.root","RECREATE","tmpf");
  tmpf->cd();
  Write();
  tmpf->Close();
}

void snExtract::Write()
{
  for(int i=_iZero;i<_iNum+_iZero;i++)
  {
    for(int j=_jZero;j<_jNum+_jZero;j++)
    {
      rawPix[i][j].Write();
      //printf("_%02d_%02d : %9d\n",i,j,pixels[i][j].GetN());
    }
  }
}

void snExtract::Fit()
{
  cout<<"   =======Fitting======="<<endl;
  for(int i=_iZero;i<_iNum+_iZero;i++)
  {
    for(int j=_jZero;j<_jNum+_jZero;j++)
    {
      //rawPix[i][j].FindLine();
      //_Ks->Fill(rawPix[i][j].k);
      //_Bs->Fill(rawPix[i][j].b);
      rawPix[i][j].ReFill();
      rawPix[i][j].ReFill2();
    }
  }
}

bool snExtract::Select(double x,double y)
{
  if(x<250||y<250) return false;
  if(x>SNmaxCh||y>SNmaxCh) return false;
  if(y>1.2*x+300) return false;
  if(x>1.2*y+300) return false;
  return true;
}