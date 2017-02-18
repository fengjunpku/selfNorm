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
  _dname = tele + dssd;
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
}

void snExtract::Init()
{
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
  Weight();
  char fname[20];
  sprintf(fname,"pix_%s_%04d.root",_dname.Data(),_runNum);
  Save(fname);
  Output();
  Clear();
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

void snExtract::Save(const char* tmpfile)
{
  TFile *tmpf = new TFile(tmpfile,"RECREATE","tmpf");
  tmpf->cd();
  Write();
  tmpf->Close();
}

void snExtract::Clear()
{
  for(int i=_iZero;i<_iNum+_iZero;i++)
    for(int j=_jZero;j<_jNum+_jZero;j++)
      rawPix[i][j].Clear();
}

void snExtract::Fit()
{
  cout<<"   =======Fitting======="<<endl;
  for(int i=_iZero;i<_iNum+_iZero;i++)
  {
    for(int j=_jZero;j<_jNum+_jZero;j++)
    {
      rawPix[i][j].Fit();
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

void snExtract::Weight()
{
  cout<<"   =======Weight======="<<endl;
  int refBack = _jZero+_jNum/2;
  cout<<" *Front:"<<endl;
  for(int q=_iZero;q<_iNum+_iZero;q++)
  {
    
    W_f[q][0] = rawPix[q][refBack].w0;
    W_f[q][1] = rawPix[q][refBack].w1;
    P_f[q][0] = W_f[q][0] * rawPix[q][refBack].p0;
    P_f[q][1] = W_f[q][1] * rawPix[q][refBack].p1;
    e_f[q][0] = 0;
    e_f[q][1] = 0;
    fNode[q] = 0.1;
    for(int i=_iZero;i<_iNum+_iZero;i++)
    {
      if(i==q) continue;
      for(int j=_jZero;j<_jNum+_jZero;j++)
      {
        if(j==refBack) continue;
        if(rawPix[q][j].Drop||rawPix[i][j].Drop||rawPix[i][refBack].Drop) continue;
        fNode[q] += 1;
        Double_t K1 = rawPix[q][j].p1;
        Double_t B1 = rawPix[q][j].p0;
        Double_t ek1 = rawPix[q][j].e1;
        Double_t eb1 = rawPix[q][j].e0;
        
        Double_t K2 = rawPix[i][j].p1;
        Double_t B2 = rawPix[i][j].p0;
        Double_t ek2 = rawPix[i][j].e1;
        Double_t eb2 = rawPix[i][j].e0;

        Double_t K3 = rawPix[i][refBack].p1;
        Double_t B3 = rawPix[i][refBack].p0;
        Double_t ek3 = rawPix[i][refBack].e1;
        Double_t eb3 = rawPix[i][refBack].e0;
        //Double_t Wtmp = W1*W2*W3/(W1*W2+W2*W3+W3*W1);
        Double_t Ktmp = K1*K3/K2;
        Double_t Wktp = 1./(pow2(K3*ek1/K2)+pow2(K1*ek3/K2)+pow2(K1*K3*ek2/K2/K2));
        Double_t Btmp = B3+(B1-B2)*K3/K2;
        Double_t Wbtp = 1./(eb3*eb3+pow2((B1-B2)*ek3/K2)+pow2(K3*eb1/K2)+pow2(K3*eb2/K2)+pow2((B1-B2)*K3*ek2/K2/K2));
        P_f[q][1] += Wktp*Ktmp;
        P_f[q][0] += Wbtp*Btmp;
        W_f[q][1] += Wktp;
        W_f[q][0] += Wbtp;
      }
    }
    if(W_f[q][1]>0)
    {
      P_f[q][1] /= W_f[q][1];
      e_f[q][1] = 1./TMath::Sqrt( W_f[q][1]);
    }
    if(W_f[q][0]>0)
    {
      P_f[q][0] /= W_f[q][0];
      e_f[q][0] = 1./TMath::Sqrt( W_f[q][0]);
    }
    printf("front ch %02d;k: %8.6f +/- %8.6f, b: %8.4f +/- %8.4f, n: %10.1f\n",q,P_f[q][1],e_f[q][1],P_f[q][0],e_f[q][0],fNode[q]);
  }
  //========================================
  cout<<" *Back:"<<endl;
  for(int p=_jZero;p<_jNum+_jZero;p++)
  {
    W_b[p][0] = 0;
    W_b[p][1] = 0;
    P_b[p][0] = 0;
    P_b[p][1] = 0;
    e_b[p][0] = 0;
    e_b[p][1] = 0;
    bNode[p] = 0.1;
    for(int i=_iZero;i<_iNum+_iZero;i++)
    {
      if(p==refBack) 
      {
        P_b[p][1] = 1;
        break;
      }
      if(rawPix[i][p].Drop||rawPix[i][refBack].Drop) continue;
      bNode[p] += 1;
      Double_t K1 = rawPix[i][p].p1;
      Double_t B1 = rawPix[i][p].p0;
      Double_t ek1 = rawPix[i][p].e1;
      Double_t eb1 = rawPix[i][p].e0;
      
      Double_t K2 = rawPix[i][refBack].p1;
      Double_t B2 = rawPix[i][refBack].p0;
      Double_t ek2 = rawPix[i][refBack].e1;
      Double_t eb2 = rawPix[i][refBack].e0;
      //Double_t Wtmp = W1*W2/(W1+W2);
      Double_t Ktmp = K2/K1;
      Double_t Wktp = 1./(pow2(ek2/K1)+pow2(K2*ek1/K1/K1));
      Double_t Btmp = B2-B1*K2/K1;
      Double_t Wbtp = 1./(pow2(eb2)+pow2(eb1*K2/K1)+pow2(B1*ek2/K1)+pow2(B1*K2*ek1/K1/K1));
      P_b[p][1] += Wktp*Ktmp;
      P_b[p][0] += Wbtp*Btmp;
      W_b[p][1] += Wktp;
      W_b[p][0] += Wbtp;
    }
    if(W_b[p][0]>0)
    {
      P_b[p][0] /= W_b[p][0];
      e_b[p][0] = 1./TMath::Sqrt( W_b[p][0]);
    }
    if(W_b[p][1]>0)
    {
      P_b[p][1] /= W_b[p][1];
      e_b[p][1] = 1./TMath::Sqrt( W_b[p][1]);
    }
    printf("back  ch %02d;k: %8.6f +/- %8.6f, b: %8.4f +/- %8.4f, n: %10.1f\n",p,P_b[p][1],e_b[p][1],P_b[p][0],e_b[p][0],bNode[p]);
  }
}

void snExtract::Output()
{
  FILE *out;
  char buff[20];
  sprintf(buff,"sn_%sf_%04d.txt",_dname.Data(),_runNum);
  out = fopen(buff,"w+");
  if(!out) MiaoError("scExtract::Output, Can't open output file !");
  for(int i=_iZero;i<_iNum+_iZero;i++)
    fprintf(out, "    %f     %f     %f     %f\n",P_f[i][0],e_f[i][0],P_f[i][1],e_f[i][1]);
  fclose(out);
  //---
  sprintf(buff,"sn_%sb_%04d.txt",_dname.Data(),_runNum);
  out = fopen(buff,"w+");
  if(!out) MiaoError("scExtract::Output, Can't open output file !");
  for(int j=_jZero;j<_jNum+_jZero;j++)
    fprintf(out, "    %f     %f     %f     %f\n",P_b[j][0],e_b[j][0],P_b[j][1],e_b[j][1]);
  fclose(out);
}