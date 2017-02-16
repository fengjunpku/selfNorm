#include "scRawPixel.hh"

scRawPixel::scRawPixel()
{
  fch = bch = -1;
  Drop = false;
  k = b = w = 1;
}

scRawPixel::~scRawPixel()
{
  if(Pixel) delete Pixel;
}

void scRawPixel::Init(int frontCh,int backCh)
{
  if(frontCh<0||frontCh>SCDiNum||backCh<0||backCh>SCDjNum)
  {
    MiaoError("scRawPixel::scRawPixel, bad parameters !");
    exit(0);
  }
  fch = frontCh;
  bch = backCh;
  char NameBuff[20];
  sprintf(NameBuff,"rawPix_%02d_%02d",fch,bch);
  Pixel = new TGraph();
  Pixel->SetName(NameBuff);
  Pixel->SetDrawOption("AP*");
  sprintf(NameBuff,"fg_%02d_%02d",fch,bch);
  fg = new TGraph();
  fg->SetName(NameBuff);
}

void scRawPixel::Fill(double x,double y)
{
  Pixel->SetPoint(GetN(),x,y);
}

void scRawPixel::Write()
{
  if(!Pixel) {MiaoError("scRawPixel::Write(), Graph is null !");exit(0);}
  //printf("%02d_%02d : Slope : %4f+-%4f; Intercept : %4f+-%4f\n",fch,bch,Slope,ErSlope,Intercept,ErIntercept);
  Pixel->Write();
  fg->Write();
}

int scRawPixel::GetN()
{
  if(!Pixel) {MiaoError("scRawPixel::GetN(), Graph is null !");exit(0);}
  return Pixel->GetN();
}

void scRawPixel::Clear()
{
  Pixel->Clear();
}

void scRawPixel::FindLine()
{
  TF1 *line = new TF1("line","pol1",0,SCDmaxCh);
  int pNum = GetN();
  if(pNum<5) 
  {
    k=1;
    b=w=0;
    line->SetParameters(b,k);
    line->SetLineColor(kBlue);
    Pixel->GetListOfFunctions()->Add(line);
    Drop = true;
    return;
  }
  char NameBuff[20];
  sprintf(NameBuff,"htPix_%02d_%02d",fch,bch);
  TH2F hough(NameBuff,NameBuff,800,-65,-25,SCDmaxCh,-0.5*SCDmaxCh,0.5*SCDmaxCh);
  TF1 ht("ht",htCurve,-180,180,2);
  Double_t MaxDeg = hough.GetXaxis()->GetXmax();
  Double_t MinDeg = hough.GetXaxis()->GetXmin();
  Int_t nX = hough.GetXaxis()->GetNbins();
  Double_t stepX = (MaxDeg-MinDeg)/nX;
  Double_t MaxHei = hough.GetYaxis()->GetXmax();
  Double_t MinHei = hough.GetYaxis()->GetXmin();
  Int_t nY = hough.GetYaxis()->GetNbins();
  Double_t stepY = (MaxHei-MinHei)/nY;
  for(int i=0;i<pNum;i++)
  {
    Double_t x,y;
    Pixel->GetPoint(i,x,y);
    ht.SetParameters(x,y);
    for(int iht=0;iht<nX;iht++)
    {
      Double_t xht = iht*stepX+MinDeg;
      hough.Fill(xht,ht.Eval(xht));
    }
  }
  Int_t xx,yy,zz;
  hough.GetMaximumBin(xx,yy,zz);
  Double_t th = xx*stepX+MinDeg;
  Double_t hh = yy*stepY+MinHei;
  b = hh/TMath::Sin(th*deg);
  k = -1./TMath::Tan(th*deg);
  w = hough.GetMaximum()/stepX/stepY;
  if(hough.GetMaximum()<3) {k=1;b=w=0;Drop=true;}
  line->SetParameters(b,k);
  line->SetLineColor(kBlue);
  Pixel->GetListOfFunctions()->Add(line);
  printf("  * find %02d,%02d : k= %10.4f, b= %8.2f, w= %8.0f\n",fch,bch,k,b,w);
  char titleBuff[40];
  sprintf(titleBuff,"P%02d_%02d_k_%f_b_%f_w_%f",fch,bch,k,b,w);
  Pixel->SetTitle(titleBuff);
}

void scRawPixel::ReFill()
{
  int pNum = GetN();
  int fi = 0;
  for(int i=0;i<pNum;i++)
  {
    Double_t x,y;
    Pixel->GetPoint(i,x,y);
    Double_t d = TMath::Abs(y-k*x-b);
    if(d<5) fg->SetPoint(fi++,x,y);
  }
  fg->Fit("pol1","Q");
  double p0,p1,e0,e1,chi2;
  int ndf;
  p0 = fg->GetFunction("pol1")->GetParameter(0);
  p1 = fg->GetFunction("pol1")->GetParameter(1);
  e0 = fg->GetFunction("pol1")->GetParError(0);
  e1 = fg->GetFunction("pol1")->GetParError(1);
  chi2 = fg->GetFunction("pol1")->GetChisquare();
  ndf = fg->GetFunction("pol1")->GetNDF();
  printf("  - fit  %02d,%02d : k= %10.4f, b= %8.2f, chi2= %8.0f\n",fch,bch,p1,p0,chi2);
  char titleBuff[40];
  sprintf(titleBuff,"F%02d_%02d_k_%f_b_%f_c2_%f",fch,bch,p1,p0,chi2);
  fg->SetTitle(titleBuff);
}
