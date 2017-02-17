#include "snPixel.hh"

snPixel::snPixel()
{
  fch = bch = -1;
  Drop = false;
  k = b = w = 1;
}

snPixel::~snPixel()
{
  if(Pixel) delete Pixel;
}

void snPixel::Init(int frontCh,int backCh)
{
  if(frontCh<0||frontCh>32||backCh<0||backCh>32)
  {
    MiaoError("snPixel::snPixel, bad parameters !");
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
  sprintf(NameBuff,"fg2_%02d_%02d",fch,bch);
  fg2 = new TGraph();
  fg2->SetName(NameBuff);
  n2k = 0;
}

void snPixel::Fill(double x,double y)
{
  Pixel->SetPoint(GetN(),x,y);
}

void snPixel::Write()
{
  if(!Pixel) {MiaoError("snPixel::Write(), Graph is null !");exit(0);}
  //printf("%02d_%02d : Slope : %4f+-%4f; Intercept : %4f+-%4f\n",fch,bch,Slope,ErSlope,Intercept,ErIntercept);
  Pixel->Write();
  fg->Write();
  fg2->Write();
}

int snPixel::GetN()
{
  if(!Pixel) {MiaoError("snPixel::GetN(), Graph is null !");exit(0);}
  return Pixel->GetN();
}

void snPixel::Clear()
{
  Pixel->Clear();
}

void snPixel::ReFill()
{
  int pNum = GetN();
  int fi = 0;
  for(int i=0;i<pNum;i++)
  {
    Double_t x,y;
    Pixel->GetPoint(i,x,y);
    //Double_t d = TMath::Abs(y-k*x-b);
    //Double_t limitD = x*0.05;
    //if(d<limitD) fg->SetPoint(fi++,x,y);
    if(x<2000)
    {
      if(n2k<20)
      {
        fg->SetPoint(fi++,x,y);
        n2k++;
      }
    }
    else
      fg->SetPoint(fi++,x,y);
  }
  double chi2;
  int ndf;
  if(fg->GetN()>4)
  {
    fg->Fit("pol1","Q");
    p0 = fg->GetFunction("pol1")->GetParameter(0);
    p1 = fg->GetFunction("pol1")->GetParameter(1);
    e0 = fg->GetFunction("pol1")->GetParError(0);
    e1 = fg->GetFunction("pol1")->GetParError(1);
    chi2 = fg->GetFunction("pol1")->GetChisquare();
    ndf = fg->GetFunction("pol1")->GetNDF();
  }
  else
  {
    p0 = b;
    p1 = k;
    e0 = 100;
    e1 = 10;
    chi2 = 100;
  }
  printf("  - fit  %02d,%02d : k= %10.4f, b= %8.2f, chi2= %8.0f\n",fch,bch,p1,p0,chi2);
  char titleBuff[40];
  sprintf(titleBuff,"F%02d_%02d_k_%f_b_%f_c2_%f",fch,bch,p1,p0,chi2);
  fg->SetTitle(titleBuff);
}

void snPixel::ReFill2()
{
  int pNum = GetN();
  int fi = 0;
  for(int i=0;i<pNum;i++)
  {
    Double_t x,y;
    Pixel->GetPoint(i,x,y);
    Double_t d = TMath::Abs(y-p1*x-p0);
    Double_t limitD = x*0.03;
    n2k = 0;
    if(d<limitD) 
    {
      if(x<2000)
      {
        if(n2k<20)
        {
          fg2->SetPoint(fi++,x,y);
          n2k++;
        }
      }
      else
        fg2->SetPoint(fi++,x,y);
    }
  }
  double chi2;
  int ndf;
  if(fg2->GetN()>4)
  {
    fg2->Fit("pol1","Q");
    p0 = fg2->GetFunction("pol1")->GetParameter(0);
    p1 = fg2->GetFunction("pol1")->GetParameter(1);
    e0 = fg2->GetFunction("pol1")->GetParError(0);
    e1 = fg2->GetFunction("pol1")->GetParError(1);
    chi2 = fg2->GetFunction("pol1")->GetChisquare();
    ndf = fg2->GetFunction("pol1")->GetNDF();
  }
  else
  {
    p0 = b;
    p1 = k;
    e0 = 100;
    e1 = 10;
    chi2 = 100;
  }
  printf("  - fit2 %02d,%02d : k= %10.4f, b= %8.2f, chi2= %8.0f\n",fch,bch,p1,p0,chi2);
  char titleBuff[40];
  sprintf(titleBuff,"F2%02d_%02d_k_%f_b_%f_c2_%f",fch,bch,p1,p0,chi2);
  fg2->SetTitle(titleBuff);
}
