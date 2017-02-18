#include "snPixel.hh"

snPixel::snPixel()
{
  fch = bch = -1;
  p0 = 0;
  p1 = 1;
  e0 = 1000;
  e1 = 10;
  w0 = 0;
  w1 = 0;
  Drop = true;
}

snPixel::~snPixel()
{
  if(rg) delete rg;
  if(fg) delete fg;
  if(fg2) delete fg2;
}

void snPixel::Init(int frontCh,int backCh)
{
  if(frontCh<0||frontCh>32||backCh<0||backCh>32)
    MiaoError("snPixel::snPixel, bad parameters !");
  fch = frontCh;
  bch = backCh;
  char NameBuff[20];
  sprintf(NameBuff,"raw_%02d_%02d",fch,bch);
  rg = new TGraph();
  rg->SetName(NameBuff);
  rg->SetDrawOption("AP*");
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
  rg->SetPoint(GetN(),x,y);
}

void snPixel::Write()
{
  if(!rg) {MiaoError("snPixel::Write(), Graph is null !");exit(0);}
  //printf("%02d_%02d : Slope : %4f+-%4f; Intercept : %4f+-%4f\n",fch,bch,Slope,ErSlope,Intercept,ErIntercept);
  rg->Write();
  fg->Write();
  fg2->Write();
}

int snPixel::GetN()
{
  if(!rg) {MiaoError("snPixel::GetN(), Graph is null !");exit(0);}
  return rg->GetN();
}

void snPixel::Clear()
{
  rg->Clear();
  fg->Clear();
  fg2->Clear();
}

void snPixel::Fit()
{
  FitG(rg);
  printf(" -- fit %02d,%02d  ------  p0  --------  e0  -------  p1  -------  e1  -- \n",fch,bch);
  printf("                     %4.6f +/- %4.6f,|  %8.4f +/- %8.4f\n",p1,e1,p0,e0);
  char titleBuff[40];
  sprintf(titleBuff,"F0_%02d,%02d_k_%f_b_%f",fch,bch,p1,p0);
  rg->SetTitle(titleBuff);
}

void snPixel::ReFill()
{
  int pNum = GetN();
  int fi = 0;
  for(int i=0;i<pNum;i++)
  {
    Double_t x,y;
    rg->GetPoint(i,x,y);
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
  FitG(fg);
  //printf("  - fit  %02d,%02d : k= %10.4f, b= %8.2f, chi2= %8.0f\n",fch,bch,p1,p0,chi2);
  printf("                     %4.6f +/- %4.6f,|  %8.4f +/- %8.4f\n",p1,e1,p0,e0);
  char titleBuff[40];
  sprintf(titleBuff,"F1_%02d,%02d_k_%f_b_%f",fch,bch,p1,p0);
  fg->SetTitle(titleBuff);
}

void snPixel::ReFill2()
{
  int pNum = GetN();
  int fi = 0;
  for(int i=0;i<pNum;i++)
  {
    Double_t x,y;
    rg->GetPoint(i,x,y);
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
  Drop = !FitG(fg2);
  //printf("  - fit2 %02d,%02d : k= %10.4f, b= %8.2f, chi2= %8.0f\n",fch,bch,p1,p0,chi2);
  printf("                     %4.6f +/- %4.6f,|  %8.4f +/- %8.4f\n",p1,e1,p0,e0);
  char titleBuff[40];
  sprintf(titleBuff,"F2_%02d,%02d_k_%f_b_%f",fch,bch,p1,p0);
  fg2->SetTitle(titleBuff);
}

bool snPixel::FitG(TGraph *g)
{
  if(g->GetN()>3)
  {
    g->Fit("pol1","Q");
    p0 = g->GetFunction("pol1")->GetParameter(0);
    p1 = g->GetFunction("pol1")->GetParameter(1);
    e0 = g->GetFunction("pol1")->GetParError(0);
    e1 = g->GetFunction("pol1")->GetParError(1);
    w0 = 1./e0/e0;
    w1 = 1./e1/e1;
    return true;
  }
  return false;
}