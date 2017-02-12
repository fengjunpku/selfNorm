#include "scExtract.hh"

scExtract::scExtract(int runNo)
{
  char filename[20];
  runNum = runNo;
  sprintf(filename,"data%04d.root",runNo);
  TString rootfile = SCDpath + filename;
  inputFile = TFile::Open(rootfile.Data());
  if(!inputFile || !inputFile->IsOpen())
    MiaoError("Can not Open the input file! ");
  Ks = new TStatistic("Sta_K");
  Bs = new TStatistic("Sta_B");
  Init();
}

scExtract::~scExtract()
{
  inputFile->Close();
  delete inputFile;
  delete Ks;
  delete Bs;
}

void scExtract::Init()
{
  inputFile->GetObject(SCDtreeName.Data(),dtree);
  dtree->SetBranchStatus("*",0);
  dtree->SetBranchStatus(SCDadcName.Data(),1);
  dtree->SetBranchAddress(SCDadcName.Data(),adc);
  for(int i=0;i<SCDiNum;i++)
    for(int j=0;j<SCDjNum;j++)
      rawPix[i][j].Init(i,j);
}

void scExtract::Load(int face,int back)
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
    for(int i=SCDiZero;i<SCDiNum;i++)
    {
      for(int j=SCDjZero;j<SCDjNum;j++)
      {
        double x = adc[face][i];
        double y = adc[back][j];
        if(Select(x,y))
          rawPix[i][j].Fill(x,y);
      }
    }
    }
  Fit();
  Weight();
  char buff[20];
  sprintf(buff,"pixels_data%04d_%02d_%02d.root",runNum,face,back);
  Save(buff);
  sprintf(buff,"pars_data%04d_%02d_%02d.txt",runNum,face,back);
  Output(buff);
  Clear();
}

void scExtract::Write()
{
  for(int i=0;i<SCDiNum;i++)
  {
    for(int j=0;j<SCDjNum;j++)
    {
      rawPix[i][j].Write();
      //printf("_%02d_%02d : %9d\n",i,j,pixels[i][j].GetN());
    }
  }
}

void scExtract::Save(const char* tmpfile)
{
  TFile *tmpf = new TFile(tmpfile,"RECREATE","tmpf");
  tmpf->cd();
  Write();
  tmpf->Close();
}

bool scExtract::Select(double x,double y)
{
  if(x<200||y<200) return false;
  if(x>SCDmaxCh||y>SCDmaxCh) return false;
  if(y>1.2*x+300) return false;
  if(x>1.2*y+300) return false;
  return true;
}

void scExtract::Fit()
{
  cout<<"   =======Fitting======="<<endl;
  for(int i=0;i<SCDiNum;i++)
    for(int j=0;j<SCDjNum;j++)
    {
      rawPix[i][j].FindLine();
      Ks->Fill(rawPix[i][j].k);
      Bs->Fill(rawPix[i][j].b);
    }
}

void scExtract::Clear()
{
  for(int i=0;i<SCDiNum;i++)
    for(int j=0;j<SCDjNum;j++)
    {  
      rawPix[i][j].Clear();
      rawPix[i][j].Init(i,j);
    }
}

void scExtract::Weight()
{
  cout<<"   =======Weight======="<<endl;
  int refBack = SCDjNum/2;
  cout<<" *Front:"<<endl;
  for(int q=0;q<SCDiNum;q++)
  {
    Wf[q] = rawPix[q][refBack].w * 1000;//1000 by default
    Kf[q] = Wf[q] * rawPix[q][refBack].k;
    Bf[q] = Wf[q] * rawPix[q][refBack].b;
    fNode[q] = 0.1;
    for(int i=0;i<SCDiNum;i++)
    {
      if(i==q) continue;
      for(int j=0;j<SCDjNum;j++)
      {
        if(j==refBack) continue;
        if(rawPix[q][j].Drop||rawPix[i][j].Drop||rawPix[i][refBack].Drop) continue;
        fNode[q] += 1;
        Double_t K1 = rawPix[q][j].k;
        Double_t B1 = rawPix[q][j].b;
        Double_t W1 = rawPix[q][j].w;
        Double_t K2 = rawPix[i][j].k;
        Double_t B2 = rawPix[i][j].b;
        Double_t W2 = rawPix[i][j].w;
        Double_t K3 = rawPix[i][refBack].k;
        Double_t B3 = rawPix[i][refBack].b;
        Double_t W3 = rawPix[i][refBack].w;
        Double_t Wtmp = W1*W2*W3/(W1*W2+W2*W3+W3*W1);
        Double_t Ktmp = K1*K3/K2;
        Double_t Btmp = B3+(B1-B2)*K3/K2;
        Kf[q] += Wtmp*Ktmp;
        Bf[q] += Wtmp*Btmp;
        Wf[q] += Wtmp;
      }
    }
    if(Wf[q]>0)
    {
      Kf[q] /= Wf[q];
      Bf[q] /= Wf[q];
    }
    printf("front ch %02d;k: %8.4f, b: %8.2f, w: %16.4f, n: %10.1f\n",q,Kf[q],Bf[q],Wf[q],fNode[q]);
  }
  cout<<" *Back:"<<endl;
  for(int p=0;p<SCDjNum;p++)
  {
    Wb[p] = 1000;  //default 1000
    Kb[p] = Wb[p]; //w*1
    Bb[p] = 0;     //w*0
    bNode[p] = 0.1;
    for(int i=0;i<SCDiNum;i++)
    {
      if(p==refBack) break;
      if(rawPix[i][p].Drop||rawPix[i][refBack].Drop) continue;
      bNode[p] += 1;
      Double_t K1 = rawPix[i][p].k;
      Double_t B1 = rawPix[i][p].b;
      Double_t W1 = rawPix[i][p].w;
      Double_t K2 = rawPix[i][refBack].k;
      Double_t B2 = rawPix[i][refBack].b;
      Double_t W2 = rawPix[i][refBack].w;
      Double_t Wtmp = W1*W2/(W1+W2);
      Double_t Ktmp = K2/K1;
      Double_t Btmp = B2-B1*K2/K1;
      Kb[p] += Wtmp*Ktmp;
      Bb[p] += Wtmp*Btmp;
      Wb[p] += Wtmp;
    }
    if(Wb[p]>0)
    {
      Kb[p] /= Wb[p];
      Bb[p] /= Wb[p];
    }
    printf("back ch %02d;k: %8.4f, b: %8.2f, w: %16.4f, n: %10.1f\n",p,Kb[p],Bb[p],Wb[p],bNode[p]);
  }
}

void scExtract::Output(const char* outfile)
{
  FILE *out;
  out = fopen(outfile,"w+");
  if(!out) {MiaoError("scExtract::Output, Can't open output file !");exit(0);}
  for(int i=0;i<SCDiNum;i++)
    fprintf(out, "f_%02d_ref: %10.4f %8.2f %16.4f %10.1f\n",i,Kf[i],Bf[i],Wf[i],fNode[i]);
  for(int j=0;j<SCDjNum;j++)
    fprintf(out, "b_%02d_ref: %10.4f %8.2f %16.4f %10.1f\n",j,Kb[j],Bb[j],Wb[j],bNode[j]);
  fclose(out);
}
