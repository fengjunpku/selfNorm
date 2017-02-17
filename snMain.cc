#include <iostream>
#include <string>
#include "stdio.h"
#include "stdlib.h"

#include "JunErrors.hh"
#include "snExtract.hh"
#include "TStopwatch.h"

using namespace std;

int main(int argc,char** argv)
{
  TStopwatch watch;
  if(argc != 2 && argc != 4)
    MiaoError("Bad input! Please give a run num!");
  int runNum = atoi(argv[1]);
  TString tele(argv[2]);
  TString dssd(argv[3]);
  snExtract extr(runNum,tele,dssd);
  extr.Load();
  printf("CPU_Time: %f, RealTime: %f\n",watch.CpuTime(),watch.RealTime());
}
