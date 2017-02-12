#include <iostream>
#include <string>
#include "stdio.h"
#include "stdlib.h"

#include "JunErrors.hh"
#include "selCali_def.hh"
#include "scRawPixel.hh"
#include "scExtract.hh"

#include "TStopwatch.h"

using namespace std;
using namespace selCaliDef;

int main(int argc,char** argv)
{
  TStopwatch watch;
  if(argc != 2 && argc != 4)
    MiaoError("Bad input! Please give a run num!");
  int runNum = atoi(argv[1]);
  int fn,bn;
  if(argc == 2)
  {
   fn = SCDfaceNum;
   bn = SCDbackNum;
  }
  else
  {
    fn = atoi(argv[2]);
    bn = atoi(argv[3]);
  }
  scExtract extr(runNum);
  extr.Load(fn,bn);
  extr.Load(bn,fn);
  printf("CPU_Time: %f, RealTime: %f\n",watch.CpuTime(),watch.RealTime());
}
