#include <iostream>
#include <string>

#include "include/ncollFunctions_5TeV.h"

int relativeWeights()
{
  const int nHiBins = 6;
  const int hiBinsLow[nHiBins] = {0, 20, 60, 100, 100, 140};
  const int hiBinsHi[nHiBins] = {20, 60, 100, 200, 140, 200};

  double ncollSum[nHiBins]= {0, 0, 0, 0, 0, 0};

  for(int i = 0; i < nHiBins; ++i){
    
    for(int j = hiBinsLow[i]; j < hiBinsHi[i]; ++j){
      ncollSum[i] += findNcoll(j);
    }
  }

  for(int i = 0; i < nHiBins; ++i){
    std::cout << hiBinsLow[i] << "-" << hiBinsHi[i] << ": " << ncollSum[i] << std::endl;
  }


  return 0;
}

int main()
{
  relativeWeights();
  return 0;
}
