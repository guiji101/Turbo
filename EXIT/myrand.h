#ifndef MYRAND_H
#define MYRAND_H

#include <vector>

using std::vector;

void setSeed(unsigned int seed);

unsigned long int randInt();

unsigned int randBit();

double random();

double rayrnd();

double normrnd(double mu,double sigma);

long int randInRange(long int a, long int b);//rand int in [a,b)

void shuffle(vector<int>& vec);

#endif