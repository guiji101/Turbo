#include "myrand.h"
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>///////////////////
#include <cassert>

#define _MAX_ 2147483647L
#define PI 3.1415926535898
unsigned long int _SEED_=1;

void setSeed(unsigned int seed)
{
	_SEED_=seed;
}

unsigned long int randInt()
{
	//unsigned int _MAX_=2147483647;
	unsigned long int a=16807;
	static unsigned long int z=_SEED_;
	return z=(z*a)%_MAX_;
}


unsigned int randBit()
{
	return randInt()%2;
}

double random()
{
	//unsigned int M=2147483647;
	return double(randInt())/_MAX_;
}

long int randInRange(long int a, long int b)
{
	assert(b > a);
	double r = random();
	return a + floor((b-a) * r);
}

double rayrnd()
{
	return sqrt((-2)*log(random()));
}

double normrnd(double mu,double sigma)
{
	double ray=sqrt((-2)*log(random()));
	double fai=(2*random()-1)*PI;
	return mu+sigma*ray*cos(fai);
}

void shuffle(vector<int>& vec)
{
	int len = vec.size();
	for (int i = 0; i < len - 1; i++)
	{
		int r = randInRange(i, len);
		int tmp = vec[i];
		vec[i] = vec[r];
		vec[r] = tmp;
	}
}