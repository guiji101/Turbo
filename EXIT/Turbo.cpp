#include "Turbo.h"
#include <fstream>
#include <cmath>
using std::ofstream;

TurboCode::TurboCode(char *filename1, char *filename2, int k)
:_conv1(filename1,k,true),_conv2(filename2,k,false)
{
	_k = k;
	_interleaver_len = 10000;
	designInterleaver(_interleaveMap,_k,10);//////////////////////////
	_deMap.resize(_interleaveMap.size());
	inverseMap(_interleaveMap,_deMap);
}

void TurboCode::encode(const vector<int>& info, vector<int>& u, vector<int>& p1, vector<int>& p2)
{
	_conv1.encode(info,u,p1);
	vector<int> info_pai(info.size());
	interleave(info, info_pai,_interleaveMap);
	_conv2.encode(info_pai,p2);
}

void TurboCode::decode(vector<double>& Lu, vector<double>& Lp1, vector<double>& Lp2, vector<int>& Lout)
{
	_conv1.clear();
	_conv2.clear();
	_conv1.BCJRinit();
	_conv2.BCJRinit();
	
	vector<double> L21(_k);
	vector<double> L12(_k);
	vector<double> Lu_pai(_k);

	int i;
	for (i = 0; i < _k; i++)
	{
		L21[i] = 0;
		L12[i] = 0;
	}
	vector<double> Lu1(Lu.begin(), Lu.begin()+_k);
	interleave(Lu1,Lu_pai,_interleaveMap);

	int ITER=50;

	_conv2.setFirstIter();
	for (i = 0; i < ITER; i++)
	{
		interleave(L21,_deMap);
		_conv1.BCJR1(Lu,Lp1,L21,L12);
		interleave(L12,_interleaveMap);
		_conv2.BCJR2(Lu_pai,Lp2,L12,L21);

	}
	interleave(L21,_deMap);
	double tmp;
	Lout.resize(_k);
	for (int i = 0; i < _k; i++)
	{
		tmp = Lu[i] + L21[i];
		if (tmp > 0)
		{
			Lout[i] = 0;
		}
		else
		{
			Lout[i] = 1;
		}
	}
}