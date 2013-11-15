#ifndef TURBO_H
#define TURBO_H

#include "ConvCode.h"
#include "simulation.h"
#include <vector>
using std::vector;

class TurboCode
{
public:
	TurboCode(char* filename1,char* filename2,int k);
	void encode(const vector<int>& info, vector<int>& u, vector<int>& p1, vector<int>& p2);
	void decode(vector<double>& Lu, vector<double>& Lp1, vector<double>& Lp2, vector<int>& Lout);
private:
	int			_k;//length of information
	int			_interleaver_len;
	ConvCode	_conv1;
	ConvCode	_conv2;

	vector<int> _interleaveMap;
	vector<int> _deMap;
};
#endif