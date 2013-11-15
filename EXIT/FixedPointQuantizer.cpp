#include "FixedPointQuantizer.h"

FixedPointQuantizer::FixedPointQuantizer()
{
	init(1,1);
}
FixedPointQuantizer::FixedPointQuantizer(int n, int nf)
{
	init(n, nf);
}

void FixedPointQuantizer::init(int n, int nf)
{
	_n = n;
	_nf = nf;
	if (_n > 30)
	{
		_n = 30;
	}
	if (_nf > 20)
	{
		_nf = 20;
	}
	_ni = _n - _nf;
	if (_ni < 0)
	{
		_ni = 0;
		_n = _nf;
	}
	_powN = 1 << _nf;
	_upperLimit = ( 1<<(_ni-1) ) - ( 1.0 / (1<<_nf) );//2^(ni-1) - 1 / 2^nf;
}
/*
double FixedPointQuantizer::quantify(double val)
{
	double ret = int(val * _powN) / double(_powN);
	if (ret > _upperLimit)
	{
		return _upperLimit;
	}
	else if (ret < -_upperLimit)
	{
		return -_upperLimit;
	}
	else
	{
		return ret;
	}
}
*/
void FixedPointQuantizer::quantifyVec(vector<double>& vec)
{
	for (vector<double>::iterator it = vec.begin(); it != vec.end(); it++)
	{
		*it = quantify(*it);
	}
}