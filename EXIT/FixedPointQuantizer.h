#ifndef FIXED_POINT_QUANTIZER
#define FIXED_POINT_QUANTIZER

#include <vector>
using std::vector;

class FixedPointQuantizer
{
public:
	FixedPointQuantizer();
	FixedPointQuantizer(int n, int nf);
	void init(int n, int nf);
	double quantify(double val)
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
	void quantifyVec(vector<double>& vec);
private:
	int _n;//number of total bits
	int	_ni;//number of bits for integral part(include 1 bit for sign)
	int	_nf;//number of bits for fractional part
	double _powN;//2^_nf
	double _upperLimit;//max number can be represented
};
#endif