/***********************
Description:rate 1/2 convolutional code 
************************/

#ifndef CONVCODE_H
#define CONVCODE_H

#include "FixedPointQuantizer.h"
#include <vector>
#include <cassert>
#include <iostream>
using std::vector;

const double INF	= 1e6;

const int ARRY_SZ	= 5000;
const double DELTA	= 0.01;
const double FACTOR	= 1/DELTA;

class Link
{
	int		_sl;
	int		_sr;
	int		_input;
	double		_output[2];
	friend class ConvCode;
};

class State
{
	int		_next[2];//next state
	int		_forward[2];//forward state
	friend class ConvCode;
};

class ConvCode
{
public:
	ConvCode(char* filename, int k, bool dotailbitting);
	~ConvCode();

	void initialize();	//build trellis
	void clear();		//clear _alpha, _beta, _gama
	void encode(const vector<int>& input,vector<int>& u, vector<int>& p);// do tailbitting
	void encode(const vector<int>& input, vector<int>& p);
	void BCJR1(vector<double>& Lu, vector<double>& Lp, vector<double>& La
		,vector<double>& Lout);
	void BCJR2(vector<double>& Lu, vector<double>& Lp, vector<double>& La
		,vector<double>& Lout);
	void BCJRinit();

	void setFirstIter() { _isFirstIter = true; }
	//void clearFirstIter() { _isFirstIter = false; }
private:
	//enum{NINPUT=2,NSTATE=4};
	static const int NINPUT	=2;//number of all possible inputs

	int MU;		//number of registers (number of tail bits)
	int NSTATE;	//number of states in the trellis (=2^MU)
	
	int			_k;//number of information bits, tail bits are not included

	bool		_dotailbitting;

	State*	_states;//len=NSTATE
	Link*	_links;//len=2*NSTATE
	int**	_nxtState;//NINPUT * NSTATE
	int**	_output;//NINPUT * NSTATE
	int*	_tailBit;//len=NSTATE

	double**	_alpha;// LEN * NSTATE
	double**	_beta;
	double**	_gama;

	bool _isFirstIter;

	FixedPointQuantizer _rQuantizer;
	FixedPointQuantizer _LcQuantizer;
	FixedPointQuantizer _LeQuantizer;
	FixedPointQuantizer _stateQuantizer;
	FixedPointQuantizer _gammaQuantizer;

	double& alpha(const int i,const int s)
	{
		assert( 0<=i && i<=_k+MU );
		assert( 0<=s && s<NSTATE );
		return _alpha[i][s];
	}
	double& beta(const int i,const int s)
	{
		assert( 0<=i && i<=_k+MU );
		assert( 0<=s && s<NSTATE );
		return _beta[i][s];
	}
	double& gama(const int i,const int sl,const int sr)
	{
		
		assert( 0<=i && i<_k+MU );
		assert( 0<=sl && sl<NSTATE );
		assert( 0<=sr && sr<NSTATE );
		return _gama[i][ sl*NSTATE + sr ];
	}

	
	void caclGama(const vector<double>& Lu, const vector<double>& Lp, const vector<double>& La);
	void caclAlpha();
	void caclBeta();
	void caclL(vector<double>& Lout);
	void caclLe(const vector<double>& Lp,vector<double>& Lout);
	
	
	static double plus(const double x,const double y)
	{
		if ( x <= -INF )
		{
			return -INF;
		} 
		else
		{
			return x+y;
		}
	}
	void buildS();
	double maxAsterisk(double x,double y);//max*(x,y)
	enum {_array_sz=ARRY_SZ};
	double _s[_array_sz];

	friend class TurboCode;
};

#endif