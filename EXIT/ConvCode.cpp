#include "ConvCode.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>

#define DO_QUANTIFICATION
const double EPSILON = 1e-6;
inline int getLow(int val)
{
	return val & 1u;
}

inline int getHigh(int val)
{
	return ( (val>>1) & 1u );
}

void ConvCode::initialize()
{
	int i;
	for (i = 0; i < NSTATE; i++)
	{
		_states[i]._next[0]	= _nxtState[0][i];
		_states[i]._next[1]	= _nxtState[1][i];
		_states[ _nxtState[0][i] ]._forward[0] = i;
		_states[ _nxtState[1][i] ]._forward[1] = i;
	}
	for (i = 0; i < NSTATE; i++)
	{
		_links[i]._sl			= i;
		_links[i]._sr			= _nxtState[0][i];
		_links[i]._input		= 0;
		_links[i]._output[0]	= 1.0;//system bits
		_links[i]._output[1]	= 1.0-2*_output[0][i];//parity check bits
	}
	for (i = 0; i < NSTATE; i++)
	{
		_links[i+NSTATE]._sl		= i;
		_links[i+NSTATE]._sr		= _nxtState[1][i];
		_links[i+NSTATE]._input		= 1;
		_links[i+NSTATE]._output[0]	= -1.0;//system bits
		_links[i+NSTATE]._output[1]	= 1.0-2*_output[1][i];//parity check bits
	}

	int ub=0;
	if (_dotailbitting)
	{
		ub = _k+MU;
	}
	else
	{
		ub = _k;
	}
	
	_alpha = new double* [ub+1];
	for (i = 0; i < ub+1; i++)
	{
		_alpha[i] = new double [NSTATE];
	}
	
	_beta = new double* [ub+1];
	for (i = 0; i < ub+1; i++)
	{
		_beta[i] = new double [NSTATE];
	}
	_gama = new double* [ub];
	for (i = 0; i < ub; i++)
	{
		_gama[i] = new double [NSTATE*NSTATE];
	}
}

void ConvCode::clear()
{
	int i;
	int ub = 0;
	if (_dotailbitting)
	{
		ub = _k+MU;
	}
	else
	{
		ub = _k;
	}
	for (i = 0; i < ub+1; i++)
	{
		memset(_alpha[i],0,sizeof(double)*NSTATE);
	}
	for (i = 0; i < ub+1; i++)
	{
		memset(_beta[i],0,sizeof(double)*NSTATE);
	}
	for (i = 0; i < ub; i++)
	{
		memset(_gama[i],0,sizeof(double)*NSTATE*NSTATE);
	}
}

ConvCode::ConvCode(char* filename, int k, bool dotailbitting)
{
	
	std::ifstream in(filename);
	if (!in)
	{
		std::cout<<"can not open "<<filename<<"\n";
		system("PAUSE");
		exit(0);
	}
	
	in >> MU;
	NSTATE = 1 << MU;
	buildS();
	_k				= k;
	_dotailbitting	= dotailbitting;

	//allocate memory
	_states		= new State [NSTATE];
	_links		= new Link [2*NSTATE];
	_tailBit	= new int [NSTATE];

	int i,j;
	
	_nxtState = new int* [NINPUT];
	for (i = 0; i < NINPUT; i++)
	{
		_nxtState[i] = new int [NSTATE];
	}

	_output = new int* [NINPUT];
	for (i = 0; i < NINPUT; i++)
	{
		_output[i] = new int [NSTATE];
	}

	//read trellis data
	for (i = 0; i < NINPUT; i++)
	{
		for (j = 0; j < NSTATE; j++)
		{
			in >> _nxtState[i][j];
		}
	}
	for (i = 0; i < NINPUT; i++)
	{
		for (j = 0; j < NSTATE; j++)
		{
			in >> _output[i][j];
		}
	}
	for (i = 0; i < NSTATE; i++)
	{
		in>>_tailBit[i];
	}

	initialize();

	clear();
	
#ifdef DO_QUANTIFICATION
	std::ifstream quantizerConfig("quantizer_config.txt");
	int n, nf;
	quantizerConfig>>n>>nf;
	_rQuantizer.init(n, nf);
	quantizerConfig>>n>>nf;
	_LcQuantizer.init(n, nf);
	quantizerConfig>>n>>nf;
	_LeQuantizer.init(n, nf);
	quantizerConfig>>n>>nf;
	_stateQuantizer.init(n, nf);
	quantizerConfig>>n>>nf;
	_gammaQuantizer.init(n, nf);
#endif
}

ConvCode::~ConvCode()
{
	int i;

	delete [] _states;
	delete [] _links;
	delete [] _tailBit;

	for (i = 0; i < NINPUT; i++)
	{
		delete [] _nxtState[i];
	}
	delete [] _nxtState;

	for (i = 0; i < NINPUT; i++)
	{
		delete [] _output[i];
	}
	delete [] _output;

	int ub = 0;
	if (_dotailbitting)
	{
		ub = _k+MU;
	}
	else
	{
		ub = _k;
	}
	for (i = 0; i < ub+1; i++)
	{
		delete [] _alpha[i];
	}
	delete [] _alpha;
	for (i = 0; i < ub+1; i++)
	{
		delete [] _beta[i];
	}
	delete [] _beta;
	for (i = 0; i < ub; i++)//////////////////////
	{
		delete [] _gama[i];
	}
	delete [] _gama;
}

void ConvCode::encode(const vector<int>& input,vector<int>& u, vector<int>& p)
{
	assert( input.size() == _k );
	if (_dotailbitting)
	{
		u.resize(_k+MU);
		p.resize(_k+MU);
	}
	else
	{
		u.resize(_k);
		p.resize(_k);
	}
	int i;
	int curState = 0;
	int po = 0;

	for (i = 0; i < _k; i++)
	{
		u[po]		= input[i];
		p[po]		= _output[ input[i] ][curState];
		curState	= _nxtState[ input[i] ][curState];
		po++;
	}
	//do tailbitting
	
	int tail	= _tailBit[curState];
	int tailbit	= 0;
	for(i = 0; i < MU; i++)
	{
		tailbit = tail%2;
		tail /= 2;
		u[po]		= tailbit;
		p[po]		= _output[ tailbit ][curState];
		curState	= _nxtState[ tailbit ][curState];
		po++;
	}
	assert( curState == 0 );///////////////////////////////////
	
}

void ConvCode::encode(const vector<int>& input, vector<int>& p)
{//no tailbitting
	//assert( input.size() == _k );

	p.resize(_k);

	int i;
	int curState = 0;
	int po = 0;

	for (i = 0; i < _k; i++)
	{
		p[po]		= _output[ input[i] ][curState];
		curState	= _nxtState[ input[i] ][curState];
		po++;
	}
}


void ConvCode::BCJRinit()
{
	int i;
	alpha(0,0) = 0;
	for (i = 1; i < NSTATE; i++)
	{
		alpha(0,i) = -INF;
	}

	if (_dotailbitting)
	{
		beta(_k+MU,0) = 0;
		for (i = 1; i < NSTATE; i++)
		{
			beta(_k+MU,i) = -INF;
		}
	}
	else
	{
		for (i = 0; i < NSTATE; i++)
		{
			beta(_k, i) = 1.0 / NSTATE;
		}
	}
}

void ConvCode::caclGama(const vector<double>& Lu, const vector<double>& Lp, const vector<double>& La)
{
	int i,j;
	int sl,sr;
	for (i = 0; i < _k; i++)
	{
		for (j = 0; j < 2*NSTATE; j++)
		{
			sl = _links[j]._sl;
			sr = _links[j]._sr;
#ifdef DO_QUANTIFICATION
			double d1 = _gammaQuantizer.quantify(La[i]*_links[j]._output[0]);
			double d2 = _gammaQuantizer.quantify(Lu[i]*_links[j]._output[0]);
			double d3 = _gammaQuantizer.quantify(Lp[i]*_links[j]._output[1]);
			double d4 = _gammaQuantizer.quantify(d1 + d2);
			double d5 = _gammaQuantizer.quantify(d4 + d3);
			gama(i,sl,sr) = _gammaQuantizer.quantify(d5 * 0.5);
#else	
			gama(i,sl,sr) =
				( La[i]*_links[j]._output[0] + Lu[i]*_links[j]._output[0] + Lp[i]*_links[j]._output[1] ) *0.5 ;
#endif		
		}
	}
	//for tail bits, no Lin for them
	//gama=u*La/2 + (Lc/2)r*v  
	//    and (Lc/2)r*v  (for tail bits)
	if (_dotailbitting)
	{
		for (i = _k; i < _k+MU; i++)
		{
			for (j = 0; j < 2*NSTATE; j++)
			{
				sl = _links[j]._sl;
				sr = _links[j]._sr;
#ifdef DO_QUANTIFICATION
				double d1 = _gammaQuantizer.quantify(Lu[i]*_links[j]._output[0]);
				double d2 = _gammaQuantizer.quantify(Lp[i]*_links[j]._output[1]);
				double d3 = _gammaQuantizer.quantify(d1 + d2);
				gama(i,sl,sr) = _gammaQuantizer.quantify(d3 * 0.5);
#else
			gama(i,sl,sr) =( Lu[i]*_links[j]._output[0] + Lp[i]*_links[j]._output[1] ) * 0.5;
#endif
				
			}
		}
	}
	
}
void ConvCode::caclAlpha()
{
	int i,j;
	double tmpd1,tmpd2;
	int sl,sr;
	int ub=0;
	if (_dotailbitting)
	{
		ub = _k+MU;
	}
	else
	{
		ub = _k;
	}
	for (i = 1; i <= ub; i++)
	{
		for (j = 0; j < NSTATE; j++)
		{
			sr		= j;
			
			sl		= _states[j]._forward[0];
#ifdef DO_QUANTIFICATION
			tmpd1	= plus( alpha(i-1,sl) , _stateQuantizer.quantify(gama(i-1,sl,sr)) ) ;
#else
			tmpd1	= plus( alpha(i-1,sl) , gama(i-1,sl,sr) ) ;
#endif
			
			
			sl		= _states[j]._forward[1];
#ifdef DO_QUANTIFICATION
			tmpd2	= plus( alpha(i-1,sl) , _stateQuantizer.quantify(gama(i-1,sl,sr)) ) ;
#else
			tmpd2	= plus( alpha(i-1,sl) , gama(i-1,sl,sr) ) ;
#endif
				
			alpha(i,j) = maxAsterisk(tmpd1,tmpd2);

		}
		//normalize
	
		double minVal = alpha(i,0);
		for (j = 1; j < NSTATE; j++)
		{
			if (alpha(i,j) < minVal)
			{
				minVal = alpha(i,j);
			}
		}
		if (minVal < -INF + EPSILON)
		{
			continue;
		}
		for (j = 0; j < NSTATE; j++)
		{
			alpha(i,j) -= minVal;
#ifdef DO_QUANTIFICATION
			alpha(i,j) = _stateQuantizer.quantify(alpha(i,j));	
#endif
		}
	}
}

void ConvCode::caclBeta()
{
	int i,j;
	double tmpd1,tmpd2;
	int sl,sr;
	int ub = 0;
	if (_dotailbitting)
	{
		ub = _k+MU;
	} 
	else
	{
		ub = _k;
	}
	for (i = ub-1; i >= 0; i--)
	{
		for (j = 0; j < NSTATE; j++)
		{
			sl		= j;
			
			sr		= _states[j]._next[0];
#ifdef DO_QUANTIFICATION
			tmpd1	= plus( beta(i+1,sr), _stateQuantizer.quantify(gama(i,sl,sr)) ) ;
#else
			tmpd1	= plus( beta(i+1,sr),gama(i,sl,sr) ) ;
#endif
			
			
			sr		= _states[j]._next[1];
#ifdef DO_QUANTIFICATION
			tmpd2	= plus( beta(i+1,sr), _stateQuantizer.quantify(gama(i,sl,sr)) ) ;
#else
			tmpd2	= plus( beta(i+1,sr),gama(i,sl,sr) ) ;
#endif
			
			beta(i,j) = maxAsterisk(tmpd1,tmpd2);
		}
		//normalize
		double minVal = beta(i,0);
		for (j = 1; j < NSTATE; j++)
		{
			if (beta(i,j) < minVal)
			{
				minVal = beta(i,j);
			}
		}
		if (minVal < -INF + EPSILON)
		{
			continue;
		}
		for (j = 0; j < NSTATE; j++)
		{
			beta(i,j) -= minVal;
#ifdef QUANTIFICATION
			beta(i,j) = _stateQuantizer.quantify(beta(i,j));
#endif
		}
	}
}



void ConvCode::caclL(vector<double>& Lout)
{
	int i,j;
	double pos,neg;
	double sum;
	int sl,sr;
	for (i = 0; i < _k; i++)
	{
		pos = -INF;
		neg = -INF;

		for (j = 0; j < NSTATE; j++)
		{
			sl	= _links[j]._sl;
			sr	= _links[j]._sr;
#ifdef DO_QUANTIFICATION
			double d1 = _LeQuantizer.quantify(beta(i+1,sr));
			double d2 = _LeQuantizer.quantify(gama(i,sl,sr));
			double d3 = _LeQuantizer.quantify(alpha(i,sl)); 
			double d4 = _LeQuantizer.quantify(d1 + d2);
			sum = _LeQuantizer.quantify(d4 + d3);
			pos = _LeQuantizer.quantify(maxAsterisk(pos,sum));
#else
			sum	= beta(i+1,sr) + gama(i,sl,sr) + alpha(i,sl);
			pos = maxAsterisk(pos,sum);
#endif
		}

		for (j = NSTATE; j < 2*NSTATE; j++)
		{
			sl	= _links[j]._sl;
			sr	= _links[j]._sr;
#ifdef DO_QUANTIFICATION
			double d1 = _LeQuantizer.quantify(beta(i+1,sr));
			double d2 = _LeQuantizer.quantify(gama(i,sl,sr));
			double d3 = _LeQuantizer.quantify(alpha(i,sl));
			double d4 = _LeQuantizer.quantify(d1 + d2);
			sum = _LeQuantizer.quantify(d4 + d3);
			neg = _LeQuantizer.quantify(maxAsterisk(neg,sum));
#else
			sum	= beta(i+1,sr) + gama(i,sl,sr) + alpha(i,sl);
			neg = maxAsterisk(neg,sum);
#endif
		}
#ifdef DO_QUANTIFICATION
		Lout[i] = _LeQuantizer.quantify(pos - neg);
#else
		Lout[i] = pos - neg;
#endif
	}	
}
void ConvCode::caclLe(const vector<double>& Lp,vector<double>& Lout)
{
	int i,j;
	double pos,neg;
	double sum;
	int sl,sr;
	//double halfLc=1.0/(sigma*sigma);
	for (i = 0; i < _k; i++)
	{
		pos = -INF;
		neg = -INF;

		for (j = 0; j < NSTATE; j++)
		{
			sl	= _links[j]._sl;
			sr	= _links[j]._sr;
#ifdef DO_QUANTIFICATION
			double d1 = _LeQuantizer.quantify(beta(i+1,sr));
			double d2 = _LeQuantizer.quantify(0.5*Lp[i]*_links[j]._output[1]);
			double d3 = _LeQuantizer.quantify(alpha(i,sl));
			double d4 = _LeQuantizer.quantify(d1 + d2);
			sum = _LeQuantizer.quantify(d4 + d3);
			pos = _LeQuantizer.quantify(maxAsterisk(pos,sum));
#else
			sum	= beta(i+1,sr) + 0.5*Lp[i]*_links[j]._output[1] + alpha(i,sl);
			pos = maxAsterisk(pos,sum);
#endif
		}

		for (j = NSTATE; j < 2*NSTATE; j++)
		{
			sl	= _links[j]._sl;
			sr	= _links[j]._sr;
#ifdef DO_QUANTIFICATION
			double d1 = _LeQuantizer.quantify(beta(i+1,sr));
			double d2 = _LeQuantizer.quantify(0.5*Lp[i]*_links[j]._output[1]);
			double d3 = _LeQuantizer.quantify(alpha(i,sl));
			double d4 = _LeQuantizer.quantify(d1 + d2);
			sum = _LeQuantizer.quantify(d4 + d3);
			neg = _LeQuantizer.quantify(maxAsterisk(neg,sum));
#else
			sum	= beta(i+1,sr) + 0.5*Lp[i]*_links[j]._output[1] + alpha(i,sl);
			neg = maxAsterisk(neg,sum);
#endif
		}

#ifdef DO_QUANTIFICATION
		Lout[i] = _LeQuantizer.quantify(pos - neg);
#else
		Lout[i] = pos - neg;
#endif
	}
}
void ConvCode::BCJR1(vector<double>& Lu, vector<double>& Lp, vector<double>& La,
					vector<double>& Lout)
{

	Lout.resize(La.size());

#ifdef DO_QUANTIFICATION
	_rQuantizer.quantifyVec(Lu);
	_rQuantizer.quantifyVec(Lp);
	_LeQuantizer.quantifyVec(La);
#endif

	caclGama(Lu,Lp,La);
	caclAlpha();
	caclBeta();
	//caclL(Lout);
	caclLe(Lp,Lout);
}

void ConvCode::BCJR2(vector<double>& Lu, vector<double>& Lp, vector<double>& La,
					 vector<double>& Lout)
{
	/*
	if (_isFirstIter)
	{
		for (int i = 0; i < NSTATE; i++)
		{
			beta(_k,i) = alpha(_k,i);
		}
		_isFirstIter = false;
	}*/
	
	Lout.resize(La.size());
#ifdef DO_QUANTIFICATION
	_rQuantizer.quantifyVec(Lu);
	_rQuantizer.quantifyVec(Lp);
	_LeQuantizer.quantifyVec(La);
#endif

	caclGama(Lu,Lp,La);
	caclAlpha();
	caclBeta();
	//caclL(Lout);
	caclLe(Lp,Lout);
}



void ConvCode::buildS()
{
	int i;
	double x;
	for (i = 0; i < _array_sz; i++)
	{
		x		= double(i) * DELTA;
		_s[i]	= log( 1+exp(-x) );
	}
}

double ConvCode::maxAsterisk(double x,double y)
{
	if ( x <= -INF)
	{
		return y;
	}
	if ( y <= -INF)
	{
		return x;
	}
	double m	= x>y ? x : y;
	int n		= int( fabs(x-y) * FACTOR );/////////?
	//double a	= double(647)*DELTA;
	double s	= 0;
	if ( n < _array_sz )
	{
		s = _s[n];
	}
#ifdef DO_QUANTIFICATION
	s = _stateQuantizer.quantify(s);////////improve
	return _stateQuantizer.quantify(m+s);//////////improve
#else 
	return m + s;
#endif
}
