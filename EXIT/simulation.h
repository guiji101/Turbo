#ifndef SIMULATION_H
#define SIMULATION_H

#include "myrand.h"
#include <vector>
using std::vector;

void codeGen(vector<int>& vec);

void modulate(const vector<int>& ivec,vector<double>& fvec);

void awgn(const vector<double>& input,vector<double>& output,double sigma);

void dispCode(const vector<int>& vec,std::ostream& out);

//void dispCode(const vector<int>& vec);

int codeComp(const vector<int>& c1, const vector<int>& c2);

void designInterleaver(vector<int>& map,int len,int s);

void inverseMap(const vector<int>& map,vector<int>& demap);

template<typename T>
void interleave(const vector<T>& input,vector<T>& output,vector<int>& map)
{
	int interleaveLen=map.size();
	int inputLen=input.size();
	assert( output.size() == inputLen );
	//assert( inputLen%interleaveLen == 0 );
	int i,j;
	int T=inputLen/interleaveLen;
	for(i = 0; i < T; i++)
	{
		for (j = 0; j < interleaveLen; j++)
		{
			output[ i*interleaveLen + map[j] ] 
				= input[  i*interleaveLen + j ] ;
		}
	}
}

template<typename T>
void interleave(vector<T>& input,vector<int>& map)
{
	vector<T> tmpvec(input);
// 	for (int i = 0; i < input.size(); i++)
// 	{
// 		input[ map[i] ] = tmpvec[i];
// 	}
	interleave(input,tmpvec,map);
	for (int i = 0; i < input.size(); i++)
	{
		input[i] = tmpvec[i];
	}
}

#endif