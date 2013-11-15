#include "simulation.h"
#include <cassert>
#include <iostream>
#include <fstream>
void codeGen(vector<int>& vec)
{
	for (int i=0;i<vec.size();i++)
	{
		vec[i]=randBit();
		//vec[i] = 0;
	}
}

void modulate(const vector<int>& ivec,vector<double>& fvec)
{
	for (int i=0;i<ivec.size();i++)
	{
		if (ivec[i]%2)
		{//1
			fvec[i]=-1.0;
		} 
		else
		{//0
			fvec[i]=1.0;
		}
	}
}


void awgn(const vector<double>& input,vector<double>& output,double sigma)
{
	assert(input.size()==output.size());

	for (int i=0;i<input.size();i++)
	{
		output[i]=input[i]+normrnd(0,sigma);
	}
}



void dispCode(const vector<int>& vec,std::ostream& out)
{
	out<< "<" ;
	for (int i=0;i<vec.size();i++)
	{
		out<<vec[i]<<" ";
	}
	out<< ">" <<std::endl;
}

// void dispCode(const vector<int>& vec)
// {
// 	out<< "<" ;
// 	for (int i=0;i<vec.size();i++)
// 	{
// 		out<<vec[i]<<" ";
// 	}
// 	out<< ">" <<std::endl;
// }

int codeComp(const vector<int>& c1, const vector<int>& c2)
{
	assert( c1.size() == c2.size() );
	int cnt=0;
	for (int i=0;i<c1.size();i++)
	{
		if ( c1[i] != c2[i] )
		{
			cnt++;
		}
	}
	return cnt;
}

void designInterleaver(vector<int>& map,int len,int s)
{
	map.resize(len);

	int i;
	for (i = 0; i < len; i++)
	{
		map[i]=i;
	}
	shuffle(map);
}

void inverseMap(const vector<int>& map,vector<int>& demap)
{
	assert( map.size() == demap.size() );
	for (int i = 0; i < map.size(); i++)
	{
		demap[ map[i] ] = i;
	}
}