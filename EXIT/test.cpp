#include "myrand.h"
#include "simulation.h"
#include "ConvCode.h"
#include "Turbo.h"

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdio>

using namespace std;

void toLLR(vector<double>& received, vector<double>& LLR, double sigma_ch);

int main()
{
	int len=10000;
	int v = 3;

	vector<int> infoBits(len);
	vector<int> codedU(len+v);
	vector<int> codedP1(len+v);
	vector<int> codedP2(len);
	vector<double> modulatedU(len+v);
	vector<double> modulatedP1(len+v);
	vector<double> modulatedP2(len);
	vector<double> receivedU(len+v);
	vector<double> receivedP1(len+v);
	vector<double> receivedP2(len);
	vector<int> decoded(len);
	vector<double> Lu(len+v);
	vector<double> Lp1(len+v);
	vector<double> Lp2(len);
	

	TurboCode turbo("(13,15)trellis_data.txt","(13,15)trellis_data.txt",len);
	FILE * fp = fopen("result.txt","w");

	for (double EbN0_indB = -0.3; EbN0_indB <= 1.0; EbN0_indB += 0.1)
	{
		double EbN0 = pow(10.0,EbN0_indB/10);
		double sigma_ch = sqrt(1.5/EbN0);

		double blks = 0;
		double biterr = 0;
		double blkerr = 0;
		int cnt = 0;

		while (1)
		{
			codeGen(infoBits);
			turbo.encode(infoBits,codedU,codedP1,codedP2);
			modulate(codedU,modulatedU);
			modulate(codedP1,modulatedP1);
			modulate(codedP2,modulatedP2);
			awgn(modulatedU,receivedU,sigma_ch);
			awgn(modulatedP1,receivedP1,sigma_ch);
			awgn(modulatedP2,receivedP2,sigma_ch);
			toLLR(receivedU, Lu, sigma_ch);
			toLLR(receivedP1, Lp1, sigma_ch);
			toLLR(receivedP2, Lp2, sigma_ch);

			turbo.decode(Lu,Lp1,Lp2,decoded);

			cnt = 0;
			for (int i = 0; i < len; i++)
			{
				if (decoded[i] != infoBits[i])
				{
					cnt++;
				}
			}

			biterr += cnt;
			blks += 1;
			if (cnt > 0)
			{
				blkerr += 1.0;
			}

			printf("EbN0=%lf\tBLKS=%lf\tBER=%lf\tFER=%lf\n",EbN0_indB,blks,(biterr/(len*blks)),blkerr/blks);
			if (blkerr >= 100 || blks >= 10000)
			{
				fprintf(fp,"EbN0=%lf\tBLKS=%lf\tBER=%lf\tFER=%lf\n",EbN0_indB,blks,(biterr/(len*blks)),blkerr/blks);
				fflush(fp);
				break;
			}
		}
	}
	fclose(fp);
	system("PAUSE");
	return 0;
}

void toLLR(vector<double>& received, vector<double>& LLR, double sigma_ch)
{
	LLR.resize(received.size());
	int sz = received.size();
	double factor = 2.0 / (sigma_ch * sigma_ch);
	for (int i = 0; i < sz; i++)
	{
		LLR[i] = received[i] * factor;
	}
}