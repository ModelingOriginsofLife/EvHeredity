#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// 3 input columns, 3 output columns
// In the output column, if the index is >=SDEPTH/2 then output 1 else output 0. Should be fully general
#define CYCLEN 2
#define SDEPTH 4
#define SBITS 2
#define REDUNDANCY 4
#define EXCESS 900

#define SLEN (CYCLEN*SDEPTH*2*SBITS*REDUNDANCY + EXCESS)

#define HIGH_ERROR_RATE 1e-3
#define ERROR_RATE 1e-3

char xorKey[SLEN];

void InitXOR()
{
	for (int i=0;i<SLEN;i++)
		xorKey[i] = rand()%2;
}

class FSM
{
	public:
		int links[CYCLEN*SDEPTH*2];
		
		void fromBitstring(char *bstr);
		void runOnBitstring(char *bstr);
};

void FSM::fromBitstring(char *bstr)
{
	int idx;
	
	for (int cyc=0;cyc<CYCLEN;cyc++)
	{
		for (int depth=0;depth<SDEPTH;depth++)
		{
			for (int iftype=0;iftype<2;iftype++)
			{
				int val = 0;
				
				for (int sbit = 0;sbit<SBITS;sbit++)
				{
					int idx = REDUNDANCY*(sbit + SBITS*iftype + SBITS*2*depth + SBITS*2*SDEPTH*cyc);
					
					val += (1<<sbit) * ( bstr[idx] ^ xorKey[idx]);
				}
				
				links[iftype + 2*depth + 2*SDEPTH*cyc] = val;
			}
		}
	}
}

void FSM::runOnBitstring(char *bstr)
{
	int x,y;
	int bitidx = 0;
	
	x = 0; y = 0;
	
	while (bitidx < SLEN)
	{
		int nexty = links[bstr[bitidx] + 2*y + 2*SDEPTH*x];
		
		if (x>=CYCLEN/2)
		{
			bstr[bitidx] = (y>=SDEPTH/2);
			if (rand()%10000000<10000000*ERROR_RATE)
				bstr[bitidx] = rand()%2;
		}
		
		y = nexty;
		x++; bitidx++;
		
		if (x==CYCLEN/2)
		{
			bitidx -= CYCLEN/2;
		}
		
		if (x>=CYCLEN)
		{
			x = 0;
		}
	}	
}

int main(int argc, char **argv)
{
	char bstr[SLEN];
	
	for (int i=0;i<SLEN;i++)
		bstr[i] = rand()%2;
		
	int iter=0;
	double error = 0;
	while (iter<50000000)
	{
		FSM F;
		char etest[EXCESS/REDUNDANCY];
		
		for (int i=0;i<EXCESS/REDUNDANCY;i++)
		{
			etest[i] = rand()%2;
			for (int j=0;j<REDUNDANCY;j++)
			{
				bstr[CYCLEN*SDEPTH*2*SBITS*REDUNDANCY + i*REDUNDANCY + j] = etest[i];
			}
		}
		
		F.fromBitstring(bstr);
		
		for (int bitidx=0;bitidx<SLEN;bitidx++)
		{
			if (rand()%10000000<10000000*HIGH_ERROR_RATE)
				bstr[bitidx] = rand()%2;
		}
		
		F.runOnBitstring(bstr);		
		
		for (int i=0;i<EXCESS/REDUNDANCY;i++)
		{
			for (int j=0;j<REDUNDANCY;j++)
			{
				error += fabs(bstr[CYCLEN*SDEPTH*2*SBITS*REDUNDANCY + i*REDUNDANCY + j] - etest[i]);
			}
		}
		
		iter++;
		
		if (iter%100==0)
		{
			FILE *f = fopen("error.txt","a");
			fprintf(f,"%d %.9g\n",iter,error/(double)(100.0*EXCESS));
			fclose(f);
			error = 0;
			/*FILE *f = fopen("log.txt","a");
			
			for (int i=0;i<SLEN;i++)
				fprintf(f,"%d ", bstr[i]);
			fprintf(f,"\n");
			fclose(f);*/
		}
	}
}
