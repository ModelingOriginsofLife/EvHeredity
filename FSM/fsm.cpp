#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>

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
		void runOnBitstring(char *bstr, char *n1, char *n2);
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
				
				if (iftype==1)
				{
					if (links[1 + 2*depth + 2*SDEPTH*cyc] ==
						links[0 + 2*depth + 2*SDEPTH*cyc])
						links[1 + 2*depth + 2*SDEPTH*cyc] = (val+1)%SDEPTH;
				}
			}
		}
	}	
}

void FSM::runOnBitstring(char *bstr, char *n1, char *n2)
{
	int x,y;
	int bitidx = 0;
	int nidx = 0;
	
	x = 0; y = 0;
	
	while (bitidx < SLEN)
	{
		int bitval = 0;
		int sum = bstr[bitidx];
		
		if (sum>1) bitval=1; else bitval = 0;
		
		int nexty = links[bitval + 2*y + 2*SDEPTH*x];
		
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
			nidx++;
		}
/*		else if (x==2*CYCLEN/3)
		{
			bitidx -= CYCLEN/3;
			nidx++;
		}		*/
		if (x>=CYCLEN)
		{
			x = 0;
			nidx = 0;
		}
	}	
}

vector<char*> Pop;
vector<char*> Proxy;

int main(int argc, char **argv)
{
	for (int j=0;j<1000;j++)
	{
		char *bstr = (char*)malloc(SLEN);
		for (int i=0;i<SLEN;i++)
			bstr[i] = rand()%2;
		Pop.push_back(bstr);

		bstr = (char*)malloc(SLEN);
		for (int i=0;i<SLEN;i++)
			bstr[i] = rand()%2;
		
		Proxy.push_back(bstr);
	}
			
	int iter=0;
	double error = 0;
	while (iter<50000000)
	{
		char etest[EXCESS/REDUNDANCY];
		for (int i=0;i<EXCESS/REDUNDANCY;i++)
		{
			etest[i] = rand()%2;
		}
		
		for (int p=0;p<Pop.size();p++)
		{		
			for (int i=0;i<EXCESS/REDUNDANCY;i++)
			{
				for (int j=0;j<REDUNDANCY;j++)
				{
					Pop[p][CYCLEN*SDEPTH*2*SBITS*REDUNDANCY + i*REDUNDANCY + j] = etest[i];
				}
			}					
			for (int bitidx=0;bitidx<SLEN;bitidx++)
			{
				if (rand()%10000000<10000000*HIGH_ERROR_RATE)
					Pop[p][bitidx] = rand()%2;
			}
		}
		
		for (int p=0;p<Pop.size();p++)
		{
			int j,k;
			FSM F;
			char mStr[SLEN];
			
			j = p-1; if (j<0) j+=Pop.size();
			k = p+1; if (k>=Pop.size()) k-=Pop.size();
			
			for (int i=0;i<SLEN;i++)
			{
				mStr[i] = (Pop[j][i] + Pop[p][i] + Pop[k][i] > 1);
			}
			F.fromBitstring(mStr);
			memcpy(Proxy[p],Pop[p],SLEN);
			F.runOnBitstring(Proxy[p],Pop[j],Pop[k]);
		}
			
		for (int p=0;p<Pop.size();p++)
		{
			memcpy(Pop[p],Proxy[p],SLEN);
			for (int i=0;i<EXCESS/REDUNDANCY;i++)
			{
				for (int j=0;j<REDUNDANCY;j++)
				{
					error += fabs(Pop[p][CYCLEN*SDEPTH*2*SBITS*REDUNDANCY + i*REDUNDANCY + j] - etest[i]);
				}
			}
		}
		
/*		for (int i=0;i<Pop.size()/4;i++)
		{
			int p1=rand()%Pop.size();
			int p2=rand()%Pop.size();
			
			if (p2 != p1)
				memcpy(Pop[p2], Pop[p1], SLEN);			
		}*/
		iter++;
		
		//if (iter%100==0)
		{
			FILE *f = fopen("error.txt","a");
			fprintf(f,"%d %.9g\n",iter,error/(double)(Pop.size()*EXCESS));
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
