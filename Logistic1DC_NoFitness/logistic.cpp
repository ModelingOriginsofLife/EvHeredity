#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <SDL/SDL.h>
#include <SDL/SDL_image.h>

#include "fblib.h"
#include "imgload.h"
#include "inputs.h"

#include "globals.h"

class Site
{
	public:
		float logis[NVAR], targ[NVAR];
		int bits[NVAR];
		float fixed[2];
		float W, persist;
};

float logistic(float x)
{
	return LOGISTIC_R * x * (1.0f - x);
}

Site *Grid;

void iterateLogistic(Site *Grid, int idx)
{   
	if (idx<NSITES*NVAR)
	{   
		int site = idx/NVAR;
		int var = idx%NVAR;
		
		for (int i=0;i<NITER;i++)
		{
			Grid[site].logis[var] = logistic( Grid[site].logis[var] );
		}
		
		float d1 = abs(Grid[site].logis[var] - Grid[site].fixed[0]);
		float d2 = abs(Grid[site].logis[var] - Grid[site].fixed[1]);
		
		Grid[site].bits[var] = (d2<d1);
	}
}

void getFixed(Site *Grid, int idx)
{
	if (idx<NSITES*3)
	{   
		int fid = idx/NSITES;
		int sid = idx%NSITES;
		int result = 0;
		
		int bidx = fid*NBITS;
		float dist = 0;
				
		for (int i=0;i<NBITS;i++)
		{
			result += (1<<(NBITS-i-1))*(Grid[sid].bits[bidx] ^ rseq[bidx]);
			dist += (Grid[sid].bits[bidx] != (TARGET>>i)%2);
			if (fid==3)
			{
				Grid[sid].bits[bidx] = (TARGET>>i)%2;
			}
			bidx++;
		}
		
		if (fid<2)
		{
			Grid[sid].fixed[fid] = ((float)result)/(float)(1<<NBITS);
		}
		else
		{
			Grid[sid].persist = fmaxf(0.0f,((1.0f-dist/(float)NBITS)));
			Grid[sid].W = WBIAS * (0.54 - FITNESS_THRESHOLD);//WBIAS * ( (1.0-dist/(float)NBITS) - 0.5 );
		}
	}
}

void Contract(Site *Grid, int idx)
{
	if (idx<NSITES*NVAR)
	{   
		int site = idx/NVAR;
		int var = idx%NVAR;

		// Determine target point based on neighbors
		Grid[site].targ[var] = Grid[site].fixed[ Grid[site].bits[var] ];
		float norm = 1.0f;
		
		for (int j=0;j<2;j++)
		{
			int k=2*j-1 + site;
			
			if (k<0) k += NSITES;
			if (k>=NSITES) k -= NSITES;
			
			norm += Grid[k].W;
			Grid[site].targ[var] += Grid[k].fixed[ Grid[k].bits[var] ] * Grid[k].W;			
		}

		Grid[site].targ[var] /= norm;

		Grid[site].logis[var] += CRATE*(Grid[site].targ[var] - Grid[site].logis[var]); 
		if (Grid[site].logis[var]<0.0001) Grid[site].logis[var] = 0.0001;
		if (Grid[site].logis[var]>0.9999) Grid[site].logis[var] = 0.9999;
	}
}

void Iterate()
{
/*	int block_size = BSIZE;
	int n_blocks1 = NSITES*NVAR/block_size + (NSITES*NVAR%block_size == 0 ? 0 : 1);  
	int n_blocks2 = (NSITES*3)/block_size + ((NSITES*3)%block_size == 0 ? 0 : 1);  
	
	iterateLogistic <<< n_blocks1, block_size >>> (Grid);
	getFixed <<< n_blocks2, block_size >>> (Grid);
	Contract <<< n_blocks1, block_size >>> (Grid);*/
	
	for (int i=0;i<NSITES*NVAR;i++)
	{
		iterateLogistic(Grid,i);
	}
	
	for (int i=0;i<NSITES*3;i++)
	{
		getFixed(Grid,i);
	}

	for (int i=0;i<NSITES*NVAR;i++)
	{
		Contract(Grid,i);
	}
}

void Init()
{
	Grid=(Site*)malloc(sizeof(Site)*NSITES);
	
	for (int i=0;i<NSITES;i++)
	{
		Grid[i].fixed[0] = 0.25;
		Grid[i].fixed[1] = 0.75;
		for (int j=0;j<NVAR;j++)
			Grid[i].logis[j] = (rand()%1000000)/1000000.0;			
		//Grid[i].W = (rand()%100001)/100000.0;
	}
	
	
	XRes = WIDTH; YRes = HEIGHT; Bpp = 3;
	InitSDL();
	ScreenBuf=(unsigned char*)malloc(XRes*YRes*Bpp);
}

int iter = 0;

void Render()
{
	int y = iter%HEIGHT;
	
	for (int x=0;x<NSITES;x++)
	{
		int r,g,b;
		float w = Grid[x].persist; //(1.0/(1.0-FITNESS_THRESHOLD))*(Grid[x].W/WBIAS);
				
		g=192-92*(w-0.5)*(w-0.5)*4;
		r=255*(1.0-w);
		b=255*w;
		
		if (w != w) printf("NaN!\n");
		//if (w<0) r=g=b=255;
		
		if (r<0) r=0; if (r>255) r=255;
		if (g<0) g=0; if (g>255) g=255;
		if (b<0) b=0; if (b>255) b=255;
		
		for (int y2=0;y2<CELLSIZE;y2++)
		{
			for (int x2=0;x2<CELLSIZE;x2++)
			{
				int xm = CELLSIZE*x+x2;
				int ym = CELLSIZE*y+y2;
				ScreenBuf[(xm+ym*XRes)*Bpp]=b;
				ScreenBuf[(xm+ym*XRes)*Bpp+1]=g;
				ScreenBuf[(xm+ym*XRes)*Bpp+2]=r;
			}
		}
	}
}

int main(int argc, char **argv)
{
	Img I;
	
	srand(time(NULL));
	Init();
	
	I.Width = WIDTH; I.Height = HEIGHT; I.Image=ScreenBuf;
	
	while (1)
	{
		int Ch=ReadKey();
		
		if (Ch=='q') return 0;

		Iterate();
		Render();
		BlitBuf(ScreenBuf,0,0,WIDTH,HEIGHT);
		iter++;
		
		double fbar = 0;
		
		for (int i=0;i<NSITES;i++)
		{
			fbar += Grid[i].persist;//(Grid[i].W/WBIAS+0.5);
		}
		fbar /= (double)NSITES;
		
		if (iter%100==0)
			printf("%.6g\n",fbar);
			
		//~ if (iter == HEIGHT)
		//~ {
			//~ PNMSave("tmp.pnm",I);
			//~ return 0;
		//~ }
	}
}
