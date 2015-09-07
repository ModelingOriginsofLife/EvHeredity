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

class Vector
{
	public:
		float x,y;
};

class Site
{
	public:
		Vector logis[NVAR], targ[NVAR];
		int bits[NVAR];
		float fixed[4];
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
		int site = idx/(NVAR);
		int var = idx%(NVAR);
				
		for (int i=0;i<NITER;i++)
		{
			float nx = 1 - ALPHA*pow(Grid[site].logis[var].x,2)+Grid[site].logis[var].y;
			float ny = BETA*Grid[site].logis[var].x;
			
			if (pow(nx,2)+pow(ny,2)>=3*3)
			{
				nx = (rand()%2000001-1000000.0)/1000000.0;
				ny = (rand()%2000001-1000000.0)/1000000.0;
			}
			
			Grid[site].logis[var].x = nx;//0.2*(rand()%200001-100000.0)/100000.0;
			Grid[site].logis[var].y = ny;//0.2*(rand()%200001-100000.0)/100000.0;
		}
		
		float d1 = pow(Grid[site].logis[var].x - Grid[site].fixed[0],2) + pow(Grid[site].logis[var].y - Grid[site].fixed[1],2);
		float d2 = pow(Grid[site].logis[var].x - Grid[site].fixed[2],2) + pow(Grid[site].logis[var].y - Grid[site].fixed[3],2);
		
		Grid[site].bits[var] = (d2<d1);
	}
}

void getFixed(Site *Grid, int idx)
{
	if (idx<NSITES*5)
	{   
		int fid = idx/NSITES;
		int sid = idx%NSITES;
		int result = 0;
		
		int bidx = fid*NBITS;
		float dist = 0;
				
		for (int i=0;i<NBITS;i++)
		{
			result += (1<<(NBITS-i-1))*(Grid[sid].bits[bidx] ^ ((i+fid)%7 == 0));
			dist += (Grid[sid].bits[bidx] != (TARGET>>i)%2);
			if (fid==4)
			{
				Grid[sid].bits[bidx] = (TARGET>>i)%2;
			}
			bidx++;
		}
		
		if (fid==4)
		{
			Grid[sid].persist = fmaxf(0.0f,((1.0f-dist/(float)NBITS)));
			Grid[sid].W =  WBIAS * (NEIGHBOR_W - FITNESS_THRESHOLD);//WBIAS * ( (1.0-dist/(float)NBITS) - 0.5 ); 
		}
		else
		{			
			Grid[sid].fixed[fid] = 2.0*( ((float)result)/(float)(1<<NBITS) - 0.5);
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
		Grid[site].targ[var].x = Grid[site].fixed[ 0 + 2*Grid[site].bits[var] ];
		Grid[site].targ[var].y = Grid[site].fixed[ 1 + 2*Grid[site].bits[var] ];
		float norm = 1.0f;
		
		for (int j=0;j<2;j++)
		{
			int k=2*j-1 + site;
			
			if (k<0) k += NSITES;
			if (k>=NSITES) k -= NSITES;
			
			norm += Grid[k].W;
			Grid[site].targ[var].x += Grid[k].fixed[ 0 + 2*Grid[k].bits[var] ] * Grid[k].W;
			Grid[site].targ[var].y += Grid[k].fixed[ 1 + 2*Grid[k].bits[var] ] * Grid[k].W;
		}

		Grid[site].targ[var].x /= norm;	Grid[site].targ[var].y /= norm;

		Grid[site].logis[var].x += CRATE*(Grid[site].targ[var].x - Grid[site].logis[var].x); 
		//~ if (Grid[site].logis[var].x<-2) Grid[site].logis[var].x = -2;
		//~ if (Grid[site].logis[var].x>2) Grid[site].logis[var].x = 2;

		Grid[site].logis[var].y += CRATE*(Grid[site].targ[var].y - Grid[site].logis[var].y); 
		//~ if (Grid[site].logis[var].y<-2) Grid[site].logis[var].y = -2;
		//~ if (Grid[site].logis[var].y>2) Grid[site].logis[var].y = 2;
		
		int oob = 0;
		
		if (pow(Grid[site].logis[var].x,2)+pow(Grid[site].logis[var].y,2) >= 2*2)
		{
			Grid[site].logis[var].x = (rand()%2000001-1000000.0)/1000000.0;
			Grid[site].logis[var].y = (rand()%2000001-1000000.0)/1000000.0;
		}
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
	
	for (int i=0;i<NSITES*5;i++)
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
		Grid[i].fixed[0] = 0.63-0.15;
		Grid[i].fixed[1] = 0.189-0.1;
		Grid[i].fixed[2] = 0.63+0.15;
		Grid[i].fixed[3] = 0.189+0.1;
		for (int j=0;j<NVAR;j++)
		{
			Grid[i].logis[j].x = 2*( (rand()%1000000)/1000000.0 - 0.5);
			Grid[i].logis[j].y = 2*( (rand()%1000000)/1000000.0 - 0.5);
		}
		//Grid[i].W = (rand()%100001)/100000.0;
	}
	
	
	XRes = WIDTH; YRes = HEIGHT; Bpp = 3;
	InitSDL();
	ScreenBuf=(unsigned char*)malloc(XRes*YRes*Bpp);
}

int iter = 0;

void RenderStats()
{
	int y = iter%(HEIGHT/CELLSIZE);
	
	for (int x=0;x<NVAR;x++)
	{
		int r,g,b;
		double w = 0;
		
		for (int x2=0;x2<NSITES;x2++)
		{
			w += (Grid[x2].bits[x] == 1)/(double)NSITES;
		}
		
		w = 4*w*(1-w);
		
		g=192-92*(w-0.5)*(w-0.5)*4;
		r=255*(1.0-w);
		b=255*w;

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

void RenderState()
{
	int y = iter%(HEIGHT/CELLSIZE);
	
	for (int x=0;x<NSITES;x++)
	{
		int r,g,b;
		float w = Grid[x].persist; //(1.0/(1.0-FITNESS_THRESHOLD))*(Grid[x].W/WBIAS);
		//float w = (Grid[x].logis[0].x+1.0)/2.0;
		
		w = (Grid[x].fixed[0]+1)/2.0;
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
				int xm = CELLSIZE*x+x2 + WIDTH1;
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
	
	int t=0;
	
	while (1)
	{
		int Ch=ReadKey();
		
		t++;
		//CRATE = 0.25*(-cos(2*M_PI*t/1000.0)+1);
		//~ CRATE += 1e-3;
		//~ if (CRATE > 1.0) CRATE = 1.0;
		printf("%.6g %.6g %.6g %.6g: %d %d %d %d\n",Grid[0].fixed[0], Grid[0].fixed[1],Grid[0].fixed[2], Grid[0].fixed[3], Grid[0].bits[0],Grid[0].bits[1],Grid[0].bits[2],Grid[0].bits[3]);
		if (Ch=='q') return 0;

		Iterate();
		RenderStats();
		RenderState();
		BlitBuf(ScreenBuf,0,0,WIDTH,HEIGHT);
		iter++;
		
		double fbar = 0;
		
		for (int i=0;i<NSITES;i++)
		{
			fbar += Grid[i].persist;//(Grid[i].W/WBIAS+0.5);
		}
		fbar /= (double)NSITES;
		
		//~ if (iter%100==0)
			//~ printf("%.6g\n",fbar);
			
		//~ if (iter == HEIGHT)
		//~ {
			//~ PNMSave("tmp.pnm",I);
			//~ return 0;
		//~ }
	}
}
