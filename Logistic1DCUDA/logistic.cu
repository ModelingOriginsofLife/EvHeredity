#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <cuda_runtime_api.h>
#include <curand.h>
#include <curand_kernel.h>

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
		float W;
};

__device__ float logistic(float x)
{
	return fmaxf(0.0001f,fminf(0.9999f,LOGISTIC_R * x * (1.0f - x)));
	//~ return x;
}

Site *hostGrid, *devGrid;

__global__ void iterateLogistic(Site *devGrid)
{
	long int idx = blockIdx.x*blockDim.x + threadIdx.x;
   
	if (idx<NSITES*NVAR)
	{   
		int site = idx/NVAR;
		int var = idx%NVAR;
		
		for (int i=0;i<NITER;i++)
		{
			devGrid[site].logis[var] = logistic( devGrid[site].logis[var] );
		}
		
		float d1 = abs(devGrid[site].logis[var] - devGrid[site].fixed[0]);
		float d2 = abs(devGrid[site].logis[var] - devGrid[site].fixed[1]);
		
		devGrid[site].bits[var] = (d2<d1);
	}
}

__global__ void getFixed(Site *devGrid)
{
	long int idx = blockIdx.x*blockDim.x + threadIdx.x;
	
	int rseq[NVAR] =
	{
		1,0,1,
		1,0,0,
		1,1,0,
		1,0,0,
		1,0,0,
		1,0,0,
		1,0,0,
		1,0,0,
		1,0,0,
		1,0,0,
		1,0,0,
		1,0,0
	};
	
	if (idx<NSITES*3)
	{   
		int fid = idx/NSITES;
		int sid = idx%NSITES;
		int result = 0;
		
		int bidx = fid*NBITS;
		float dist = 0;
				
		for (int i=0;i<NBITS;i++)
		{
			result += (1<<i)*(devGrid[sid].bits[bidx] ^ rseq[bidx]);
			dist += (devGrid[sid].bits[bidx] != (TARGET>>i)%2);
			bidx++;
		}
		
		if (fid<2)
		{
			devGrid[sid].fixed[fid] = ((float)result)/(float)(1<<NBITS);
		}
		else
		{
			devGrid[sid].W = WBIAS * fmaxf(0.0f,((1.0f-dist/(float)NBITS) - 0.5f));//WBIAS * ( (1.0-dist/(float)NBITS) - 0.5 );
		}
	}
}

__global__ void Contract(Site *devGrid)
{
	long int idx = blockIdx.x*blockDim.x + threadIdx.x;
   
	if (idx<NSITES*NVAR)
	{   
		int site = idx/NVAR;
		int var = idx%NVAR;

		// Determine target point based on neighbors
		devGrid[site].targ[var] = devGrid[site].fixed[ devGrid[site].bits[var] ];
		float norm = 1.0f;
		
		for (int j=0;j<2;j++)
		{
			int k=2*j-1 + site;
			
			if (k<0) k += NSITES;
			if (k>=NSITES) k -= NSITES;
			
			norm += devGrid[k].W;
			devGrid[site].targ[var] += devGrid[k].fixed[ devGrid[k].bits[var] ] * devGrid[k].W;			
		}

		devGrid[site].targ[var] /= norm;

		//~ if (devGrid[site].W > devGrid[(site+1)%NSITES].W)
			//~ devGrid[site].targ[var] = devGrid[site].fixed[ devGrid[site].bits[var] ];
		//~ else
			//~ devGrid[site].targ[var] = devGrid[(site+1)%NSITES].fixed[ devGrid[(site+1)%NSITES].bits[var] ];
		devGrid[site].logis[var] += CRATE*(devGrid[site].targ[var] - devGrid[site].logis[var]);
	}
}

void Iterate()
{
	int block_size = BSIZE;
	int n_blocks1 = NSITES*NVAR/block_size + (NSITES*NVAR%block_size == 0 ? 0 : 1);  
	int n_blocks2 = (NSITES*3)/block_size + ((NSITES*3)%block_size == 0 ? 0 : 1);  
	
	iterateLogistic <<< n_blocks1, block_size >>> (devGrid);
	getFixed <<< n_blocks2, block_size >>> (devGrid);
	Contract <<< n_blocks1, block_size >>> (devGrid);
	
	cudaMemcpy(hostGrid, devGrid, sizeof(Site)*NSITES, cudaMemcpyDeviceToHost);
}

void Init()
{
	hostGrid=(Site*)malloc(sizeof(Site)*NSITES);
	cudaMalloc((void**)&devGrid, NSITES*sizeof(Site));
	
	for (int i=0;i<NSITES;i++)
	{
		hostGrid[i].fixed[0] = 0.25;
		hostGrid[i].fixed[1] = 0.75;
		for (int j=0;j<NVAR;j++)
			hostGrid[i].logis[j] = (rand()%1000000)/1000000.0;			
		//hostGrid[i].W = (rand()%100001)/100000.0;
	}
	
	cudaMemcpy(devGrid, hostGrid, sizeof(Site)*NSITES, cudaMemcpyHostToDevice);
	
	XRes = WIDTH; YRes = HEIGHT; Bpp = 3;
	InitSDL();
	ScreenBuf=(unsigned char*)malloc(XRes*YRes*Bpp);
	memset(ScreenBuf,0,XRes*YRes*Bpp);
}

int iter = 0;

void Render()
{
	int y = iter%NSITES;
	
	for (int x=0;x<NSITES;x++)
	{
		int r,g,b;
		float w = (hostGrid[x].W/WBIAS+0.5);
		
		g=192-92*w;
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
	Init();
	
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
			fbar += (hostGrid[i].W/WBIAS+0.5);
		}
		fbar /= (double)NSITES;
		
		if (iter%100==0)
			printf("%.6g\n",fbar);
	}
}
