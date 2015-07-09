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

#include "fblib.h"
#include "imgload.h"

#define BASES 200
#define DIFF 0.1

double *Interaction;
double BETA, FILTER;

int XR=256, YR=256;

/* Very simple model:
 * - Bases form dimers in solution according to their interaction energy Interaction[i,j]
 * - Periodic washing removes all spare monomers, replaces them with a uniform distribution
 * - Dimers are broken down to monomers, fixed count are sampled for the next cycle
 * 
 * In bonding phase: 
 * P(a=b) ~ Na * Nb * exp(-I(a,b)/T)
 * 
 * Nx(t+1) = (1/B)*sum(Ny(t)) + sum(Nyz(t) * ( delta(y=x) + delta(z=x) ))
 */

/* Make a random matrix of interaction energies */
void generateInteractions()
{
	Interaction = (double*)malloc(BASES*BASES*sizeof(double));
	
	for (int j=0;j<BASES;j++)
	{
		for (int i=0;i<BASES;i++)
		{
			if (j<i)
				Interaction[i+j*BASES] = BETA * (rand()%2000001-1000000.0)/1000000.0;
			else
				Interaction[i+j*BASES] = Interaction[j+i*BASES];
			if (i==j) Interaction[i+j*BASES] = 1.0;
		}
	}
}

double *conc, *conc2, *cPairs, *cPairs2;
double red[BASES], green[BASES], blue[BASES];

void Init()
{
	conc = (double*)malloc(sizeof(double)*BASES*XR*YR);
	conc2 = (double*)malloc(sizeof(double)*BASES*XR*YR);
	cPairs = (double*)malloc(sizeof(double)*BASES*BASES);
	cPairs2 = (double*)malloc(sizeof(double)*BASES*BASES);
	
	for (int i=0;i<BASES*XR*YR;i++)
		conc[i] = 1.0/(double)BASES;
		
	/* To visualize this, we're just going to assign each compound a unique RGB color and then calculate the average color. 
	 * This isn't guaranteed to be unique, but it tends to be sufficient to visually distinguish contiguous regions with
	 * the same composition */
	for (int i=0;i<BASES;i++)
	{
		red[i] = (rand()%1000001)/1000000.0;
		green[i] = (rand()%1000001)/1000000.0;
		blue[i] = (rand()%1000001)/1000000.0;
	}
}

void Iterate()
{
	double norm = 0;
	memset(cPairs,0,BASES*BASES*sizeof(double));
	memset(cPairs2,0,BASES*BASES*sizeof(double));
	
	for (int y=0;y<YR;y++)
	for (int x=0;x<XR;x++)
	{
		norm = 0;
		memset(cPairs,0,BASES*BASES*sizeof(double));
		memset(cPairs2,0,BASES*BASES*sizeof(double));
		int coord = (x+y*XR)*BASES;
		
		/* Calculate the fraction of dimers of type i,j in equilibrium */
		for (int i=0;i<BASES;i++)
		{
			for (int j=0;j<=i;j++)
			{
				cPairs[i+j*BASES] = conc[i+coord]*conc[j+coord]*exp(-Interaction[i+j*BASES]);
				norm += cPairs[i+j*BASES];
			}
		}

		for (int i=0;i<BASES;i++)
		{
			for (int j=0;j<=i;j++)
			{
				cPairs[i+j*BASES] *= (0.5 - FILTER/2.0)/norm; // Normalize so that a fixed percentage of the number of molecules are inside of dimers
			}
		}
		
		/* Wash the supernatant with a random mixture of 10 different bases */
		int baseList[BASES];
		memset(baseList,0,sizeof(int)*BASES);
		
		for (int i=0;i<10;i++)
			baseList[rand()%BASES]++;
		
		for (int i=0;i<BASES;i++)
		{
			conc[i+coord] = baseList[i] * FILTER/(double)10.0;
		}
		
		/* Add the contents of the dimers back in as monomers */
		for (int i=0;i<BASES;i++)
		{
			for (int j=0;j<=i;j++)
			{
				conc[i+coord] += cPairs[i+j*BASES];
				conc[j+coord] += cPairs[i+j*BASES];
			}
		}
	}
	
	/* Diffusion */
	for (int y=0;y<YR;y++)
	for (int x=0;x<XR;x++)
	{
		int xm,ym;
		int coord = (x+y*XR)*BASES;
		
		for (int i=0;i<BASES;i++)
		{
			conc2[i+coord] = conc[i+coord];
			
			xm = x-1; ym = y; if (xm<0) xm+=XR;
			conc2[i+coord] += DIFF*( conc[i+(xm+ym*XR)*BASES] - conc[i+coord] );
			xm = x+1; ym = y; if (xm>=XR) xm-=XR;
			conc2[i+coord] += DIFF*( conc[i+(xm+ym*XR)*BASES] - conc[i+coord] );
			xm = x; ym = y-1; if (ym<0) ym+=YR;
			conc2[i+coord] += DIFF*( conc[i+(xm+ym*XR)*BASES] - conc[i+coord] );
			xm = x; ym = y+1; if (ym>=YR) ym-=YR;
			conc2[i+coord] += DIFF*( conc[i+(xm+ym*XR)*BASES] - conc[i+coord] );
		}
	}
	
	memcpy(conc,conc2,XR*YR*BASES*sizeof(double));	
}

void Render()
{
	for (int y=0;y<YR;y++)
	for (int x=0;x<XR;x++)
	{
		double rf=0,gf=0,bf=0;
		int r,g,b;
		
		for (int i=0;i<BASES;i++)
		{
			rf += conc[i+(x+y*XR)*BASES]*red[i];
			gf += conc[i+(x+y*XR)*BASES]*green[i];
			bf += conc[i+(x+y*XR)*BASES]*blue[i];
		}
		
		r=255*rf; g=255*gf; b=255*bf;
		if (r<0) r=0; if (r>255) r=255;
		if (g<0) g=0; if (g>255) g=255;
		if (b<0) b=0; if (b>255) b=255;
		
		ScreenBuf[(x+y*XR)*Bpp]=b;
		ScreenBuf[(x+y*XR)*Bpp+1]=g;
		ScreenBuf[(x+y*XR)*Bpp+2]=r;
	}
}

int main(int argc, char **argv)
{
	char Str[512];
	Img I;
	
	Bpp=3;
	ScreenBuf=(unsigned char*)malloc(XR*YR*Bpp);
	I.Image = ScreenBuf; I.Width=XR; I.Height=YR;
	
	srand(12345);
	Init();
	
	int t=0;
	FILE *f;
	
	BETA = 5.0;
	FILTER = 0.2;
	
	for (int i=0;i<BASES*XR*YR;i++)
		conc[i] = 0;
	
	/* Because the total concentration of non-water is treated as conserved, we have to initialize this with something.
	 * In this case, we pick a random base for each grid site and set its concentration to 100% at that location */
	for (int x=0;x<XR*YR;x++)
		conc[(rand()%BASES)+x*BASES] = 1.0;
	
	generateInteractions();
		
	t=0;
	while (t<4000)
	{		
		Iterate();
		
		t++;
		
		printf("%d\n",t);
		Render();
		PNMSave("tmp.pnm",I);
		sprintf(Str,"pnmtopng -force tmp.pnm > frames/%.6d.png",t);
		system(Str);
	}
}
