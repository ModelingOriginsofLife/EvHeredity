// Evolving population

#define TARGET (1+2+4+8+16+32)*(1+2)
#define NBITS 24
#define NITER 6
#define NCONTRACT 1
#define LOGISTIC_R 3.7f
double CRATE = 0.75f;
#define WBIAS 2.5f
#define NEIGHBOR_W 0.54

#define ALPHA 1.4
#define BETA 0.3

#define FITNESS_THRESHOLD 0.5

#define NSITES 512
#define NVAR (5*NBITS)

#define CELLSIZE 4
#define WIDTH1 (CELLSIZE*NVAR)
#define WIDTH2 (CELLSIZE*NSITES)
#define WIDTH (WIDTH1+WIDTH2)
#define HEIGHT (CELLSIZE*NVAR*2)

#define BSIZE 128

int rseq[3*36] =
{
	0,1,0,
	1,1,1,
	0,1,1,
	0,1,1,
	1,0,1,
	1,1,1,
	0,1,1,
	1,0,1,
	1,0,1,
	0,0,0,
	0,1,0,
	0,0,1,
	1,1,1,
	1,0,1,
	1,1,0,
	1,0,0,
	1,0,0,
	1,0,1,
	0,0,1,
	1,1,0,
	0,0,1,
	1,0,0,
	0,1,0,
	1,1,1,
	0,0,1,
	1,0,1,
	1,1,0,
	1,0,0,
	1,1,0,
	1,0,1,
	1,1,0,
	0,1,0,
	0,0,1,
	1,0,0,
	0,1,1,
	1,0,1


};
	
