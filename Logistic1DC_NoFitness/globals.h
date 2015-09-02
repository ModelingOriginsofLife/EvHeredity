// Evolving population

#define TARGET (1+2+4+8+16+32+64+128)
#define NBITS 16
#define NITER 6
#define NCONTRACT 1
#define LOGISTIC_R 3.7f
#define CRATE 0.9f
#define WBIAS 2.5f

#define FITNESS_THRESHOLD 0.5

#define NSITES 2048
#define NVAR (3*NBITS)

#define CELLSIZE 1
#define WIDTH (CELLSIZE*NSITES)
#define HEIGHT (CELLSIZE*NSITES/2)

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
	
