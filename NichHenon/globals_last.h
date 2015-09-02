// Evolving population

#define TARGET 1+2+4+8+16+32
#define NBITS 12
#define NITER 4
#define NCONTRACT 1
#define LOGISTIC_R 3.7f
#define CRATE 0.8f
#define WBIAS 2.5f

#define FITNESS_THRESHOLD 0.65

#define NSITES 1024
#define NVAR (3*NBITS)

#define CELLSIZE 2
#define WIDTH (CELLSIZE*NSITES)
#define HEIGHT (CELLSIZE*NSITES/4)

#define BSIZE 128

int rseq[NVAR] =
{
	0,0,1,
	1,0,1,
	1,1,0,
	0,0,0,
	1,0,0,
	1,0,1,
	0,0,0,
	1,1,0,
	0,0,1,
	1,0,0,
	0,1,0,
	1,0,1
};
	
