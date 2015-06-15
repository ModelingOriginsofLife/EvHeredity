#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <deque>
#include <vector>

using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <unordered_map>

#define KSCALE 1e-1
#define SPLITRATIO 0.25
#define TARGETSIZE 8
#define EQUILIBRIUM 0

double rateConstants[16];

unordered_map<string,int> lastFrame, lastFrame2;

char reverseCase(char c)
{
	switch(c)
	{
		case 'A': return 'a';
		case 'B': return 'b';
		case 'C': return 'C'; // These are their own inverse
		case 'D': return 'D'; // These are their own inverse
		case 'E': return 'E'; // These are their own inverse
		case 'F': return 'F'; // These are their own inverse
		case 'G': return 'g';
		case 'H': return 'h';
		case 'I': return 'i';
		case 'J': return 'j';
		case 'a': return 'A';
		case 'b': return 'B';
		case 'c': return 'C'; // These are their own inverse
		case 'd': return 'D'; // These are their own inverse
		case 'e': return 'E'; // These are their own inverse
		case 'f': return 'F'; // These are their own inverse
		case 'g': return 'G';
		case 'h': return 'H';
		case 'i': return 'I';
		case 'j': return 'J';
	}
}

double getPairRate(char c1, char c2)
{
	int leftidx, rightidx;
	// A: 0=1
	// B: 2=3
	// a: 1=0
	// b: 3=2
	// C,c: 0=0
	// D,d: 1=1
	// E,e: 2=2
	// F,f: 3=3
	// G: 0=2
	// H: 0=3
	// I: 1=2
	// J: 1=3
	// g: 2=0
	// h: 3=0
	// i: 2=1
	// j: 3=1
	
	switch(c1)
	{
		case 'A': leftidx=1; break;
		case 'B': leftidx=3; break;
		case 'a': leftidx=0; break;
		case 'b': leftidx=2; break;
		case 'C':
		case 'c': leftidx=0; break;
		case 'D':
		case 'd': leftidx=1; break;
		case 'E':
		case 'e': leftidx=2; break;
		case 'F':
		case 'f': leftidx=3; break;
		
		case 'G': leftidx=2; break;
		case 'H': leftidx=3; break;
		case 'I': leftidx=2; break;
		case 'J': leftidx=3; break;
		case 'g': leftidx=0; break;
		case 'h': leftidx=0; break;
		case 'i': leftidx=1; break;
		case 'j': leftidx=1; break;
	}

	switch(c2)
	{
		case 'A': rightidx=0; break;
		case 'B': rightidx=2; break;
		case 'a': rightidx=1; break;
		case 'b': rightidx=3; break;
		case 'C':
		case 'c': rightidx=0; break;
		case 'D':
		case 'd': rightidx=1; break;
		case 'E':
		case 'e': rightidx=2; break;
		case 'F':
		case 'f': rightidx=3; break;
		case 'g': rightidx=2; break;
		case 'h': rightidx=3; break;
		case 'i': rightidx=2; break;
		case 'j': rightidx=3; break;
		case 'G': rightidx=0; break;
		case 'H': rightidx=0; break;
		case 'I': rightidx=1; break;
		case 'J': rightidx=1; break;
	}
	
	return rateConstants[leftidx + rightidx*4];
}

class Polymer
{
	public:
		deque<char> seq;
		
		bool attemptMerge(Polymer &reactant);
		bool attemptSplit(Polymer *P1, Polymer *P2);
		deque<char> Canonicalize();
		void reverse();
};

void Polymer::reverse()
{
	deque<char> newSeq;
	
	for (int i=seq.size()-1;i>=0;i--)
	{
		newSeq.push_back(reverseCase(seq[i]));
	}
	
	seq.clear();
	seq = newSeq;
}

bool Polymer::attemptMerge(Polymer &reactant)
{
	double k = getPairRate(seq[seq.size()-1], reactant.seq[0]);
	
	if (rand()%10000000<10000000.0*k*KSCALE*SPLITRATIO)
	{
		for (int i=0;i<reactant.seq.size();i++)
			seq.push_back(reactant.seq[i]);
		
		return true;
	}
	
	return false;
}

bool Polymer::attemptSplit(Polymer *P1, Polymer *P2)
{
	double ktotal = 0;
	
	if (seq.size()<=1) return false;
	
	for (int i=1;i<seq.size();i++)
	{
		double k = getPairRate(seq[i-1],seq[i]);
		
		ktotal += k * KSCALE * SPLITRATIO;
	}
	
	double r = (rand()%10000000)/10000000.0;
	
	if (r>ktotal) return false;
	
	int i=1;
	do
	{
		double k = getPairRate(seq[i-1],seq[i]);
		
		r -= k * KSCALE * SPLITRATIO;
		i++;
	} while ((r>=0)&&(i<seq.size()));
	
	i--;
	
	P1->seq.clear(); P2->seq.clear();
	
	for (int j=0;j<i;j++)
	{
		P1->seq.push_back(seq[j]);
	}
	for (int j=i;j<seq.size();j++)
	{
		P2->seq.push_back(seq[j]);
	}
	
	return true;
}

void setRates()
{
	for (int i=0;i<16;i++)
		rateConstants[i] = pow(2,-8*(rand()%1000001)/1000000.0);
}

vector<Polymer> Bath;
vector<Polymer> Sample;

void Crack()
{
	// Attempt Split
	int i = rand()%Bath.size();
	Polymer N1,N2;
	if (Bath[i].attemptSplit(&N1,&N2))
	{
		if ((N1.seq.size()>0)&&(N2.seq.size()>0))
		{
			Bath[i].seq.clear();
			Bath[i] = N1;
			Bath.push_back(N2);
		}
	}
}

void Iterate(int target, int crackAttempts)
{
	if (!Bath.size()) return;
	
	int i=rand()%Bath.size();
	int j=rand()%Bath.size();
	
	// Attempt Join
	if (i!=j)
	{
		if (Bath[i].attemptMerge(Bath[j]))
		{
			Bath[j].seq.clear();
			Bath.erase(Bath.begin()+j);
		}
	}
	
	// Attempt Split
	for (int k=0;k<crackAttempts;k++)
		Crack();

	// Flip
	i = rand()%Bath.size();
	Bath[i].reverse();
	
	// Extract
	if (!EQUILIBRIUM)
	{
		i = rand()%Bath.size();
		if (Bath[i].seq.size() == target)
		{
			Sample.push_back(Bath[i]);
			Bath.erase(Bath.begin()+i);
		}
	}
}

void initBath()
{
	for (int i=0;i<600;i++)
	{
		Polymer P;
		char letter = reverseCase('A'+rand()%2); // Initial reverse case to map away self-inverses
		if (rand()%2 == 0) letter = reverseCase(letter);
		
		P.seq.push_back(letter);		
		Bath.push_back(P);
	}
}

void addMass()
{
	int mass=0;
	for (int i=0;i<Bath.size();i++)
		mass += Bath[i].seq.size();
	
	for (int i=0;i<600-mass;i++)
	{
		Polymer P;
		char letter = reverseCase('A'+rand()%2); // Initial reverse case to map away self-inverses
		if (rand()%2 == 0) letter = reverseCase(letter);
		
		P.seq.push_back(letter);		
		Bath.push_back(P);
	}
}

deque<char> Polymer::Canonicalize()
{
	deque<char> newSeq;
	
	for (int i=seq.size()-1;i>=0;i--)
	{
		int j=seq.size()-i-1;
		
		newSeq.push_back(reverseCase(seq[i]));
	}
	
	for (int i=0;i<seq.size();i++)
	{
		if (seq[i] > newSeq[i])
		{
			newSeq.clear();
			return seq;
		}
		else if (seq[i] < newSeq[i])
			return newSeq;
	}
	
	newSeq.clear();
	return seq;
}

int main(int argc, char **argv)
{
	srand(12345);
	setRates();
	initBath();
	
	int t=0;
	int recur = 0, recur2 = 0;
	
	while (1)
	{
		if (EQUILIBRIUM)
		{
			for (int i=0;i<2000;i++)
				Iterate(0,1);
			t++;
		
			if (t%40==0)
			{
				FILE *f=fopen("equilibrium.txt","a");
				for (int i=0;i<Bath.size();i++)
				{
					if (Bath[i].seq.size() == TARGETSIZE)
					{
						deque<char> canon = Bath[i].Canonicalize();
						for (int j=0;j<canon.size();j++)
						{
							fprintf(f,"%c",canon[j]);
						}
						canon.clear();
						fprintf(f,"\n");					
					}
				}
				fclose(f);
			}
		}
		else
		{
			while (Bath.size()>1)
			{
				Iterate(TARGETSIZE,1);
			}
			
			Bath.clear();
			Bath = Sample;
			Sample.clear();
			
			char Str[512];
			int oldsize = Bath.size();
			sprintf(Str,"data/%.6d.txt",t);
			FILE *f=fopen(Str,"wb");
			for (int i=0;i<Bath.size();i++)
			{
				deque<char> canon = Bath[i].Canonicalize();
				string s = "";
				for (int j=0;j<canon.size();j++)
				{
					fprintf(f,"%c",canon[j]);
					s = s + canon[j];
				}
				canon.clear();
				fprintf(f,"\n");
				
				if (lastFrame[s] == 1)
				{
					recur++;
					lastFrame[s] = 2;
				}
			}
			fclose(f);
			
			lastFrame.clear();
			for (int i=0;i<Bath.size();i++)
			{
				deque<char> canon = Bath[i].Canonicalize();
				string s = "";
				for (int j=0;j<canon.size();j++)
				{
					s = s + canon[j];
				}
				lastFrame[s] = 1;
			}
						
			while (Bath.size()>1)
			{
				Iterate(4,2);
			}
				
			Bath.clear();
			Bath = Sample;
			Sample.clear();

			sprintf(Str,"data2/%.6d.txt",t);
			f=fopen(Str,"wb");
			for (int i=0;i<Bath.size();i++)
			{
				deque<char> canon = Bath[i].Canonicalize();
				string s = "";
				for (int j=0;j<canon.size();j++)
				{
					fprintf(f,"%c",canon[j]);
					s = s + canon[j];
				}
				canon.clear();
				fprintf(f,"\n");
				
				if (lastFrame2[s] == 1)
				{
					recur2++;
					lastFrame2[s] = 2;
				}
			}
			fclose(f);
			
			lastFrame2.clear();
			for (int i=0;i<Bath.size();i++)
			{
				deque<char> canon = Bath[i].Canonicalize();
				string s = "";
				for (int j=0;j<canon.size();j++)
				{
					s = s + canon[j];
				}
				lastFrame2[s] = 1;
			}
			
			if (t%20==0)
			{
				f=fopen("log.txt","a");
				fprintf(f,"%d %d %d %d %d\n",t,recur,recur2, oldsize*20, Bath.size()*20);
				fclose(f);
				recur = 0;
				recur2 = 0;
			}
			// Add mass back in for lost monomers
			addMass();
			
			t++;
		}
	}
}
