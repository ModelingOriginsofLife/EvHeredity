#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <unordered_map>

unordered_map<string,int> cHash;

class Entry
{
	public:
		string name;
		int value;
};

bool compareEntry(Entry A, Entry B)
{
	return (A.value > B.value);
}

int main(int argc, char **argv)
{
	unordered_map<string,int>::iterator it;
	
	FILE *f=fopen(argv[1],"rb");
	char Str[512];
	
	while (fscanf(f,"%s\n",Str)!=EOF)
	{
		cHash[Str] += 1;
	}
	
	fclose(f);
	
	vector<Entry> cList;
	
	for (it=cHash.begin();it!=cHash.end();++it)
	{
		Entry E;
		E.name = it->first;
		E.value = it->second;
		
		cList.push_back(E);
	}
	
	sort(cList.begin(), cList.end(), compareEntry);
	
	f=fopen(argv[2],"wb");
	
	for (int i=0;i<cList.size();i++)
	{
		fprintf(f,"%s %d\n",cList[i].name.c_str(), cList[i].value);
	}
	fclose(f);
}
