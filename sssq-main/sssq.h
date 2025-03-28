#ifndef _sssq
#define _sssq
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>
#include <map>
#include <fstream>
#include "more/ore.h"
using namespace std;

#define BORDER 12
struct bentry {
	double maxX;
	double maxY;
	double minX;
	double minY;
};

class sssq {
public:
	vector<pair<double, double> > dataset;
	vector<bentry> border;
	//vector<vector<int> > regions; 
	vector<pair<double, double> > vertex;
	//map<pair<int,int>,vector<int> > edges;
	//vector<vector<int> > adjv;	//adjacent vertex	
public:
	sssq() {};
	int Process(char*,char*);
};
#endif
