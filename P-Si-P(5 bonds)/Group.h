#pragma once

#include <iostream>
#include <random>
#include "CONSTANT.h"
#include "Tool.h"

class Group
{
private:
	int orderofO;
	T directCoorO[3], cartesianCoorO[3];
	T directCoorSi[3], cartesianCoorSi[3];
	T directCoorP[3], cartesianCoorP[3];
	T allCoorO[27][3];
	int orderP_O_Si[3], orderSifromP_O_Si[3], orderP_O_P[2];
	int orderAllBondedOwithP[5], orderAllBondedOwithSi[4];
public:
	Group();
	~Group();

public:
	void Graph();
	void allPossibleCoordinates(const T (&coor)[3], T (&allCoor)[27][3]);//for direct coordinate
	void groupfromGraph();
	void cutBond(int graph[][129]);
	int generator_int(int begin, int end);
	
};

