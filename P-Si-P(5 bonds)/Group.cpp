#include "Group.h"



Group::Group()
{
	//std::cout << "Group constructor" << std::endl;
}


Group::~Group()
{
	//std::cout << "Group deconstructor" << std::endl;
}


/*find all bonded atoms*/
void Group::Graph()
{
	//call constructor
	Tool groupTool;
	
	// in this code, O is the first atom(need to be changed when use different POSCAR)
	orderofO = 0;

	//get each real O atom coordinate x, y, x from atomCoordinate
	for (int i = startPoint[orderofO]; i < atomQuantity[orderofO] * 3; i += 3)// i is NOT the order of O
	{	
		directCoorO[0] = atomCoordinate[i];//get x
		directCoorO[1] = atomCoordinate[i + 1];//get y
		directCoorO[2] = atomCoordinate[i + 2];//get z

		//get all possible coordinates of this O atom
		allPossibleCoordinates(directCoorO, allCoorO);


		//calcualte distance between Si and O;
		for (int j = startPoint[1]; j < (atomQuantity[0] + atomQuantity[1]) * 3; j += 3)//get Si coordinates(j is NOT equal to the order of Si)
		{
			directCoorSi[0] = atomCoordinate[j];//get x
			directCoorSi[1] = atomCoordinate[j + 1];//get y
			directCoorSi[2] = atomCoordinate[j + 2];//get z

			//convert direct to cartesian for both Si and O coordinates
			groupTool.directtoCart(directCoorSi, cartesianCoorSi);//for Si

			for (int k = 0; k < 27; ++k)//find shortest distance between Si and 27 possible O 
			{
				
				directCoorO[0] = allCoorO[k][0];//get O coordinate of x
				directCoorO[1] = allCoorO[k][1];//get O coordinate of y
				directCoorO[2] = allCoorO[k][2];//get O coordinate of z
				groupTool.directtoCart(directCoorO, cartesianCoorO);//for O

				T distance = groupTool.distanceinCart(cartesianCoorO, cartesianCoorSi);//get distance

				//determine whether they bond together
				int temp = 0;
				//1.56 need to be changed when use different POSCAR 
				if (distance < 2.00)//the bond length between Si and O is around 1.55A
				{
					//[(j - startPoint[1]) / 3] is the order of Si, [i / 3] is the order of O.
					graph[(j - startPoint[1]) / 3][i / 3] = 1;//if Si and O bonded together, set value = 1.
					temp++;
					
				}
				if(temp == 4)
					break;

			}
			

		}

		//calcualte distance between P and O;
		for (int j = startPoint[2]; j < totalAtomQuantity * 3; j += 3)//get P coordinates(j is NOT equal to the order of P)
		{
			directCoorP[0] = atomCoordinate[j];//get x
			directCoorP[1] = atomCoordinate[j + 1];//get y
			directCoorP[2] = atomCoordinate[j + 2];//get z

			//convert direct to cartesian for both P and O coordinates
			groupTool.directtoCart(directCoorP, cartesianCoorP);//for Si

			for (int k = 0; k < 27; ++k)//find shortest distance between P and 27 possible O 
			{

				directCoorO[0] = allCoorO[k][0];//get O coordinate of x
				directCoorO[1] = allCoorO[k][1];//get O coordinate of y
				directCoorO[2] = allCoorO[k][2];//get O coordinate of z
				groupTool.directtoCart(directCoorO, cartesianCoorO);//for O

				T distance = groupTool.distanceinCart(cartesianCoorO, cartesianCoorP);//get distance

				//determine whether they bond together
				//1.56 need to be changed when use different POSCAR 
				if (distance < 2.00)//the bond length between P and O is around 1.55A
				{
					//[(j - startPoint[2]) / 3] is the order of P, [i / 3] is the order of O.
					graph[((j - startPoint[2]) / 3)+atomQuantity[1]][i / 3] = 1;//if Si and O bonded together, set value = 1.
				}

			}


		}

	}
}


/*get all possible coordinates for one atom*/
void Group::allPossibleCoordinates(const T (&coor)[3], T (&allCoor)[27][3])
{
	int loop = 0;
	for (int i = -1; i <= 1; i++)
	{
		for (int j = -1; j <= 1; j++)
		{
			for (int k = -1; k <= 1; k++)
			{
				allCoor[loop][0] = coor[0] + i;//x
				allCoor[loop][1] = coor[1] + j;//y
				allCoor[loop][2] = coor[2] + k;//z
				++loop;
			}
		}
	}
}


//search bonded group and store orders in each array
void Group::groupfromGraph()
{
	//group Si and store O's order in array
	for (int i = 0; i < atomQuantity[1]; ++i)
	{
		int n = 0;
		for (int j = 0; j < atomQuantity[0]; ++j) {
			if (graph[i][j] == 1 && n < 4) {
				//store O's order
				groupSi[i][n] = j;
				++n;
			}
		}
	}
	//group P and store O's order in array
	for (int i = atomQuantity[1]; i < (atomQuantity[1]+atomQuantity[2]); ++i)
	{
		int l = 0;
		for (int j = 0; j < atomQuantity[0]; ++j) {
			if (graph[i][j] == 1 && l < 5) {
				//store O's order
				groupP[i-atomQuantity[1]][l] = j;
				++l;
			}
		}
	}

	//group O and store Si's order in array
	for (int i = 0; i < atomQuantity[0]; ++i)
	{
		int m = 0;
		for (int j = 0; j < (atomQuantity[1]+atomQuantity[2]); ++j)
		{
			if (graph[j][i] == 1 && m < 2)
			{
				groupO[i][m] = j;
				++m;
			}
		}
	}
};


//cut bonds
void Group::cutBond(int graph[][129])
{
	//choose a P atom randomly
	int center_order = generator_int(0, atomQuantity[2] - 1);//need to know the quantity of P
	int P_order = center_order + atomQuantity[1];//real P order in graph

	//find a random Si atom bounded to same O atom as P
	//get all order of O atoms bounded to the P atom
	int a = 0;
	int b = 0;
	for (int i = 0; i < atomQuantity[0]; ++i)//get all O in graph
	{
		if (graph[P_order][i] == 1 && (a+b)<5)//center_order+atomQuantity[1] is the order of P atom
		{
			for (int j = 0; j < (totalAtomQuantity - atomQuantity[0]); ++j)
			{
				if (graph[j][i] == 1 && j < atomQuantity[1])//means this O is bonded to Si
				{
					orderP_O_Si[a] = i;
					orderSifromP_O_Si[a] = j;
					++a;
				}
				else if (graph[j][i] == 1 && j >= atomQuantity[1] && j!=P_order)
				{
					orderP_O_P[b] = i;
					++b;
				}
			}
		}
	}

	int O_number = generator_int(0, 2);//choose an O atom at random from orderP_O_Si[3]
	int O_order = orderP_O_Si[O_number];//real O order in graph
	orderAllBondedOwithP[0] = O_order;
	int Si_order = orderSifromP_O_Si[O_number];//real Si order in graphi


	std::cout << "P before cut: " << P_order + 1 << "    ";
	for (int i = 0; i < atomQuantity[0]; ++i)
	{
		if (graph[P_order][i] == 1)
			std::cout << i + 1 << "    ";
	}
	std::cout << std::endl;
	std::cout << "Si before cut: " << Si_order + 1 << "    ";
	for (int i = 0; i < atomQuantity[0]; ++i)
	{
		if (graph[Si_order][i] == 1)
			std::cout << i + 1 << "    ";
	}
	std::cout << std::endl;



	//First exchange
	//exchange positions between choosen P and Si atoms
	//get start point for P and Si in atomcoordinates
	int startPointforPinCoordinates = startPoint[2] + center_order * 3;
	int startPointforSiinCoordinates = startPoint[1] + Si_order * 3;

	        std::cout << "Si begin :"<<"    ";
        for(int i=0; i<3;++i)
                std::cout << atomCoordinate[startPointforSiinCoordinates + i]<<"          ";
        std::cout << std::endl;
        std::cout << "P begin: "<<"    ";
        for(int i=0; i<3;++i)
                std::cout << atomCoordinate[startPointforPinCoordinates + i]<<"          ";

	//exchange positions
	T *temp = new T[3];
	//store Si position into temp
	for (int i = 0; i < 3; ++i)
		temp[i] = atomCoordinate[startPointforSiinCoordinates + i];
	//store P position into Si coordinates
	for (int i = 0; i < 3; ++i)
		atomCoordinate[startPointforSiinCoordinates + i] = atomCoordinate[startPointforPinCoordinates + i];
	//store temp to P coordiantes
	for (int i = 0; i < 3; ++i)
		atomCoordinate[startPointforPinCoordinates + i] = temp[i];
	delete[] temp;

	//Second exchange or cut
	//determine bonded O with Si and P
	//randomly choose an O atom except the O atom bonded with both P and Si chosen.
	int cutO = generator_int(0, 2);//choose one O atom to be cut from the two chooses bonded with Si
	while (orderP_O_Si[cutO] == O_order)
		cutO = generator_int(0, 2);
	orderAllBondedOwithP[1] = orderP_O_Si[cutO];

	//find all O atoms bonded with Si and store three orders into new P neighbors
	int p = 2;
	for (int i = 0; i < atomQuantity[0]; ++i)
	{
		if (graph[Si_order][i] == 1 && i != O_order)
		{
			orderAllBondedOwithP[p] = i;
			++p;
		}
	}
	
	for(int i=1; i<3; ++i)
		orderAllBondedOwithSi[i] =  orderP_O_P[i-1];
	orderAllBondedOwithSi[0] = O_order;
	for(int i=0; i<3; ++i)
	{
		if(orderP_O_Si[i] != O_order && orderP_O_Si[i] != orderAllBondedOwithP[1] )
			orderAllBondedOwithSi[3] = orderP_O_Si[i];
	}

	//reswt to 0 
	for(int i=0; i< atomQuantity[0]; ++i)
	{
		graph[Si_order][i] = 0;
		graph[P_order][i] = 0;
	}

	//reset bonded O with P and Si
	for(int i=0; i<4; ++i)
		graph[Si_order][orderAllBondedOwithSi[i]] = 1;
	for(int i=0;i<5;++i )
		graph[P_order][orderAllBondedOwithP[i]] = 1;
	
	std::cout << "Si final :"<<"    ";
	for(int i=0; i<3;++i)
		std::cout << atomCoordinate[startPointforSiinCoordinates + i]<<"          ";
	std::cout << std::endl;
	std::cout << "P final: "<<"    ";	
	for(int i=0; i<3;++i)
                std::cout << atomCoordinate[startPointforPinCoordinates + i]<<"          ";


	std::cout << "P after cut: " << P_order + 1 << "    ";
	for (int i = 0; i < atomQuantity[0]; ++i)
	{
		if (graph[P_order][i] == 1)
			std::cout << i + 1 << "    ";
	}
	std::cout << std::endl;
	std::cout << "Si after cut: " << Si_order + 1 << "    ";
	for (int i = 0; i < atomQuantity[0]; ++i)
	{
		if (graph[Si_order][i] == 1)
			std::cout << i + 1 << "    ";
	}
	std::cout << std::endl;

};



//generate a integer in range (begin, end)
int Group::generator_int(int begin, int end)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(begin, end);
	return dis(gen);
}



