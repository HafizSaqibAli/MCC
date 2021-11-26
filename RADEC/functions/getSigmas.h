/* LJA/LJB to Sigma Code
**
** This function calculates sigmas from Lennard-Jones A and B coefficients. 
**
*/

#ifndef GET_SIGMAS
#define GET_SIGMAS

#include <math.h>

vector<double> getSigmas(string topFile, vector<atom> atoms){

	//Basic IO variables
	double temp;
	int i = 0, j = 0;
	bool flag = true;
	ifstream TopologyIn;
	streampos ACoeffStart = findACoeffStart(topFile);
	streampos BCoeffStart = findBCoeffStart(topFile);

	vector<double> ACoeffs;
	vector<double> BCoeffs;

	vector<int> NBPIndex = getNonBondedParmIndex(topFile);

	int noTypes = 0;
	for(int i = 0; i < atoms.size(); i++){
		if(atoms[i].getType() > noTypes){
			noTypes = atoms[i].getType();
		}
	}

	//Open solute topology file
	TopologyIn.open(topFile.c_str());

	//Find and read LJ A Coefficients
	TopologyIn.seekg(ACoeffStart);
	while(TopologyIn){
		string str;
		getline(TopologyIn, str);
		if(str[0] == '%'){
			break;
		}
		stringstream ss(str);
		while(ss >> temp){
			ACoeffs.push_back(temp);
		}
	}

	//Find and read LJ B Coefficients
	TopologyIn.seekg(BCoeffStart);
	while(TopologyIn){
		string str;
		getline(TopologyIn, str);
		if(str[0] == '%'){
			break;
		}
		stringstream ss(str);
		while(ss >> temp){
			BCoeffs.push_back(temp);
		}
	}

	vector<double> sigmas;
	for(i = 0; i < noTypes; i++){
		for(j = 0; j < noTypes; j++){
			sigmas.push_back(0);
		}
	}

	for(i = 0; i < noTypes; i++){
		for(j = 0; j < noTypes; j++){
			int index = -1;
			if(i <= j){
				index = NBPIndex[((i)*noTypes)+j];
			}
			if(i > j){
				index = NBPIndex[((j)*noTypes)+i];
			}

			if(index > 0 && index <= sigmas.size()){
				if(ACoeffs[index-1] == 0 || BCoeffs[index-1] == 0){
					sigmas[((i)*noTypes)+j] = 0;
				}
				else{
					sigmas[((i)*noTypes)+j] = (pow((ACoeffs[index-1] / BCoeffs[index-1]), 1.0/6.0));
				}
			}
		}
	}

	//Close input file
	TopologyIn.close();

	return sigmas;
}

#endif
