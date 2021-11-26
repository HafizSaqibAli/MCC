/* Molecule functions, Dr Jonathan Higham
** Created 13/02/2014
** Modified 25/02/2014
**
** This collection of functions acts on the Force-Torque Covariance 
** matricies of the molecule class.
**
*/

#ifndef FTCM_OPERATIONS
#define FTCM_OPERATIONS

void updateFTCMatrices(vector<molecule> &molecules);
void solveFTCMatrices(vector<molecule> &molecules);
void divideFTCMatricies(vector<molecule> &molecules, double d);
void averageFTCMatricies(vector<molecule> molecules, vector<molecule> &moleculeTypes);

void updateFTCMatrices(vector<molecule> &molecules){

	for (int i = 0; i < molecules.size(); i++){
		molecules[i].updateFTCMatrix();
	}
}

void solveFTCMatrices(vector<molecule> &molecules){
	for (int i = 0; i < molecules.size(); i++){
		molecules[i].solveFTCMatrix();
	}
}

void divideFTCMatricies(vector<molecule> &molecules, double d){

	for (int i = 0; i < molecules.size(); i++){
		molecules[i].divideFTCMatrix(d);
	}
}

void averageFTCMatricies(vector<molecule> molecules, vector<molecule> &moleculeTypes){

	for (int t = 0; t < moleculeTypes.size(); t++){
		moleculeTypes[t].createFTCMatrix();
	}

	vector<double> occurence;
	for (int t = 0; t < moleculeTypes.size(); t++){
		occurence.push_back(0);
	}

	for (int i = 0; i < molecules.size(); i++){
		for (int t = 0; t < moleculeTypes.size(); t++){
			if (molecules[i] == moleculeTypes[t]){
				occurence[t]++;
				moleculeTypes[t].addFTCMatrix(molecules[i]);
			}
		}
	}

	for (int t = 0; t < moleculeTypes.size(); t++){
		if (occurence[t] > 0){
			moleculeTypes[t].divideFTCMatrix(occurence[t]);
		}
	}
}

void removeNonBondedInteractionsFromFTCMatricies(vector<molecule> &molecules){

	for (int i = 0; i < molecules.size(); i++){
		molecules[i].removeNonBondedInteractionsFTCMatrix();
	}
}

#endif
