/* Fragment RAD, Dr Jonathan Higham
** Created 22/08/2015
** Modified 26/07/2016
**
** This function applies the RAD criteria to a vector of fragments (usually heavy atom + hydrogen groups)
** and determines which fragments share a binary correlation.
**
*/

#ifndef RAD_MOLECULE_SHELLS
#define RAD_MOLECULE_SHELLS


void RADMoleculeShells(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells, float* boxSize);
void RADMoleculeShells(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells, float* boxSize, molecule* moleculei, molecule* moleculej);

void RADMoleculeShells(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells, float* boxSize){
	RADMoleculeShells(molecules, moleculeShells, boxSize, nullptr, nullptr);
}

void RADMoleculeShells(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells, float* boxSize, molecule* moleculei, molecule* moleculej){

	for (int i = 0; i < molecules.size(); i++){

		moleculeShell tempShell;
		moleculeShells[i] = tempShell;
	}

	vector< vector<double> > distances(molecules.size(), vector<double>(molecules.size(), 0));
	
	//Loop over each fragment
	for (int i = 0; i < molecules.size(); i++){

		//Loop over each fragment again
		for (int j = 0; j < molecules.size(); j++){

			double r = 0;

			//Add the squared distance in x,y,z
			double vecx = molecules[i].getCOMPositionx() - molecules[j].getCOMPositionx();
			double vecy = molecules[i].getCOMPositiony() - molecules[j].getCOMPositiony();
			double vecz = molecules[i].getCOMPositionz() - molecules[j].getCOMPositionz();

			//Update anything that goes over a box boundary
			if (vecx > (boxSize[0] / 2)){
				vecx = vecx - boxSize[0];
			}
			if (vecx < -(boxSize[0] / 2)){
				vecx = vecx + boxSize[0];
			}
			if (vecy > (boxSize[1] / 2)){
				vecy = vecy - boxSize[1];
			}
			if (vecy < -(boxSize[1] / 2)){
				vecy = vecy + boxSize[1];
			}
			if (vecz > (boxSize[2] / 2)){
				vecz = vecz - boxSize[2];
			}
			if (vecz < -(boxSize[2] / 2)){
				vecz = vecz + boxSize[2];
			}

			//Update distances
			distances[i][j] = sqrt((vecx*vecx) + (vecy*vecy) + (vecz*vecz));
			distances[j][i] = distances[i][j];

		}
	}

	for (int i = 0; i < molecules.size(); i++){

		vector<bool> moleculeNotAlreadyInShell(molecules.size(), true);

		bool inShell = true;
		vector<int> neighboursj;

		bool checki = false;
		if(moleculei == nullptr){
			checki = true;
		}
		else if(molecules[i].getType() == moleculei->getType()){
			checki = true;
		}

		while (inShell && checki){

			int j12;
			double r12 = 0;

			//Loop over each fragment again
			for (int j = 0; j < molecules.size(); j++){

				if ((distances[i][j] < r12 || r12 == 0) && moleculeNotAlreadyInShell[j] && j != i){
					j12 = j;
					r12 = distances[i][j];
				}
			}
				
			neighboursj.push_back(j12);
			moleculeNotAlreadyInShell[j12] = false;
			double a12 = 1 / pow(r12, 2);
			double vecx12 = molecules[i].getCOMPositionx() - molecules[j12].getCOMPositionx();
			double vecy12 = molecules[i].getCOMPositiony() - molecules[j12].getCOMPositiony();
			double vecz12 = molecules[i].getCOMPositionz() - molecules[j12].getCOMPositionz();

			//Update anything that goes over a box boundary
			if (vecx12 > (boxSize[0] / 2))vecx12 = vecx12 - boxSize[0]; if (vecx12 < -(boxSize[0] / 2))vecx12 = vecx12 + boxSize[0];
			if (vecy12 > (boxSize[1] / 2))vecy12 = vecy12 - boxSize[1]; if (vecy12 < -(boxSize[1] / 2))vecy12 = vecy12 + boxSize[1];
			if (vecz12 > (boxSize[2] / 2))vecz12 = vecz12 - boxSize[2]; if (vecz12 < -(boxSize[2] / 2))vecz12 = vecz12 + boxSize[2];

			bool checkj = false;
			if(moleculej == nullptr){
				checkj = true;
			}
			else if(molecules[j12].getType() == moleculej->getType()){
				checkj = true;
			}

			bool inShell2 = true;
			if(checkj){
				for (int d = 0; d < neighboursj.size() - 1; d++){
					//Distance between atom and probe
					int j13 = neighboursj[d];
					double r13 = distances[i][j13];
					double a13 = 1 / pow(r13, 2);

					double vecx13 = molecules[i].getCOMPositionx() - molecules[j13].getCOMPositionx();
					double vecy13 = molecules[i].getCOMPositiony() - molecules[j13].getCOMPositiony();
					double vecz13 = molecules[i].getCOMPositionz() - molecules[j13].getCOMPositionz();
					if (vecx13 >(boxSize[0] / 2))vecx13 = vecx13 - boxSize[0]; if (vecx13 < -(boxSize[0] / 2))vecx13 = vecx13 + boxSize[0];
					if (vecy13 > (boxSize[1] / 2))vecy13 = vecy13 - boxSize[1]; if (vecy13 < -(boxSize[1] / 2))vecy13 = vecy13 + boxSize[1];
					if (vecz13 > (boxSize[2] / 2))vecz13 = vecz13 - boxSize[2]; if (vecz13 < -(boxSize[2] / 2))vecz13 = vecz13 + boxSize[2];

					double r13dotr12 = (vecx12*vecx13) + (vecy12*vecy13) + (vecz12*vecz13);
					double cosTheta = (r13dotr12 / (r13*r12));

					if (a13*cosTheta > a12 && i != j13){
						inShell = false;
						moleculeNotAlreadyInShell[j12] = true;
					}
					if (a13*cosTheta > a12 && i == j13){
						inShell2 = false;
						moleculeNotAlreadyInShell[j12] = true;
					}
				}
				if (inShell && inShell2 && (i != j12)){
					moleculeShells[i].addMolecule(&molecules[j12]);
				}
			}
		}
	}
}//End of function

#endif