/* Force Parsing Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 14/02/2014
**
** This function updates a vector of each different molecule + shell size combination in
** a vector of molecules.
**
*/

#ifndef UPDATE_SHELL_SIZE_TYPES
#define UPDATE_SHELL_SIZE_TYPES

void updateMoleculeShellSizeTypes(vector<molecule> &molecules, vector<molecule> &moleculeShellTypes);
void updateMoleculeShellSizeTypes(vector<molecule> &molecules, vector<molecule> &moleculeShellTypes, molecule* moleculei);

void updateMoleculeShellSizeTypes(vector<molecule> &molecules, vector<molecule> &moleculeShellTypes){
	updateMoleculeShellSizeTypes(molecules, moleculeShellTypes, nullptr);
}

void updateMoleculeShellSizeTypes(vector<molecule> &molecules, vector<molecule> &moleculeShellTypes, molecule* moleculei){

	for (int i = 0; i < molecules.size(); i++){

		bool checki = false;
		if(moleculei == nullptr){
			checki = true;
		}
		else if(molecules[i].getType() == moleculei->getType()){
			checki = true;
		}

		if(checki){	

			bool inTypes = false;

			for (int t = 0; t < moleculeShellTypes.size(); t++){

				for (int j = 0; j < molecules[i].getSize(); j++){
					if (molecules[i].getShell(j)->getSize() == moleculeShellTypes[t].getShell(j)->getSize()){
						inTypes = true;
						moleculeShellTypes[t].incrementOccurence();
					}
				}
			}
			if (inTypes == false){
				moleculeShellTypes.push_back(molecules[i]);
				moleculeShellTypes[moleculeShellTypes.size() - 1].incrementOccurence();
			}
		}
	}
}

#endif
