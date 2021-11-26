/* Molecule sorting Code, Dr Jonathan Higham
**
** This function runs through a vector of the molecule class
** and sorts it based on some metric
**
*/

#ifndef SORT_MOLECULE
#define SORT_MOLECULE

void sortMoleculesBySize(vector<molecule> &molecules);
void sortMoleculesByShellSize(vector<molecule> &molecules);

void sortMoleculesBySize(vector<molecule> &molecules){


	for(int j = 0; j < molecules.size(); j++){
		for(int i = 0; i < molecules.size()-1; i++){
	
			//Find the total shell sizes for each molecule
			int iShellSize = 0;

				iShellSize = molecules[i].getSize();

			int jShellSize = 0;

				jShellSize = molecules[i+1].getSize();


			if(iShellSize > jShellSize){
				molecule tempMolecule;
				tempMolecule = molecules[i];
				molecules[i] = molecules[i+1];
				molecules[i+1] = tempMolecule;
			}
		}
	}
}

void sortMoleculesByShellSize(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells){


	if(molecules.size() == moleculeShells.size()){
		for(int j = 0; j < molecules.size(); j++){
			for(int i = 0; i < molecules.size()-1; i++){
		
				//Find the total shell sizes for each molecule
				int iShellSize = moleculeShells[i].getSize();
				int jShellSize = moleculeShells[i+1].getSize();

				if(iShellSize > jShellSize){
					molecule tempMolecule;
					tempMolecule = molecules[i];
					molecules[i] = molecules[i+1];
					molecules[i+1] = tempMolecule;
				
					moleculeShell tempMoleculeShell;
					tempMoleculeShell = moleculeShells[i];
					moleculeShells[i] = moleculeShells[i+1];
					moleculeShells[i+1] = tempMoleculeShell;
				}
			}
		}
	}
}

#endif
