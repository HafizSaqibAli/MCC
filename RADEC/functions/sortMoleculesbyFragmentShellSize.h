/* Molecule sorting Code, Dr Jonathan Higham
**
** This function runs through a vector of the molecule class
** and sorts it based on some metric
**
*/

#ifndef SORT_MOLECULES_BY_FRAGMENT_SHELL_SIZE
#define SORT_MOLECULES_BY_FRAGMENT_SHELL_SIZE

void sortMoleculesByFragmentShellSize(vector<molecule> &molecules);

void sortMoleculesByFragmentShellSize(vector<molecule> &molecules){


	for(int j = 0; j < molecules.size(); j++){
		for(int i = 0; i < molecules.size()-1; i++){
	
			//Find the total shell sizes for each molecule
			int iShellSize = 0;
			for(int s1 = 0; s1 < molecules[i].getSize(); s1++){
				iShellSize = iShellSize + molecules[i].getFragmentShell(s1)->getSize();
			}
			int jShellSize = 0;
			for(int s2 = 0; s2 < molecules[i+1].getSize(); s2++){
				jShellSize = jShellSize + molecules[i+1].getFragmentShell(s2)->getSize();
			}

			if(iShellSize > jShellSize){
				molecule tempMolecule;
				tempMolecule = molecules[i];
				molecules[i] = molecules[i+1];
				molecules[i+1] = tempMolecule;
			}
		}
	}
}


#endif
