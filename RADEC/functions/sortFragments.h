/* Molecule sorting Code, Dr Jonathan Higham
**
** This function runs through a vector of the fragment class
** and sorts it based on some metric
**
*/

#ifndef SORT_FRAGMENTS
#define SORT_FRAGMENTS

void sortFragmentsByShellSize(vector<fragment> &fragments, vector<shell> &fragmentShells);

void sortFragmentsByShellSize(vector<fragment> &fragments, vector<shell> &fragmentShells){


	for(int j = 0; j < fragments.size(); j++){
		for(int i = 0; i < fragments.size()-1; i++){
	
			//Find the total shell sizes for each molecule
			int iShellSize = fragmentShells[i].getSize();
			int jShellSize = fragmentShells[i+1].getSize();

			if(iShellSize > jShellSize){
				fragment tempFragment;
				tempFragment = fragments[i];
				fragments[i] = fragments[i+1];
				fragments[i+1] = tempFragment;

				shell tempShell;
				tempShell = fragmentShells[i];
				fragmentShells[i] = fragmentShells[i+1];
				fragmentShells[i+1] = tempShell;
			}
		}
	}
}

#endif