/* Force Parsing Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 14/02/2014
**
** This function updates a vector of each different molecule + shell combination in
** a vector of molecules.
**
*/

#ifndef CREATE_MOLECULE_SHELLS
#define CREATE_MOLECULE_SHELLS

void createMoleculeShells(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells);

void createMoleculeShells(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells){

	for (int i = 0; i < molecules.size(); i++){

		moleculeShell tempShell;
		moleculeShells.push_back(tempShell);
	}
}

#endif