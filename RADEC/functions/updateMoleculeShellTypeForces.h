/* Force Parsing Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 14/02/2014
**
** This function updates a vector of each different molecule + shell combination in
** a vector of molecules.
**
*/

#ifndef UPDATE_MOLECULE_SHELL_TYPE_FORCES
#define UPDATE_MOLECULE_SHELL_TYPE_FORCES

void updateMoleculeShellTypeForces(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells, vector<molecule> &moleculeShellMoleculeTypes, vector<moleculeShell> &moleculeShellShellTypes);
void updateMoleculeShellTypeForces(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells, vector<molecule> &moleculeShellMoleculeTypes, vector<moleculeShell> &moleculeShellShellTypes, molecule* moleculei);

void updateMoleculeShellTypeForces(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells, vector<molecule> &moleculeShellMoleculeTypes, vector<moleculeShell> &moleculeShellShellTypes){
	updateMoleculeShellTypeForces(molecules, moleculeShells, moleculeShellMoleculeTypes, moleculeShellShellTypes, nullptr);
}

void updateMoleculeShellTypeForces(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells, vector<molecule> &moleculeShellMoleculeTypes, vector<moleculeShell> &moleculeShellShellTypes, molecule* moleculei){

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
		
			for (int t = 0; t < moleculeShellMoleculeTypes.size(); t++){

				if (molecules[i] == moleculeShellMoleculeTypes[t] && moleculeShells[i] == moleculeShellShellTypes[t]){

					inTypes = true;
					moleculeShellMoleculeTypes[t].incrementOccurence();
					moleculeShellMoleculeTypes[t].addForceAverage(molecules[i].getForceAverage());
					moleculeShellMoleculeTypes[t].addForceSquared(molecules[i].getForceSquared());
				}
			}

			if (inTypes == false){
				moleculeShellMoleculeTypes.push_back(molecules[i]);
				moleculeShellShellTypes.push_back(moleculeShells[i]);

				moleculeShellMoleculeTypes[moleculeShellMoleculeTypes.size() - 1].incrementOccurence();
				moleculeShellMoleculeTypes[moleculeShellMoleculeTypes.size() - 1].addForceAverage(molecules[i].getForceAverage());
				moleculeShellMoleculeTypes[moleculeShellMoleculeTypes.size() - 1].addForceSquared(molecules[i].getForceSquared());
			}
		}
	}
}

#endif
