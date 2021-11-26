/* Force Parsing Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 14/02/2014
**
** This function updates a vector of each different molecule + shell combination in
** a vector of molecules.
**
*/

#ifndef UPDATE_SHELL_TYPE_FORCES
#define UPDATE_SHELL_TYPE_FORCES

void updateFragmentShellTypeForces(vector<molecule> &molecules, vector<fragment> &fragmentShellFragmentTypes, vector<shell> &fragmentShellShellTypes);
void updateFragmentShellTypeForces(vector<molecule> &molecules, vector<fragment> &fragmentShellFragmentTypes, vector<shell> &fragmentShellShellTypes, molecule* moleculei);

void updateFragmentShellTypeForces(vector<molecule> &molecules, vector<fragment> &fragmentShellFragmentTypes, vector<shell> &fragmentShellShellTypes){
	updateFragmentShellTypeForces(molecules, fragmentShellFragmentTypes, fragmentShellShellTypes, nullptr);
}

void updateFragmentShellTypeForces(vector<molecule> &molecules, vector<fragment> &fragmentShellFragmentTypes, vector<shell> &fragmentShellShellTypes, molecule* moleculei){

	for (int i = 0; i < molecules.size(); i++){
		bool checki = false;
		if(moleculei == nullptr){
			checki = true;
		}
		else if(molecules[i].getType() == moleculei->getType()){
			checki = true;
		}
		
		if(checki){			
			
			for (int j = 0; j < molecules[i].getSize(); j++){
		
				bool inTypes = false;
				for (int t = 0; t < fragmentShellFragmentTypes.size(); t++){

					if (molecules[i].returnFragment(j) == fragmentShellFragmentTypes[t] && molecules[i].returnShell(j) == fragmentShellShellTypes[t]){
							inTypes=true;
							fragmentShellFragmentTypes[t].incrementOccurence();
							fragmentShellFragmentTypes[t].addForceAverage(molecules[i].getForceAverage());
							fragmentShellFragmentTypes[t].addForceSquared(molecules[i].getForceSquared());
					}
				}

				if (inTypes == false){
					fragmentShellShellTypes.push_back(molecules[i].returnShell(j));
					fragmentShellFragmentTypes.push_back(molecules[i].returnFragment(j));
					fragmentShellFragmentTypes[fragmentShellFragmentTypes.size() - 1].addForceAverage(molecules[i].getForceAverage());
					fragmentShellFragmentTypes[fragmentShellFragmentTypes.size() - 1].addForceSquared(molecules[i].getForceSquared());
					fragmentShellFragmentTypes[fragmentShellFragmentTypes.size() - 1].incrementOccurence();
				}
			}
		}
	}
}

#endif