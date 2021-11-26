/* Molecule Typing Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 14/02/2014
**
** This collection of functions runs through a vector of the molecule class
** and either returns a vector containing one of each type, or gives them 
** a type equal to their position within that vector.
**
*/

#ifndef MOLECULE_TYPES
#define MOLECULE_TYPES

void setMoleculeTypes(vector<molecule> &molecules);
vector <molecule> getMoleculeTypes(vector<molecule> &molecules);

void setMoleculeTypes(vector<molecule> &molecules){

	vector<molecule> moleculeTypes;

	for(int i = 0; i < molecules.size(); i++){
		bool inTypes = false;		
		for(int j = 0; j < moleculeTypes.size(); j++){
	
			if(molecules[i] == moleculeTypes[j]){
				inTypes = true;
				molecules[i].setType(j);
			}
		}
		if(inTypes == false){
			molecules[i].setType(moleculeTypes.size());
			moleculeTypes.push_back(molecules[i]);
		}
	}
}

vector <molecule> getMoleculeTypes(vector<molecule> &molecules){

	vector<molecule> moleculeTypes;

	for(int i = 0; i < molecules.size(); i++){
		bool inTypes = false;
		for(int j = 0; j < moleculeTypes.size(); j=j+1){
	
			if(molecules[i] == moleculeTypes[j]){
				inTypes = true;
			}
		}
		if(inTypes == false){
			moleculeTypes.push_back(molecules[i]);
		}
	}

	return moleculeTypes;
}


#endif
