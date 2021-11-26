/* Force Parsing Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 14/02/2014
**
** This collection of functions creates and acts on an array of "fragments", 
** each of which is an object containing a vector of atoms
** 
**
*/

#ifndef FRAGMENT_TYPES
#define FRAGMENT_TYPES

void setFragmentTypes(vector<fragment> &fragments);
vector <fragment> getFragmentTypes(vector<fragment> &fragments);

void setFragmentTypes(vector<fragment> &fragments){

	vector<fragment> fragmentTypes;

	for(int i = 0; i < fragments.size(); i++){
		bool inTypes = false;		
		for(int j = 0; j < fragmentTypes.size(); j++){
	
			//Add hydrogens to each fragment
			if(fragments[i] == fragmentTypes[j]){
				inTypes = true;
				fragments[i].setType(j);
			}
		}
		if(inTypes == false){
			fragments[i].setType(fragmentTypes.size());
			fragmentTypes.push_back(fragments[i]);
		}
	}
}

vector <fragment> getFragmentTypes(vector<fragment> &fragments){

	vector<fragment> fragmentTypes;

	for(int i = 0; i < fragments.size(); i++){
		bool inTypes = false;
		for(int j = 0; j < fragmentTypes.size(); j=j+1){
	
			//Add hydrogens to each fragment
			if(fragments[i] == fragmentTypes[j]){
				inTypes = true;
			}
		}
		if(inTypes == false){
			fragmentTypes.push_back(fragments[i]);
		}
	}

	return fragmentTypes;
}


#endif
