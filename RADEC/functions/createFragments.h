/* Force Parsing Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 14/02/2014
**
** This collection of functions creates and acts on an array of "fragments", 
** each of which is an object containing a vector of atoms
** 
**
*/

#ifndef CREATE_FRAGMENTS
#define CREATE_FRAGMENTS

void createFragments(string topFile, vector<atom> &atoms, vector<fragment> &fragments);
vector <fragment> getFragmentCopy(vector <fragment> inputFragments);

void createFragments(string topFile, vector<atom> &atoms, vector<fragment> &fragments){

	char A; int I;
	bool flag = true;
	int temp;

	vector<int> HBonds;
	vector<int> coBonds;
	vector<bool> ignorefragment;
	fragment tempfragment;


	//Open solute topology file	
	ifstream TopologyIn;
	TopologyIn.open(topFile.c_str());

	//Find the start of the bond declarations
	streampos HBondStart = findHBondStart(topFile);
	TopologyIn.seekg(HBondStart);

	int lineCheck = 0;	
	while(TopologyIn){
		//Read in the two atoms in each bond

		string str;
		getline(TopologyIn, str);
		if(str[0] == '%'){
			break;
		}
		stringstream ss(str);

		while(ss >> temp){
			if(lineCheck < 2){
				HBonds.push_back(temp);
				lineCheck++;
			}
			else if(lineCheck == 2){
				lineCheck = 0;
			}
		}
	}

	for(int i = 0; i < HBonds.size(); i++){
		HBonds[i] = (abs(HBonds[i])/3);
	}

	for(int i = 0; i < atoms.size(); i++){
		//Construct the central atom of each fragment
		fragment tempfragment;
		tempfragment.addAtom(&atoms[i]);
		fragments.push_back(tempfragment);
		//Ignore the fragment unless it has a mass (to deal with point charges and etc.)
		if(atoms[i].getMass() > 0){
			ignorefragment.push_back(0);
		}
		else{
			ignorefragment.push_back(1);
		}
	}

	int l = 0;
	for(int i = 0; i < fragments.size()-l; i++){
		for(int j = 0; j < HBonds.size()-(2*l); j=j+2){
	
			//Add hydrogens to each fragment

			if(atoms[HBonds[j+1]].getName() == "H"){
				if(fragments[i].getAtom(0)->getIndex() == HBonds[j] && !ignorefragment[i]){
					fragments[i].addAtom(&atoms[HBonds[j+1]]);
					fragments[i].getAtom(fragments[i].getSize()-1)->setIndex(HBonds[j+1]);
					ignorefragment[HBonds[j+1]] = 1; //Ignore fragment after addition
				}
			}
			else if(atoms[HBonds[j]].getName() == "H"){
				if(fragments[i].getAtom(0)->getIndex() == HBonds[j+1] && !ignorefragment[i]){
					fragments[i].addAtom(&atoms[HBonds[j]]);
					fragments[i].getAtom(fragments[i].getSize()-1)->setIndex(HBonds[j]);
					ignorefragment[HBonds[j]] = 1; //Ignore fragment after addition
				}
			}
		}
	}

	for(int i = fragments.size()-1; i >= 0; i--){
		if(ignorefragment[i] == 1){
			fragments.erase(fragments.begin()+i);
		}
	}

	for(int i = 0; i < fragments.size(); i++){
		fragments[i].orderAtomsUsingMass();
	}
}

vector <fragment> getFragmentCopy(vector <fragment> inputFragments){
	vector <fragment> outputFragments;
	for(int i = 0; i < inputFragments.size(); i++){
		outputFragments.push_back(inputFragments[i]);
	}
	return outputFragments;
}

#endif
