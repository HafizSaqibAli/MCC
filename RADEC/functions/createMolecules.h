/* 
** Created 13/02/2014
** Modified 15/07/2016
**
** Function to add fragments together to form molecules based on an AMBER topology file
**
*/

#ifndef CREATE_MOLECULES
#define CREATE_MOLECULES

void createMolecules(string topFile, vector<fragment> &fragments, vector<molecule> &molecules);

void createMolecules(string topFile, vector<fragment> &fragments, vector<molecule> &molecules){

	vector<int> coBonded1;
	vector<int> coBonded2;
	vector<int> centralFragment;
	vector<bool> isCentralFragment;

	for(int i = 0; i < fragments.size(); i++){
		//Construct the central atom of each fragment
		molecule tempMolecule;
		tempMolecule.addFragment(&fragments[i]);
		molecules.push_back(tempMolecule);

		centralFragment.push_back(-1);
		isCentralFragment.push_back(false);
	}

	//Open solute topology file
	ifstream TopologyIn;
	TopologyIn.open(topFile.c_str());

	//Find the start of the covalent bond declarations
	streampos coBondStart = findCoBondStart(topFile);
	TopologyIn.seekg(coBondStart);

	int temp;
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
			if(lineCheck == 0){
				coBonded1.push_back(temp);
				lineCheck++;
			}
			else if(lineCheck == 1){
				coBonded2.push_back(temp);
				lineCheck++;
			}
			else if(lineCheck == 2){
				lineCheck = 0;
			}
		}
	}

	if(coBonded1.size() > 0){

		for(int c = 0; c < coBonded1.size(); c++){
			coBonded1[c] = (abs(coBonded1[c])/3);
			coBonded2[c] = (abs(coBonded2[c])/3); 
		}

		for(int c = 0; c < coBonded1.size(); c++){

			int mIndex1 = -1;
			int mIndex2 = -1;

			for(int i = 0; i < molecules.size(); i++){
				if(molecules[i].getFragment(0)->getAtom(0)->getIndex() == coBonded1[c]){
					mIndex1 = i;
				}
				if(molecules[i].getFragment(0)->getAtom(0)->getIndex() == coBonded2[c]){
					mIndex2 = i;
				}
			}

			if (mIndex1 != -1 && mIndex2 != -1){

				//If neither is a central fragment nor connected to a central fragment
				if (!isCentralFragment[mIndex1] && !isCentralFragment[mIndex2] &&
					centralFragment[mIndex1] == -1 && centralFragment[mIndex2] == -1){
					molecules[mIndex2].addFragment(molecules[mIndex1].getFragment(0));
					centralFragment[mIndex1] = mIndex2;
					isCentralFragment[mIndex2] = true;
				}
				//If 1 is central fragment and 2 isn't, and 2 is not currently assigned to a central fragment
				else if (isCentralFragment[mIndex1] && !isCentralFragment[mIndex2] &&
					centralFragment[mIndex2] == -1){
					molecules[mIndex1].addFragment(molecules[mIndex2].getFragment(0));
					centralFragment[mIndex2] = mIndex1;
				}
				//If 2 is central fragment and 1 isn't, and 1 is not currently assigned to a central fragment
				else if (isCentralFragment[mIndex2] && !isCentralFragment[mIndex1] &&
					centralFragment[mIndex1] == -1){
					molecules[mIndex2].addFragment(molecules[mIndex1].getFragment(0));
					centralFragment[mIndex1] = mIndex2;
				}
				//If neither is a central fragment, but 1 is connected to a central fragment
				else if (!isCentralFragment[mIndex1] && !isCentralFragment[mIndex2] &&
					centralFragment[mIndex1] > -1 && centralFragment[mIndex2] == -1){
					molecules[centralFragment[mIndex1]].addFragment(molecules[mIndex2].getFragment(0));
					centralFragment[mIndex2] = centralFragment[mIndex1];
				}
				//If neither is a central fragment, but 2 is connected to a central fragment
				else if (!isCentralFragment[mIndex1] && !isCentralFragment[mIndex2] &&
					centralFragment[mIndex1] == -1 && centralFragment[mIndex2] > -1){
					molecules[centralFragment[mIndex2]].addFragment(molecules[mIndex1].getFragment(0));
					centralFragment[mIndex1] = centralFragment[mIndex2];
				}

				//If both fragments are central fragments
				else if (isCentralFragment[mIndex1] && isCentralFragment[mIndex2]){
					molecules[mIndex2].addMolecule(molecules[mIndex1]);
					centralFragment[mIndex1] = mIndex2;
					isCentralFragment[mIndex1] = false;
				}
				//If one fragment is a central fragment and the other is connected to a different one
				else if (isCentralFragment[mIndex1] && centralFragment[mIndex2] > -1 &&
					centralFragment[mIndex2] != mIndex1){
					molecules[centralFragment[mIndex2]].addMolecule(molecules[mIndex1]);
					centralFragment[mIndex1] = mIndex2;
					isCentralFragment[mIndex1] = false;
				}
				else if (isCentralFragment[mIndex2] && centralFragment[mIndex1] > -1 &&
					centralFragment[mIndex1] != mIndex2){
					molecules[centralFragment[mIndex1]].addMolecule(molecules[mIndex2]);
					centralFragment[mIndex2] = mIndex1;
					isCentralFragment[mIndex1] = false;
				}

				//If both fragments are linked to each other AND different central fragments
				else if (centralFragment[mIndex1] > -1 && centralFragment[mIndex2] > -1 &&
					centralFragment[mIndex1] != centralFragment[mIndex2]){

					molecules[centralFragment[mIndex1]].addMolecule(molecules[centralFragment[mIndex2]]);

					centralFragment[centralFragment[mIndex2]] = centralFragment[mIndex1];
					isCentralFragment[centralFragment[mIndex2]] = false;
	
					int centralFragmentToBeRemoved = centralFragment[mIndex2];
					for(int i = 0; i < molecules.size(); i++){
						if(centralFragment[i] == centralFragmentToBeRemoved ){
							centralFragment[i] = centralFragment[mIndex1];
						}
					}
				}

				//If both fragments are already on the same molecule due to other bonds
				else if (centralFragment[mIndex1] > -1 && centralFragment[mIndex2] > -1 &&
					centralFragment[mIndex1] == centralFragment[mIndex2]){
				}
				//If both fragments are already on the same molecule due to other bonds and 1 is the central fragment
				else if (centralFragment[mIndex1] == -1 && centralFragment[mIndex2] == mIndex1){
				}
				//If both fragments are already on the same molecule due to other bonds and 2 is the central fragment
				else if (centralFragment[mIndex2] == -1 && centralFragment[mIndex1] == mIndex2){
				}
				//Catch case just in case the topology is in an odd amber format
				else{
					cerr << "Weirdly ordered topology file (4), please re-generate " << endl;
					exit(0);
				}

			}
		}

		for(int i = molecules.size()-1; i >= 0; i--){
			if(!isCentralFragment[i] && centralFragment[i] > -1){
				molecules.erase(molecules.begin()+i);
				centralFragment.erase(centralFragment.begin()+i);
			}
		}

		//Loop over everything and add each fragment to its own bonded list first
		for(int i = 0; i < molecules.size(); i++){
			for(int j = 0; j < molecules[i].getSize(); j++){
				molecules[i].getFragment(j)->addBondedFragment(j);
			}
		}

		//Run through the cobonded list a final time and link covalently bonded fragments 
		for(int c = 0; c < coBonded1.size(); c++){
			//Loop over molecules and fragments
			for(int i = 0; i < molecules.size(); i++){
				for(int j = 0; j < molecules[i].getSize(); j++){					
					//Find the first atom in the molecule list
					if(molecules[i].getFragment(j)->getAtom(0)->getIndex() == coBonded1[c]){
						//Loop over other fragments in the molecule
						for(int j2 = 0; j2 < molecules[i].getSize(); j2++){
							//Find the second atom in the molecule list
							if(molecules[i].getFragment(j2)->getAtom(0)->getIndex() == coBonded2[c]){
								molecules[i].getFragment(j)->addBondedFragment(j2);
								molecules[i].getFragment(j2)->addBondedFragment(j);
							}
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < molecules.size(); i++){
		molecules[i].createFTCMatrix();
	}
}

#endif
