#ifndef fragment_CLASS
#define fragment_CLASS

class fragment{
	private:
		//Basic declarations for each fragment
		int fragmentSize;
		int fragmentType;
		std::string fragmentName;
		double fragmentMass;
		double fragmentForcex;
		double fragmentForcey;
		double fragmentForcez;
		double fragmentTorquex;
		double fragmentTorquey;
		double fragmentTorquez;

		//Useful bonus variables for averaging
		int occurrence;
		double fragmentForceAverage;
		double fragmentForceSquared;
		
		//Vector of pointers to each atom in the fragment
		std::vector<atom*> atoms;

		//Vector of integers listing the position of covalently bonded fragments in the same molecule
		std::vector<int> bondedFragments;
				
	public:
		//Default constructor
		fragment(){fragmentSize = 0; fragmentType = 0; fragmentMass = 0; fragmentName = "";
				fragmentForcex = 0; fragmentForcey = 0; fragmentForcez = 0; 
				fragmentTorquex = 0; fragmentTorquey = 0; fragmentTorquez = 0;
				occurrence = 0; fragmentForceAverage = 0; fragmentForceSquared = 0;}

		//Overloaded member operators
		bool operator==(fragment test) const;
		bool operator!=(fragment test) const;
		
		//Pushback functions to add an atom
		void addAtom(atom* A){
			atoms.push_back(A);
			fragmentName.append(A->getName());
			fragmentMass = fragmentMass + A->getMass();
			fragmentSize = fragmentSize + 1;
		}
		//Function to delete an atom. Expensive, avoid.
		void delAtom(int i){
			fragmentSize = fragmentSize - 1;
			fragmentMass = fragmentMass - atoms[i]->getMass();
			atoms.erase(atoms.begin()+i);
		}
		
		//Add a covelently bonded fragment
		void addBondedFragment(int F){
			bondedFragments.push_back(F);
		}

		//Delete a covelently bonded fragment
		void delBondedFragment(int i){
			bondedFragments.erase(bondedFragments.begin()+i);
		}

		//Basic fragment get functions
		int getSize(){return fragmentSize;}
		int getType(){return fragmentType;}
		std::string getName(){return fragmentName;}
		double getMass(){return fragmentMass;}
		double getForcex(){return fragmentForcex;}
		double getForcey(){return fragmentForcey;}
		double getForcez(){return fragmentForcez;}
		double getTorquex(){return fragmentTorquex;}
		double getTorquey(){return fragmentTorquey;}
		double getTorquez(){return fragmentTorquez;}
		double getForceAverage(){return fragmentForceAverage;}
		double getForceSquared(){return fragmentForceSquared;}
		int getOccurrence(){ return occurrence; }

		//The heavy atom defines the position of the fragment
		double getPositionx(){return atoms[0]->getPositionx();}
		double getPositiony(){return atoms[0]->getPositiony();}
		double getPositionz(){return atoms[0]->getPositionz();}

		//Atom get functions
		atom* getAtom(int i){return atoms[i];}
		atom returnAtom(int i){return *atoms[i];}

		//Bonded fragment get functions
		int getBondedFragment(int i){return bondedFragments[i];}
		int getNumBondedFragments(){return bondedFragments.size();}

		//Basic fragment set functions
		void setType(int f){ fragmentType = f;}
		void setMass(double f){fragmentMass = f;}
		void setForcex(double f){fragmentForcex = f;}
		void setForcey(double f){fragmentForcey = f;}
		void setForcez(double f){fragmentForcez = f;}
		void setTorquex(double t){fragmentTorquex = t;}
		void setTorquey(double t){fragmentTorquey = t;}
		void setTorquez(double t){fragmentTorquez = t;}

		//Functions to increment properties
		void addMass(double f){fragmentMass = fragmentMass + f;}
		void addForcex(double f){fragmentForcex = fragmentForcex + f;}
		void addForcey(double f){fragmentForcey = fragmentForcey + f;}
		void addForcez(double f){fragmentForcez = fragmentForcez + f;}
		void addTorquex(double t){fragmentTorquex = fragmentTorquex + t;}
		void addTorquey(double t){ fragmentTorquey = fragmentTorquey + t; }
		void addTorquez(double t){ fragmentTorquez = fragmentTorquez + t; }
		void setForceAverage(double f){fragmentForceAverage = f;}
		void setForceSquared(double f){fragmentForceSquared = f;}

		void addForceAverage(double f){fragmentForceAverage = fragmentForceAverage+ f;}
		void addForceSquared(double f){fragmentForceSquared = fragmentForceSquared+ f;}
		void incrementOccurence(){ occurrence = occurrence + 1; }

		void orderAtomsUsingIndex();
		void orderAtomsUsingMass();
};

void fragment::orderAtomsUsingIndex(){
	for(int b = 0; b < this->getSize(); b++){
		for(int a = 1; a < this->getSize(); a++){
			if(atoms[a]->getIndex() < atoms[a-1]->getIndex()){
				atom *tempAtom1 = atoms[a];
				atoms[a] = atoms[a-1];
				atoms[a-1] = tempAtom1;
			}
			fragmentName = "";
			for(int i = 0; i < atoms.size(); i++){
				fragmentName.append(atoms[i]->getName());
			}
		}
	}
}

void fragment::orderAtomsUsingMass(){
	for(int b = 0; b < this->getSize(); b++){
		for(int a = 1; a < this->getSize(); a++){
			if(atoms[a]->getMass() > atoms[a-1]->getMass()){
				atom *tempAtom1 = atoms[a];
				atoms[a] = atoms[a-1];
				atoms[a-1] = tempAtom1;
			}
			fragmentName = "";
			for(int i = 0; i < atoms.size(); i++){
				fragmentName.append(atoms[i]->getName());
			}
		}
	}
}	

bool fragment::operator==(fragment test)const{

	bool isTheSame = true;
	if (fragmentSize != test.getSize()){
		isTheSame = false;
	}

	for (int i = 0; i < fragmentSize && isTheSame; i++){
		if (*atoms[i] != test.returnAtom(i)){
			isTheSame = false;
		}
	}

	return isTheSame;
}

bool fragment::operator!=(fragment test)const{
	
	bool isTheSame = true;
	if (fragmentSize != test.getSize()){
		isTheSame = false;
	}

	for (int i = 0; i < fragmentSize && isTheSame; i++){
		if (*atoms[i] != test.returnAtom(i)){
			isTheSame = false;
		}
	}

	return !isTheSame;
}

#endif
