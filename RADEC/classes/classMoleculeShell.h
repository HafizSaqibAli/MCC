#ifndef MOLECULE_SHELL_CLASS
#define MOLECULE_SHELL_CLASS

class moleculeShell{
	private:
		//Basic declarations for a shell
		int shellSize;
		std::vector<molecule*> molecules;		

	public:
		//Default constructor
		moleculeShell(){ shellSize = 0; }

		//Overloaded member operators
		bool operator==(moleculeShell test) const;
		bool operator!=(moleculeShell test) const;
		
		//Pushback function to add a molecule
		void addMolecule(molecule* F){
			molecules.push_back(F);
			shellSize = shellSize + 1;
		}

		//Function to delete a molecule. Expensive, avoid.
		void delMolecule(int i){
			shellSize = shellSize - 1;
			molecules.erase(molecules.begin() + i);
		}

		void clearShell(){
			molecules.clear();
			shellSize = 0;
		}
		
		//Basic molecule get functions
		int getSize(){return shellSize;}
		molecule* getMolecule(int i){ return molecules[i]; }
		molecule returnMolecule(int i){ return *molecules[i]; }

		void ordermoleculesUsingMass();
		void ordermoleculesUsingType();
};

void moleculeShell::ordermoleculesUsingMass(){
	for (int b = 0; b < this->getSize(); b++){
		for (int a = 1; a < this->getSize(); a++){
			if (molecules[a]->getMass() < molecules[a - 1]->getMass()){
				molecule *tempMolecule = molecules[a];
				molecules[a] = molecules[a - 1];
				molecules[a - 1] = tempMolecule;
			}
		}
	}
}

void moleculeShell::ordermoleculesUsingType(){
	for (int b = 0; b < this->getSize(); b++){
		for (int a = 1; a < this->getSize(); a++){
			if (molecules[a]->getType() < molecules[a - 1]->getType()){
				molecule *tempMolecule = molecules[a];
				molecules[a] = molecules[a - 1];
				molecules[a - 1] = tempMolecule;
			}
		}
	}
}

bool moleculeShell::operator==(moleculeShell test)const{

	bool isTheSame = true;
	if (shellSize != test.getSize()){
		isTheSame = false;
	}

	for (int i = 0; i < shellSize && isTheSame; i++){
		if (*molecules[i] != test.returnMolecule(i)){
			isTheSame = false;
		}
	}

	return isTheSame;
}

bool moleculeShell::operator!=(moleculeShell test)const{
	
	bool isTheSame = true;
	if (shellSize != test.getSize()){
		isTheSame = false;
	}

	for (int i = 0; i < shellSize && isTheSame; i++){
		if (*molecules[i] != test.returnMolecule(i)){
			isTheSame = false;
		}
	}

	return !isTheSame;
}

#endif