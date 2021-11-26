#ifndef shell_CLASS
#define shell_CLASS

class shell{
	private:
		//Basic declarations for a shell
		int shellSize;
		std::vector<fragment*> fragments;		

	public:
		//Default constructor
		shell(){ shellSize = 0; }

		//Overloaded member operators
		bool operator==(shell test) const;
		bool operator!=(shell test) const;
		
		//Pushback function to add a fragment
		void addFragment(fragment* F){
			fragments.push_back(F);
			shellSize = shellSize + 1;
		}

		//Function to delete a fragment. Expensive, avoid.
		void delFragment(int i){
			shellSize = shellSize - 1;
			fragments.erase(fragments.begin() + i);
		}

		void clearShell(){
			fragments.clear();
			shellSize = 0;
		}
		
		//Basic fragment get functions
		int getSize(){return shellSize;}
		fragment* getFragment(int i){ return fragments[i]; }
		fragment returnFragment(int i){ return *fragments[i]; }

		void orderfragmentsUsingMass();
		void orderfragmentsUsingType();
};

void shell::orderfragmentsUsingMass(){
	for (int b = 0; b < this->getSize(); b++){
		for (int a = 1; a < this->getSize(); a++){
			if (fragments[a]->getMass() < fragments[a - 1]->getMass()){
				fragment *tempFragment = fragments[a];
				fragments[a] = fragments[a - 1];
				fragments[a - 1] = tempFragment;
			}
		}
	}
}

void shell::orderfragmentsUsingType(){
	for (int b = 0; b < this->getSize(); b++){
		for (int a = 1; a < this->getSize(); a++){
			if (fragments[a]->getType() < fragments[a - 1]->getType()){
				fragment *tempFragment = fragments[a];
				fragments[a] = fragments[a - 1];
				fragments[a - 1] = tempFragment;
			}
		}
	}
}

bool shell::operator==(shell test)const{

	bool isTheSame = true;
	if (shellSize != test.getSize()){
		isTheSame = false;
	}

	for (int i = 0; i < shellSize && isTheSame; i++){
		if (*fragments[i] != test.returnFragment(i)){
			isTheSame = false;
		}
	}

	return isTheSame;
}

bool shell::operator!=(shell test)const{
	
	bool isTheSame = true;
	if (shellSize != test.getSize()){
		isTheSame = false;
	}

	for (int i = 0; i < shellSize && isTheSame; i++){
		if (*fragments[i] != test.returnFragment(i)){
			isTheSame = false;
		}
	}

	return !isTheSame;
}

#endif