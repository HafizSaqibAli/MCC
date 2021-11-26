#ifndef ATOM_CLASS
#define ATOM_CLASS

class atom{
	private:
		//Declarations for all atom variables; type, charge, force and position
		int index;
		int type;
		std::string name;
		double mass;
		double charge;
		double xForce;
		double yForce;
		double zForce;
		double xPosition;
		double yPosition;
		double zPosition;

	public:
		//Default constructor
		atom(){index = 0; type = 'X'; mass = 0; charge = 0;
			xPosition = yPosition = zPosition = xForce = yForce = zForce = 0;}

		bool operator==(atom test) const;
		bool operator!=(atom test) const;

		//Get functions
		int getIndex(){return index;}
		int getType(){return type;}
		std::string getName(){return name;}
		double getMass(){return mass;}
		double getCharge(){return charge;}
		double getForcex(){return xForce;}
		double getForcey(){return yForce;}
		double getForcez(){return zForce;}
		double getPositionx(){return xPosition;}
		double getPositiony(){return yPosition;}
		double getPositionz(){return zPosition;}

		//Set functions
		void setIndex(int i){index = i;}
		void setType(int A){type = A;}
		void setName(std::string A){name = A;}
		void setMass(double m){mass = m;}
		void setCharge(double c){charge = c;}
		void setForcex(double xF){xForce = xF;}
		void setForcey(double yF){yForce = yF;}
		void setForcez(double zF){zForce = zF;}
		void setPositionx(double xF){xPosition = xF;}
		void setPositiony(double yF){yPosition = yF;}
		void setPositionz(double zF){zPosition = zF;}

		//Destructor
		~atom(){};
};

bool atom::operator==(atom test)const{

	if (type == test.getType()){
		return true;
	}
	else{
		return false;
	}
}

bool atom::operator!=(atom test)const{

	if (type != test.getType()){
		return true;
	}
	else{
		return false;
	}
}


#endif
