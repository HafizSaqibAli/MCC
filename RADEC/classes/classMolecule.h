#ifndef molecule_CLASS
#define molecule_CLASS

#include "removeEigenRow.h"
#include "removeEigenCol.h"

class molecule{
	private:
		//Basic declarations for each fragment
		int numAtoms;
		int moleculeSize;
		int moleculeType;
		double moleculeMass;
		std::string moleculeName;
		double moleculeForcex;
		double moleculeForcey;
		double moleculeForcez;
		double moleculeTorquex;
		double moleculeTorquey;
		double moleculeTorquez;

		//An array of pointers to each fragment in the molecule
		std::vector<fragment*> fragments;
		//An array of the shells of each fragment
		std::vector<shell> fragmentShells;

		//Force Covariance Matrix
		Eigen::MatrixXd forceCovariance;
		Eigen::EigenSolver<Eigen::MatrixXd> forceEigenvalues;
		Eigen::MatrixXd torqueCovariance;
		Eigen::EigenSolver<Eigen::MatrixXd> torqueEigenvalues;

		Eigen::MatrixXd forceTorqueCovariance;
		Eigen::EigenSolver<Eigen::MatrixXd> eigenvalues;
		Eigen::MatrixXd wholeMoleculeForceTorqueCovariance;
		Eigen::EigenSolver<Eigen::MatrixXd> wholeMoleculeEigenvalues;

		//Useful bonus variables for averaging
		int occurrence;
		double moleculeForceAverage;
		double moleculeForceSquared;

	public:
		//Default constructor
		molecule(){ numAtoms = 0;  moleculeSize = 0; moleculeMass = 0; moleculeName = "";
			    occurrence = 0; moleculeForceAverage = 0; moleculeForceSquared = 0; }
		//Overloaded member operators
		bool operator==(molecule test) const;
		bool operator!=(molecule test) const;

		//FCM member functions
		void createFTCMatrix();
		void updateFTCMatrix();
		void addFTCMatrix(molecule tempMolecule);
		void divideFTCMatrix(double d);
		void multiplyFTCMatrix(double d);
		void removeZerosFTCMatrix();
		void massweightFTCMatrix();
		void removeNonBondedInteractionsFTCMatrix();
		void solveFTCMatrix();

		//Function to add an atom
		void addFragment(fragment* frag){
			fragments.push_back(frag);
			moleculeName.append(frag->getName());
			moleculeSize = moleculeSize + 1;
			moleculeMass = moleculeMass + frag->getMass();
			numAtoms = numAtoms + frag->getSize();

			shell tempShell;
			fragmentShells.push_back(tempShell);
		}
		//Function to add molecule
		void addMolecule(molecule mol){
			for(int i = 0; i < mol.getSize(); i++){
				fragments.push_back(mol.getFragment(i));
				fragmentShells.push_back(mol.returnFragmentShell(i));

				moleculeSize = moleculeSize + 1;
				moleculeMass = moleculeMass + mol.getFragment(i)->getMass();
				moleculeName.append(mol.getFragment(i)->getName());
				numAtoms = numAtoms + mol.getFragment(i)->getSize();
			}
		}

		void delFragment(int i){
			numAtoms = numAtoms - fragments[i]->getSize();
			moleculeMass = moleculeMass - fragments[i]->getMass();
			moleculeSize = moleculeSize - 1;

			fragments.erase(fragments.begin() + i);
			fragmentShells.erase(fragmentShells.begin() + i);
		}

		fragment *getFragment(int i){return fragments[i];}
		fragment returnFragment(int i){ return *fragments[i]; }
		shell* getFragmentShell(int i){ return &fragmentShells[i]; }
		shell returnFragmentShell(int i){ return fragmentShells[i]; }

		Eigen::MatrixXd getFCMatrix(){ return forceCovariance; }
		Eigen::MatrixXd getTCMatrix(){ return torqueCovariance; }
		Eigen::MatrixXd getFTCMatrix(){ return forceTorqueCovariance; }
		Eigen::EigenSolver<Eigen::MatrixXd> getFEigenvalues(){ return forceEigenvalues; }
		Eigen::EigenSolver<Eigen::MatrixXd> getTEigenvalues(){ return torqueEigenvalues; }
		Eigen::EigenSolver<Eigen::MatrixXd> getFTEigenvalues(){ return eigenvalues; }

		Eigen::MatrixXd getWMFTCMatrix(){ return wholeMoleculeForceTorqueCovariance; }
		Eigen::EigenSolver<Eigen::MatrixXd> getWMEigenvalues(){ return wholeMoleculeEigenvalues; }

		//Basic molecule get functions
		int getNumAtoms(){ return numAtoms; }
		int getSize(){ return moleculeSize; }
		int getType(){return moleculeType;}
		double getMass(){ return moleculeMass; }
		std::string getName(){ return moleculeName; }
		double getForcex(){return moleculeForcex;}
		double getForcey(){return moleculeForcey;}
		double getForcez(){return moleculeForcez;}
		double getTorquex(){return moleculeTorquex;}
		double getTorquey(){return moleculeTorquey;}
		double getTorquez(){return moleculeTorquez;}
		double getForceAverage(){return moleculeForceAverage;}
		double getForceSquared(){return moleculeForceSquared;}

		int getOccurrence(){ return occurrence; }
		double getFTCMElement(int i, int j){ return forceTorqueCovariance(i, j);}
		double getWMFTCMElement(int i, int j){ return wholeMoleculeForceTorqueCovariance(i, j);}

		//Basic molecule set functions
		void setType(int t){ moleculeType = t;}
		void setMass(double m){moleculeMass = m;}
		void setForcex(double f){moleculeForcex = f;}
		void setForcey(double f){moleculeForcey = f;}
		void setForcez(double f){moleculeForcez = f;}
		void setTorquex(double t){moleculeTorquex = t;}
		void setTorquey(double t){moleculeTorquey = t;}
		void setTorquez(double t){moleculeTorquez = t;}
		void setForceAverage(double f){moleculeForceAverage = f;}
		void setForceSquared(double f){moleculeForceSquared = f;}

		//Functions to increment properties
		void addMass(double m){ moleculeMass = moleculeMass + m; }
		void addForcex(double f){moleculeForcex = moleculeForcex + f;}
		void addForcey(double f){moleculeForcey = moleculeForcey + f;}
		void addForcez(double f){moleculeForcez = moleculeForcez + f;}
		void addTorquex(double t){moleculeTorquex = moleculeTorquex + t;}
		void addTorquey(double t){moleculeTorquey = moleculeTorquey + t;}
		void addTorquez(double t){moleculeTorquez = moleculeTorquez + t;}
		void addForceAverage(double f){moleculeForceAverage = moleculeForceAverage + f;}
		void addForceSquared(double f){moleculeForceSquared = moleculeForceSquared + f;}
		void incrementOccurence(){ occurrence = occurrence + 1; }

		double getCOMPositionx(){
			double ComX = 0;
			double massNorm = 1 / moleculeMass;

			for(int j = 0; j < fragments.size(); j++){

				for (int k = 0; k < fragments[j]->getSize(); k++){//Loop atoms
					double massFrac = fragments[j]->getAtom(k)->getMass()*massNorm;
					ComX = ComX + fragments[j]->getAtom(k)->getPositionx()*massFrac;
				}
			}
			return ComX;
		}
		double getCOMPositiony(){
			double ComY = 0;
			double massNorm = 1 / moleculeMass;

			for(int j = 0; j < fragments.size(); j++){

				for (int k = 0; k < fragments[j]->getSize(); k++){//Loop atoms
					double massFrac = fragments[j]->getAtom(k)->getMass()*massNorm;
					ComY = ComY + fragments[j]->getAtom(k)->getPositiony()*massFrac;
				}
			}
			return ComY;
		}
		double getCOMPositionz(){
			double ComZ = 0;
			double massNorm = 1 / moleculeMass;

			for(int j = 0; j < fragments.size(); j++){

				for (int k = 0; k < fragments[j]->getSize(); k++){//Loop atoms
					double massFrac = fragments[j]->getAtom(k)->getMass()*massNorm;
					ComZ = ComZ + fragments[j]->getAtom(k)->getPositionz()*massFrac;
				}
			}
			return ComZ;
		}
};

bool molecule::operator==(molecule test)const{

	bool isTheSame = true;
	if (moleculeSize != test.getSize()){
		isTheSame = false;
	}

	for (int i = 0; i < moleculeSize && isTheSame; i++){
		if (*fragments[i] != test.returnFragment(i)){
			isTheSame = false;
		}
	}

	return isTheSame;
}

bool molecule::operator!=(molecule test)const{

	bool isTheSame = true;
	if (moleculeSize != test.getSize()){
		isTheSame = false;
	}

	for (int i = 0; i < moleculeSize && isTheSame; i++){
		if (*fragments[i] != test.returnFragment(i)){
			isTheSame = false;
		}
	}

	return !isTheSame;
}

//Eigen functions for force-torque covariance matrix
void molecule::createFTCMatrix(){
	forceTorqueCovariance.resize(fragments.size() * 6, fragments.size() * 6);
	for (int i = 0; i < fragments.size()*6; i++){
		for (int j = 0; j < fragments.size()*6; j++){
			forceTorqueCovariance(i, j) = 0;
		}
	}
	forceCovariance.resize(fragments.size() * 3, fragments.size() * 3);
	for (int i = 0; i < fragments.size()*3; i++){
		for (int j = 0; j < fragments.size()*3; j++){
			forceCovariance(i, j) = 0;
		}
	}
	torqueCovariance.resize(fragments.size() * 3, fragments.size() * 3);
	for (int i = 0; i < fragments.size()*3; i++){
		for (int j = 0; j < fragments.size()*3; j++){
			torqueCovariance(i, j) = 0;
		}
	}
	wholeMoleculeForceTorqueCovariance.resize(6, 6);
	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 6; j++){
			wholeMoleculeForceTorqueCovariance(i, j) = 0;
		}
	}
}

void molecule::updateFTCMatrix(){

	//Loop over all fragments for i
	int i = 0;
	for (int i1 = 0; i1 < moleculeSize; i1++){

		std::vector <double> iforceTorque;
		iforceTorque.push_back(fragments[i1]->getForcex());
		iforceTorque.push_back(fragments[i1]->getForcey());
		iforceTorque.push_back(fragments[i1]->getForcez());
		iforceTorque.push_back(fragments[i1]->getTorquex());
		iforceTorque.push_back(fragments[i1]->getTorquey());
		iforceTorque.push_back(fragments[i1]->getTorquez());

		//Loop over all fragments for j
		int j = 0;
		for (int j1 = 0; j1 < moleculeSize; j1++){

			std::vector <double> jforceTorque;
			jforceTorque.push_back(fragments[j1]->getForcex());
			jforceTorque.push_back(fragments[j1]->getForcey());
			jforceTorque.push_back(fragments[j1]->getForcez());
			jforceTorque.push_back(fragments[j1]->getTorquex());
			jforceTorque.push_back(fragments[j1]->getTorquey());
			jforceTorque.push_back(fragments[j1]->getTorquez());

			for(int FCi = 0; FCi < 6; FCi++){
				for(int FCj = 0; FCj < 6; FCj++){

					forceTorqueCovariance(i + FCi, j + FCj) = forceTorqueCovariance(i + FCi, j + FCj)
						+ (iforceTorque[FCi] * jforceTorque[FCj]);
				}
			}

			j = j + 6;
		}
		i = i + 6;
	}

	std::vector <double> forceTorque(6,0.0);

	forceTorque[0] = moleculeForcex;
	forceTorque[1] = moleculeForcey;
	forceTorque[2] = moleculeForcez;
	forceTorque[3] = moleculeTorquex;
	forceTorque[4] = moleculeTorquey;
	forceTorque[5] = moleculeTorquez;

	for(int i1 = 0; i1 < 6; i1++){
		for(int j1 = 0; j1 < 6; j1++){
			wholeMoleculeForceTorqueCovariance(i1, j1) = wholeMoleculeForceTorqueCovariance(i1, j1) + forceTorque[i1]*forceTorque[j1];
		}
	}
}

void molecule::addFTCMatrix(molecule tempMolecule){

	for (int i = 0; i < fragments.size()*6; i++){
		for (int j = 0; j < fragments.size()*6; j++){
			forceTorqueCovariance(i, j) = forceTorqueCovariance(i, j) + tempMolecule.getFTCMElement(i, j);
		}
	}
	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 6; j++){
			wholeMoleculeForceTorqueCovariance(i, j) = wholeMoleculeForceTorqueCovariance(i, j) + tempMolecule.wholeMoleculeForceTorqueCovariance(i, j);
		}
	}
}

void molecule::divideFTCMatrix(double d){

	for (int i = 0; i < fragments.size()*6; i++){
		for (int j = 0; j < fragments.size()*6; j++){
			forceTorqueCovariance(i, j) = forceTorqueCovariance(i, j) / d;
		}
	}
	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 6; j++){
			wholeMoleculeForceTorqueCovariance(i, j) = wholeMoleculeForceTorqueCovariance(i, j) / d;
		}
	}
}

void molecule::multiplyFTCMatrix(double d){

	for (int i = 0; i < fragments.size()*6; i++){
		for (int j = 0; j < fragments.size()*6; j++){
			forceTorqueCovariance(i, j) = forceTorqueCovariance(i, j) * d;
		}
	}
	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 6; j++){
			wholeMoleculeForceTorqueCovariance(i, j) = wholeMoleculeForceTorqueCovariance(i, j) * d;
		}
	}
}

void molecule::solveFTCMatrix(){


	//Loop over all fragments for i
	int i = 0;
	for (int i1 = 0; i1 < moleculeSize; i1++){

		//Loop over all fragments for j
		int j = 0;
		for (int j1 = 0; j1 < moleculeSize; j1++){

			for(int FCi = 0; FCi < 3; FCi++){
				for(int FCj = 0; FCj < 3; FCj++){

					forceCovariance(i + FCi, j + FCj) = forceTorqueCovariance((2*i) + FCi, (2*j) + FCj);
					torqueCovariance(i + FCi, j + FCj) = forceTorqueCovariance((2*i) + FCi + 3, (2*j) + FCj + 3);
				}
			}

			j = j + 3;
		}
		i = i + 3;
	}

	for (int i = 0; i < torqueCovariance.rows(); i++){
		bool zeros = true;
		for (int j = 0; j < torqueCovariance.cols(); j++){
			if(torqueCovariance(i, j) != 0){
				zeros = false;
			}
		}
		if(zeros){
			removeRow(torqueCovariance, i);
			i--;
		}
	}

	for (int j = 0; j < torqueCovariance.cols(); j++){
		bool zeros = true;
		for (int i = 0; i < torqueCovariance.rows(); i++){
			if(torqueCovariance(i, j) != 0){
				zeros = false;
			}
		}
		if(zeros){
			removeCol(torqueCovariance, j);
			j--;
		}
	}

	for (int i = 0; i < forceCovariance.rows(); i++){
		bool zeros = true;
		for (int j = 0; j < forceCovariance.cols(); j++){
			if(forceCovariance(i, j) != 0){
				zeros = false;
			}
		}
		if(zeros){
			removeRow(forceCovariance, i);
			i--;
		}
	}

	for (int j = 0; j < forceCovariance.cols(); j++){
		bool zeros = true;
		for (int i = 0; i < forceCovariance.rows(); i++){
			if(forceCovariance(i, j) != 0){
				zeros = false;
			}
		}
		if(zeros){
			removeCol(forceCovariance, j);
			j--;
		}
	}

	for (int i = 0; i < forceTorqueCovariance.rows(); i++){
		bool zeros = true;
		for (int j = 0; j < forceTorqueCovariance.cols(); j++){
			if(forceTorqueCovariance(i, j) != 0){
				zeros = false;
			}
		}
		if(zeros){
			removeRow(forceTorqueCovariance, i);
			i--;
		}
	}

	for (int j = 0; j < forceTorqueCovariance.cols(); j++){
		bool zeros = true;
		for (int i = 0; i < forceTorqueCovariance.rows(); i++){
			if(forceTorqueCovariance(i, j) != 0){
				zeros = false;
			}
		}
		if(zeros){
			removeCol(forceTorqueCovariance, j);
			j--;
		}
	}

	for (int i = 0; i < wholeMoleculeForceTorqueCovariance.rows(); i++){
		bool zeros = true;
		for (int j = 0; j < wholeMoleculeForceTorqueCovariance.cols(); j++){
			if(wholeMoleculeForceTorqueCovariance(i, j) != 0){
				zeros = false;
			}
		}
		if(zeros){
			removeRow(wholeMoleculeForceTorqueCovariance, i);
			i--;
		}
	}

	for (int j = 0; j < wholeMoleculeForceTorqueCovariance.cols(); j++){
		bool zeros = true;
		for (int i = 0; i < wholeMoleculeForceTorqueCovariance.rows(); i++){
			if(wholeMoleculeForceTorqueCovariance(i, j) != 0){
				zeros = false;
			}
		}
		if(zeros){
			removeCol(wholeMoleculeForceTorqueCovariance, j);
			j--;
		}
	}

	if(forceTorqueCovariance.rows() > 0 && forceTorqueCovariance.rows() > 0){
		eigenvalues.compute(forceTorqueCovariance, true);
	}
	else{
		Eigen::MatrixXd tempMatrix = Eigen::MatrixXd::Ones(1,1);
		eigenvalues.compute(tempMatrix, true);
	}
	if(forceCovariance.rows() > 0 && forceCovariance.rows() > 0){
		forceEigenvalues.compute(forceCovariance, true);
	}
	else{
		Eigen::MatrixXd tempMatrix = Eigen::MatrixXd::Ones(1,1);
		forceEigenvalues.compute(tempMatrix, true);
	}
	if(torqueCovariance.rows() > 0 && torqueCovariance.rows() > 0){
		torqueEigenvalues.compute(torqueCovariance, true);
	}
	else{
		Eigen::MatrixXd tempMatrix = Eigen::MatrixXd::Ones(1,1);
		torqueEigenvalues.compute(tempMatrix, true);
	}
	if(wholeMoleculeForceTorqueCovariance.rows() > 0 && wholeMoleculeForceTorqueCovariance.rows() > 0){
		wholeMoleculeEigenvalues.compute(wholeMoleculeForceTorqueCovariance, true);
	}
	else{
		Eigen::MatrixXd tempMatrix = Eigen::MatrixXd::Ones(1,1);
		wholeMoleculeEigenvalues.compute(tempMatrix, true);
	}
}

void molecule::removeNonBondedInteractionsFTCMatrix(){

	//Loop over all fragments for i
	int i = 0;
	for (int i1 = 0; i1 < moleculeSize; i1++){

		if(fragments[i1]->getNumBondedFragments() > 0){
			//Loop over all fragments for j
			int j = 0;
			for (int j1 = 0; j1 < moleculeSize; j1++){

				bool inBondedList = false;
				for(int k = 0; k < fragments[i1]->getNumBondedFragments(); k++){
					if(fragments[i1]->getBondedFragment(k) == j1){
						inBondedList = true;
					}
				}

				if(inBondedList == false){
					for(int FCi = 0; FCi < 6; FCi++){
						for(int FCj = 0; FCj < 6; FCj++){

							forceTorqueCovariance(i + FCi, j + FCj) = 0;
						}
					}
				}
				j = j + 6;
			}
			i = i + 6;
		}
	}
}



#endif



/*bool shellEqualTo(molecule tempMolecule){

	bool isTheSame = true;
	if (moleculeSize != test.getSize()){
		isTheSame = false;
	}

	for (int i = 0; i < moleculeSize && isTheSame; i++){
		if (*fragments[i] != test.returnFragment(i)){
			isTheSame = false;
		}
	}

	return !isTheSame;

	return true;
}
*/