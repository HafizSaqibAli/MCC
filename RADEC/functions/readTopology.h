/* Topology Parsing Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 15/06/2016
**
** This collection of functions acts on an AMBER topology 
** file to extract any relevant data about the each atom
**
*/

#ifndef READ_TOPOLOGY
#define READ_TOPOLOGY

using namespace std;

streampos findAtomStart(string topFile);
streampos findNameStart(string topFile);
streampos findMassStart(string topFile);
streampos findChargeStart(string topFile);
streampos findHBondStart(string topFile);
streampos findCoBondStart(string topFile);
streampos findACoeffStart(string topFile);
streampos findBCoeffStart(string topFile);
streampos findNonbondedParmStart(string topFile);

int getNumAtoms(string topFile);
vector<atom> getAtoms(string topFile);
void readNames(vector<atom> &atoms);
void readCharges(vector<atom> &atoms);
void readMasses(vector<atom> &atoms);
vector<int> getNonBondedParmIndex(string topFile);

streampos findAtomStart(string topFile){

	//Basic IO variables
	char A[3];
	bool flag = true, flag2 = false, flag3 = false;
	ifstream TopologyIn;

	//Open solute topology file and check for errors
	TopologyIn.open(topFile.c_str());
	if(!TopologyIn){
		cerr << "Cannot open the atom topology file: " << topFile << endl; 
		exit(0);
	}

	//Check for the start of the atom declarations
	while(flag){

		A[0] = A[1];
		A[1] = A[2];
		TopologyIn >> noskipws >> A[2];

		//Searches for %FLAG ATOM_NAME
		if(A[0] == 'E' && A[1] == '_' && A[2] == 'I'){ 
			flag2 = true;
		}
		//Finds the ) of the following format declaration
		if(flag2){
			if(A[2] == ')'){
				flag3 = true;
			}
		}
		//And finally finds the end of the format lines
		if(flag3){
			if(A[2] == '\n'){
				flag = false;
			}
		}
	}
	//At this point the file is now pointing to the start of the atom declaration array
	streampos atomStart = TopologyIn.tellg();

	//Close input file and return
	TopologyIn.close();
	return(atomStart);
}

streampos findMassStart(string topFile){

	//Basic IO variables
	char A[3];
	streampos massStart;
	bool flag = true, flag2 = false, flag3 = false;
	ifstream TopologyIn;

	//Open solute topology file and check for errors
	TopologyIn.open(topFile.c_str());
	if(!TopologyIn){
		cerr << "Cannot open the charge topolgy file." << topFile << endl; 
		exit(1);
	}

	//Check for the start of the atom declarations
	while(flag){
		A[0] = A[1];		
		A[1] = A[2];
		TopologyIn >> noskipws >> A[2];

		//Searches for %FLAG ATOM_NAME
		if(A[0] == 'M' && A[1] == 'A' && A[2] == 'S'){ 
			flag2 = true;
		}
		//Finds the ) of the following format declaration
		if(flag2){
			if(A[2] == ')'){
				flag3 = true;
			}
		}
		//And finally finds the end of the format lines
		if(flag3){
			if(A[2] == '\n'){
				flag = false;
			}
		}
	}
	//At this point the file is now pointing to the start of the atom declaration array
	massStart = TopologyIn.tellg();

	//Close input file and return
	TopologyIn.close();
	return(massStart);
}

streampos findNameStart(string topFile){

	//Basic IO variables
	char A[3];
	bool flag = true, flag2 = false, flag3 = false;
	ifstream TopologyIn;

	//Open solute topology file and check for errors
	TopologyIn.open(topFile.c_str());
	if(!TopologyIn){
		cerr << "Cannot open the atom topology file: " << topFile << endl; 
		exit(0);
	}

	//Check for the start of the atom declarations
	while(flag){

		A[0] = A[1];
		A[1] = A[2];
		TopologyIn >> noskipws >> A[2];

		//Searches for %FLAG ATOM_NAME
		if(A[0] == 'A' && A[1] == 'M' && A[2] == 'E'){ 
			flag2 = true;
		}
		//Finds the ) of the following format declaration
		if(flag2){
			if(A[2] == ')'){
				flag3 = true;
			}
		}
		//And finally finds the end of the format lines
		if(flag3){
			if(A[2] == '\n'){
				flag = false;
			}
		}
	}
	//At this point the file is now pointing to the start of the atom declaration array
	streampos atomStart = TopologyIn.tellg();

	//Close input file and return
	TopologyIn.close();
	return(atomStart);
}

streampos findChargeStart(string topFile){

	//Basic IO variables
	char A[3];
	streampos chargeStart;
	bool flag = true, flag2 = false, flag3 = false;
	ifstream TopologyIn;

	//Open solute topology file and check for errors
	TopologyIn.open(topFile.c_str());
	if(!TopologyIn){
		cerr << "Cannot open the charge topolgy file." << topFile << endl; 
		exit(1);
	}

	//Check for the start of the atom declarations
	while(flag){
		A[0] = A[1];		
		A[1] = A[2];
		TopologyIn >> noskipws >> A[2];

		//Searches for %FLAG ATOM_NAME
		if(A[0] == 'R' && A[1] == 'G' && A[2] == 'E'){ 
			flag2 = true;
		}
		//Finds the ) of the following format declaration
		if(flag2){
			if(A[2] == ')'){
				flag3 = true;
			}
		}
		//And finally finds the end of the format lines
		if(flag3){
			if(A[2] == '\n'){
				flag = false;
			}
		}
	}
	//At this point the file is now pointing to the start of the atom declaration array
	chargeStart = TopologyIn.tellg();

	//Close input file and return
	TopologyIn.close();
	return(chargeStart);
}

streampos findHBondStart(string topFile){

	//Basic IO variables
	char A[4];
	streampos bondStart;
	bool flag = true, flag2 = false, flag3 = false;
	ifstream TopologyIn;

	//Open solute topology file and check for errors
	TopologyIn.open(topFile.c_str());
	if(!TopologyIn){
		cerr << "Cannot open the HBond topology file." << topFile << endl; 
		exit(1);
	}

	//Check for the start of the atom declarations
	while(flag){
		A[0] = A[1];		
		A[1] = A[2];		
		A[2] = A[3];
		TopologyIn >> noskipws >> A[3];

		//Searches for %FLAG ATOM_NAME
		if(A[0] == 'D' && A[1] == 'S' && A[2] == '_' && A[3] == 'I'){ 
			flag2 = true;
		}
		//Finds the ) of the following format declaration
		if(flag2){
			if(A[3] == ')'){
				flag3 = true;
			}
		}
		//And finally finds the end of the format lines
		if(flag3){
			if(A[3] == '\n'){
				flag = false;
			}
		}
	}
	//At this point the file is now pointing to the start of the atom declaration array
	bondStart = TopologyIn.tellg();

	//Close input file and return
	TopologyIn.close();
	return(bondStart);
}

streampos findCoBondStart(string topFile){

	//Basic IO variables
	char A[4];
	streampos bondStart;
	bool flag = true, flag2 = false, flag3 = false;
	ifstream TopologyIn;

	//Open solute topology file and check for errors
	TopologyIn.open(topFile.c_str());
	if(!TopologyIn){
		cerr << "Cannot open the coBond topology file." << topFile << endl; 
		exit(1);
	}

	//Check for the start of the atom declarations
	while(flag){
		A[0] = A[1];		
		A[1] = A[2];		
		A[2] = A[3];
		TopologyIn >> noskipws >> A[3];

		//Searches for %FLAG ATOM_NAME
		if(A[0] == 'D' && A[1] == 'S' && A[2] == '_' && A[3] == 'W'){ 
			flag2 = true;
		}
		//Finds the ) of the following format declaration
		if(flag2){
			if(A[3] == ')'){
				flag3 = true;
			}
		}
		//And finally finds the end of the format lines
		if(flag3){
			if(A[3] == '\n'){
				flag = false;
			}
		}
	}
	//At this point the file is now pointing to the start of the atom declaration array
	bondStart = TopologyIn.tellg();

	//Close input file and return
	TopologyIn.close();
	return(bondStart);
}

streampos findACoeffStart(string topFile){

	//Basic IO variables
	char A[3];
	streampos ACoeffStart;
	bool flag = true, flag2 = false, flag3 = false;
	ifstream TopologyIn;

	//Open solute topology file and check for errors
	TopologyIn.open(topFile.c_str());
	if(!TopologyIn){
		cerr << "Cannot open the charge topolgy file." << topFile << endl; 
		exit(1);
	}

	//Check for the start of the atom declarations
	while(flag){
		A[0] = A[1];		
		A[1] = A[2];
		TopologyIn >> noskipws >> A[2];

		//Searches for %FLAG ATOM_NAME
		if(A[0] == 'S' && A[1] == '_' && A[2] == 'A'){ 
			flag2 = true;
		}
		//Finds the ) of the following format declaration
		if(flag2){
			if(A[2] == ')'){
				flag3 = true;
			}
		}
		//And finally finds the end of the format lines
		if(flag3){
			if(A[2] == '\n'){
				flag = false;
			}
		}
	}
	//At this point the file is now pointing to the start of the LJ_ACoeff declaration array
	ACoeffStart = TopologyIn.tellg();

	//Close input file and return
	TopologyIn.close();
	return(ACoeffStart);
}

streampos findBCoeffStart(string topFile){

	//Basic IO variables
	char A[3];
	streampos BCoeffStart;
	bool flag = true, flag2 = false, flag3 = false;
	ifstream TopologyIn;

	//Open solute topology file and check for errors
	TopologyIn.open(topFile.c_str());
	if(!TopologyIn){
		cerr << "Cannot open the charge topolgy file." << topFile << endl; 
		exit(1);
	}

	//Check for the start of the atom declarations
	while(flag){
		A[0] = A[1];		
		A[1] = A[2];
		TopologyIn >> noskipws >> A[2];

		//Searches for %FLAG ATOM_NAME
		if(A[0] == 'S' && A[1] == '_' && A[2] == 'B'){ 
			flag2 = true;
		}
		//Finds the ) of the following format declaration
		if(flag2){
			if(A[2] == ')'){
				flag3 = true;
			}
		}
		//And finally finds the end of the format lines
		if(flag3){
			if(A[2] == '\n'){
				flag = false;
			}
		}
	}
	//At this point the file is now pointing to the start of the LJ_BCoeff declaration array
	BCoeffStart = TopologyIn.tellg();

	//Close input file and return
	TopologyIn.close();
	return(BCoeffStart);
}

streampos findNonbondedParmStart(string topFile){

	//Basic IO variables
	char A[3];
	streampos NBPStart;
	bool flag = true, flag2 = false, flag3 = false;
	ifstream TopologyIn;

	//Open solute topology file and check for errors
	TopologyIn.open(topFile.c_str());
	if(!TopologyIn){
		cerr << "Cannot open the charge topolgy file." << topFile << endl; 
		exit(1);
	}

	//Check for the start of the atom declarations
	while(flag){
		A[0] = A[1];		
		A[1] = A[2];
		TopologyIn >> noskipws >> A[2];

		//Searches for %FLAG ATOM_NAME
		if(A[0] == 'M' && A[1] == '_' && A[2] == 'I'){ 
			flag2 = true;
		}
		//Finds the ) of the following format declaration
		if(flag2){
			if(A[2] == ')'){
				flag3 = true;
			}
		}
		//And finally finds the end of the format lines
		if(flag3){
			if(A[2] == '\n'){
				flag = false;
			}
		}
	}
	//At this point the file is now pointing to the start of the LJ_BCoeff declaration array
	NBPStart = TopologyIn.tellg();

	//Close input file and return
	TopologyIn.close();
	return(NBPStart);
}

int getNumAtoms(string topFile){

	//Number of atoms
	int nAtoms;
	//Basic IO variables
	char A[3];
	bool flag = true, flag2 = false, flag3 = false;
	ifstream TopologyIn;


	//Open solute topology file and check for errors
	TopologyIn.open(topFile.c_str());
	if(!TopologyIn){
		cerr << "Cannot open the atom topology file: " << topFile << endl; 
		exit(0);
	}

	//Check for the start of the atom declarations
	while(flag){

		A[0] = A[1];
		A[1] = A[2];
		TopologyIn >> noskipws >> A[2];

		//Searches for %FLAG ATOM_NAME
		if(A[0] == ' ' && A[1] == 'P' && A[2] == 'O'){ 
			flag2 = true;
		}
		//Finds the ) of the following format declaration
		if(flag2){
			if(A[2] == ')'){
				flag3 = true;
			}
		}
		//And finally finds the end of the format lines
		if(flag3){
			if(A[2] == '\n'){
				flag = false;
			}
		}
	}

	string atomPointers;
	getline(TopologyIn, atomPointers);
	std::istringstream atomPointerString(atomPointers);

	atomPointerString >> nAtoms;

	return nAtoms;
}

vector<atom> getAtoms(string topFile){

	int t;
	bool flag = true;
	//Find the number of atoms
	int nAtoms = getNumAtoms(topFile);

	vector<atom> atoms;

	//Open solute topology file
	ifstream TopologyIn;
	TopologyIn.open(topFile.c_str());
	//Get to the start of the atom declaration
	streampos atomStart = findAtomStart(topFile);
	TopologyIn.seekg(atomStart);

	for(int i = 0; i < nAtoms; i++){

		//Read in an integer
		TopologyIn >> t;

		//Set atom type into an array
		atom tempAtom;
		tempAtom.setType(t);
		tempAtom.setIndex(i);
		atoms.push_back(tempAtom);
	}
	return atoms;
}

void readNames(string topFile, vector<atom> &atoms){

	//Basic IO variables
	std::string temp;
	int i = 0, j = 0;
	bool flag = true;
	ifstream TopologyIn;
	streampos nameStart = findNameStart(topFile);

	//Open solute topology file
	TopologyIn.open(topFile.c_str());

	//Find the start of the charge declaration
	TopologyIn.seekg(nameStart);
		
	while(TopologyIn){
		
		string str;
		getline(TopologyIn, str);
		if(str[0] == '%'){
			break;
		}
		stringstream ss(str);

		while(ss >> temp){
			temp.erase(std::remove_if(temp.begin(), temp.end(), ::isdigit), temp.end());
			if(temp[0] == 'H' && isupper(temp[1])){
				temp = "H";
			}
			atoms[i].setName(temp);
			i++;
			if(i == atoms.size()){
				break;
			}
		}
	}
	//Close input file
	TopologyIn.close();
}

void readMasses(string topFile, vector<atom> &atoms){

	//Basic IO variables
	double temp;
	int i = 0, j = 0;
	bool flag = true;
	ifstream TopologyIn;
	streampos massStart = findMassStart(topFile);

	//Open solute topology file
	TopologyIn.open(topFile.c_str());

	//Find the start of the charge declaration
	TopologyIn.seekg(massStart);
		
	while(TopologyIn){
		
		string str;
		getline(TopologyIn, str);
		if(str[0] == '%'){
			break;
		}
		stringstream ss(str);

		while(ss >> temp){
			atoms[i].setMass(temp);
			i++;
			if(i == atoms.size()){
				break;
			}
		}
	}

	//Close input file
	TopologyIn.close();
}

void readCharges(string topFile, vector<atom> &atoms){

	//Basic IO variables
	double temp;
	int i = 0, j = 0;
	bool flag = true;
	ifstream TopologyIn;
	streampos chargeStart = findChargeStart(topFile);

	//Open solute topology file
	TopologyIn.open(topFile.c_str());

	//Find the start of the charge declaration
	TopologyIn.seekg(chargeStart);
		
	while(TopologyIn){
		
		string str;
		getline(TopologyIn, str);
		if(str[0] == '%'){
			break;
		}
		stringstream ss(str);

		while(ss >> temp){
			atoms[i].setCharge(temp);
			i++;
			if(i == atoms.size()){
				break;
			}
		}
	}

	//Close input file
	TopologyIn.close();
}

vector<int> getNonBondedParmIndex(string topFile){

	//Basic IO variables
	int temp;
	ifstream TopologyIn;
	streampos NBPStart = findNonbondedParmStart(topFile);

	//Open solute topology file
	TopologyIn.open(topFile.c_str());

	//Find the start of the charge declaration
	TopologyIn.seekg(NBPStart);
	
	vector<int> NBPIndex;

	while(TopologyIn){
		
		string str;
		getline(TopologyIn, str);
		if(str[0] == '%'){
			break;
		}
		stringstream ss(str);
		
		while(ss >> temp){
			NBPIndex.push_back(temp);
		}
	}

	//Close input file
	TopologyIn.close();

	return NBPIndex;
}

#endif
