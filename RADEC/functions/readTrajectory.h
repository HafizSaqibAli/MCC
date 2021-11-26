/* Coordinate Parsing Code, Dr Jonathan Higham
** Created 12/02/2014
** Modified 13/02/2014
**
** Opens an AMBER coordinate file specified by a string,
** at a specific frame specificed as number of bytes since the start of the file,
** and then reads xyz coordinate data into a vector of atom objects. 
**
*/

#ifndef READ_TRAJECTORY
#define READ_TRAJECTORY

using namespace std;

streampos initialFrameStart(string crdFile);
streampos frameRead(string crdFile, streampos frameStart, vector<molecule> &molecules, float *boxSize);

streampos initialFrameStart(string crdFile){
	
	streampos frameStart;
	ifstream moleculesIn;

	//Open solute topology file and check for errors
	moleculesIn.open(crdFile.c_str());
	if(!moleculesIn){
		cerr << "Cannot open the trajectory file.\n" << endl; 
		exit(0);
	}

	//Skip nameline
	string tempString;
	std::getline(moleculesIn, tempString);

	//Set stream position and return
	frameStart = moleculesIn.tellg();
	return(frameStart);
}

streampos frameRead(string crdFile, streampos frameStart, vector<atom> &atoms, float *boxSize){

	double x, y, z;
	ifstream moleculesIn;

	//Open solute topology file and check for errors
	moleculesIn.open(crdFile.c_str());
	if(!moleculesIn){
		cerr << "Cannot open an input file.\n" << endl; 
		exit(0);
	}
	moleculesIn.seekg(frameStart);

	//Read positions into the fragment vector
	for (int i = 0; i < atoms.size(); i++){
		moleculesIn >> x;
		moleculesIn >> y;
		moleculesIn >> z;
		atoms[i].setPositionx(x);
		atoms[i].setPositiony(y);
		atoms[i].setPositionz(z);
	}

	//Read in box size
	moleculesIn >> boxSize[0] >> boxSize[1] >> boxSize[2];
	//Move frameStart and return
	frameStart = moleculesIn.tellg();
	return(frameStart);
}

#endif
