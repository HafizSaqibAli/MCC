/* Force Parsing Code, Dr Jonathan Higham
** Created 12/02/2014
** Modified 13/02/2014
**
** Opens an AMBER coordinate file specified by a string,
** at a specific frame specificed as number of bytes since the start of the file,
** and then reads xyz force data into a vector of atom objects.
**
*/

#ifndef FORCE_PARSING
#define FORCE_PARSING

using namespace std;

streampos initialForceStart(string frcFile);
streampos forceRead(string frcFile, streampos frameStart, vector<atom> &atoms, float *boxSize);

streampos initialForceStart(string frcFile){
	
	streampos frameStart;
	ifstream ForcesIn;

	//Open input file and check for errors
	ForcesIn.open(frcFile.c_str());
	if(!ForcesIn){
		cerr << "Cannot open the forces file." << endl; 
		exit(0);
	}

	string tempString;
	std::getline(ForcesIn, tempString);
	frameStart = ForcesIn.tellg();

	ForcesIn.close();

	return(frameStart);
}

streampos forceRead(string frcFile, streampos frameStart, vector<atom> &atoms, float *boxSize){

	double x, y, z;
	ifstream ForcesIn;

	//Open solute topology file and check for errors
	ForcesIn.open(frcFile.c_str());
	if(!ForcesIn){
		cerr << "Cannot open an input file.\n" << endl; 
		exit(0);
	}

	//Seek to the start of the current frame
	ForcesIn.seekg(frameStart);

	//Read positions into the molecule vector
	for (int i = 0; i < atoms.size(); i++){
		ForcesIn >> x;
		ForcesIn >> y;
		ForcesIn >> z;
		atoms[i].setForcex(x);
		atoms[i].setForcey(y);
		atoms[i].setForcez(z);
	}

	//Move forceStart to the end of the current frame
	frameStart = ForcesIn.tellg();
	return(frameStart);
}

#endif
