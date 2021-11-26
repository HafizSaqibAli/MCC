/* Force Parsing Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 14/02/2014
**
** This collection of functions acts on vectors of the molecule class.
**
*/

#ifndef OUTPUT_MOLECULES
#define OUTPUT_MOLECULES

void outputMolecules(vector<molecule> &molecules, string outFile);
void outputMolecules(vector<molecule> &molecules, string outFile, molecule *moleculeType);
void outputAverageForces(vector<molecule> &molecules, string outFile);

void outputMolecules(vector<molecule> &molecules, string outFile){
	outputMolecules(molecules, outFile, nullptr);
}

void outputMolecules(vector<molecule> &molecules, string outFile, molecule *moleculeType){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		bool check = false;

		if(!moleculeType){
			check = true;
		}
		else if(molecules[i].getType() == moleculeType->getType()){
			check = true;
		}

		if(check){
			for (int j = 0; j < molecules[i].getSize(); j++){
				for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
					out << molecules[i].getFragment(j)->getAtom(k)->getName();
				}
				out << ", ";
			}
			out << endl;
		}
	}

}

void outputPositions(vector<molecule> &molecules, float *boxSize, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	out << fixed;
	out << endl;

	int line = 0;
	for (int i = 0; i < molecules.size(); i++){
		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				if(line == 10){out << endl; line = 0;}
				out << setw(8) << setprecision(3) << molecules[i].getFragment(j)->getAtom(k)->getPositionx();
				line++;
				if(line == 10){out << endl; line = 0;}
				out << setw(8) << setprecision(3) << molecules[i].getFragment(j)->getAtom(k)->getPositiony();
				line++;
				if(line == 10){out << endl; line = 0;}
				out << setw(8) << setprecision(3) << molecules[i].getFragment(j)->getAtom(k)->getPositionz();
				line++;
			}
		}
	}
	
	out << endl;
	out << setw(8) << setprecision(3) << 2*boxSize[0] << " " << 2*boxSize[1] << " " << 2*boxSize[2];
	out << scientific;
}

#endif
