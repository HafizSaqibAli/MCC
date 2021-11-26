/* Shell Output Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 25/05/2017
**
** Output the shells of fragments or molecules.
**
*/

#ifndef OUTPUT_SHELLS
#define OUTPUT_SHELLS

void outputMoleculeShells(vector<molecule> &molecules, string outFile);
void outputMoleculeShells(vector<molecule> &molecules, string outFile, molecule *moleculeType);

void outputMoleculeShells(vector<molecule> &molecules, string outFile){
	outputMoleculeShells(molecules, outFile, nullptr);
}

void outputMoleculeShells(vector<molecule> &molecules, string outFile, molecule *moleculeType){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	sortMoleculesByShellSize(molecules);

	int line = 0;
	for (int i = 0; i < molecules.size(); i++){

		bool check = false;

		if(!moleculeType){
			check = true;
			for (int j = 0; j < molecules[i].getSize(); j++){
				for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
					out << molecules[i].getFragment(j)->getType(); //getAtom(k)->getName();
					out << " ";
				}
			}
		}
		else if(molecules[i].getType() == moleculeType->getType()){
			check = true;
		}
		
		if(check){
			for (int j = 0; j < molecules[i].getSize(); j++){
				for(int k = 0; k < molecules[i].getShell(j)->getSize(); k++){
					out << molecules[i].getShell(j)->getFragment(k)->getType();
				}
			}
			out << " " << molecules[i].getOccurrence() << endl;
		}
	}
}

#endif