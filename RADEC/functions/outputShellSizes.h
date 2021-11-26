/* Force Parsing Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 25/05/2017
**
** Output the shell sizes of fragments or molecules.
**
*/

#ifndef OUTPUT_SHELL_SIZES
#define OUTPUT_SHELL_SIZES

void outputMoleculeShellSizes(vector<molecule> &molecules, string outFile);
void outputMoleculeShellSizes(vector<molecule> &molecules, string outFile, molecule *moleculeType);

void outputMoleculeShellSizes(vector<molecule> &molecules, string outFile){
	outputMoleculeShellSizes(molecules, outFile, nullptr);
}
void outputMoleculeShellSizes(vector<molecule> &molecules, string outFile, molecule *moleculeType){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	sortMoleculesByShellSize(molecules);

	int line = 0;
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
				out << molecules[i].getShell(j)->getSize() << " ";
			}
			out << " " << molecules[i].getOccurrence() << endl;
		}
	}
}

#endif
