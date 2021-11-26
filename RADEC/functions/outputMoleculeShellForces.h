/* Shell Output Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 25/05/2017
**
** Output the shells of fragments or molecules.
**
*/

#ifndef OUTPUT_MOLECULE_SHELL_FORCES
#define OUTPUT_MOLECULE_SHELL_FORCES

void outputMoleculeShellForces(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells, string outFile);
void outputMoleculeShellForces(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells, string outFile, molecule *moleculeType);

void outputMoleculeShellForces(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells, string outFile){
	outputMoleculeShellForces(molecules, moleculeShells, outFile, nullptr);
}

void outputMoleculeShellForces(vector<molecule> &molecules, vector<moleculeShell> &moleculeShells, string outFile, molecule *moleculeType){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	sortMoleculesByShellSize(molecules, moleculeShells);

	for (int i = 0; i < molecules.size(); i++){

		bool check = false;

		if (!moleculeType){
			check = true;
		}
		else if (molecules[i].getType() == moleculeType->getType()){
			check = true;
		}

		if (check){
			out << molecules[i].getType() << " ";

			for (int k = 0; k < moleculeShells[i].getSize(); k++){
				out << moleculeShells[i].getMolecule(k)->getType();
			}		
			
			double FAvg = molecules[i].getForceAverage() / ((double)molecules[i].getOccurrence());
			double FSq = molecules[i].getForceSquared() / ((double)molecules[i].getOccurrence());
			out << " " << molecules[i].getOccurrence() << " " << FAvg << " " << FSq << endl;
		}
	}
}

#endif