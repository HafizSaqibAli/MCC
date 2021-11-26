/* Shell Output Code, Dr Jonathan Higham
** Created 13/02/2014
** Modified 25/05/2017
**
** Output the shells of fragments or molecules.
**
*/

#ifndef OUTPUT_SHELL_FORCES
#define OUTPUT_SHELL_FORCES

void outputFragmentShellForces(vector<fragment> &fragments, vector<shell> &fragmentShells, string outFile);
void outputFragmentShellForces(vector<fragment> &fragments, vector<shell> &fragmentShells, string outFile, fragment *fragmentType);

void outputFragmentShellForces(vector<fragment> &fragments, vector<shell> &fragmentShells, string outFile){
	outputFragmentShellForces(fragments, fragmentShells, outFile, nullptr);
}

void outputFragmentShellForces(vector<fragment> &fragments, vector<shell> &fragmentShells, string outFile, fragment *fragmentType){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	sortFragmentsByShellSize(fragments, fragmentShells);

	for (int i = 0; i < fragments.size(); i++){

		bool check = false;

		if (!fragmentType){
			check = true;
		}
		else if (fragments[i].getType() == fragmentType->getType()){
			check = true;
		}

		if (check){
			out << fragments[i].getType();
			out << " ";

			for (int j = 0; j < fragmentShells[i].getSize(); j++){
					out << fragmentShells[i].getFragment(j)->getType();
			}
			
			double FAvg = fragments[i].getForceAverage() / ((double)fragments[i].getOccurrence());
			double FSq = fragments[i].getForceSquared() / ((double)fragments[i].getOccurrence());
			out << " " << fragments[i].getOccurrence() << " " << FAvg << " " << FSq << endl;
		}
	}
}

#endif