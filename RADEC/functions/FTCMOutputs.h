/* Molecule functions, Dr Jonathan Higham
** Created 13/02/2014
** Modified 14/02/2014
**
** This collection of functions acts outputs the Force-Torque Covariance 
** matricies, eigenvalues and eigenvectors of the molecule class.
**
*/
#ifndef FCTM_OUTPUTS
#define FCTM_OUTPUTS

void outputFCMatrix(vector<molecule> &molecules, string outFile);
void outputTCMatrix(vector<molecule> &molecules, string outFile);
void outputFTCMatrix(vector<molecule> &molecules, string outFile);
void outputWMFTCMatrix(vector<molecule> &molecules, string outFile);
void outputEigenvalues(vector<molecule> &molecules, string outFile);
void outputWMEigenvalues(vector<molecule> &molecules, string outFile);
void outputEigenvectors(vector<molecule> &molecules, string outFile);
void outputWMEigenvectors(vector<molecule> &molecules, string outFile);

void outputFCMatrices(vector<molecule> &molecules, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				out << molecules[i].getFragment(j)->getAtom(k)->getName();
			}
			out << ", ";
		}
		out << endl;

		out << molecules[i].getFCMatrix() << endl << endl;
	}
}

void outputTCMatrices(vector<molecule> &molecules, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				out << molecules[i].getFragment(j)->getAtom(k)->getName();
			}
			out << ", ";
		}
		out << endl;

		out << molecules[i].getTCMatrix() << endl << endl;
	}
}


void outputFTCMatrices(vector<molecule> &molecules, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				out << molecules[i].getFragment(j)->getAtom(k)->getName();
			}
			out << ", ";
		}
		out << endl;

		out << molecules[i].getFTCMatrix() << endl << endl;
	}
}

void outputWMFTCMatrices(vector<molecule> &molecules, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				out << molecules[i].getFragment(j)->getAtom(k)->getName();
			}
			out << ", ";
		}
		out << endl;

		out << molecules[i].getWMFTCMatrix() << endl << endl;
	}
}


void outputFEigenvalues(vector<molecule> &molecules, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				out << molecules[i].getFragment(j)->getAtom(k)->getName();
			}
			out << ", ";
		}
		out << endl;

		Eigen::EigenSolver<Eigen::MatrixXd> ev = molecules[i].getFEigenvalues();
		out << ev.eigenvalues().real() << endl << endl;
	}
}

void outputTEigenvalues(vector<molecule> &molecules, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				out << molecules[i].getFragment(j)->getAtom(k)->getName();
			}
			out << ", ";
		}
		out << endl;

		Eigen::EigenSolver<Eigen::MatrixXd> ev = molecules[i].getTEigenvalues();
		out << ev.eigenvalues().real() << endl << endl;
	}
}

void outputFTEigenvalues(vector<molecule> &molecules, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				out << molecules[i].getFragment(j)->getAtom(k)->getName();
			}
			out << ", ";
		}
		out << endl;

		Eigen::EigenSolver<Eigen::MatrixXd> ev = molecules[i].getFTEigenvalues();
		out << ev.eigenvalues().real() << endl << endl;
	}
}

void outputWMEigenvalues(vector<molecule> &molecules, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				out << molecules[i].getFragment(j)->getAtom(k)->getName();
			}
			out << ", ";
		}
		out << endl;

		Eigen::EigenSolver<Eigen::MatrixXd> ev = molecules[i].getWMEigenvalues();
		out << ev.eigenvalues().real() << endl << endl;
	}
}

void outputFEigenvectors(vector<molecule> &molecules, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				out << molecules[i].getFragment(j)->getAtom(k)->getName();
			}
			out << ", ";
		}
		out << endl;

		Eigen::EigenSolver<Eigen::MatrixXd> ev = molecules[i].getFEigenvalues();
		out << std::fixed <<  std::setprecision(2) << ev.eigenvectors().real() << endl << endl;
	}
}

void outputTEigenvectors(vector<molecule> &molecules, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				out << molecules[i].getFragment(j)->getAtom(k)->getName();
			}
			out << ", ";
		}
		out << endl;

		Eigen::EigenSolver<Eigen::MatrixXd> ev = molecules[i].getTEigenvalues();
		out << std::fixed <<  std::setprecision(2) << ev.eigenvectors().real() << endl << endl;
	}
}

void outputFTEigenvectors(vector<molecule> &molecules, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				out << molecules[i].getFragment(j)->getAtom(k)->getName();
			}
			out << ", ";
		}
		out << endl;

		Eigen::EigenSolver<Eigen::MatrixXd> ev = molecules[i].getFTEigenvalues();
		out << std::fixed <<  std::setprecision(2) << ev.eigenvectors().real() << endl << endl;
	}
}

void outputWMEigenvectors(vector<molecule> &molecules, string outFile){

	ofstream out;
	//Open output file
	out.open(outFile.c_str());

	for (int i = 0; i < molecules.size(); i++){

		for (int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
				out << molecules[i].getFragment(j)->getAtom(k)->getName();
			}
			out << ", ";
		}
		out << endl;

		Eigen::EigenSolver<Eigen::MatrixXd> ev = molecules[i].getWMEigenvalues();
		out << std::fixed <<  std::setprecision(2) << ev.eigenvectors().real() << endl << endl;
	}
}

#endif
