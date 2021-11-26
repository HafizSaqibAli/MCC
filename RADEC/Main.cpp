#include<cmath>
#include<cctype>
#include<vector>
#include<cstdlib>
#include<string>
#include<fstream>
#include<sstream>
#include<iostream>
#include<iomanip>

#include"eigen/Dense"

#include"classes/classAtom.h"
#include"classes/classFragment.h"
#include"classes/classShell.h"
#include"classes/classMolecule.h"
#include"classes/classMoleculeShell.h"

#include"functions.h"

using namespace std;

int main(int argc, char *argv[]){

	cout << "Reading data from topology file." << endl;

	//Number of frames
	const int tMax = atoi(argv[1]);
	const int freq = atoi(argv[2]);

	//Input filenames
	string crdInFile = "trajectory.crd";
	string frcInFile = "forces.frc";
	string topInFile = "molecules.top";

	//Output filenames
	string outFileShells = "Shells.txt";
	string outFileWMShells = "WMShells.txt";
	string outFileFC = "FCMatrix.txt";
	string outFileTC = "TCMatrix.txt";
	string outFileFTC = "FTCMatrix.txt";
	string outFileWMFC = "WMFCMatrix.txt";
	string outFileFEVal = "FEigenValues.txt";
	string outFileTEVal = "TEigenValues.txt";
	string outFileFTEVal = "FTEigenValues.txt";
	string outFileWMEVal = "WMEigenValues.txt";
	string outFileFEVec = "FEigenVectors.txt";
	string outFileTEVec = "TEigenVectors.txt";
	string outFileFTEVec = "FTEigenVectors.txt";
	string outFileWMEVec = "WMEigenVectors.txt";
	string outFilePositions = "positions.crd";

	//Read a list of atoms and their properties
	vector<atom> atoms;
	atoms = getAtoms(topInFile);
	readNames(topInFile, atoms);
	readMasses(topInFile, atoms);
	readCharges(topInFile, atoms);

	//Read add the atoms to fragments
	vector<fragment> fragments;
	createFragments(topInFile, atoms, fragments);
	setFragmentTypes(fragments);

	//Create another vector with one of each fragment type for averaging etc.
	vector<fragment> fragmentTypes;
	fragmentTypes = getFragmentTypes(fragments);

	//Finally, compose fragments into molecules
	vector<molecule> molecules;
	createMolecules(topInFile, fragments, molecules);
	cout << "Molecules created from topology file." << endl;

	setMoleculeTypes(molecules);
	vector<moleculeShell> moleculeShells;
	createMoleculeShells(molecules, moleculeShells);
	cout << "Molecule shells initialised." << endl;

	//Create another vector with one of each molecule type for averaging etc.
	vector<molecule> moleculeTypes;
	moleculeTypes = getMoleculeTypes(molecules);

	//A vector separating fragments based on their shell type
	vector<fragment> RADFragmentTypes;
	vector<shell> RADFragmentShellTypes;

	//A vector separating molecules based on their shell type
	vector<molecule> RADMoleculeTypes;
	vector<moleculeShell> RADMoleculeShellTypes;

	//Read the first frame and set time
	int t = 1;
	float boxSize[3];

	cout << "Commencing analysis of force and coordinate data." << endl;
	streampos frameStart = initialFrameStart(crdInFile);
	streampos forceStart = initialForceStart(frcInFile);
	frameStart = frameRead(crdInFile, frameStart, atoms, boxSize);
	forceStart = forceRead(frcInFile, forceStart, atoms, boxSize);

	int dfreq = freq;
	int totalFramesRun = 0;

	//Loop over the trajectory to calculate cosTheta shells
	while ((t < tMax + 1) && (frameStart != -1)){

		if(dfreq == freq){

			//Calculate RAD molecule shells and update
			RADMoleculeShells(molecules, moleculeShells, boxSize, nullptr, nullptr);
			localiseForces(molecules);
			updateFTCMatrices(molecules);
			updateMoleculeShellTypeForces(molecules, moleculeShells, RADMoleculeTypes, RADMoleculeShellTypes, nullptr);

			//calculate RAD fragment shells
			RADFragmentShells(molecules, boxSize, nullptr, nullptr);
			updateFragmentShellTypeForces(molecules, RADFragmentTypes, RADFragmentShellTypes, nullptr);

			//Print out a progress message
			cout << "shell calculation, frame " << t << "/" << tMax << " complete." << endl;

			//Update time and frame position
			frameStart = frameRead(crdInFile, frameStart, atoms, boxSize);
			forceStart = forceRead(frcInFile, forceStart, atoms, boxSize);

			dfreq = 1;
			totalFramesRun = totalFramesRun + 1;
		}
		else{	
			//Update time and frame position
			frameStart = frameRead(crdInFile, frameStart, atoms, boxSize);
			forceStart = forceRead(frcInFile, forceStart, atoms, boxSize);
			cout << "shell calculation, frame " << t << "/" << tMax << " skipped." << endl;
			dfreq=dfreq+1;
		}
		t = t + 1;
	}
	//Print out a message if tmax is larger than the input file
	if(t < tMax){
		cout << "Reached end of trajectory file, skipping remaining steps" << endl;
		cout << "Note: Localised positions will not be output" << endl;
	}

	cout << "Calculations complete, outputting results. Total frames used = " << totalFramesRun << endl;

	//Output RAD Shells
	outputFragmentShellForces(RADFragmentTypes, RADFragmentShellTypes, outFileShells, nullptr);
	outputMoleculeShellForces(RADMoleculeTypes, RADMoleculeShellTypes, outFileWMShells, nullptr);

	//Average all FTC matrices for each type, normalise by time, then solve
	averageFTCMatricies(molecules, moleculeTypes);
	divideFTCMatricies(moleculeTypes, (double)totalFramesRun);
//	removeNonBondedInteractionsFromFTCMatricies(moleculeTypes);
	solveFTCMatrices(moleculeTypes);

	//Output matrices and eigenvalues/vectors
	outputFCMatrices(moleculeTypes, outFileFC);
	outputTCMatrices(moleculeTypes, outFileTC);
	outputFTCMatrices(moleculeTypes, outFileFTC);
	outputWMFTCMatrices(moleculeTypes, outFileWMFC);

	outputFEigenvalues(moleculeTypes, outFileFEVal);
	outputTEigenvalues(moleculeTypes, outFileTEVal);
	outputFTEigenvalues(moleculeTypes, outFileFTEVal);
	outputWMEigenvalues(moleculeTypes, outFileWMEVal);

	outputFEigenvectors(moleculeTypes, outFileFEVec);
	outputTEigenvectors(moleculeTypes, outFileTEVec);
	outputFTEigenvectors(moleculeTypes, outFileFTEVec);
	outputWMEigenvectors(moleculeTypes, outFileWMEVec);

	localisePositions(molecules);
	outputPositions(molecules, boxSize, outFilePositions);

	cout << "Done." << endl;

	return 0;
}
