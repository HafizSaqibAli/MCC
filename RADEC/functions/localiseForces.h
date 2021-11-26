/* Force Calculation  Code, Dr Jonathan Higham
** Modified 19/05/2017
**
** This function acts on a vector of molecule 
** objects to calculate the average forces and torques 
** on each fragment within each molecule
**
*/

#ifndef LOCALISE_FORCES
#define LOCALISE_FORCES

#include <math.h>

void localiseForces(vector<molecule> &molecules);
void localisePositions(vector<molecule> &molecules);
void localiseForces(vector<molecule> &molecules, bool localisePositions);

void localiseForces(vector<molecule> &molecules){
	localiseForces(molecules, false);
}

void localisePositions(vector<molecule> &molecules){
	localiseForces(molecules, true);
}

void localiseForces(vector<molecule> &molecules, bool localisePositions){

	//Loop over every molecule
	for (int i = 0; i < molecules.size(); i++){

		//Total mass for each molecule
		double massTotal = molecules[i].getMass();
		double massNorm = 1 / massTotal;
	
		//Whole molecule coordinate axes
		double WMPrincipalXaxisx, WMPrincipalXaxisy, WMPrincipalXaxisz;
		double WMPrincipalYaxisx, WMPrincipalYaxisy, WMPrincipalYaxisz;
		double WMPrincipalZaxisx, WMPrincipalZaxisy, WMPrincipalZaxisz;
		double xaxisMI, yaxisMI, zaxisMI;

		//Calculate the centre of mass for each molecule
		double ComX = 0;
		double ComY = 0;
		double ComZ = 0;

		for(int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){//Loop atoms
				double massFrac = molecules[i].getFragment(j)->getAtom(k)->getMass()*massNorm;
				ComX = ComX + molecules[i].getFragment(j)->getAtom(k)->getPositionx()*massFrac;
				ComY = ComY + molecules[i].getFragment(j)->getAtom(k)->getPositiony()*massFrac;
				ComZ = ComZ + molecules[i].getFragment(j)->getAtom(k)->getPositionz()*massFrac;
			}
		}

		//Calculate the total force on the molecule in the system frame
		double sysFrameForcex = 0;
		double sysFrameForcey = 0;
		double sysFrameForcez = 0;

		for(int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){//Loop atoms
				sysFrameForcex = sysFrameForcex + (molecules[i].getFragment(j)->getAtom(k)->getForcex()*4184*1e10)/(2*6.02214e23);
				sysFrameForcey = sysFrameForcey + (molecules[i].getFragment(j)->getAtom(k)->getForcey()*4184*1e10)/(2*6.02214e23);
				sysFrameForcez = sysFrameForcez + (molecules[i].getFragment(j)->getAtom(k)->getForcez()*4184*1e10)/(2*6.02214e23);
			}
		}

		//Calculate Moment of Inertia Tensor
		Eigen::Matrix3d MITensor;
		MITensor << 0, 0, 0, 0, 0, 0, 0, 0, 0;
	
		for(int j = 0; j < molecules[i].getSize(); j++){
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){//Loop atoms		

				double x = (molecules[i].getFragment(j)->getAtom(k)->getPositionx() - ComX)/1e10;
				double y = (molecules[i].getFragment(j)->getAtom(k)->getPositiony() - ComY)/1e10;
				double z = (molecules[i].getFragment(j)->getAtom(k)->getPositionz() - ComZ)/1e10;
				double m = (molecules[i].getFragment(j)->getAtom(k)->getMass())/(6.02214e23*1e3);

				MITensor(0, 0) = MITensor(0, 0) + ((y*y + z*z)*m);
				MITensor(0, 1) = MITensor(0, 1) - ((x*y)*m);
				MITensor(0, 2) = MITensor(0, 2) - ((x*z)*m);

				MITensor(1, 0) = MITensor(1, 0) - ((x*y)*m);
				MITensor(1, 1) = MITensor(1, 1) + ((z*z + x*x)*m);
				MITensor(1, 2) = MITensor(1, 2) - ((y*z)*m);

				MITensor(2, 0) = MITensor(2, 0) - ((x*z)*m);
				MITensor(2, 1) = MITensor(2, 1) - ((y*z)*m);
				MITensor(2, 2) = MITensor(2, 2) + ((x*x + y*y)*m);
			}
		}

		//Solve MITensor and calculate eigenvectors
		Eigen::EigenSolver<Eigen::Matrix3d> MISolver;
		MISolver.compute(MITensor, true);

		double maxEigenvalue = abs(MISolver.eigenvalues()[0].real());
		if(MISolver.eigenvalues()[1].real() > maxEigenvalue)maxEigenvalue = MISolver.eigenvalues()[1].real();
		if(MISolver.eigenvalues()[2].real() > maxEigenvalue)maxEigenvalue = MISolver.eigenvalues()[2].real();

		double minEigenvalue = abs(MISolver.eigenvalues()[0].real());
		if(MISolver.eigenvalues()[1].real() < minEigenvalue)minEigenvalue = MISolver.eigenvalues()[1].real();
		if(MISolver.eigenvalues()[2].real() < minEigenvalue)minEigenvalue = MISolver.eigenvalues()[2].real();

		for(int e = 0; e < 3; e++){
			if(MISolver.eigenvalues()[e].real() == maxEigenvalue){
				WMPrincipalXaxisx = MISolver.eigenvectors().col(e)[0].real();
				WMPrincipalXaxisy = MISolver.eigenvectors().col(e)[1].real();
				WMPrincipalXaxisz = MISolver.eigenvectors().col(e)[2].real();
				xaxisMI = MISolver.eigenvalues()[e].real();
			}

			else if(MISolver.eigenvalues()[e].real() == minEigenvalue){
				WMPrincipalZaxisx = MISolver.eigenvectors().col(e)[0].real();
				WMPrincipalZaxisy = MISolver.eigenvectors().col(e)[1].real();
				WMPrincipalZaxisz = MISolver.eigenvectors().col(e)[2].real();
				zaxisMI = MISolver.eigenvalues()[e].real();
			}

			else{
				WMPrincipalYaxisx = MISolver.eigenvectors().col(e)[0].real();
				WMPrincipalYaxisy = MISolver.eigenvectors().col(e)[1].real();
				WMPrincipalYaxisz = MISolver.eigenvectors().col(e)[2].real();
				yaxisMI = MISolver.eigenvalues()[e].real();
			}
		}

		double xcrossyx = WMPrincipalXaxisy*WMPrincipalYaxisz - WMPrincipalXaxisz*WMPrincipalYaxisy;
		double xcrossyy  = WMPrincipalXaxisz*WMPrincipalYaxisx - WMPrincipalXaxisx*WMPrincipalYaxisz;
		double xcrossyz  = WMPrincipalXaxisx*WMPrincipalYaxisy - WMPrincipalXaxisy*WMPrincipalYaxisx;

		double xcrossydotz = xcrossyx*WMPrincipalZaxisx + xcrossyy*WMPrincipalZaxisy + xcrossyz*WMPrincipalZaxisz;

		if(xcrossydotz < 0){
			WMPrincipalZaxisx = -WMPrincipalZaxisx;
			WMPrincipalZaxisy = -WMPrincipalZaxisy;
			WMPrincipalZaxisz = -WMPrincipalZaxisz;
		}

		//Don't localise single-atom molecules
		if (molecules[i].getSize() == 1 && molecules[i].getFragment(0)->getSize() == 1){
			molecules[i].getFragment(0)->setForcex(molecules[i].getFragment(0)->getAtom(0)->getForcex()
				/ ((sqrt(massTotal)*6.02214e23) / (sqrt(6.02214e23*1e3)*0.5*4184e10)));
			molecules[i].getFragment(0)->setForcey(molecules[i].getFragment(0)->getAtom(0)->getForcey()
				/ ((sqrt(massTotal)*6.02214e23) / (sqrt(6.02214e23*1e3)*0.5*4184e10)));
			molecules[i].getFragment(0)->setForcez(molecules[i].getFragment(0)->getAtom(0)->getForcez()
				/ ((sqrt(massTotal)*6.02214e23) / (sqrt(6.02214e23*1e3)*0.5*4184e10)));

			molecules[i].setForcex(molecules[i].getFragment(0)->getForcex());
			molecules[i].setForcey(molecules[i].getFragment(0)->getForcey());
			molecules[i].setForcez(molecules[i].getFragment(0)->getForcez());

			molecules[i].setTorquex(0);
			molecules[i].setTorquey(0);
			molecules[i].setTorquez(0);

			molecules[i].setForceAverage((abs(molecules[i].getForcex()) + abs(molecules[i].getForcey()) + abs(molecules[i].getForcez())) / 3);
			molecules[i].setForceSquared(((molecules[i].getForcex()*molecules[i].getForcex()) + (molecules[i].getForcey()*molecules[i].getForcey()) + (molecules[i].getForcez()*molecules[i].getForcez())) / 3);
		}//End of single atom molecule condition

		//Localise single fragment molecules into one local frame
		else if (molecules[i].getSize() == 1 && molecules[i].getFragment(0)->getSize() > 1){

			//Get the first fragment only
			double j = 0;
			
			molecules[i].getFragment(j)->setForcex(0);
			molecules[i].getFragment(j)->setForcey(0);
			molecules[i].getFragment(j)->setForcez(0);

			//Rotate molecule/fragment forces into molecule/fragment coordinate frame
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){

				double forcex = (molecules[i].getFragment(j)->getAtom(k)->getForcex() * 4184 * 1e10) / (2 * 6.02214e23);
				double forcey = (molecules[i].getFragment(j)->getAtom(k)->getForcey() * 4184 * 1e10) / (2 * 6.02214e23);
				double forcez = (molecules[i].getFragment(j)->getAtom(k)->getForcez() * 4184 * 1e10) / (2 * 6.02214e23);

				double newForcex = WMPrincipalXaxisx*forcex
					+ WMPrincipalXaxisy*forcey
					+ WMPrincipalXaxisz*forcez;
				double newForcey = WMPrincipalYaxisx*forcex
					+ WMPrincipalYaxisy*forcey
					+ WMPrincipalYaxisz*forcez;
				double newForcez = WMPrincipalZaxisx*forcex
					+ WMPrincipalZaxisy*forcey
					+ WMPrincipalZaxisz*forcez;

				molecules[i].getFragment(j)->addForcex(newForcex / sqrt(massTotal / (6.02214e23*1e3)));
				molecules[i].getFragment(j)->addForcey(newForcey / sqrt(massTotal / (6.02214e23*1e3)));
				molecules[i].getFragment(j)->addForcez(newForcez / sqrt(massTotal / (6.02214e23*1e3)));

				molecules[i].setForcex(molecules[i].getFragment(j)->getForcex());
				molecules[i].setForcey(molecules[i].getFragment(j)->getForcey());
				molecules[i].setForcez(molecules[i].getFragment(j)->getForcez());
			}

			molecules[i].setForceAverage((abs(molecules[i].getForcex()) + abs(molecules[i].getForcey()) + abs(molecules[i].getForcez())) / 3);
			molecules[i].setForceSquared(((molecules[i].getForcex()*molecules[i].getForcex()) + (molecules[i].getForcey()*molecules[i].getForcey()) + (molecules[i].getForcez()*molecules[i].getForcez())) / 3);

			molecules[i].getFragment(j)->setTorquex(0);
			molecules[i].getFragment(j)->setTorquey(0);
			molecules[i].getFragment(j)->setTorquez(0);

			//Rotate molecule/fragment torques into molecule/fragment coordinate frame
			for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){

				double newx = (WMPrincipalXaxisx*(molecules[i].getFragment(j)->getAtom(k)->getPositionx() - ComX)
					+ WMPrincipalXaxisy*(molecules[i].getFragment(j)->getAtom(k)->getPositiony() - ComY)
					+ WMPrincipalXaxisz*(molecules[i].getFragment(j)->getAtom(k)->getPositionz() - ComZ));
				double newy = (WMPrincipalYaxisx*(molecules[i].getFragment(j)->getAtom(k)->getPositionx() - ComX)
					+ WMPrincipalYaxisy*(molecules[i].getFragment(j)->getAtom(k)->getPositiony() - ComY)
					+ WMPrincipalYaxisz*(molecules[i].getFragment(j)->getAtom(k)->getPositionz() - ComZ));
				double newz = (WMPrincipalZaxisx*(molecules[i].getFragment(j)->getAtom(k)->getPositionx() - ComX)
					+ WMPrincipalZaxisy*(molecules[i].getFragment(j)->getAtom(k)->getPositiony() - ComY)
					+ WMPrincipalZaxisz*(molecules[i].getFragment(j)->getAtom(k)->getPositionz() - ComZ));

				double forcex = (molecules[i].getFragment(j)->getAtom(k)->getForcex() * 4184 * 1e10) / (2 * 6.02214e23);
				double forcey = (molecules[i].getFragment(j)->getAtom(k)->getForcey() * 4184 * 1e10) / (2 * 6.02214e23);
				double forcez = (molecules[i].getFragment(j)->getAtom(k)->getForcez() * 4184 * 1e10) / (2 * 6.02214e23);

				double newForcex = WMPrincipalXaxisx*forcex
					+ WMPrincipalXaxisy*forcey
					+ WMPrincipalXaxisz*forcez;
				double newForcey = WMPrincipalYaxisx*forcex
					+ WMPrincipalYaxisy*forcey
					+ WMPrincipalYaxisz*forcez;
				double newForcez = WMPrincipalZaxisx*forcex
					+ WMPrincipalZaxisy*forcey
					+ WMPrincipalZaxisz*forcez;

				//Torque calculation
				double newTorquex = (newy*newForcez - newz*newForcey) / 1e10;
				double newTorquey = (newz*newForcex - newx*newForcez) / 1e10;
				double newTorquez = (newx*newForcey - newy*newForcex) / 1e10;

				if (xaxisMI != 0)molecules[i].getFragment(j)->addTorquex(newTorquex / sqrt(xaxisMI));
				if (yaxisMI != 0)molecules[i].getFragment(j)->addTorquey(newTorquey / sqrt(yaxisMI));
				if (zaxisMI != 0)molecules[i].getFragment(j)->addTorquez(newTorquez / sqrt(zaxisMI));

				molecules[i].setTorquex(molecules[i].getFragment(j)->getTorquex());
				molecules[i].setTorquey(molecules[i].getFragment(j)->getTorquey());
				molecules[i].setTorquez(molecules[i].getFragment(j)->getTorquez());
			}
		}//End of single fragment molecule condition

		else if (molecules[i].getSize() > 1){

			molecules[i].setForcex(0);
			molecules[i].setForcey(0);
			molecules[i].setForcez(0);

			//Rotate molecule forces into molecule coordinate frame
			for(int j = 0; j < molecules[i].getSize(); j++){
				for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){

					double forcex = (molecules[i].getFragment(j)->getAtom(k)->getForcex()*4184*1e10)/(2*6.02214e23);
					double forcey = (molecules[i].getFragment(j)->getAtom(k)->getForcey()*4184*1e10)/(2*6.02214e23);
					double forcez = (molecules[i].getFragment(j)->getAtom(k)->getForcez()*4184*1e10)/(2*6.02214e23);

					forcex = forcex - 1000000*sysFrameForcex;
					forcex = forcey - 1000000*sysFrameForcey;
					forcex = forcez - 1000000*sysFrameForcez;

					double newForcex = WMPrincipalXaxisx*forcex
							+ WMPrincipalXaxisy*forcey
							+ WMPrincipalXaxisz*forcez;
					double newForcey = WMPrincipalYaxisx*forcex
							+ WMPrincipalYaxisy*forcey
							+ WMPrincipalYaxisz*forcez;
					double newForcez = WMPrincipalZaxisx*forcex
							+ WMPrincipalZaxisy*forcey
							+ WMPrincipalZaxisz*forcez;

					molecules[i].addForcex(newForcex / sqrt(massTotal/(6.02214e23*1e3)));
					molecules[i].addForcey(newForcey / sqrt(massTotal/(6.02214e23*1e3)));
					molecules[i].addForcez(newForcez / sqrt(massTotal/(6.02214e23*1e3)));
				}
			}
			
			molecules[i].setForceAverage((abs(molecules[i].getForcex()) + abs(molecules[i].getForcey()) + abs(molecules[i].getForcez())) / 3);
			molecules[i].setForceSquared(((molecules[i].getForcex()*molecules[i].getForcex()) + (molecules[i].getForcey()*molecules[i].getForcey()) + (molecules[i].getForcez()*molecules[i].getForcez())) / 3);

			molecules[i].setTorquex(0);
			molecules[i].setTorquey(0);
			molecules[i].setTorquez(0);

			//Rotate molecule torques into molecule coordinate frame
			for(int j = 0; j < molecules[i].getSize(); j++){
				for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){
		
					double newx = (WMPrincipalXaxisx*(molecules[i].getFragment(j)->getAtom(k)->getPositionx() - ComX)
							+ WMPrincipalXaxisy*(molecules[i].getFragment(j)->getAtom(k)->getPositiony() - ComY)
							+ WMPrincipalXaxisz*(molecules[i].getFragment(j)->getAtom(k)->getPositionz() - ComZ));
					double newy = (WMPrincipalYaxisx*(molecules[i].getFragment(j)->getAtom(k)->getPositionx() - ComX)
							+ WMPrincipalYaxisy*(molecules[i].getFragment(j)->getAtom(k)->getPositiony() - ComY)
							+ WMPrincipalYaxisz*(molecules[i].getFragment(j)->getAtom(k)->getPositionz() - ComZ));
					double newz = (WMPrincipalZaxisx*(molecules[i].getFragment(j)->getAtom(k)->getPositionx() - ComX)
							+ WMPrincipalZaxisy*(molecules[i].getFragment(j)->getAtom(k)->getPositiony() - ComY)
							+ WMPrincipalZaxisz*(molecules[i].getFragment(j)->getAtom(k)->getPositionz() - ComZ));

					double forcex = (molecules[i].getFragment(j)->getAtom(k)->getForcex()*4184*1e10)/(2*6.02214e23);
					double forcey = (molecules[i].getFragment(j)->getAtom(k)->getForcey()*4184*1e10)/(2*6.02214e23);
					double forcez = (molecules[i].getFragment(j)->getAtom(k)->getForcez()*4184*1e10)/(2*6.02214e23);

					double newForcex = WMPrincipalXaxisx*forcex
							+ WMPrincipalXaxisy*forcey
							+ WMPrincipalXaxisz*forcez;
					double newForcey = WMPrincipalYaxisx*forcex
							+ WMPrincipalYaxisy*forcey
							+ WMPrincipalYaxisz*forcez;
					double newForcez = WMPrincipalZaxisx*forcex
							+ WMPrincipalZaxisy*forcey
							+ WMPrincipalZaxisz*forcez;

					//Torque calculation
					double newTorquex = (newy*newForcez - newz*newForcey)/1e10;
					double newTorquey = (newz*newForcex - newx*newForcez)/1e10;
					double newTorquez = (newx*newForcey - newy*newForcex)/1e10;

					if(xaxisMI != 0)molecules[i].addTorquex(newTorquex / sqrt(xaxisMI));
					if(yaxisMI != 0)molecules[i].addTorquey(newTorquey / sqrt(yaxisMI));
					if(zaxisMI != 0)molecules[i].addTorquez(newTorquez / sqrt(zaxisMI));
				}
			}


			//Loop over each fragment in the molecule
			for(int j = 0; j < molecules[i].getSize(); j++){

				double fragmentMassTotal = molecules[i].getFragment(j)->getMass();
				double fragmentMassNorm = 1 / fragmentMassTotal;

				//Fragment coordinate axes
				double FTxaxisx, FTxaxisy, FTxaxisz;
				double FTyaxisx, FTyaxisy, FTyaxisz;
				double FTzaxisx, FTzaxisy, FTzaxisz;
				
				//Calculate the forces on each fragment in the fragment MI frame
				molecules[i].getFragment(j)->setForcex(0);
				molecules[i].getFragment(j)->setForcey(0);
				molecules[i].getFragment(j)->setForcez(0);
			
				//Set fragment forces in molecule coordinate frame
				for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){

					double forcex = (molecules[i].getFragment(j)->getAtom(k)->getForcex()*4184*1e10)/(2*6.02214e23);
					double forcey = (molecules[i].getFragment(j)->getAtom(k)->getForcey()*4184*1e10)/(2*6.02214e23);
					double forcez = (molecules[i].getFragment(j)->getAtom(k)->getForcez()*4184*1e10)/(2*6.02214e23);

					double newForcex = WMPrincipalXaxisx*forcex
							+ WMPrincipalXaxisy*forcey
							+ WMPrincipalXaxisz*forcez;
					double newForcey = WMPrincipalYaxisx*forcex
							+ WMPrincipalYaxisy*forcey
							+ WMPrincipalYaxisz*forcez;
					double newForcez = WMPrincipalZaxisx*forcex
							+ WMPrincipalZaxisy*forcey
							+ WMPrincipalZaxisz*forcez;

					molecules[i].getFragment(j)->addForcex(newForcex / sqrt(fragmentMassTotal / (6.02214e23*1e3)));
					molecules[i].getFragment(j)->addForcey(newForcey / sqrt(fragmentMassTotal / (6.02214e23*1e3)));
					molecules[i].getFragment(j)->addForcez(newForcez / sqrt(fragmentMassTotal / (6.02214e23*1e3)));
				}

				//Torques on current fragment
				molecules[i].getFragment(j)->setTorquex(0);
				molecules[i].getFragment(j)->setTorquey(0);
				molecules[i].getFragment(j)->setTorquez(0);

				//Calculate the torques on each molecule with <= 1 fragment
				if(molecules[i].getFragment(j)->getSize() == 0 || molecules[i].getFragment(j)->getSize() == 1){
					//Do nothing, zero torque
				}

				//Calculate the torques on each fragment with >1 atom
				else if(molecules[i].getFragment(j)->getSize() > 0){

					//Calculate the centre of mass (heavy atom position) for each fragment
					double fragmentComX = molecules[i].getFragment(j)->getAtom(0)->getPositionx();
					double fragmentComY = molecules[i].getFragment(j)->getAtom(0)->getPositiony();
					double fragmentComZ = molecules[i].getFragment(j)->getAtom(0)->getPositionz();
	
					//Torque axes for fragments with exactly one other bonded fragment (2 bonded fragments including themselves), and at least one Hydrogen
					if (molecules[i].getFragment(j)->getNumBondedFragments() == 2){
						int k = molecules[i].getFragment(j)->getBondedFragment(1);

						FTxaxisx = (molecules[i].getFragment(k)->getAtom(0)->getPositionx() - fragmentComX) / 1e10;
						FTxaxisy = (molecules[i].getFragment(k)->getAtom(0)->getPositiony() - fragmentComY) / 1e10;
						FTxaxisz = (molecules[i].getFragment(k)->getAtom(0)->getPositionz() - fragmentComZ) / 1e10;

						//Calculate XH vector
						double XHvectorx = molecules[i].getFragment(j)->getAtom(0)->getPositionx() -
							molecules[i].getFragment(j)->getAtom(1)->getPositionx();
						double XHvectory = molecules[i].getFragment(j)->getAtom(0)->getPositiony() -
							molecules[i].getFragment(j)->getAtom(1)->getPositiony();
						double XHvectorz = molecules[i].getFragment(j)->getAtom(0)->getPositionz() -
							molecules[i].getFragment(j)->getAtom(1)->getPositionz();

						//Torque y axis is the cross product of the covalent bond and the first hydrogen bond (i.e. perp to plane of XH)
						FTyaxisx = XHvectory*FTxaxisz - XHvectorz*FTxaxisy;
						FTyaxisy = XHvectorz*FTxaxisx - XHvectorx*FTxaxisz;
						FTyaxisz = XHvectorx*FTxaxisy - XHvectory*FTxaxisx;

						//Torque z axis is the cross product of other two torque axes
						FTzaxisx = FTyaxisy*FTxaxisz - FTyaxisz*FTxaxisy;
						FTzaxisy = FTyaxisz*FTxaxisx - FTyaxisx*FTxaxisz;
						FTzaxisz = FTyaxisx*FTxaxisy - FTyaxisy*FTxaxisx;
					}
					
					//Torque axes for fragments with more than one other bonded fragment (>2 bonded fragments including themselves)
					if (molecules[i].getFragment(j)->getNumBondedFragments() > 2){

						FTxaxisx = 0;
						FTxaxisy = 0;
						FTxaxisz = 0;
						//First axis is the average of all covalent bonds
						for (int k2 = 1; k2 < molecules[i].getFragment(j)->getNumBondedFragments(); k2++){
							int k = molecules[i].getFragment(j)->getBondedFragment(k2);

							FTxaxisx = FTxaxisx + (molecules[i].getFragment(k)->getAtom(0)->getPositionx() - fragmentComX) / 1e10;
							FTxaxisy = FTxaxisy + (molecules[i].getFragment(k)->getAtom(0)->getPositiony() - fragmentComY) / 1e10;
							FTxaxisz = FTxaxisz + (molecules[i].getFragment(k)->getAtom(0)->getPositionz() - fragmentComZ) / 1e10;
						}

						int k1 = molecules[i].getFragment(j)->getBondedFragment(1);
						int k2 = molecules[i].getFragment(j)->getBondedFragment(2);

						//Calculate X1X2 vector
						double X1X2vectorx = molecules[i].getFragment(k1)->getAtom(0)->getPositionx() -
							molecules[i].getFragment(k2)->getAtom(0)->getPositionx();
						double X1X2vectory = molecules[i].getFragment(k1)->getAtom(0)->getPositiony() -
							molecules[i].getFragment(k2)->getAtom(0)->getPositiony();
						double X1X2vectorz = molecules[i].getFragment(k1)->getAtom(0)->getPositionz() -
							molecules[i].getFragment(k2)->getAtom(0)->getPositionz();

						//Torque y axis is the cross product of the average covalent bond and the first covalent bond (i.e. perp to plane of heavy atoms X1 and X2)
						FTyaxisx = X1X2vectory*FTxaxisz - X1X2vectorz*FTxaxisy;
						FTyaxisy = X1X2vectorz*FTxaxisx - X1X2vectorx*FTxaxisz;
						FTyaxisz = X1X2vectorx*FTxaxisy - X1X2vectory*FTxaxisx;

						//Torque z axis is the cross product of other two torque axes
						FTzaxisx = FTyaxisy*FTxaxisz - FTyaxisz*FTxaxisy;
						FTzaxisy = FTyaxisz*FTxaxisx - FTyaxisx*FTxaxisz;
						FTzaxisz = FTyaxisx*FTxaxisy - FTyaxisy*FTxaxisx;
					}

					double modFTx = FTxaxisx*FTxaxisx + FTxaxisy*FTxaxisy + FTxaxisz*FTxaxisz;
					double modFTy = FTyaxisx*FTyaxisx + FTyaxisy*FTyaxisy + FTyaxisz*FTyaxisz;
					double modFTz = FTzaxisx*FTzaxisx + FTzaxisy*FTzaxisy + FTzaxisz*FTzaxisz;

					FTxaxisx = FTxaxisx / sqrt(modFTx);
					FTxaxisy = FTxaxisy / sqrt(modFTx);
					FTxaxisz = FTxaxisz / sqrt(modFTx);

					FTyaxisx = FTyaxisx / sqrt(modFTy);
					FTyaxisy = FTyaxisy / sqrt(modFTy);
					FTyaxisz = FTyaxisz / sqrt(modFTy);

					FTzaxisx = FTzaxisx / sqrt(modFTz);
					FTzaxisy = FTzaxisy / sqrt(modFTz);
					FTzaxisz = FTzaxisz / sqrt(modFTz);

					//Take the dot product of the first principle axis with the CoM->atom1 vector
					double x = (molecules[i].getFragment(j)->getAtom(0)->getPositionx() - fragmentComX) / 1e10;
					double y = (molecules[i].getFragment(j)->getAtom(0)->getPositiony() - fragmentComY) / 1e10;
					double z = (molecules[i].getFragment(j)->getAtom(0)->getPositionz() - fragmentComZ) / 1e10;

					//Flip wach principle axis if it is pointing out of the fragment
					double dotProd = (x*FTxaxisx) + (y*FTxaxisy) + (z*FTxaxisz);
					if (dotProd < 0){
						FTxaxisx = -FTxaxisx;
						FTxaxisy = -FTxaxisy;
						FTxaxisz = -FTxaxisz;
					}
					dotProd = (x*FTyaxisx) + (y*FTyaxisy) + (z*FTyaxisz);
					if (dotProd < 0){
						FTyaxisx = -FTyaxisx;
						FTyaxisy = -FTyaxisy;
						FTyaxisz = -FTyaxisz;
					}
					dotProd = (x*FTzaxisx) + (y*FTzaxisy) + (z*FTzaxisz);
					if (dotProd < 0){
						FTzaxisx = -FTzaxisx;
						FTzaxisy = -FTzaxisy;
						FTzaxisz = -FTzaxisz;
					}

					//Calculate the moment of inertia in the fragment coordinate system
					double xaxisMI = 0;
					double yaxisMI = 0;
					double zaxisMI = 0;

					for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){

						double dx = (molecules[i].getFragment(j)->getAtom(k)->getPositionx() - fragmentComX) / 1e10;
						double dy = (molecules[i].getFragment(j)->getAtom(k)->getPositiony() - fragmentComY) / 1e10;
						double dz = (molecules[i].getFragment(j)->getAtom(k)->getPositionz() - fragmentComZ) / 1e10;

						double FTaxisdxx = FTxaxisy*dz - FTxaxisz*dy;
						double FTaxisdxy = FTxaxisz*dx - FTxaxisx*dz;
						double FTaxisdxz = FTxaxisx*dy - FTxaxisy*dx;

						double FTaxisdyx = FTyaxisy*dz - FTyaxisz*dy;
						double FTaxisdyy = FTyaxisz*dx - FTyaxisx*dz;
						double FTaxisdyz = FTyaxisx*dy - FTyaxisy*dx;

						double FTaxisdzx = FTzaxisy*dz - FTzaxisz*dy;
						double FTaxisdzy = FTzaxisz*dx - FTzaxisx*dz;
						double FTaxisdzz = FTzaxisx*dy - FTzaxisy*dx;

						double daxisx = (sqrt(FTaxisdxx*FTaxisdxx + FTaxisdxy*FTaxisdxy + FTaxisdxz*FTaxisdxz));
						double daxisy = (sqrt(FTaxisdyx*FTaxisdyx + FTaxisdyy*FTaxisdyy + FTaxisdyz*FTaxisdyz));
						double daxisz = (sqrt(FTaxisdzx*FTaxisdzx + FTaxisdzy*FTaxisdzy + FTaxisdzz*FTaxisdzz));

						xaxisMI = xaxisMI + daxisx*daxisx*molecules[i].getFragment(j)->getAtom(k)->getMass() / (6.02214e23*1e3);
						yaxisMI = yaxisMI + daxisy*daxisy*molecules[i].getFragment(j)->getAtom(k)->getMass() / (6.02214e23*1e3);
						zaxisMI = zaxisMI + daxisz*daxisz*molecules[i].getFragment(j)->getAtom(k)->getMass() / (6.02214e23*1e3);

					}
					
					//Remove the lowest MI degree of freedom if the fragment has a single hydrogen
					if (molecules[i].getFragment(j)->getSize() == 2){

						if (xaxisMI < yaxisMI && xaxisMI < zaxisMI){
							xaxisMI = 0;
						}
						if (yaxisMI < zaxisMI && yaxisMI < xaxisMI){
							yaxisMI = 0;
						}
						if (zaxisMI < xaxisMI && zaxisMI < yaxisMI){
							zaxisMI = 0;
						}
					}

//					xaxisMI = xaxisMI/fragmentMassTotal;
//					yaxisMI = yaxisMI/fragmentMassTotal;
//					zaxisMI = zaxisMI/fragmentMassTotal;

					//Rotate fragment into local coordinate frame and add up torques
					for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){

						double newx = (FTxaxisx*(molecules[i].getFragment(j)->getAtom(k)->getPositionx() - fragmentComX)
									+ FTxaxisy*(molecules[i].getFragment(j)->getAtom(k)->getPositiony() - fragmentComY)
									+ FTxaxisz*(molecules[i].getFragment(j)->getAtom(k)->getPositionz() - fragmentComZ));
						double newy = (FTyaxisx*(molecules[i].getFragment(j)->getAtom(k)->getPositionx() - fragmentComX)
									+ FTyaxisy*(molecules[i].getFragment(j)->getAtom(k)->getPositiony() - fragmentComY)
									+ FTyaxisz*(molecules[i].getFragment(j)->getAtom(k)->getPositionz() - fragmentComZ));
						double newz = (FTzaxisx*(molecules[i].getFragment(j)->getAtom(k)->getPositionx() - fragmentComX)
									+ FTzaxisy*(molecules[i].getFragment(j)->getAtom(k)->getPositiony() - fragmentComY)
									+ FTzaxisz*(molecules[i].getFragment(j)->getAtom(k)->getPositionz() - fragmentComZ));

						double forcex = (molecules[i].getFragment(j)->getAtom(k)->getForcex()*4184*1e10)/(2*6.02214e23);
						double forcey = (molecules[i].getFragment(j)->getAtom(k)->getForcey()*4184*1e10)/(2*6.02214e23);
						double forcez = (molecules[i].getFragment(j)->getAtom(k)->getForcez()*4184*1e10)/(2*6.02214e23);

						double newForcex = FTxaxisx*forcex
								+ FTxaxisy*forcey
								+ FTxaxisz*forcez;
						double newForcey = FTyaxisx*forcex
								+ FTyaxisy*forcey
								+ FTyaxisz*forcez;
						double newForcez = FTzaxisx*forcex
								+ FTzaxisy*forcey
								+ FTzaxisz*forcez;

						//Torque calculation
						double newTorquex = (newy*newForcez - newz*newForcey)/1e10;
						double newTorquey = (newz*newForcex - newx*newForcez)/1e10;
						double newTorquez = (newx*newForcey - newy*newForcex)/1e10;

						if(xaxisMI != 0)molecules[i].getFragment(j)->addTorquex(newTorquex / sqrt(xaxisMI));
						if(yaxisMI != 0)molecules[i].getFragment(j)->addTorquey(newTorquey / sqrt(yaxisMI));
						if(zaxisMI != 0)molecules[i].getFragment(j)->addTorquez(newTorquez / sqrt(zaxisMI));
					}
				}//End of condition for fragments with >1 atom
			}//End of loop over fragment in multifragment molecule
		}///End of multiple fragment molecule condition

		//The following localises the positions in the molecule frame. 
		//It serves no purpose other than as a test/debug environment for any frame issues
		if ((molecules[i].getSize() > 1 || molecules[i].getFragment(0)->getSize() > 1) && localisePositions){
			for(int j = 0; j < molecules[i].getSize(); j++){				//Loop over every fragment
				for (int k = 0; k < molecules[i].getFragment(j)->getSize(); k++){	//Loop over every atom
			
					double newx = (WMPrincipalXaxisx*(molecules[i].getFragment(j)->getAtom(k)->getPositionx() - ComX)
							+ WMPrincipalXaxisy*(molecules[i].getFragment(j)->getAtom(k)->getPositiony() - ComY)
							+ WMPrincipalXaxisz*(molecules[i].getFragment(j)->getAtom(k)->getPositionz() - ComZ));
					double newy = (WMPrincipalYaxisx*(molecules[i].getFragment(j)->getAtom(k)->getPositionx() - ComX)
							+ WMPrincipalYaxisy*(molecules[i].getFragment(j)->getAtom(k)->getPositiony() - ComY)
							+ WMPrincipalYaxisz*(molecules[i].getFragment(j)->getAtom(k)->getPositionz() - ComZ));
					double newz = (WMPrincipalZaxisx*(molecules[i].getFragment(j)->getAtom(k)->getPositionx() - ComX)
							+ WMPrincipalZaxisy*(molecules[i].getFragment(j)->getAtom(k)->getPositiony() - ComY)
							+ WMPrincipalZaxisz*(molecules[i].getFragment(j)->getAtom(k)->getPositionz() - ComZ));

					molecules[i].getFragment(j)->getAtom(k)->setPositionx(newx);
					molecules[i].getFragment(j)->getAtom(k)->setPositiony(newy);
					molecules[i].getFragment(j)->getAtom(k)->setPositionz(newz);
				}
			}
		}
	}//End of molecule loop
}//end of function

#endif