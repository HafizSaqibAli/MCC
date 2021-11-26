/* Fragment RAD, Dr Jonathan Higham
** Created 22/08/2015
** Modified 26/07/2016
**
** This function applies the RAD criteria to a vector of fragments (usually heavy atom + hydrogen groups)
** and determines which fragments share a binary correlation.
**
*/

#ifndef CUTOFF_SHELLS
#define CUTOFF_SHELLS

void cutoffShells(vector<molecule> &molecules, double cutoff, float* boxSize, fragment* fragmenti, fragment* fragmentj){

	for (int i = 0; i < molecules.size(); i++){
		for (int j = 0; j < molecules[i].getSize(); j++){
			molecules[i].getShell(j)->clearShell();
		}
	}

	int totalFragments = 0;
	for (int i = 0; i < molecules.size(); i++){
		for (int j = 0; j < molecules[i].getSize(); j++){
			totalFragments = totalFragments + 1;
		}
	}

	vector< vector<double> > distances(totalFragments, vector<double>(totalFragments, 0));

	//Loop over each fragment
	int i = 0;
	for (int mi = 0; mi < molecules.size(); mi++){
		for (int fi = 0; fi < molecules[mi].getSize(); fi++){

			//Loop over each fragment again
			int j = 0;
			for (int mj = 0; mj < molecules.size(); mj++){
				for (int fj = 0; fj < molecules[mj].getSize(); fj++){

					double r = 0;

					//Add the squared distance in x,y,z
					double vecx = molecules[mi].getFragment(fi)->getPositionx() - molecules[mj].getFragment(fj)->getPositionx();
					double vecy = molecules[mi].getFragment(fi)->getPositiony() - molecules[mj].getFragment(fj)->getPositiony();
					double vecz = molecules[mi].getFragment(fi)->getPositionz() - molecules[mj].getFragment(fj)->getPositionz();

					//Update anything that goes over a box boundary
					if (vecx >(boxSize[0] / 2)){
						vecx = vecx - boxSize[0];
					}
					if (vecx < -(boxSize[0] / 2)){
						vecx = vecx + boxSize[0];
					}
					if (vecy >(boxSize[1] / 2)){
						vecy = vecy - boxSize[1];
					}
					if (vecy < -(boxSize[1] / 2)){
						vecy = vecy + boxSize[1];
					}
					if (vecz >(boxSize[2] / 2)){
						vecz = vecz - boxSize[2];
					}
					if (vecz < -(boxSize[2] / 2)){
						vecz = vecz + boxSize[2];
					}

					//Update distances
					distances[i][j] = sqrt((vecx*vecx) + (vecy*vecy) + (vecz*vecz));
					distances[j][i] = distances[i][j];

					j = j + 1;
				}
			}
			i = i + 1;
		}
	}

	i = 0;
	for (int mi = 0; mi < molecules.size(); mi++){
		for (int fi = 0; fi < molecules[mi].getSize(); fi++){

			bool inShell = true;
			vector<int> neighboursm;
			vector<int> neighboursf;
			vector<int> neighboursj;
			vector<bool> fragmentNotAlreadyInShell(totalFragments, true);

			bool checki = false;
			if(!fragmenti){
				checki = true;
			}
			else if(molecules[mi].getFragment(fi)->getType() == fragmenti->getType()){
				checki = true;
			}

			while (inShell && checki){

				int j12, m12, f12;
				double r12 = 0;

				//Loop over each fragment again
				int j = 0;
				for (int mj = 0; mj < molecules.size(); mj++){
					for (int fj = 0; fj < molecules[mj].getSize(); fj++){

						bool checkj = false;
						if(!fragmentj){
							checkj = true;
						}
						else if(molecules[mj].getFragment(fj)->getType() == fragmentj->getType()){
							checkj = true;
						}

						if ((distances[i][j] < r12 || r12 == 0) && fragmentNotAlreadyInShell[j] && j != i && checkj){
							m12 = mj;
							f12 = fj;
							j12 = j;
							r12 = distances[i][j];
						}
						j = j + 1;
					}
				}

				if (r12 >= cutoff){
					inShell = false;
				}

				if (inShell){
					fragmentNotAlreadyInShell[j12] = false;
					molecules[mi].getShell(fi)->addFragment(molecules[m12].getFragment(f12));
				}
			}
			i = i + 1;
		}
	}

	for (int mi = 0; mi < molecules.size(); mi++){
		for (int fi = 0; fi < molecules[mi].getSize(); fi++){
			molecules[mi].getShell(fi)->orderfragmentsUsingMass();
		}
	}
}//End of function

#endif
