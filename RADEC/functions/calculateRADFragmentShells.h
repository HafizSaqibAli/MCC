/* Fragment RAD, Dr Jonathan Higham
** Created 22/08/2015
** Modified 26/07/2016
**
** This function applies the RAD criteria to a vector of fragments (usually heavy atom + hydrogen groups)
** and determines which fragments share a binary correlation.
**
*/

#ifndef RAD_FRAGMENT_SHELLS
#define RAD_FRAGMENT_SHELLS


void RADFragmentShells(vector<molecule> &molecules, float* boxSize);
void RADFragmentShells(vector<molecule> &molecules, float* boxSize, fragment* fragmenti, fragment* fragmentj);

void RADFragmentShells(vector<molecule> &molecules, float* boxSize){
	RADFragmentShells(molecules, boxSize, nullptr, nullptr);
}

void RADFragmentShells(vector<molecule> &molecules, float* boxSize, fragment* fragmenti, fragment* fragmentj){

	for (int i = 0; i < molecules.size(); i++){
		for (int j = 0; j < molecules[i].getSize(); j++){
			molecules[i].getFragmentShell(j)->clearShell();
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
					if (vecx > (boxSize[0] / 2)){
						vecx = vecx - boxSize[0];
					}
					if (vecx < -(boxSize[0] / 2)){
						vecx = vecx + boxSize[0];
					}
					if (vecy > (boxSize[1] / 2)){
						vecy = vecy - boxSize[1];
					}
					if (vecy < -(boxSize[1] / 2)){
						vecy = vecy + boxSize[1];
					}
					if (vecz > (boxSize[2] / 2)){
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
			if(fragmenti == nullptr){
				checki = true;
			}
			else if(molecules[mi].getFragment(fi)->getType() == fragmenti->getType()){
				checki = true;
			}

			while (inShell && checki){

				int j12, m12, f12, mb;
				double r12 = 0;

				//Loop over each fragment again
				int j = 0;
				for (int mj = 0; mj < molecules.size(); mj++){
					for (int fj = 0; fj < molecules[mj].getSize(); fj++){

						if ((distances[i][j] < r12 || r12 == 0) && fragmentNotAlreadyInShell[j] && j != i){
							m12 = mj;
							f12 = fj;
							j12 = j;
							r12 = distances[i][j];
						}
						j = j + 1;
					}
				}

				neighboursm.push_back(m12);
				neighboursf.push_back(f12);
				neighboursj.push_back(j12);
				fragmentNotAlreadyInShell[j12] = false;
				double a12 = 1 / pow(r12, 2);
				double vecx12 = molecules[mi].getFragment(fi)->getPositionx() - molecules[m12].getFragment(f12)->getPositionx();
				double vecy12 = molecules[mi].getFragment(fi)->getPositiony() - molecules[m12].getFragment(f12)->getPositiony();
				double vecz12 = molecules[mi].getFragment(fi)->getPositionz() - molecules[m12].getFragment(f12)->getPositionz();

				//Update anything that goes over a box boundary
				if (vecx12 > (boxSize[0] / 2))vecx12 = vecx12 - boxSize[0]; if (vecx12 < -(boxSize[0] / 2))vecx12 = vecx12 + boxSize[0];
				if (vecy12 > (boxSize[1] / 2))vecy12 = vecy12 - boxSize[1]; if (vecy12 < -(boxSize[1] / 2))vecy12 = vecy12 + boxSize[1];
				if (vecz12 > (boxSize[2] / 2))vecz12 = vecz12 - boxSize[2]; if (vecz12 < -(boxSize[2] / 2))vecz12 = vecz12 + boxSize[2];

				bool checkj = false;
				if(fragmentj == nullptr){
					checkj = true;
				}
				else if(molecules[m12].getFragment(f12)->getType() == fragmentj->getType()){
					checkj = true;
				}
				
				if(checkj){
					for (int d = 0; d < neighboursj.size() - 1; d++){
						//Distance between atom and probe
						int m13 = neighboursm[d];
						int f13 = neighboursf[d];
						int j13 = neighboursj[d];
						double r13 = distances[i][j13];
						double a13 = 1 / pow(r13, 2);

						double vecx13 = molecules[mi].getFragment(fi)->getPositionx() - molecules[m13].getFragment(f13)->getPositionx();
						double vecy13 = molecules[mi].getFragment(fi)->getPositiony() - molecules[m13].getFragment(f13)->getPositiony();
						double vecz13 = molecules[mi].getFragment(fi)->getPositionz() - molecules[m13].getFragment(f13)->getPositionz();
						if (vecx13 >(boxSize[0] / 2))vecx13 = vecx13 - boxSize[0]; if (vecx13 < -(boxSize[0] / 2))vecx13 = vecx13 + boxSize[0];
						if (vecy13 > (boxSize[1] / 2))vecy13 = vecy13 - boxSize[1]; if (vecy13 < -(boxSize[1] / 2))vecy13 = vecy13 + boxSize[1];
						if (vecz13 > (boxSize[2] / 2))vecz13 = vecz13 - boxSize[2]; if (vecz13 < -(boxSize[2] / 2))vecz13 = vecz13 + boxSize[2];

						double r13dotr12 = (vecx12*vecx13) + (vecy12*vecy13) + (vecz12*vecz13);
						double cosTheta = (r13dotr12 / (r13*r12));

						if (a13*cosTheta > a12){
							inShell = false;
							mb = m13;
						}
					}
					
					if (inShell && mi != m12){
					//	cout << "test5 " << endl;
					//	cout << mi << " " << molecules[mi].getName() << endl;
					//	cout << mi << " " << fi << " " << molecules[mi].getFragment(fi)->getName() << endl;
					//	cout << m12 << " " << molecules[m12].getName() << endl;
					//	cout << m12 << " " << f12 << " " << molecules[m12].getFragment(f12)->getName() << endl;

						molecules[mi].getFragmentShell(fi)->addFragment(molecules[m12].getFragment(f12));
					//	cout << "test end" << endl;
					}
					if(mi == mb){
						inShell = true;
					}
				}
			}
			i = i + 1;
		}
	}

	for (int mi = 0; mi < molecules.size(); mi++){
		for (int fi = 0; fi < molecules[mi].getSize(); fi++){
			molecules[mi].getFragmentShell(fi)->orderfragmentsUsingType();
		}
	}
}//End of function

#endif
