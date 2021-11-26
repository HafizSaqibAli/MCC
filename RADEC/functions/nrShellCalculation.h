/* g(r) cutoff Code, Dr Jonathan Higham
** Modified 15/06/2016
**
** This function acts on a vector of fragment 
** objects to calculate the radial distribution
** function of those objects
**
*/


#ifndef NR_SHELL_CALCULATION
#define NR_SHELL_CALCULATION

using namespace std;

void nrShellUpdate(vector<molecule> molecules, vector <int> &nr, float* boxSize, double dr, int maxBin);
void nrShellUpdate(vector<molecule> molecules, vector <int> &nr, float* boxSize, double dr, int maxBin, fragment* fragmenti, fragment* fragmentj);

void nrShellUpdate(vector<molecule> molecules, vector <int> &nr, float* boxSize, double dr, int maxBin){
	nrShellUpdate(molecules, nr, boxSize, dr, maxBin, nullptr, nullptr);
}

void nrShellUpdate(vector<molecule> molecules, vector <int> &nr, float* boxSize, double dr, int maxBin, fragment* fragmenti, fragment* fragmentj){

	//Loop over each molecule
	for(int mi = 0; mi < molecules.size(); mi++){
		for(int fi = 0; fi < molecules[mi].getSize(); fi++){
			bool checki = false;
			if(!fragmenti){
				checki = true;
			}
			else if(molecules[mi].getFragment(fi)->getType() == fragmenti->getType()){
				checki = true;
			}

			if(checki){
				//Loop over each shell molecule
				for(int s = 0; s < molecules[mi].getShell(fi)->getSize(); s++){

					bool checkj = false;
					if(!fragmentj){
						checkj = true;
					}
					else if(molecules[mi].getShell(fi)->getFragment(s)->getType() == fragmentj->getType()){
						checkj = true;
					}

					if(checkj){
						//Add the squared distance in x,y,z
						double vecx = molecules[mi].getFragment(fi)->getPositionx() - molecules[mi].getShell(fi)->getFragment(s)->getPositionx();
						double vecy = molecules[mi].getFragment(fi)->getPositiony() - molecules[mi].getShell(fi)->getFragment(s)->getPositiony();
						double vecz = molecules[mi].getFragment(fi)->getPositionz() - molecules[mi].getShell(fi)->getFragment(s)->getPositionz();

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

						//Calculate distance
						double r = sqrt((vecx*vecx) + (vecy*vecy) + (vecz*vecz));
						//Find bin
						int bin = ceil(r/dr);
						//Increment bin
						if(bin < nr.size()){
							nr[bin]++;
						}
					}
				}
			}//End of comparison fragment loop
		}
	}//End of main fragment loop
}

#endif
