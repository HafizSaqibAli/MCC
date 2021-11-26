/* g(r) cutoff Code, Dr Jonathan Higham
** Modified 15/06/2016
**
** This function acts on a vector of fragment 
** objects to calculate the radial distribution
** function of those objects
**
*/


#ifndef NR_CALCULATION
#define NR_CALCULATION

using namespace std;

void nrUpdate(vector<fragment> fragments, vector <int> &nr, float* boxSize, double dr, int maxBin, fragment* fragmenti, fragment* fragmentj);

void nrUpdate(vector<fragment> fragments, vector <int> &nr, float* boxSize, double dr, int maxBin, fragment* fragmenti, fragment* fragmentj){

	//Loop over each fragment
	for(int i = 0; i < fragments.size(); i++){
		//Loop over each comparison fragment
		for(int j = 0; j < fragments.size(); j++){

			if(i != j){

				bool checki = false;
				bool checkj = false;

				if(!fragmenti){
					checki = true;
				}
				else if(fragments[i].getType() == fragmenti->getType()){
					checki = true;
				}
				
				if(!fragmentj){
					checkj = true;
				}
				else if(fragments[j].getType() == fragmentj->getType()){
					checkj = true;
				}

				if(checki && checkj){
					//Add the squared distance in x,y,z
					double vecx = fragments[i].getPositionx() - fragments[j].getPositionx();
					double vecy = fragments[i].getPositiony() - fragments[j].getPositiony();
					double vecz = fragments[i].getPositionz() - fragments[j].getPositionz();

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
	}//End of main fragment loop
}

#endif
