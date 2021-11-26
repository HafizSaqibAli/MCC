/* g(r) cutoff Code, Dr Jonathan Higham
** Modified 15/06/2016
**
** This function acts on a vector of fragment 
** objects to calculate the radial distribution
** function of those objects
**
*/


#ifndef GR_CALCULATION
#define GR_CALCULATION

using namespace std;

vector <double> grCalculation(vector <int> nr, int totalFragments, int tMax, float* boxSize, double dr, int maxBin);
int getFirstMinima(vector <double> gr);

vector <double> grCalculation(vector <int> nr, int totalFragments, int tMax, float* boxSize, double dr, int maxBin){

	//Some constants for normalising the g(r) curve
	const double pi = 3.14159;
	const double p = (totalFragments/(boxSize[0]*boxSize[1]*boxSize[2]));

	vector <double> gr(maxBin, 0.0);

	for(int i = 1; i < maxBin; i++){
		//Normalise for frames, fragments and double counting
		gr[i] = (double)nr[i]/((double)tMax*totalFragments);
		//Normalise by average density
		double r = (double)i * dr;
		gr[i] = gr[i]/(4*pi*r*r*dr*p);
	}
	
	return gr;
}

int getFirstMinima(vector <double> gr){

	//Loop over g(r) to find the first minima
	int grMinBin = 0; 
	bool grFlag = true;
	int gri = 2;

	while(gr[gri] <= 0){
		gri++;
	}

	for(int i = gri; i < gr.size()-1; i++){
		if(gr[i-1] >= gr[i] && gr[i-1] >= gr[i] && gr[i] <= gr[i+1] && grFlag){
			grMinBin = i; 
			grFlag = false;
		}
	}

	return grMinBin;
}

#endif
