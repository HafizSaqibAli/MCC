/* g(r) cutoff Code, Dr Jonathan Higham
** Modified 15/06/2016
**
** This function acts on a vector of fragment 
** objects to calculate the radial distribution
** function of those objects
**
*/


#ifndef SMOOTH_CURVE
#define SMOOTH_CURVE

using namespace std;

void smoothCurve(vector <double> &curve);

void smoothCurve(vector <double> &curve){

	for(int i = 1; i < curve.size()-1; i++){
		//Normalise for frames, fragments and double counting
		curve[i] = (curve[i-1] + curve[i] + curve[i+1]) / 3;
	}
}

#endif
