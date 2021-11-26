/* 
**
** Output some vectors. Pretty simple.
**
*/

#ifndef OUTPUT_GRCURVE
#define OUTPUT_GRCURVE

void outputgrCurve(vector<double> vec, string vecOutFile);
void outputgrCurveWithCutoff(vector<double> vec, string vecOutFile, int cutoff);

void outputgrCurve(vector<double> vec, string vecOutFile){
	
	ofstream out;
	out.open(vecOutFile.c_str());

	for(int i = 0; i < vec.size(); i++){
		out << vec[i] << endl;
	}
}


void outputgrCurveWithCutoff(vector<double> vec, string vecOutFile, int cutoff){
	
	ofstream out;
	out.open(vecOutFile.c_str());

	for(int i = 0; i < vec.size(); i++){
		if(i <= cutoff){
			out << vec[i] << " " << vec[i] << endl;
		}
		else{
			out << vec[i] << " 0" << endl;
		}
	}
}


#endif
