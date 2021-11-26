/* 
**
** Output some vectors. Pretty simple.
**
*/

#ifndef OUTPUT_VECTOR
#define OUTPUT_VECTOR

void outputVector(vector<int> vec, string vecOutFile);
void outputVector(vector<double> vec, string vecOutFile);
void appendVector(vector<double> vec, string vecOutFile);

void outputVector(vector<int> vec, string vecOutFile){
	
	ofstream out;
	out.open(vecOutFile.c_str());

	for(int i = 0; i < vec.size(); i++){
		out << vec[i] << endl;
	}
}

void outputVector(vector<int> vec, string vecOutFile, int start, int end){
	
	ofstream out;
	out.open(vecOutFile.c_str());

	for(int i = start; i < vec.size() && i < end; i++){
		out << vec[i] << endl;
	}
}

void outputVector(vector<double> vec, string vecOutFile){
	
	ofstream out;
	out.open(vecOutFile.c_str());

	for(int i = 0; i < vec.size(); i++){
		out << vec[i] << endl;
	}
}

void outputVector(vector<double> vec, string vecOutFile, int start, int end){
	
	ofstream out;
	out.open(vecOutFile.c_str());

	for(int i = start; i < vec.size() && i < end; i++){
		out << vec[i] << endl;
	}
}

void appendVector(vector<double> vec, string vecOutFile){
	
	ofstream out;
	out.open(vecOutFile.c_str(), std::ios_base::app);

	for(int i = 0; i < vec.size(); i++){
		out << vec[i] << " ";
	}
	out << endl;
}

#endif
