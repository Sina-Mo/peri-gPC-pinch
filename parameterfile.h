#ifndef PARAMETER_FILE_INCLUDE
#define PARAMETER_FILE_INCLUDE

// Struct that holds all input parameters from a file
struct ParameterFile {

  // List expected parameters and type here
  double L1;
  int N1;
  double Let;
  double tdiameter;
  double tR2;
  double pamp;
  double freq;
  double speed; 
  double sigma;
  int numPts;
  int Nendtarget;
  int Nendspring;
  double kappa1;

	// Debugging function to print the values of parameters
  //	void print(){
	  //		cout <<"L1 = " <<L1 <<endl;
	  //	cout <<"N1 = " <<N1 <<endl;
	  //	cout <<"Let = " <<Let <<endl;
	  //	cout <<"diameter = " <<diameter <<endl;
	  //	}

	// Set parameter values given input file name
	void set(string & param_file){

	        ifstream fsin(param_file.c_str());
		string tmp_line;

		while (getline(fsin, tmp_line)){
			istringstream ssin(tmp_line.substr(tmp_line.find("=") + 1 ));
			if (tmp_line.find("L1") != -1){ ssin >> L1; }
			else if (tmp_line.find("N1") != -1){ ssin >> N1; }
			else if (tmp_line.find("Let") != -1){ ssin >> Let; }
			else if (tmp_line.find("tdiameter") != -1){ ssin >> tdiameter; }
			else if (tmp_line.find("tR2") != -1){ ssin >> tR2; }
			else if (tmp_line.find("pamp") != -1){ ssin >> pamp; }
			else if (tmp_line.find("freq") != -1){ ssin >> freq; }
			else if (tmp_line.find("speed") != -1){ ssin >> speed; }
			else if (tmp_line.find("sigma") != -1){ ssin >> sigma; }
			else if (tmp_line.find("numPts") != -1){ ssin >> numPts; }
			else if (tmp_line.find("Nendtarget") != -1){ ssin >> Nendtarget; }
			else if (tmp_line.find("Nendspring") != -1){ ssin >> Nendspring; }
			else if (tmp_line.find("kappa1") != -1){ ssin >> kappa1; }
			else {
				cout <<"Parameter '" <<tmp_line <<"' not found!" <<endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	// Declare constructor and descructor
	//ParameterFile();
	//~ParameterFile();
};
#endif //ifndef PARAMETER_FILE_INCLUDE
