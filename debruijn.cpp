/************************************************************
 * Program that generates de Bruijn sequences. 
 * This program is based on Hamiltonian Cycle Program written by Basil Vandegriend 
 * For more details, please refer to:
 * http://webdocs.cs.ualberta.ca/~joe/Theses/vandegriend.html
 * 
 * Modified by Dongbo Hu
 * dongbo@mail.med.upenn.edu
 *
 * Modified by Marcelo Mattar (09/24/2010 - present)
 * mattar@sas.upenn.edu
 * Changed pieces are commented with <MMattar> and </MMattar>
 *
 *
     ``THIS SOURCE CODE IS SUPPLIED  ``AS IS'' WITHOUT WARRANTY OF ANY KIND,
     AND ITS AUTHOR AND THE JOURNAL OF ARTIFICIAL INTELLIGENCE RESEARCH
     (JAIR) AND JAIR'S PUBLISHERS AND DISTRIBUTORS, DISCLAIM ANY AND ALL
     WARRANTIES, INCLUDING BUT NOT LIMITED TO ANY IMPLIED
     WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, AND
     ANY WARRANTIES OR NON INFRINGEMENT.  THE USER ASSUMES ALL LIABILITY AND
     RESPONSIBILITY FOR USE OF THIS SOURCE CODE, AND NEITHER THE AUTHOR NOR
     JAIR, NOR JAIR'S PUBLISHERS AND DISTRIBUTORS, WILL BE LIABLE FOR
     DAMAGES OF ANY KIND RESULTING FROM ITS USE.  Without limiting the
     generality of the foregoing, neither the author, nor JAIR, nor JAIR's
     publishers and distributors, warrant that the Source Code will be
     error-free, will operate without interruption, or will meet the needs
     of the user.''


        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions
        are met:
        1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        3. All advertising materials mentioning features or use of this
        software must display the following acknowledgement:
        "This product includes software developed by the University of
        Alberta, Edmonton."
        4. Neither the name of the University nor the names of its
	contributors may be used to endorse or promote products derived 
	from this software without specific prior written permission.

        THIS SOFTWARE IS PROVIDED BY THE CONTRIBUTORS ``AS IS'' AND ANY
        EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
        THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
        PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
        CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
        SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
        NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
        LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
        HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
        CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
        OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
        EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

        THIS SOFTWARE IS SUPPLIED WITHOUT ANY SUPPORT SERVICES.

 ************************************************************/

using namespace std;

#include <iostream>
#include <set>
#include "debruijn.h"
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fftw3.h>

SeqOpt myOpt;

struct timeval g_algstart;
int g_hit_timelimit = RUN_NORMAL;


//<MMattar> Global variables for printing some data important for debugging
int debugmode = 0;
int evalmode = 0;
int SOA = 1;
int print_flag = 1;
int starttime = -1;
//</MMattar>


void BackTrackOpt::init()
{
  initvertflag = INITVERT_RANDOM;
  degsortflag  = DEGSORT_RAND;
  pruneoptflag = HC_PRUNE_NONE;

  restart_increment = 0;
  max_nodes = 0;
}

void GraphGenOpt::init() 
{
  nvertex = 50;
  makeham = HAM_DONTCARE;

  dist = 0.255;
  dim = 2;
  dflag = CONNECT_NEAR;
  wrapflag = GRAPH_WRAP;
  mindeg = 0;

  board1 = 8;
  board2 = 8;
  move1 = 1;
  move2 = 2;

  numsubgraphs = 0;

  degconst = 0.0;
  meandeg = 3.0;

  numcycles = 2.0;

  for (int loop = 0; loop < MAXNUMADDPATHS; loop++)
    pathlengths[loop] = 0.0;

  indsetsize = 0;

  /* leave these at 0 to make parse code work properly */
  degsize = 2;
  degpercent[2] = 0.0;
  degpercent[3] = 0.0;
  degpercent[4] = 0.0;
  degpercent[5] = 0.0;
}

void HeuristicOpt::init()
{
  completeflag = COMPLETE_NORM;
  visitflag = VISIT_RAND;
  cycleextendflag = NOCYCLEEXTEND;
}

void SeqOpt::init()
{
  algorithm = ALG_NOSOLVE;
  alg_timelimit = -1;		/* no timelimit */
  graphgentype = GEN_NOGRAPH;

  /* set broadest display options */
  report_flags = (REPORT_NONE);
}

/* this function initializes/emptys an allocated graph structure
 * (creates a completely unconnected graph) */
void Graph::init()
{
  numvert = 0;
  numedges = 0;
  meandeg = 0.0;
  stddevdeg = 0.0;
  mindeg = 0;
  maxdeg = 0;

  solve = HC_NOT_FOUND;
}

// set up graph data based on input n and k
void Graph::setup(int num_label, int word_len)
{  
	  
  // set up number of vertices
  unsigned tmp = 1;
  for (int i = 0; i < word_len; i++) 
    tmp *= num_label; // tmp = num_label ^ word_len

  numvert = tmp;
  // set up degree
  degree.reserve(numvert); // reserve memory space of size=numvert for the vector degree
  for (int i = 0; i < numvert; i++) 
    degree.push_back(num_label); // sets every entry of the degree vector to be num_label, which means that every node connects to other num_label nodes

  // set up neighbors of each vertex
  nbr.reserve(numvert); // reserve memory space of size=numvert for the vector nbr
  vector<int> tmpVec;
  // keep descending order for the first vertex
  for (int j = num_label - 1; j >= 0; j--) 
    tmpVec.push_back(numvert / num_label * j); // tmpVec entries decays in num_label^(word_len-1) decrements until zero
  nbr.push_back(tmpVec);
  // Add neighbors of other vertexes randomly
  vector<int> j_indx;
  for (int i = 0; i < num_label; ++i) 
    j_indx.push_back(i); // sequence from 0 to num_label

  for (int i = 1; i < numvert; i++) {
    random_shuffle (j_indx.begin(), j_indx.end()); //shuffles j_indx
    for (unsigned m = 0; m < j_indx.size(); ++m)
      tmpVec[m] = numvert / num_label * j_indx[m] + i / num_label;
    nbr.push_back(tmpVec);
  }  
      
  /* initialize variables */
  mindeg = num_label;
  maxdeg = num_label;
  numedges = num_label * numvert / 2; // number of edges in the graph

  for (int loop = 0; loop < num_label; loop++)
    deghistogram.push_back(0);

  /* calculate mean, stddev of vertex degree */
  meandeg = num_label; 
  stddevdeg = 0;
 
  /* calculate min, max degrees and total number of edges 
   * create vertex degree histogram
   */
  for (int loop = 0; loop < numvert; loop++) {
    deghistogram[degree[loop]]++;
  }

}

// Helper function to check the neighbors
void Graph::printNeighbors()
{
  for (unsigned i = 0; i < nbr.size(); i++) {
    cout << i << ": ";
    for (unsigned j = 0; j < nbr[i].size(); j++)
      cout << nbr[i][j] << " ";
    cout << endl;
  }
  cout << endl;
}


//<MMattar> Created the function
// Initialize with zero some relevant variables
void DBSeq::init()
{
	num_label = 0;
	word_len = 0;
	num_bins = 0;
}
//</MMattar>


//<MMattar> Created the whole function
// this function initializes/emptys an allocated DeBruijn sequence structure
void DBSeq::setup(int number_of_labels, int word_length, int number_of_bins, char * NeuralModel1_filename, char * NeuralModel2_filename, char * NeuralModel3_filename, char * GuideFunction_filename)
{
	num_label = number_of_labels;
	word_len = word_length;
	num_bins = number_of_bins;
	num_neuralmodels = 0;

	FILE * nm1file;
	FILE * nm2file;
	FILE * nm3file;
	FILE * gffile;

	// initializing the guide function and seq_bin
	guide_function.reserve(pow(num_label,word_len));
	seq_bin.reserve(pow(num_label,word_len));
	for (int i=0; i<pow(num_label,word_len); i++) seq_bin.push_back(0);

	// initializing the trans_seq and bins_used vectors (to save the relevant information for creating the sequence)
	trans_seq.reserve(pow(num_label,word_len));
	trans_seq2.reserve(pow(num_label,word_len));
	trans_seq3.reserve(pow(num_label,word_len));
	bins_used.reserve(pow(num_label,word_len));

	// initializing the neural model and bin matrices with zeros
	vector<double> auxvect; // auxiliary vector for initializing the neural model table
	auxvect.reserve(num_label);
	for (int i=0; i<num_label; i++) auxvect.push_back(0.0); //fill auxvect with zeros
	nm1.reserve(num_label);
	nm2.reserve(num_label);
	nm3.reserve(num_label);
	bin.reserve(num_label);
	for (int i=0; i<num_label; i++){
		nm1.push_back(auxvect); // initialize the neural model table, using auxvect (=0) num_label times
		nm2.push_back(auxvect); // initialize the neural model table, using auxvect (=0) num_label times
		nm3.push_back(auxvect); // initialize the neural model table, using auxvect (=0) num_label times
		bin.push_back(auxvect); // initialize the bins table, using auxvect (=0) num_label times
	}

	// initializing the pref_bin matrix with zeros
	vector<int> auxvect2;
	auxvect2.reserve(2);
	auxvect2.push_back(0);
	auxvect2.push_back(0);
	vector<vector<int> > auxvect3;
	auxvect3.reserve(num_bins);
	for (int i=0; i<num_bins; i++) auxvect3.push_back(auxvect2);
	pref_bin.reserve(pow(num_label,word_len));
	for (int t=0; t<pow(num_label,word_len); t++) pref_bin.push_back(auxvect3);

	// if the user entered a nm matrix path, then read the distances
	if (NeuralModel1_filename != NULL){
		nm1file = fopen (NeuralModel1_filename , "r");
		num_neuralmodels = num_neuralmodels + 1;
		// read the distances from nmfile
		for (int i=0; i<num_label; i++){
			for (int j=0; j<num_label; j++){
				fscanf (nm1file, "%lf", &nm1[i][j]); // read the distances from external file
				if (feof(nm1file) && (i+j != num_label+num_label-2)) { // if at the end of the file, check if the number of elements is less than num_label^2
					cout << endl << "The neural model file must have a " << num_label << "x" << num_label << " matrix of floating points." << endl << endl;
					exit(1);
				}
			}
		}
		// if at the end of the file, now check if the number of elements is greater than num_label^2
		if (!feof(nm1file)) {
			if (fgetc(nm1file) != '\n'){
				cout << endl << "The neural model file must have a " << num_label << "x" << num_label << " matrix of floating points." << endl << endl;
				exit(1);
			}
		}

		// look for the greatest distance in the matrix, in order to normalize the distances
		double greatest_distance = -1.0;
		for (int i=0; i<num_label; i++){
			for (int j=0; j<num_label; j++){
				if(nm1[i][j] < 0) nm1[i][j] = -1.0; // transform all NA into -1
				if(greatest_distance < nm1[i][j]) greatest_distance = nm1[i][j]; // checks if nm1[i][j] is the greatest value so far
			}
		}
		// greatest_distance will now have the greatest value of the neural model matrix (not considering NAs), to be used in the normalization

		// normalize with greatest_distance in case it is greater than 0
		if(greatest_distance > 0.00000001){ // trick to avoid division by zero
			for (int i=0; i<num_label; i++){
				for (int j=0; j<num_label; j++){
					if(nm1[i][j]>0.00000001) nm1[i][j] /= greatest_distance; // note that we only normalize non-NA distances
				}
			}
		}
		// if greatest_distance < 0, it means that all distances are NA. Then we already have nm filled with -1's
		// at this point, nm is a normalized matrix of distances, with the maximum distance being 1 and all NA distances being -1
	}
	if (NeuralModel2_filename != NULL){
		nm2file = fopen (NeuralModel2_filename , "r");
		num_neuralmodels = num_neuralmodels + 1;

		// read the distances from nmfile
		for (int i=0; i<num_label; i++){
			for (int j=0; j<num_label; j++){
				fscanf (nm2file, "%lf", &nm2[i][j]); // read the distances from external file
				if (feof(nm2file) && (i+j != num_label+num_label-2)) { // if at the end of the file, check if the number of elements is less than num_label^2
					cout << endl << "The neural model file must have a " << num_label << "x" << num_label << " matrix of floating points." << endl << endl;
					exit(1);
				}
			}
		}
		// if at the end of the file, now check if the number of elements is greater than num_label^2
		if (!feof(nm2file)) {
			if (fgetc(nm2file) != '\n'){
				cout << endl << "The neural model file must have a " << num_label << "x" << num_label << " matrix of floating points." << endl << endl;
				exit(1);
			}
		}

		// look for the greatest distance in the matrix, in order to normalize the distances
		double greatest_distance2 = -1.0;
		for (int i=0; i<num_label; i++){
			for (int j=0; j<num_label; j++){
				if(nm2[i][j] < 0) nm2[i][j] = -1.0; // transform all NA into -1
				if(greatest_distance2 < nm2[i][j]) greatest_distance2 = nm2[i][j]; // checks if nm[i][j] is the greatest value so far
			}
		}
		// greatest_distance will now have the greatest value of the neural model matrix (not considering NAs), to be used in the normalization

		// normalize with greatest_distance in case it is greater than 0
		if(greatest_distance2 > 0.00000001){ // trick to avoid division by zero
			for (int i=0; i<num_label; i++){
				for (int j=0; j<num_label; j++){
					if(nm2[i][j]>0.00000001) nm2[i][j] /= greatest_distance2; // note that we only normalize non-NA distances
				}
			}
		}
		// if greatest_distance < 0, it means that all distances are NA. Then we already have nm filled with -1's
		// at this point, nm is a normalized matrix of distances, with the maximum distance being 1 and all NA distances being -1
	}
	if (NeuralModel3_filename != NULL){
		nm3file = fopen (NeuralModel3_filename , "r");
		num_neuralmodels = num_neuralmodels + 1;

		// read the distances from nmfile
		for (int i=0; i<num_label; i++){
			for (int j=0; j<num_label; j++){
				fscanf (nm3file, "%lf", &nm3[i][j]); // read the distances from external file
				if (feof(nm3file) && (i+j != num_label+num_label-2)) { // if at the end of the file, check if the number of elements is less than num_label^2
					cout << endl << "The neural model file must have a " << num_label << "x" << num_label << " matrix of floating points." << endl << endl;
					exit(1);
				}
			}
		}
		// if at the end of the file, now check if the number of elements is greater than num_label^2
		if (!feof(nm3file)) {
			if (fgetc(nm3file) != '\n'){
				cout << endl << "The neural model file must have a " << num_label << "x" << num_label << " matrix of floating points." << endl << endl;
				exit(1);
			}
		}

		// look for the greatest distance in the matrix, in order to normalize the distances
		double greatest_distance3 = -1.0;
		for (int i=0; i<num_label; i++){
			for (int j=0; j<num_label; j++){
				if(nm3[i][j] < 0) nm3[i][j] = -1.0; // transform all NA into -1
				if(greatest_distance3 < nm3[i][j]) greatest_distance3 = nm3[i][j]; // checks if nm[i][j] is the greatest value so far
			}
		}
		// greatest_distance will now have the greatest value of the neural model matrix (not considering NAs), to be used in the normalization

		// normalize with greatest_distance in case it is greater than 0
		if(greatest_distance3 > 0.00000001){ // trick to avoid division by zero
			for (int i=0; i<num_label; i++){
				for (int j=0; j<num_label; j++){
					if(nm3[i][j]>0.00000001) nm3[i][j] /= greatest_distance3; // note that we only normalize non-NA distances
				}
			}
		}
		// if greatest_distance < 0, it means that all distances are NA. Then we already have nm filled with -1's
		// at this point, nm is a normalized matrix of distances, with the maximum distance being 1 and all NA distances being -1
	}

	// if the user entered a guide function file path, then read the function
	if (GuideFunction_filename != NULL){
		gffile = fopen (GuideFunction_filename , "r");

		if (gffile == NULL) { // if the user entered periods and not filename
			// first, see if an HRF guide function was specified
			if (GuideFunction_filename[0]=='H' && GuideFunction_filename[1]=='R' && GuideFunction_filename[2]=='F'){
				
				if (SOA == 1) { // if at the end of the file, check if the number of elements is less than num_label^2
					cout << endl << "You must specify an SOA in order to be able to use an HRF guide function" << endl << endl;
					exit(1);
				}
				
				// Set up relevant variables for the Fourier Transform
				long seq_length = pow(num_label,word_len);
				const int N = seq_length;
				fftw_complex *in, *out;
				fftw_plan p;
				fftw_plan p2;
				in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
				out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
				
				cout << "The SOA used here is: " << SOA << endl << endl << endl;
				
				// Hard-code the EndogenousFilter
				double EndogenousFilter_array[160] = {5.55112e-17,0.00767224,0.0162726,0.025849,0.0364457,0.0481024,0.0608536,0.0747281,0.0897482,0.10593,0.123281,0.141801,0.161484,0.182312,0.204261,0.227297,0.251376,0.276448,0.30245,0.329315,0.356963,0.38531,0.414261,0.443716,0.473569,0.503705,0.534007,0.564352,0.594614,0.624665,0.654374,0.68361,0.712242,0.740139,0.767176,0.793226,0.818169,0.841889,0.864277,0.885229,0.90465,0.922451,0.938554,0.952889,0.965396,0.976025,0.984737,0.991502,0.996304,0.999136,1,0.998913,0.995897,0.990989,0.984232,0.975681,0.965398,0.953452,0.939923,0.924892,0.908452,0.890696,0.871725,0.851641,0.830549,0.808556,0.78577,0.7623,0.738253,0.713734,0.688848,0.663695,0.638375,0.61298,0.5876,0.562319,0.537218,0.51237,0.487844,0.463702,0.44,0.41679,0.394115,0.372016,0.350523,0.329666,0.309467,0.289942,0.271104,0.25296,0.235516,0.218771,0.202721,0.187361,0.172683,0.158674,0.145324,0.132616,0.120536,0.109067,0.0981926,0.0878945,0.0781552,0.0689569,0.0602817,0.0521123,0.0444312,0.0372216,0.0304668,0.0241507,0.0182571,0.0127706,0.00767572,0.00295715,-0.00140033,-0.00541187,-0.00909271,-0.0124583,-0.0155242,-0.0183063,-0.0208208,-0.0230841,-0.025113,-0.0269245,-0.0285358,-0.0299644,-0.0312277,-0.0323428,-0.0333269,-0.0341965,-0.0349674,-0.0356548,-0.0362724,-0.0368327,-0.0373465,-0.0378227,-0.0382678,-0.038686,-0.0390787,-0.0394442,-0.0397777,-0.0400708,-0.0403115,-0.0404842,-0.0405693,-0.0405433,-0.0403786,-0.0400439,-0.039504,-0.0387199,-0.0376493,-0.0362464,-0.0344627,-0.0322468,-0.0295455,-0.0263035,-0.0224645,-0.0179715,-0.0127673,-0.00679528};
				
				// Upsampling the EndogenousFilter using interpolation
				double EndogenousFilter_upsampled[320];
				for (int i=0; i<320; i++){
					if ((i%2) == 0){
						EndogenousFilter_upsampled[i] = EndogenousFilter_array[i/2];
					}
					else{
					  // fix this so that at i == 319, this doesn't read EndogenousFilter_array[160]
						EndogenousFilter_upsampled[i] = EndogenousFilter_array[(int) i/2] + (((EndogenousFilter_array[(int) i/2 + 1] - EndogenousFilter_array[(int) i/2])/2)*(i%2));
					}
				}
				
				// Downsampling the EndogenousFilter, so that is has a temporal resolution equal to SOA
				double EndogenousFilter_downsampled[seq_length];
				for (int i=0; i<seq_length; i++){
					if(i*(SOA/50) < 320){
						EndogenousFilter_downsampled[i] = EndogenousFilter_upsampled[i*(SOA/50)];
					}
					else{
						EndogenousFilter_downsampled[i] = 0;
					}
				}

				/*
				//[DEBUG]
				cout << endl << "ENDOGENOUS FILTER (downsampled):" << endl;
				for (int t=0; t<seq_length; t++){
					cout << EndogenousFilter_downsampled[t] << " ";
				}
				cout << endl << endl;
				//exit(1);
				//[/DEBUG]
				*/
				
				// copy EndogenousFilter to a format suitable for the fourier transform
				for (int i=0; i<N; i++){
					in[i][0] = EndogenousFilter_downsampled[i]; //real part
					in[i][1] = 0; // complex part
				}
				
				// Calculate the FFT of EndogenousFilter
				p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
				fftw_execute(p);
				
				// Store the results on a vector
				vector<double> EndogenousKernel (2*N,0);
				for (int i=0; i<N; i++){
					EndogenousKernel[2*i] = out[i][0]; //real part
					EndogenousKernel[2*i + 1] = out[i][1]; // complex part
				}
				
				/*
				//[DEBUG]
					cout << endl << "ENDOGENOUS KERNEL:" << endl;
					for (int t=0; t<2*N; t++){
						cout << EndogenousKernel[t] << " ";
					}
					cout << endl << endl;
					exit(1);
				//[/DEBUG]
				*/

				// In the frequency domain, calculate the maximum norm and normalize EndogenousKernel
				double maxnorm = 0;
				double currentnorm = 0;
				for (int t=0; t<N*2; t += 2){
					currentnorm = sqrt(EndogenousKernel[t]*EndogenousKernel[t] + EndogenousKernel[t+1]*EndogenousKernel[t+1]);
					if (maxnorm<currentnorm) maxnorm=currentnorm;
				}
				for (int t=0; t<N*2; t++){
					EndogenousKernel[t] = EndogenousKernel[t]/maxnorm;
				}
				
				// In the frequency domain, apply the Notch filter and make the first element equals to 1
				float FilterCutoff = 0.01;
				if ((N*FilterCutoff*SOA/1000) >= 1){
					for (int t=0; t<(N*FilterCutoff*SOA/1000); t++){ // !!! Problem with small sequences
						EndogenousKernel[2*t] = 0;   // real part
						EndogenousKernel[2*t+1] = 0; // imaginary part
						EndogenousKernel[2*N-2-t] = 0;   // real part
						EndogenousKernel[2*N-1-t] = 0; // imaginary part
					}
				}

				/*
				//[DEBUG]
					cout << endl << "ENDOGENOUS KERNEL (normalized and filtered):" << endl;
					for (int t=0; t<N*2; t++){
						cout << EndogenousKernel[t] << " ";
					}
					cout << endl << endl;
					exit(1);
				//[/DEBUG]
				*/

				// Randomize the phases of the frequency representation
				vector<double> GuideFunction_fft (N*2,0);
				//double GuideFunction_fft[N*2];
				float Norm, Phase;
				int halfN = floor(N/2);
				if (halfN*2 < N) halfN++;
				
				for (int i=0; i<halfN; i++){
					Norm = sqrt(EndogenousKernel[2*i]*EndogenousKernel[2*i] + EndogenousKernel[2*i+1]*EndogenousKernel[2*i+1]);
					Phase = ((float) rand()/RAND_MAX)*2*PI;
					GuideFunction_fft[2*i] = Norm*cos(Phase); // real part
					GuideFunction_fft[2*i+1] = Norm*sin(Phase); // imaginary part
				}
				// If N is even, the value at N/2 + 1 must have zero phase, so replace it with the original data value
				if (N%2 == 0){
					GuideFunction_fft[2*halfN] = EndogenousKernel[2*(N/2)]; // real part of central element
					GuideFunction_fft[2*halfN+1] = 0; // imaginary part of central element
					// second half of FT is complex conj of first half
					for (int i=(halfN+1); i<N; i++){
						GuideFunction_fft[2*i] = GuideFunction_fft[2*(N-i)]; // real part
						GuideFunction_fft[2*i+1] = -GuideFunction_fft[2*(N-i)+1]; // imaginary part
					}
				}
				else{
					for (int i=halfN; i<N; i++){
						GuideFunction_fft[2*i] = GuideFunction_fft[2*(N-i)]; // real part
						GuideFunction_fft[2*i+1] = -GuideFunction_fft[2*(N-i)+1]; // imaginary part
					}
				}
				GuideFunction_fft[0] = 0;
				GuideFunction_fft[1] = 0;
				
				/*
				//[DEBUG]
				cout << endl << endl << "GuideFunction_fft = [";
				for (int t=0; t<N*2; t++){
					cout << GuideFunction_fft[t] << " ";
				}
				cout << "]" << endl << endl;
				//exit(1);
				//[/DEBUG]
				*/
				
				// copy GuideFunction_fft to a format suitable for the fourier transform
				for (int i=0; i<N; i++){
					in[i][0] = GuideFunction_fft[2*i]; //real part
					in[i][1] = GuideFunction_fft[2*i + 1]; // complex part
				}	
				
				// Calculate the iFFT of GuideFunction_fft
				p2 = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
				fftw_execute(p2);
				
				for (int i=0; i<N; i++){
					guide_function[i] = out[i][0];
				}
				
				/*			
				//[DEBUG]
				cout << endl << endl << "GuideFunction = [";
				for (int t=0; t<N; t++){
					cout << guide_function[t] << " ";
				}
				cout << "]" << endl << endl;
				//exit(1);
				//[/DEBUG]
				*/

				fftw_destroy_plan(p);
				fftw_destroy_plan(p2);
				fftw_free(in); fftw_free(out);
				
			}
			// otherwise, check if user entered a range of periods for the guide function
			else{
				int Tmin = 0;
				int Tmax = 0;
				int i = 1;
				while(GuideFunction_filename[i]!=','){
					Tmin = 10*Tmin + (GuideFunction_filename[i] - '0');
					i++;
				}
				i++;
				while(GuideFunction_filename[i]!=']'){
					Tmax = 10*Tmax + (GuideFunction_filename[i] - '0');
					i++;
				}

				// Check if Tmin is really less than Tmax, otherwise, invert order
				if (Tmin>Tmax){
					int aux = Tmin;
					Tmin = Tmax;
					Tmax = aux;
				}

				for (int T=Tmin; T<=Tmax; T++){
					int randphase = (((rand()%T)+T)%T);
					for (int t=0; t<pow(num_label,word_len); t++){
						guide_function[t] = guide_function[t] + sin((2*PI/T) * (t + randphase));
					}
				}
				/*
				cout << endl << endl;
				for (int t=0; t<pow(num_label,word_len); t++){
					cout << guide_function[t] << " ";
				}
				cout << endl << endl << endl;
				exit(1);
				*/
			}
		}
		else{ // if the user entered a filename and not periods
			// read the guide_function from gffile
			for (int i=0; i<pow(num_label,word_len); i++){
				fscanf (gffile, "%lf", &guide_function[i]); // read the guide function from external file
				if (feof(gffile) && (i != pow(num_label,word_len)-1)) { // if at the end of the file, check if the number of elements is less than num_label^2
					cout << endl << "The guide function file must have " << num_label << "^" << word_len << " space separated floating points." << endl << endl;
					exit(1);
				}
			}
			// if at the end of the file, now check if the number of elements is greater than num_label^2
			if (!feof(gffile)) {
				cout << endl << "The guide function file must have " << num_label << "^" << word_len << " space separated floating points." << endl << endl;
				exit(1);
			}
		}
		// normalization
		double min_value = 1e10;
		double max_value = -1e10;
		// search for the minimum value
		for (int i=0; i<pow(num_label,word_len); i++){
			if (min_value > guide_function[i]) min_value = guide_function[i];
		}
		// shift the minimum value to 0
		for (int i=0; i<pow(num_label,word_len); i++){
			guide_function[i] = guide_function[i] - min_value;
		}
		// search for the maximum value (necessarily >= 0)
		for (int i=0; i<pow(num_label,word_len); i++){
			if (max_value < guide_function[i]) max_value = guide_function[i];
		}
		// normalize such that the maximum value is now in 1
		if (max_value != 0) {
			for (int i=0; i<pow(num_label,word_len); i++){
				guide_function[i] = guide_function[i]/max_value;
			}
		}
		
		/*
		cout << endl << endl;
		for (int t=0; t<pow(num_label,word_len); t++){
			cout << guide_function[t] << " ";
		}
		cout << endl << endl << endl;
		exit(1);
		*/
	}
}
//</MMattar>


//<MMattar> Created the whole function
// allocates all distances in bins according to the period chosen, and creates a vector with the sequence of bins to take the distances from
void binning(int num_label, int word_len, int num_bins, DBSeq* myDBSeq){

	long seq_length = pow(num_label,word_len);

	// create the sequence of bins vector (by binning guide_function)
	// think about how to increase performance of this function - TOO SLOW!
	double min_value;
	int min_index;
	int elem_per_bin = floor(seq_length/num_bins);
	int rest = seq_length - elem_per_bin*num_bins; // k^n = elem_per_bin*num_bins + rest

	for (int b=1; b<=num_bins; b++){
		// for the first "rest" bins, include an extra element
		if (b<=rest) elem_per_bin = ceil(seq_length/num_bins);
		else elem_per_bin = floor(seq_length/num_bins);
		for (int i=0; i<elem_per_bin; i++){
			min_value = 10; // any value on guide_function is less than 10 (actually, is between 0 and 1)
			min_index = -1; // any index will be positive
			for (int j=0; j<seq_length; j++){
				if((min_value > myDBSeq->guide_function[j]) && (myDBSeq->seq_bin[j] == 0)){
					min_value = myDBSeq->guide_function[j]; // looks for the minimum value
					min_index = j; // stores the index of the minimum value
				}
			}
			myDBSeq->seq_bin[min_index] = b; // saves the bin of that number in seq_bin
		}
	}

	// create the pref_bin 3-dimensional matrix (k^n x B x 2)
	for (int t=0; t<seq_length; t++){
		myDBSeq->pref_bin[t][0][0] = myDBSeq->seq_bin[t];
		myDBSeq->pref_bin[t][0][1] = myDBSeq->seq_bin[t];
		for (int b=1; b<num_bins; b++){
			myDBSeq->pref_bin[t][b][0] = myDBSeq->seq_bin[t] + b;
			myDBSeq->pref_bin[t][b][1] = myDBSeq->seq_bin[t] - b;
			if ((myDBSeq->seq_bin[t]+b) > num_bins){
				if ((myDBSeq->seq_bin[t]-b) <= 0){
					myDBSeq->pref_bin[t][b][0] = 0;
					myDBSeq->pref_bin[t][b][1] = 0;
				}
				else myDBSeq->pref_bin[t][b][0] = myDBSeq->pref_bin[t][b][1];
			}
			if ((myDBSeq->seq_bin[t]-b) <= 0){
				if ((myDBSeq->seq_bin[t]+b) > num_bins){
					myDBSeq->pref_bin[t][b][0] = 0;
					myDBSeq->pref_bin[t][b][1] = 0;
				}
				else myDBSeq->pref_bin[t][b][1] = myDBSeq->pref_bin[t][b][0];
			}
		}
	}

	// allocate each distance from the nm matrix into bins
	// initialize the auxiliar vector bin_alloc
	vector<int> bin_alloc(num_bins,0);
	for (int i=0; i<(num_label*num_label); i++){
		bin_alloc[i%num_bins] ++;
	}

	double min_dist = 100;
	int min_index_i = -1;
	int min_index_j = -1;
	int random_bin = -1;

	// first, allocate NA distances in random bins
	for (int i=0; i<num_label; i++){
		for (int j=0; j<num_label; j++){
			if((myDBSeq->nm1[i][j] == -1) && (myDBSeq->bin[i][j] == 0)){
				random_bin = (((rand()%num_bins)+num_bins)%num_bins); // pick randomly from {0,1,2,3,...,num_bins-1}
				if (bin_alloc[random_bin] == 0){ // if there are no more items in that bin (not very likely, unless there are LOTS of NA distances)
					while (bin_alloc[random_bin] == 0){
						random_bin = (((rand()%num_bins)+num_bins)%num_bins); // then keep trying until we find a non-empty bin
					}
				}
				myDBSeq->bin[i][j] = random_bin + 1;
				bin_alloc[random_bin] --;
			}
		}
	}
	// then, allocate all other distances according to bin_alloc
	for (int b=1; b<=num_bins; b++){
		while (bin_alloc[b-1] > 0){
			min_dist = 100;
			for (int i=0; i<num_label; i++){
				for (int j=0; j<num_label; j++){
					if((min_dist > myDBSeq->nm1[i][j]) && (myDBSeq->bin[i][j] == 0)){ // looks for the minimum distance
						min_dist = myDBSeq->nm1[i][j]; // saves the minimum distance
						min_index_i = i; // stores the i index of the minimum value
						min_index_j = j; // stores the j index of the minimum value
					}
				}
			}
			myDBSeq->bin[min_index_i][min_index_j] = b;
			bin_alloc[b-1] --;
		}
	}
}
//</MMattar>


//<MMattar> Created the whole function
// Print out relevant information if in debug mode
void debug_mode(int num_label, int word_len, int num_bins, DBSeq* myDBSeq){

	long seq_length = pow(num_label,word_len);

	// Print the normalized nm matrix
	cout << endl << "NORMALIZED NEURAL MODEL:" << endl;
	for (int i=0; i<num_label; i++){
		for (int j=0; j<num_label; j++){
			printf ("%lf ", myDBSeq->nm1[i][j]);
		}
		printf ("\n");
	}

	// Print the normalized nm matrix
	cout << endl << "NORMALIZED NEURAL MODEL #2:" << endl;
	for (int i=0; i<num_label; i++){
		for (int j=0; j<num_label; j++){
			printf ("%lf ", myDBSeq->nm2[i][j]);
		}
		printf ("\n");
	}

	// Print the normalized nm matrix
	cout << endl << "NORMALIZED NEURAL MODEL #3:" << endl;
	for (int i=0; i<num_label; i++){
		for (int j=0; j<num_label; j++){
			printf ("%lf ", myDBSeq->nm3[i][j]);
		}
		printf ("\n");
	}

	// prints out the normalized guide function
	cout << endl << "NORMALIZED GUIDE FUNCTION:" << endl;
	for (int t=0; t<seq_length; t++){
		cout << myDBSeq->guide_function[t] << " ";
	}
	cout << endl << endl;

	// prints out the seq_bin
	cout << endl << "BIN SEQUENCE:" << endl;
	for (int t=0; t<seq_length; t++){
		cout << myDBSeq->seq_bin[t] << " ";
	}
	cout << endl << endl;

	// prints out the pref_bin sequence
	cout << endl << "SEQUENCE OF PREFERRED BINS:" << endl;
	for (int b=0; b<num_bins; b++){
		for (int t=0; t<30; t++){
			cout << myDBSeq->pref_bin[t][b][0] << myDBSeq->pref_bin[t][b][1] << "  ";
		}
		cout << "..." << endl;
	}
	cout << endl;

	// prints out the bin matrix
	cout << endl << "BIN ALLOCATION OF DISTANCES:" << endl;
	for (int i=0; i<num_label; i++){
		for (int j=0; j<num_label; j++){
			cout << myDBSeq->bin[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	// prints out the bins that were used when constructing the sequence
	cout << endl << "BINS USED:" << endl;
	for (int t=0; t<seq_length; t++){
		cout << myDBSeq->bins_used[t] << " ";
	}
	cout << endl << endl;

	// prints out the transition sequence obtained after constructing the sequence
	cout << endl << "TRANSITION SEQUENCE #1:" << endl;
	for (int t=0; t<seq_length; t++){
		cout << myDBSeq->trans_seq[t] << " ";
	}
	cout << endl << endl;

	// prints out the transition sequence obtained after constructing the sequence
	if (myDBSeq->num_neuralmodels > 1){
		cout << endl << "TRANSITION SEQUENCE #2:" << endl;
		for (int t=0; t<seq_length; t++){
			cout << myDBSeq->trans_seq2[t] << " ";
		}
	}
	cout << endl << endl;

	// prints out the transition sequence obtained after constructing the sequence
	if (myDBSeq->num_neuralmodels > 2){
		cout << endl << "TRANSITION SEQUENCE #3:" << endl;
		for (int t=0; t<seq_length; t++){
			cout << myDBSeq->trans_seq3[t] << " ";
		}
	}
	cout << endl << endl;
}
//</MMattar>


//<MMattar> Created the whole function
// Print out information for evaluating the sequence generated
void eval_mode(int num_label, int word_len, DBSeq* myDBSeq){

	long seq_length = pow(num_label,word_len);

	// correlation calculation
	// first, center the mean of trans_seq in zero
	double expected_value = 0.0;
	double expected_value2 = 0.0;
	double expected_value3 = 0.0;
	int num_posit_dist = 0;
	int num_posit_dist2 = 0;
	int num_posit_dist3 = 0;
	for (int t=0; t<seq_length; t++){
		if (myDBSeq->trans_seq[t] != -1){ // make sure we are not considering NA values
			expected_value = expected_value + myDBSeq->trans_seq[t];
			num_posit_dist ++;
		}
		if (myDBSeq->trans_seq2[t] != -1){ // make sure we are not considering NA values
			expected_value2 = expected_value2 + myDBSeq->trans_seq2[t];
			num_posit_dist2 ++;
		}
		if (myDBSeq->trans_seq3[t] != -1){ // make sure we are not considering NA values
			expected_value3 = expected_value3 + myDBSeq->trans_seq3[t];
			num_posit_dist3 ++;
		}
	}
	expected_value = expected_value/num_posit_dist;
	expected_value2 = expected_value2/num_posit_dist2;
	expected_value3 = expected_value3/num_posit_dist3;
	for (int t=0; t<seq_length; t++){
		if (myDBSeq->trans_seq[t] != -1){
			 myDBSeq->trans_seq[t] =  myDBSeq->trans_seq[t] - expected_value;
		}
		else{
			myDBSeq->trans_seq[t] = 0; // set all NA distances to zero
		}
		if (myDBSeq->trans_seq2[t] != -1){
			 myDBSeq->trans_seq2[t] =  myDBSeq->trans_seq2[t] - expected_value2;
		}
		else{
			myDBSeq->trans_seq2[t] = 0; // set all NA distances to zero
		}
		if (myDBSeq->trans_seq3[t] != -1){
			 myDBSeq->trans_seq3[t] =  myDBSeq->trans_seq3[t] - expected_value3;
		}
		else{
			myDBSeq->trans_seq3[t] = 0; // set all NA distances to zero
		}
	}

	// then, center the mean of guide_function
	expected_value = 0.0;
	for (int t=0; t<seq_length; t++){
		expected_value = expected_value + myDBSeq->guide_function[t];
	}
	expected_value = expected_value/seq_length;
	for (int t=0; t<seq_length; t++){
		myDBSeq->guide_function[t] =  myDBSeq->guide_function[t] - expected_value;
	}
	
	//at this point, both means should be zero

	// now, calculate both standard deviations
	double var_trans_seq = 0;
	double var_trans_seq2 = 0;
	double var_trans_seq3 = 0;
	double var_guide_function = 0;
	double std_trans_seq = 0;
	double std_guide_function = 0;
	for (int t=0; t<seq_length; t++){
		var_trans_seq = var_trans_seq + (myDBSeq->trans_seq[t]*myDBSeq->trans_seq[t]);
		var_guide_function = var_guide_function + (myDBSeq->guide_function[t]*myDBSeq->guide_function[t]);
	}
	std_trans_seq = sqrt(var_trans_seq/(seq_length-1));
	std_guide_function = sqrt(var_guide_function/(seq_length-1));

	// finally, calculate the correlation between trans_seq and guide_function
	double correlation = 0;
	for (int t=0; t<seq_length; t++){
		correlation = correlation + myDBSeq->guide_function[t]*myDBSeq->trans_seq[t];
	}
	correlation = correlation/(std_trans_seq*std_guide_function*(seq_length-1));


	// Print all this information
/*
	cout << endl << "ZERO_MEAN GUIDE FUNCTION:" << endl;
	for (int t=0; t<seq_length; t++){
		cout << myDBSeq->guide_function[t] << " ";
	}
	cout << endl << endl;

	cout << endl << "ZERO_MEAN TRANSITION SEQUENCE:" << endl;
	for (int t=0; t<seq_length; t++){
		cout << myDBSeq->trans_seq[t] << " ";
	}
	cout << endl << endl;
*/
	
	
	

	//*****Detection Power calculation*****
	
	int upsampled_resolution = 50; // time step after upsampling (in ms)
	
	// Upsampling the transition sequence using interpolation
	int factor = SOA/upsampled_resolution;
	int N = seq_length * factor;
	double trans_seq_upsampled[N];
	double trans_seq2_upsampled[N];
	double trans_seq3_upsampled[N];
	for (int i=0; i<N; i++){
		if ((i%factor) == 0){
			trans_seq_upsampled[i] = myDBSeq->trans_seq[i/factor];
			trans_seq2_upsampled[i] = myDBSeq->trans_seq2[i/factor];
			trans_seq3_upsampled[i] = myDBSeq->trans_seq3[i/factor];
		}
		else{
			trans_seq_upsampled[i] = myDBSeq->trans_seq[(int) i/factor] + (((myDBSeq->trans_seq[(int) i/factor + 1] - myDBSeq->trans_seq[(int) i/factor])/factor)*(i%factor));
			trans_seq2_upsampled[i] = myDBSeq->trans_seq2[(int) i/factor] + (((myDBSeq->trans_seq2[(int) i/factor + 1] - myDBSeq->trans_seq2[(int) i/factor])/factor)*(i%factor));
			trans_seq3_upsampled[i] = myDBSeq->trans_seq3[(int) i/factor] + (((myDBSeq->trans_seq3[(int) i/factor + 1] - myDBSeq->trans_seq3[(int) i/factor])/factor)*(i%factor));
		}
	}
	
	// Converting trans_seq_upsampled to a vector
	vector<double> trans_seq_vector (trans_seq_upsampled, trans_seq_upsampled + sizeof(trans_seq_upsampled) / sizeof(double) );
	vector<double> trans_seq2_vector (trans_seq2_upsampled, trans_seq2_upsampled + sizeof(trans_seq2_upsampled) / sizeof(double) );
	vector<double> trans_seq3_vector (trans_seq3_upsampled, trans_seq3_upsampled + sizeof(trans_seq3_upsampled) / sizeof(double) );
	
/*
//[DEBUG]
		cout << endl << "TRANSITION SEQUENCE (upsampled):" << endl;
		for (int t=0; t<N; t++){
			cout << trans_seq_vector[t] << " ";
		}
		cout << endl << endl;
//[/DEBUG]	
*/
	
	// Relevant variables for the FFT calculation
	fftw_complex *in, *out;
	fftw_complex *in2, *out2;
	fftw_complex *in3, *out3;
	fftw_plan planforward;
	fftw_plan planforward2;
	fftw_plan planforward3;
	fftw_plan planbackward;
	fftw_plan planbackward2;
	fftw_plan planbackward3;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	in3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	out3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	
	
	// hard-code the EndogenousFilter array
	double EndogenousFilter_array[160] = {5.55112e-17,0.00767224,0.0162726,0.025849,0.0364457,0.0481024,0.0608536,0.0747281,0.0897482,0.10593,0.123281,0.141801,0.161484,0.182312,0.204261,0.227297,0.251376,0.276448,0.30245,0.329315,0.356963,0.38531,0.414261,0.443716,0.473569,0.503705,0.534007,0.564352,0.594614,0.624665,0.654374,0.68361,0.712242,0.740139,0.767176,0.793226,0.818169,0.841889,0.864277,0.885229,0.90465,0.922451,0.938554,0.952889,0.965396,0.976025,0.984737,0.991502,0.996304,0.999136,1,0.998913,0.995897,0.990989,0.984232,0.975681,0.965398,0.953452,0.939923,0.924892,0.908452,0.890696,0.871725,0.851641,0.830549,0.808556,0.78577,0.7623,0.738253,0.713734,0.688848,0.663695,0.638375,0.61298,0.5876,0.562319,0.537218,0.51237,0.487844,0.463702,0.44,0.41679,0.394115,0.372016,0.350523,0.329666,0.309467,0.289942,0.271104,0.25296,0.235516,0.218771,0.202721,0.187361,0.172683,0.158674,0.145324,0.132616,0.120536,0.109067,0.0981926,0.0878945,0.0781552,0.0689569,0.0602817,0.0521123,0.0444312,0.0372216,0.0304668,0.0241507,0.0182571,0.0127706,0.00767572,0.00295715,-0.00140033,-0.00541187,-0.00909271,-0.0124583,-0.0155242,-0.0183063,-0.0208208,-0.0230841,-0.025113,-0.0269245,-0.0285358,-0.0299644,-0.0312277,-0.0323428,-0.0333269,-0.0341965,-0.0349674,-0.0356548,-0.0362724,-0.0368327,-0.0373465,-0.0378227,-0.0382678,-0.038686,-0.0390787,-0.0394442,-0.0397777,-0.0400708,-0.0403115,-0.0404842,-0.0405693,-0.0405433,-0.0403786,-0.0400439,-0.039504,-0.0387199,-0.0376493,-0.0362464,-0.0344627,-0.0322468,-0.0295455,-0.0263035,-0.0224645,-0.0179715,-0.0127673,-0.00679528};
	
	// Upsampling the EndogenousFilter using interpolation
	double EndogenousFilter_upsampled[320];
	for (int i=0; i<320; i++){
		if ((i%2) == 0){
			EndogenousFilter_upsampled[i] = EndogenousFilter_array[i/2];
		}
		else{
			EndogenousFilter_upsampled[i] = EndogenousFilter_array[(int) i/2] + (((EndogenousFilter_array[(int) i/2 + 1] - EndogenousFilter_array[(int) i/2])/2)*(i%2));
		}
	}
	
	// Converting EndogenousFilter_upsampled to a vector
	vector<double> EndogenousFilter (EndogenousFilter_upsampled, EndogenousFilter_upsampled + sizeof(EndogenousFilter_upsampled) / sizeof(double) );

/*
//[DEBUG]
	cout << endl << "ENDOGENOUS FILTER (upsampled):" << endl;
	for (int t=0; t<320; t++){
		cout << EndogenousFilter[t] << " ";
	}
	cout << endl << endl;
//[/DEBUG]
*/


	// copy EndogenousFilter to a format suitable for the fourier transform
	for (int i=0; i<N; i++){
		if (i<320){
			in[i][0] = EndogenousFilter[i]; //real part
		}
		else{
			in[i][0] = 0; //real part
		}
		in[i][1] = 0; // complex part
	}
	
	// Calculate the FFT of EndogenousFilter
	planforward = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(planforward);

	// Store the results on a vector
	vector<double> EndogenousKernel (2*N,0);
	for (int i=0; i<N; i++){
		EndogenousKernel[2*i] = out[i][0]; //real part
		EndogenousKernel[2*i + 1] = out[i][1]; // complex part
	}

	// copy trans_seq_fft to a format suitable for the fourier transform
	for (int i=0; i<N; i++){
		in[i][0] = trans_seq_vector[i]; //real part
		in[i][1] = 0; // complex part

		in2[i][0] = trans_seq2_vector[i]; //real part
		in2[i][1] = 0; // complex part

		in3[i][0] = trans_seq3_vector[i]; //real part
		in3[i][1] = 0; // complex part
	}
	
	// Calculate the FFT of trans_seq_vector
	planforward = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(planforward);
	
	// Calculate the FFT of trans_seq2_vector
	planforward2 = fftw_plan_dft_1d(N, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(planforward2);
	
	// Calculate the FFT of trans_seq3_vector
	planforward3 = fftw_plan_dft_1d(N, in3, out3, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(planforward3);

	// Store the results on a vector
	vector<double> trans_seq_fft (2*N,0);
	vector<double> trans_seq2_fft (2*N,0);
	vector<double> trans_seq3_fft (2*N,0);
	for (int i=0; i<N; i++){
		trans_seq_fft[2*i] = out[i][0]; //real part
		trans_seq_fft[2*i + 1] = out[i][1]; // complex part

		trans_seq2_fft[2*i] = out2[i][0]; //real part
		trans_seq2_fft[2*i + 1] = out2[i][1]; // complex part

		trans_seq3_fft[2*i] = out3[i][0]; //real part
		trans_seq3_fft[2*i + 1] = out3[i][1]; // complex part
	}

/*
//[DEBUG]
	cout << endl << "ENDOGENOUS KERNEL:" << endl; // complex conjugate?
	for (int t=0; t<N*2; t++){
		cout << EndogenousKernel[t] << " ";
	}
	cout << endl << endl;
//[/DEBUG]
*/

	// In the frequency domain, calculate the maximum norm and normalize EndogenousKernel
	double maxnorm = 0;
	double currentnorm = 0;
	for (int t=0; t<N*2; t=t+2){
		currentnorm = sqrt(EndogenousKernel[t]*EndogenousKernel[t] + EndogenousKernel[t+1]*EndogenousKernel[t+1]);
		if (maxnorm<currentnorm) maxnorm=currentnorm;
	}
	for (int t=0; t<N*2; t++){
		EndogenousKernel[t] = EndogenousKernel[t]/maxnorm;
	}

	// In the frequency domain, apply the Notch filter
	float FilterCutoff = 0.01;
	for (int t=0; t<(N*FilterCutoff*upsampled_resolution/1000); t++){
		EndogenousKernel[2*t] = 0;   // real part
		EndogenousKernel[2*t+1] = 0; // imaginary part
	}
	
	for (int t=N; t>(N - (N*FilterCutoff*upsampled_resolution/1000)); t--){
		EndogenousKernel[2*t] = 0;   // real part
		EndogenousKernel[2*t+1] = 0; // imaginary part
	}

/*
//[DEBUG]
	cout << endl << "ENDOGENOUS KERNEL (normalized and filtered):" << endl; // complex conjugate?
	for (int t=0; t<N*2; t++){
		cout << EndogenousKernel[t] << " ";
	}
	cout << endl << endl;

	cout << endl << "TRANSITION SEQUENCE (frequency domain):" << endl; // complex conjugate?
	for (int t=0; t<N*2; t++){
		cout << trans_seq_fft[t] << " ";
	}
	cout << endl << endl;
//[/DEBUG]
*/

	// Multiply (convolve) EndogenousKernel with trans_seq
	vector<double> convolution_fft (N*2,0);
	vector<double> convolution2_fft (N*2,0);
	vector<double> convolution3_fft (N*2,0);
	for (int t=0; t<N*2; t=t+2){
		// Real part
		convolution_fft[t] = trans_seq_fft[t]*EndogenousKernel[t] - trans_seq_fft[t+1]*EndogenousKernel[t+1];
		// Imaginary part
		convolution_fft[t+1] = trans_seq_fft[t]*EndogenousKernel[t+1] + trans_seq_fft[t+1]*EndogenousKernel[t];

		// Real part
		convolution2_fft[t] = trans_seq2_fft[t]*EndogenousKernel[t] - trans_seq2_fft[t+1]*EndogenousKernel[t+1];
		// Imaginary part
		convolution2_fft[t+1] = trans_seq2_fft[t]*EndogenousKernel[t+1] + trans_seq2_fft[t+1]*EndogenousKernel[t];

		// Real part
		convolution3_fft[t] = trans_seq3_fft[t]*EndogenousKernel[t] - trans_seq3_fft[t+1]*EndogenousKernel[t+1];
		// Imaginary part
		convolution3_fft[t+1] = trans_seq3_fft[t]*EndogenousKernel[t+1] + trans_seq3_fft[t+1]*EndogenousKernel[t];
	}

/*
//[DEBUG]
	cout << endl << endl << "Convolved (frequency domain) = ["; // complex conjugate?
	for (int t=0; t<N*2; t++){
		cout << convolution_fft[t] << " ";
	}
	cout << "]" << endl << endl;
//[/DEBUG]
*/

	// copy convolution_fft to a format suitable for the fourier transform
	for (int i=0; i<N; i++){
		in[i][0] = convolution_fft[2*i]; //real part
		in[i][1] = convolution_fft[2*i + 1]; // complex part

		in2[i][0] = convolution2_fft[2*i]; //real part
		in2[i][1] = convolution2_fft[2*i + 1]; // complex part

		in3[i][0] = convolution3_fft[2*i]; //real part
		in3[i][1] = convolution3_fft[2*i + 1]; // complex part
	}
	
	// Calculate the iFFT of GuideFunction_fft
	planbackward = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(planbackward);
	
	// Calculate the iFFT of GuideFunction_fft
	planbackward2 = fftw_plan_dft_1d(N, in2, out2, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(planbackward2);

	// Calculate the iFFT of GuideFunction_fft
	planbackward3 = fftw_plan_dft_1d(N, in3, out3, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(planbackward3);
	
	vector<double> convolution (N,0);
	vector<double> convolution2 (N,0);
	vector<double> convolution3 (N,0);
	for (int i=0; i<N; i++){
		convolution[i] = out[i][0]/N;
		convolution2[i] = out2[i][0]/N;
		convolution3[i] = out3[i][0]/N;
	}
	
	fftw_destroy_plan(planforward);
	fftw_destroy_plan(planforward2);
	fftw_destroy_plan(planforward3);
	fftw_destroy_plan(planbackward);
	fftw_destroy_plan(planbackward2);
	fftw_destroy_plan(planbackward3);
	fftw_free(in); fftw_free(out);
	fftw_free(in2); fftw_free(out2);
	fftw_free(in3); fftw_free(out3);
	
/*
//[DEBUG]
	cout << endl << endl << "Convolved (time domain) = [";
	for (int t=0; t<N; t++){
		cout << convolution[t] << " ";
	}
	cout << "]" << endl << endl;
//[/DEBUG]
*/

	// Centering the mean in zero
	expected_value = 0;
	expected_value2 = 0;
	expected_value3 = 0;
	for (int t=0; t<N; t++){
		expected_value =  expected_value + convolution[t];
		expected_value2 =  expected_value2 + convolution2[t];
		expected_value3 =  expected_value3 + convolution3[t];
	}
	expected_value = expected_value/N;
	expected_value2 = expected_value2/N;
	expected_value3 = expected_value3/N;
	for (int t=0; t<N; t++){
		convolution[t] =  convolution[t] - expected_value;
		convolution2[t] =  convolution2[t] - expected_value2;
		convolution3[t] =  convolution3[t] - expected_value3;
	}

	// calculate the variance of convolved_seq
	double var_convolved = 0;
	double var_convolved2 = 0;
	double var_convolved3 = 0;
	for (int t=0; t<N; t++){
		var_convolved = var_convolved + (convolution[t]*convolution[t]);
		var_convolved2 = var_convolved2 + (convolution2[t]*convolution2[t]);
		var_convolved3 = var_convolved3 + (convolution3[t]*convolution3[t]);
	}
	var_convolved = var_convolved/N;
	var_convolved2 = var_convolved2/N;
	var_convolved3 = var_convolved3/N;

	// calculate the variance of trans_seq
	var_trans_seq = 0;
	var_trans_seq2 = 0;
	var_trans_seq3 = 0;
	for (int t=0; t<N; t++){
		var_trans_seq = var_trans_seq + (trans_seq_vector[t]*trans_seq_vector[t]);
		var_trans_seq2 = var_trans_seq2 + (trans_seq2_vector[t]*trans_seq2_vector[t]);
		var_trans_seq3 = var_trans_seq3 + (trans_seq3_vector[t]*trans_seq3_vector[t]);
	}
	var_trans_seq = var_trans_seq/(N-1);
	var_trans_seq2 = var_trans_seq2/(N-1);
	var_trans_seq3 = var_trans_seq3/(N-1);

	double detection_power = var_convolved/var_trans_seq;
	double detection_power2 = var_convolved2/var_trans_seq2;
	double detection_power3 = var_convolved3/var_trans_seq3;

/*
//[DEBUG]	
	cout << "STD DISTANCES SEQUENCE =          " << std_trans_seq << endl;
	cout << "STD GUIDE FUNCTION =              " << std_guide_function << endl;
	cout << "VARIANCE OF TRANSITION SEQUENCE:  " << var_trans_seq << endl;
	cout << "VARIANCE OF CONVOLVED SEQUENCE:   " << var_convolved << endl;
//[/DEBUG]
*/
	
	cout << endl << "CORRELATION COEFFICIENT:             " << correlation << endl;
	if (detection_power2 != detection_power2){	// if detection_power2 is nan, there's only one neural model
		cout << "DETECTION POWER:                  " << detection_power << endl;
	}
	else{										// if detection_power2 is not nan, print both detection_power and detection_power2
		cout << "DETECTION POWER #1:                  " << detection_power << endl;
		cout << "DETECTION POWER #2:                  " << detection_power2 << endl;
	}
	if (detection_power3 == detection_power3){	// if detection_power3 is not nan, print it as well
		cout << "DETECTION POWER #3:                  " << detection_power3 << endl;
	}
	cout << endl;

}
//</MMattar>


// Set up sequence options
void setOptions()
{
  myOpt.algorithm = ALG_NOPRUNE_BT;
  myOpt.report_flags = REPORT_SOLUTION | REPORT_SUMMARY;

  /* default RNG seed is randomly generated from the current time
  myOpt.rng_seed = (int) time(NULL);	
  srandom(myOpt.rng_seed);
  srand48((long) myOpt.rng_seed );
 */
  // set log strings
  myOpt.stat_str = "Statistical information:\n\n";
  char tmpStr[STRLEN];
  sprintf(tmpStr, "\nRNG Seed = %d\n", myOpt.rng_seed);
  myOpt.stat_str += string(tmpStr);
  myOpt.log_str = "Program execution information:\n\n";
}


//<MMattar> Included new inputs
// Print out usage information
void usage() 
{
  cout << endl << endl;
  cout << "Usage: debruijn [-t | -v] <k> <n>" << endl;
  cout << "         [<B> <guide function> <neural model>]" << endl;
  cout << "         [<neural model #2>] [<neural model #3>] [-eval <SOA>]" << endl;
  cout << "Generate random and \"path guided\" de Bruijn cycles" << endl;
  cout << endl;
  cout << " Default output has one line of labels, separated by commas" << endl;
  cout << "     -t: terse output (no delimiters)" << endl;
  cout << "     -v: verbose output in \"necklace\" format" << endl;
  cout << endl;
  cout << "   k: number of letters (maximum is 36)" << endl;
  cout << "   n: length of the word (i.e., level of counterbalance)" << endl;
  cout << endl;
  cout << " Specify parameters for path-guided cycle" << endl;
  cout << "   B: number of bins for available transitions (integer between 1 and k^2)" << endl;
  cout << "   guide function: can be specified in three different ways -" << endl;
  cout << "                   HRF: guide with power spectrum of the BOLD HRF  -OR-" << endl;
  cout << "                   path to .txt file containing guide function -OR-" << endl;
  cout << "                   [Tmin,Tmax]: range of the periods of sinusoids" << endl;
  cout << "                     (in units of labels) to use as a guide function" << endl;
  cout << "   neural model: path to .txt file containing neural model for optimizing" << endl;
  cout << "   neural model(#2,#3): optional; detection power only reported" << endl;
  cout << endl;
  cout << " -eval: Evaluate detection power" << endl;
  cout << "   SOA: stimulus-onset asynchrony (in ms)" << endl;
  cout << endl;
  cout << "Detailed help: http://cfn.upenn.edu/aguirre/wiki/public:de_bruijn_software" << endl;
  cout << endl;
  cout << "Reference: GK Aguirre, MG Mattar, L Magis-Weinberg. (2011)" << endl;
  cout << "de Bruijn cycles for neural decoding. NeuroImage 56: 1293-1300" << endl;
  cout << endl;
  cout << "Code by: Marcelo Mattar (mattar@sas.upenn.edu)" << endl;
  cout << "         Dongbo Hu (dongbo@mail.med.upenn.edu)" << endl;
  cout << "         Hamiltonian circuit code from Basil Vandergriend" << endl;
  cout << "           (http://webdocs.cs.ualberta.ca/~joe/Theses/vandegriend.html)" << endl;
  cout << endl;
  cout << "v1.6 -- July 09, 2012" << endl;
  cout << endl;
  
}
//</MMattar>


//<MMattar> Changed the input parameters
// main function
int main(int argc, char** argv)
{
	
  //checks if the number of arguments were correct; if not, print usage information
  if (argc < 3 || argc > 12) {
    usage();
    exit(1);
  }
  
  int argsleft = argc-1;
  for (int i=1; i<argc; i++){
	  if (strcmp(argv[i], "-debug") == 0){
		  debugmode = 1;
		  argsleft = argsleft-1;
	  }
	  if (strcmp(argv[i], "-eval") == 0){
		  if((i+1)>=argc){
			  usage();
			  exit(1);
		  }
		  else{
			  evalmode = 1;
			  SOA = atoi(argv[i+1]);
			  if (SOA%50 != 0){
				  cout << endl << "The SOA must be a multiple of 50ms." << endl << endl;
				  exit(1);
			  }
			  argsleft = argsleft-2;
		  }
	  }
	  if (strcmp(argv[i], "-h") == 0){
		  usage();
		  exit(1);
	  }
  }


  /* Three output modes:
   * print_flag = 0: terse mode, one line of labels without any delimiters
   * print_flag = 1: default mode, one line of labels, separated by commas
   * print_flag = 2: verbose mode, multiple lines like a necklace */
  int indx = 1;

  int num_bins = 0;
  char * GuideFunction_filename = NULL;
  char * NeuralModel1_filename = NULL;
  char * NeuralModel2_filename = NULL;
  char * NeuralModel3_filename = NULL;
  FILE * gffile;
  FILE * nm1file;
  FILE * nm2file;
  FILE * nm3file;

  if (strcmp(argv[1], "-t") == 0){ // terse mode
	  print_flag = 0;
	  indx++;
	  argsleft = argsleft-1;
  }
  else if (strcmp(argv[1], "-v") == 0){ // verbose mode
	  print_flag = 2;
	  indx++;
	  argsleft = argsleft-1;
  }

  if (argsleft <= 1){
	  usage();
	  exit(1);
  }

  // now, read required arguments (k, n)
  int num_label = atoi(argv[indx]); // converts argument to integer
  int word_len = atoi(argv[indx + 1]); // converts argument to integer
  argsleft = argsleft-2;

  // do not allow num_label < 2 or > 36
  if (num_label <= 1 || num_label > 36) {
    cout << endl << "Invalid number of labels: " << argv[indx] << endl << endl;
    usage();
    exit(1);
  }

  // do not allow word_len < 2
  if (word_len <= 1) {
    cout << endl << "Invalid value of word length: " << argv[indx + 1] << endl << endl;
    usage();
    exit(1);
  }

  // read number of bins, neural model and guide function
  if (argsleft >= 3){
	  num_bins = atoi(argv[indx + 2]); // converts argument to integer
	  // do not allow B < 1 or > num_label^2
	  if (num_bins < 1 || num_bins > (num_label*num_label)) {
	    cout << endl << "Invalid number of bins: " << argv[indx + 2] << endl << endl;
	    usage();
	    exit(1);
	  }

	  GuideFunction_filename = argv[indx + 3];
	  gffile = fopen (GuideFunction_filename , "r");
	  if (gffile == NULL) { // if did not load the guide function file correctly
		  if (GuideFunction_filename[0]=='[' && GuideFunction_filename[strlen(GuideFunction_filename)-1]==']'){
			  int Tmin = 0;
			  int Tmax = 0;
			  int i = 1;
			  while(GuideFunction_filename[i]!=','){
				  Tmin = 10*Tmin + (GuideFunction_filename[i] - '0');
				  i++;
			  }
			  i++;
			  while(GuideFunction_filename[i]!=']'){
				  Tmax = 10*Tmax + (GuideFunction_filename[i] - '0');
				  i++;
			  }
			  if (Tmin==0 || Tmax==0){
				  cout << endl << "Invalid path or information for guide function file: " << argv[indx + 4] << endl << endl;
				  usage();
				  exit(1);
			  }
		  }
		  else if (GuideFunction_filename[0]=='H' && GuideFunction_filename[1]=='R' && GuideFunction_filename[2]=='F'){
			  cout << endl<< "Using a Hemodynamic Response Filter as a guide function." << endl;
		  }
		  else{
			  cout << endl << "Invalid path or information for guide function file: " << argv[indx + 4] << endl << endl;
			  usage();
			  exit(1);
		  }
	  }

	  NeuralModel1_filename = argv[indx + 4];
	  nm1file = fopen (NeuralModel1_filename , "r");
	  if (nm1file == NULL) { // if did not load the neural model file correctly
		  cout << endl << "Invalid path for neural model file: " << argv[indx + 3] << endl << endl;
		  usage();
		  exit(1);
	  }
	  argsleft = argsleft-3;
	  
	  if (argsleft != 0){
		  NeuralModel2_filename = argv[indx + 5];
		  nm2file = fopen (NeuralModel2_filename , "r");
		  if (nm2file == NULL) { // if did not load the neural model file correctly
			  cout << endl << "Invalid path for neural model file: " << argv[indx + 3] << endl << endl;
			  usage();
			  exit(1);
		  }
		  argsleft = argsleft-1;
	  }
	  
	  if (argsleft != 0){
		  NeuralModel3_filename = argv[indx + 6];
		  nm3file = fopen (NeuralModel3_filename , "r");
		  if (nm3file == NULL) { // if did not load the neural model file correctly
			  cout << endl << "Invalid path for neural model file: " << argv[indx + 3] << endl << endl;
			  usage();
			  exit(1);
		  }
		  argsleft = argsleft-1;
	  }
  }
  else{
	  if (evalmode == 1){
		  cout << endl << "You must specify a neural model and a guide function in order to evaluate your sequence." << endl << endl;
		  exit(1);
	  }
  }

  if (argsleft != 0){
	  //cout << endl << "argsleft = " << argsleft << endl << endl;
	    usage();
	    exit(1);
  }
	
  // Seed the random number generator with microsecond precision
  timeval t1;
  gettimeofday(&t1, NULL);
  srand(t1.tv_usec * t1.tv_sec);

  setOptions(); // Set up sequence options
  DBSeq myDBSeq; // Creates an empty sequence
  myDBSeq.setup(num_label, word_len, num_bins, NeuralModel1_filename, NeuralModel2_filename, NeuralModel3_filename, GuideFunction_filename); // initializes the DeBruijn sequence structure
  if (num_bins!=0) {
	  binning(num_label, word_len, num_bins, &myDBSeq); // allocates all transitions in bins according to the period chosen
  }
  
  start(num_label, word_len, &myDBSeq); // function to perform the experiments

  return 0;
}  
//</MMattar>


/* function to perform the experiments - calls every other important function
 *   - report solution if desired */
//<MMattar> Included parameter myDBSeq
void start(int num_label, int word_len, DBSeq* myDBSeq)
{
  starttime = time (NULL); // start to record time
  
  Graph myGraph;
  
  myGraph.setup(num_label, word_len); // Generates a graph with num_label and word_len
  //myGraph.printNeighbors();

  int trialnum = 0;
  int graphnum = 0;

  int ret = 0;
  ExpStat myStat;

  /* start looping through experiments */
  while (1) {
    ret = test_hc_alg(num_label, word_len, &myGraph, myDBSeq, &(myStat.graph[graphnum].trial[trialnum])); // tries to run the HC algorithm
    /* save graph if timelimit reached */
    if (ret == RUN_TIMELIMIT) // if reach time limit
      cout << "test_hc_alg() time limit reached" << endl;

    /* update graph stats if necessary */
    if (myStat.graph[graphnum].trial[trialnum].result == HC_FOUND) {
      if (myStat.graph[graphnum].graphham != HC_NOT_EXIST) 
        myStat.graph[graphnum].graphham = HC_FOUND;
      else
        cout << "HC_FOUND and HC_NOT_EXIST occurred for same graph." << endl;
    }
    else if (myStat.graph[graphnum].trial[trialnum].result == HC_NOT_EXIST) {
      if (myStat.graph[graphnum].graphham != HC_FOUND)
        myStat.graph[graphnum].graphham = HC_NOT_EXIST;
      else 
        cout << "HC_FOUND and HC_NOT_EXIST occurred for same graph.\n" << endl;
    }

    trialnum++;

    /* check if finished with tests on current graph */
    if (trialnum == SeqConst::num_instance_tests) {
      /* generate mindegree2 and biconnected statistics 
       * this assumes graph data structure was never changed */
      test_graph_properties(&myGraph, &(myStat.graph[graphnum]));

      graphnum++;
      trialnum = 0;
      //if (graphnum == SeqConst::num_graph_tests)
      break;	/* exit from while loop */
    }
  }  /* end of while loop through tests */
  
  //<MMattar> Force exit, in order to avoid going back to calc_noprune_bt_alg()	
  if (debugmode==1) debug_mode(num_label, word_len, myDBSeq->num_bins, myDBSeq);
  if (evalmode==1) eval_mode(num_label, word_len, myDBSeq);
  
  exit(0); //sequence found, so that's the end.
}


/* function to test if graph has certain desired properties
 * tests if graph has min degree >= 2
 * tests if graph is biconnected */
void test_graph_properties(Graph* graph, GraphStat* graphstat)
{
  int loop;

  /* test if graph has min degree >= 2 */
  graphstat->mindeg2 = 1;
  for (loop = 0; loop < graph->numvert; loop++) {
    if (graph->degree[loop] < 2) {
      graphstat->mindeg2 = 0;
      break;
    }
  }

  /* test if graph is biconnected */
  /* if minimum degree is < 2, then not biconnected */
  if (graphstat->mindeg2 == 0) {
    graphstat->biconnected = 0;
    return;
  }

  if (check_graph_cutpoints(graph) == CUTPNT_EXIST)
    graphstat->biconnected = 0;
  else
    graphstat->biconnected = 1;
  
}  /* end of test_graph_properties() */

/* testing wrapper for hamiltonian cycle algorithms
 * this executes a single trial
 *
 * pass in graph to test, and stats to update
 * graph should not be changed at end of function
 * return RUN_TIMELIMIT if hit timelimit, otherwise return RUN_NORMAL */
//<MMattar> Included parameter myDBSeq
int test_hc_alg(int num_label, int word_len, Graph* graph, DBSeq* myDBSeq, TrialStat* trialstats)
{
  int retval = RUN_NORMAL;
  struct rusage curtime;
  int solution[graph->numvert]; //array containing the solution (length = numvert)

  /* start stat timer */
  g_hit_timelimit = RUN_NORMAL;
  long timeret = getrusage(RUSAGE_SELF,&curtime);
  g_algstart.tv_sec = curtime.ru_utime.tv_sec;
  g_algstart.tv_usec = curtime.ru_utime.tv_usec;

  ///////////////////// copied from backtrack.cpp
  int hcret = master_backtrack_alg(graph, myDBSeq, trialstats, solution);

  /* stop timing, and calculate elapsed time in seconds */
  timeret = getrusage(RUSAGE_SELF,&curtime);
  g_algstart.tv_sec = curtime.ru_utime.tv_sec - g_algstart.tv_sec;
  g_algstart.tv_usec = curtime.ru_utime.tv_usec - g_algstart.tv_usec;
  trialstats->time = g_algstart.tv_sec + ((float)g_algstart.tv_usec/1000000.0);

  trialstats->result = hcret;

  /* report if hit timelimit */
  if (g_hit_timelimit == RUN_TIMELIMIT) {
    retval = RUN_TIMELIMIT;
    cout << "Warning: algorithm reached time limit." << endl;
  }

  /* verify solution */
  if (hcret == HC_FOUND) {
    if (hc_verify_solution(graph, solution) == HC_NOT_VERIFY) 
        cout << "Error:  algorithm calculated bad solution.\n" << endl;
  }

  /* print solution if found */
  if (hcret == HC_FOUND) {
    // calculate the conversion factor to convert integer in solution to vertex number 
    int factor = 1;
    for (int i = 0; i < word_len - 1; i++)
      factor *= num_label;
    /* print solution */
    printf("\nCycle found:\n");
    printf("------------\n");

    if (print_flag < 2)
      printSeq(solution, graph->numvert, factor);
    else
      printSeq_v(solution, graph->numvert, factor, word_len);

    printf("------------\n");

    // dhu: double check the solution
    if (!chkSolution(num_label, word_len, solution))
      cout << "Invalid solution \n" << endl;
    else
      cout << "Valid sequence \n" << endl;
  }
  else
    printf("\nNoCycle \n\n");
  
  return retval;  
}

// print out sequences in terse mode
void printSeq(int solution[], int n_vert, int factor)
{
  for (int i = 0; i < n_vert; i++) {
    int foo = solution[i] / factor;
    if (foo < 10)
      cout << foo;
    else
      cout << static_cast<char>(foo - 10 + 'A');

    // print a comma to separate the labels
    if (print_flag == 1 && i != n_vert -1)
      cout << ",";
  }
  cout << endl;

}

// print out sequences in verbose mode (default)
void printSeq_v(int solution[], int n_vert, int factor, int word_len)
{
  for (int i = 0; i < n_vert; i++) {
    for (int j = 0; j < word_len; j++) {
      int indx = (i + j) % n_vert;
      int foo = solution[indx] / factor;
      if (foo < 10)
	cout << foo << " ";
      else
	cout << static_cast<char>(foo - 10 + 'A') << " ";
    }
    cout << endl;
  }
  
}

// Dongbo's function to double check solution
bool chkSolution(int num_label, int word_len, int solution[])
{
  int factor = 1;
  for (int i = 0; i < word_len - 1; i++)
    factor *= num_label;

  int num_vert = factor * num_label;
  bool flag[num_vert];
  for (int i = 0; i < num_vert; i++)
    flag[i] = false;

  set<string> sol_str;
  for (int i = 0; i < num_vert; i++) {
    if (solution[i] >= num_vert || solution[i] < 0) {
      cout << "solution[" << i << "] out of range: " << solution[i] << endl;
      return false;
    }
    if (flag[solution[i]]) {
      cout << "solution[" << i << "] already exists: " << solution[i] << endl;
      return false;
    }
      
    flag[solution[i]] = true;

    string tmpStr;
    pair<set<string>::iterator, bool> ret;
    for (int j = 0; j < word_len; j++) {
      int indx = (i + j) % num_vert;
      tmpStr += num2str(solution[indx] / factor) + " ";
    }

    ret = sol_str.insert(tmpStr);
    if (!ret.second) {
      cout << "Repetitive elements found: " << tmpStr << endl;
      return false;
    }
  }
  int sol_str_size = sol_str.size();
  if (sol_str_size != num_vert) {
    cout << "Repetitive elements found\n";
    return false;
  }

  return true;

}

// Copied from backtrack.c
//<MMattar> Included parameter myDBSeq
int master_backtrack_alg(Graph* graph, DBSeq* myDBSeq, TrialStat* trialstats, int solution[])
{
  int ret;
  int loop;
  int nodecount = 0;
  int prune = 0;

  EdgeStack edgestack;   /* for saving deleted edges */
  edgestack.stack = vector<Edge>(graph->numvert * graph->degree[0] / 2);

  vector<Path> path(graph->numvert); // creates a Path vector named path, with length equal to numvert
  vector<GraphPath> graphpath(graph->numvert); // creates a GraphPath vector named graphpath, with length equal to numvert
  int pstart, pend, plength;
  int initvert;
  int tempnum;



  /* initialize variables, setting everything to -1 */
  for (loop = 0; loop < graph->numvert; loop++) {
    graphpath[loop].pathpos = -1;
    path[loop].gvert = -1;
    path[loop].next = -1;
  }

  ///// dhu: call select_initvertex() in hamcycle.cpp
  initvert = select_initvertex(graph, myOpt.bt_alg.initvertflag); // as default, selects a random vertex to be the first one (random number between 0 and numvert)

  pstart = pend = 0;
  plength = 1;
  path[pstart].gvert = initvert;
  path[pstart].next = -1;
  graphpath[initvert].pathpos = pstart;

  /* update statistics */
  trialstats->nodes = 0;
  trialstats->edgeprune = 0;

  // dhu: call copy_graph() in graphdata.c
  Graph testgraph(*graph);

  ////// dhu: call hc_do_pruning in hamcycle.cpp
  // reduces the complexity of the graph by pruning all the branches that certainly won't be used in the circuit
  // returned value indicates if an HC is impossible, or if an obligatory HC was already found
  ret = hc_do_pruning(&testgraph, &prune, HC_PRUNE_ALL, &edgestack);

  /* update initial prune statistic */
  trialstats->initprune = prune; // the way the code is now, prune is always 0

  if (ret == HC_NOT_EXIST) {
    return(HC_NOT_EXIST);
  }
  /* if have forced HC, then just let backtrack quickly find it to get
   * the actual solution */
  prune = 0;

  /* call recursive hc-backtrack algorithm */

  ret = calc_noprune_bt_alg(&testgraph, &pstart, &pend, &plength, path, graphpath, &nodecount, myDBSeq);
  trialstats->nodes = nodecount;

  if (ret == HC_FOUND) {
    /* convert path to solution */
    for (tempnum = pstart, loop = 0; loop < graph->numvert; loop++) {
      solution[loop] = path[tempnum].gvert;
      tempnum = path[tempnum].next;
    }
  }

  return (ret); // indicates if the HC was found or not
}  

// Copied from hamcycle.c
int select_initvertex(Graph* graph, int selectflag) 
{
  int selvert = 0;

  int vlist[graph->numvert];
  int numvert;
  int maxdeg;
  int degsum;
  int loop;

  switch(selectflag)
  {
    case INITVERT_RANDOM: //That's the default being used in this code
      /* select a vertex at random */
      selvert = lrand48() % graph->numvert; // Modulo (remainder) of a very large number and numvert (selects a random vertice)
      break;
    case INITVERT_MAXDEG:
      /* select a random vertex from those vertices of maximum degree
       * first calculate the maximum degree */
      for ( maxdeg = graph->degree[0], loop = 1; loop < graph->numvert; loop++) {
        if (graph->degree[loop] > maxdeg) {
          maxdeg = graph->degree[loop];
        }
      }

      /* create list of vertices of maximum degree */
      for (numvert = 0, loop = 0; loop < graph->numvert; loop++) {
        if (graph->degree[loop] == maxdeg) {
          vlist[numvert++] = loop;
        }
      }

      /* randomly select a vertex from this list */
      selvert = vlist[ (lrand48() % numvert) ];
      break;

    case INITVERT_RANDEG:
      /* select a vertex randomly, with selection probability weighted
       * according to vertex degree */

      /* first calculate summation of vertex degrees */
      for (degsum = 0, loop = 0; loop < graph->numvert; loop++) {
        degsum += graph->degree[loop];
      }

      /* do weighted random selection */
      degsum = lrand48() % degsum;
      selvert = -1;
      do {
        selvert++;
        degsum -= graph->degree[selvert];
      } while (degsum > 0);
      
      break;

    case INITVERT_FIRST:
      selvert = 0;
      break;
  }  /* end of switch */

  return (selvert);

}

// Copied from hamcycle.c
// The goal of this function is to reduce the complexity of the graph by pruning all the branches that certainly won't be used in the circuit
// There are many different rules for prunning, each consisting of a condition that has to happen in order to achieve an HC
// The returned value is
int hc_do_pruning(Graph* graph, int *prune, int prunelevel, EdgeStack* edgestack)
{
  int done;
  int ret;
  int degmrk[graph->numvert];
  int loop, eloop;
  int curprune;
  int used[graph->numvert];  
  int newvert;


  *prune = 0;
  do {
    done = 1;
    //first, do some basic tests to check if some impeditive conditions exist
    if (prunelevel & HC_PRUNE_BASIC) {
      /* check if vertex degrees are >= 2 */
      for (loop = 0; loop < graph->numvert; loop++) {
        if (graph->degree[loop] < 2) 
          return(HC_NOT_EXIST); // if any vertex degree is less than 1, then HC not exist
      }

      /* check if # of deg 2 neighbours is <= 2 */
      for (loop = 0; loop < graph->numvert; loop++) {
        degmrk[loop] = 0;
        for(eloop = 0; eloop < graph->degree[loop]; eloop++) {
          if (graph->degree[graph->nbr[loop][eloop]] == 2)
            degmrk[loop]++;
        }
  
        if (degmrk[loop] > 2)
          return(HC_NOT_EXIST); // for each vertex, if the number of neighbors with degree=2 is more than 2, then HC not exist
      }
    }

    /* initialize used[] variable */
    for (loop = 0; loop < graph->numvert; loop++) {
      used[loop] = 0; //used is an array with length numvert; set everything to 0
    }



    /* prune edges that cannot be traversed in H.C. */
    for (loop = 0; loop < graph->numvert; loop++) {
      if (prunelevel & HC_PRUNE_BASIC) {
        /* prune extra edges of vertices with 2 degree-2 neighbours */
        /* note that this vertex also becomes degree 2 */
        if ( (degmrk[loop] == 2) && (graph->degree[loop] > 2) ) {
          eloop = 0;
          while(eloop < graph->degree[loop]) {
            newvert = graph->nbr[loop][eloop];
            if (graph->degree[newvert] != 2 ) {
              rm_edge_graph(graph, loop, newvert);
              (*prune)++;
              push_edge_to_stack(loop, newvert, edgestack);
              done = 0;		/* redo entire pruning check */
            }
            else
              eloop++;
          }
        }
      }

      if (prunelevel & HC_PRUNE_CYC) {
        /* find longest path of forced edges, and remove any cycle-creating
         * edge from this path.
         * only check each forced path once, by using the used[] array */
        if ( (graph->degree[loop] == 2) && (used[loop] == 0) ) {
          curprune = *prune;
          ret = extend_forced_path(loop, graph, used, prune, edgestack);
          if (ret == HC_NOT_EXIST) {
            return(HC_NOT_EXIST);
          }
    
          /* if HC found, indicate so and return */
          if (ret == HC_FOUND) {
            return(HC_FOUND);
          }

          if (*prune > curprune)  /* pruning took place, redo pruning check */
            done = 0;
        }

      }  /* end of HC_PRUNE_CYC if statement */

    }  /* end of loop through graph vertices */

  }  while (!done); /* end of prune-graphcheck loop */

  /* check for components */
  if (prunelevel & HC_PRUNE_CONNECT) {
    if (calc_graph_components(graph) > 1)  {
      cout << "tag 2 " << endl;
      return(HC_NOT_EXIST);
    } 
  }
    
  /* check for articulation points */
  if (prunelevel & HC_PRUNE_CUTPOINT) {
    if (check_graph_cutpoints(graph) == CUTPNT_EXIST) {
      return(HC_NOT_EXIST);
    }
  }

  return(HC_NOT_FOUND);
}  

/* function to add a (deleted) edge to an edge stack 
 * v1,v2 specifies the vertex endpoints of the edge
 *
 * if edgestack is NULL then just return
 * (NULL indicates that don't want to save stack information) */
void push_edge_to_stack(int v1, int v2, EdgeStack* edges)
{
  if (edges == NULL) 
    return;

  edges->stack[edges->pointer].v1 = v1;
  edges->stack[edges->pointer].v2 = v2;
  (edges->pointer)++;
}  

/* function that starts with 1 degree 2 vertex, and forms a forced
 * path.  Any edge between the endpoints of this forced path is removed
 * (and further path-extending is then done).  
 * also returns HC_NOT_EXIST if a forced short cycle is found
 *	(start and end vertices are the same)
 * returns HC_FOUND if a forced hamiltonian cycle is found
 * otherwise returns HC_NOT_FOUND
 *
 * prune is incremented for each edge removed
 * each deg 2 vertex in the forced path is marked as used 
 *   (using the used[] array)
 *
 * deleted edges are added to the edgestack structure
 *
 * assumptions:  curvert == degree 2 */
int extend_forced_path(int curvert, Graph* graph, int* used, int *prune, EdgeStack* edgestack)
{
  int length;

  int done;
  int startvert;	/* starting vertex of path */
  int endvert;		/* ending vertex of path */
  int oldsv;		/* previous starting vertex */
  int oldev;		/* previous ending vertex */

  int tmp;


  oldsv = oldev = curvert;

  startvert = graph->nbr[curvert][0];
  endvert = graph->nbr[curvert][1];
  length = 3;
  used[curvert] = used[startvert] = used[endvert] = 1;

  /* main extend-path loop */
  do {
    done = 1;
    /* extend the start vertex */
    while (graph->degree[startvert] == 2) {
      /* get next edge along forced path */
      if (graph->nbr[startvert][0] != oldsv)
        tmp = graph->nbr[startvert][0];
      else
        tmp = graph->nbr[startvert][1];
  
      /* if have a forced cycle, check its length */
      if (tmp == endvert) {
        if (length < graph->numvert)	/* cycle is short */
          return(HC_NOT_EXIST);
        else
          return(HC_FOUND);
      }
  
      oldsv = startvert;
      startvert = tmp;
      length++;
      used[startvert] = 1;
    }
      
    /* extend the end vertex */
    while (graph->degree[endvert] == 2) {
      /* get next edge along forced path */
      if (graph->nbr[endvert][0] != oldev)
        tmp = graph->nbr[endvert][0];
      else
        tmp = graph->nbr[endvert][1];
  
      /* if have a forced cycle, check its length */
      if (tmp == startvert) {
        if (length < graph->numvert)	/* cycle is short */
          return (HC_NOT_EXIST);
        else
          return (HC_FOUND);
      }
  
      oldev = endvert;
      endvert = tmp;
      length++;
      used[endvert] = 1;
    }
    
    /* have forced path, so try to remove an edge inbetween the endpoints 
     * but do this only if current path length is less than total number of
     * vertices*/
    if (length < graph->numvert) {
      if (rm_edge_graph(graph, startvert, endvert) == EDGE_REMOVE) {
        done = 0;	/* keep trying to extend path */
        (*prune)++;
        push_edge_to_stack(startvert, endvert, edgestack);
      }
    }

  } while (!done); /* end of main while loop */

  /* unmark start and end vertices, since they are not deg 2 */
  used[startvert] = used[endvert] = 0;

  return(HC_NOT_FOUND);
} 

/* function to verify that a specified solution is indeed a hamiltonian
 * cycle of the graph
 * solution is represented by an array of vertices
 *
 * returns HC_VERIFY if it is a solution
 * returns HC_NOT_VERIFY if it is not */
int hc_verify_solution(Graph* graph, int solution[]) 
{
  int loop;
  int vertcount[graph->numvert];
  
  /* first verify that the solution is indeed a path of length = numvert */
  for (loop = 0; loop < (graph->numvert) - 1; loop++) {
    if ( check_if_edge(graph, solution[loop], solution[loop+1]) != EDGE_EXIST)
      return(HC_NOT_VERIFY);
  }

  /* check if endpoints are connected (have a cycle) */
  if (check_if_edge(graph, solution[graph->numvert-1], solution[0]) != EDGE_EXIST)
    return(HC_NOT_VERIFY);
 
  /* initialize array */
  for (loop = 0; loop < graph->numvert; loop++) { 
    vertcount[loop] = 0;
  }

  /* count # of times each vertex appears */
  for (loop = 0; loop < graph->numvert; loop++) {
    vertcount[solution[loop]] ++;
  }

  /* check that each vertex only appeared once */
  for (loop = 0; loop < graph->numvert; loop++) {
    if (vertcount[loop] != 1)
      return (HC_NOT_VERIFY);
  }

  return(HC_VERIFY);
}  

// dhu: copied from graphdata.c
/* remove undirected edge (x,y) from the graph
 *
 * returns EDGE_NOTEXIST if edge does not exist
 * returns EDGE_REMOVE if edge was removed */
int rm_edge_graph(Graph* graph,int x, int y)
{
  int loop, loop2;

  /* remove y as a neighbour of x */
  for (loop = 0; loop < graph->degree[x]; loop++) {
    if (graph->nbr[x][loop] == y) {
      for (loop2 = loop+1; loop2 < graph->degree[x]; loop2++)
        graph->nbr[x][loop2-1] = graph->nbr[x][loop2];
      break;
    }
  }

  if (loop == graph->degree[x]) { 
    return(EDGE_NOTEXIST);
  }

  /* remove x as a neighbour of y */
  for (loop = 0; loop < graph->degree[y]; loop++) {
    if (graph->nbr[y][loop] == x) {
      for (loop2 = loop+1; loop2 < graph->degree[y]; loop2++)
        graph->nbr[y][loop2-1] = graph->nbr[y][loop2];
      break;
    }
  }

  if (loop == graph->degree[y]) { 
    cout << "Error: inconsistant edge in rm_edge_graph().\n" << endl;
  }

  (graph->degree[x])--;
  (graph->degree[y])--;

  return (EDGE_REMOVE);
}  

/// dhu: copied from graphdata.c
/* checks whether 2 vertices are joined by an edge
 * (directed edge x -> y)
 * returns EDGE_EXIST if edge is there, returns EDGE_NOTEXIST if no edge */
int check_if_edge(Graph* graph, int x, int y)
{
  for (int loop = 0; loop < graph->degree[x]; loop++) {
    if (graph->nbr[x][loop] == y)
      return (EDGE_EXIST);
  }

  return (EDGE_NOTEXIST);
}  

/// dhu: copied from graphdata.c
/* depth-first-search for component checking
 * this function recursively labels vertex 'v' and all its neighbours
 * as belonging to component 'c', using the order[] array to keep track of
 * which component each vertex belongs to. */
void component_dfs(Graph* graph, int order[], int v, int c)
{
  order[v] = c;
  for (int i = 0; i < graph->degree[v]; i++) {
    if ( order[graph->nbr[v][i]] != c )
      component_dfs(graph, order, graph->nbr[v][i], c);
  }
} 

/// dhu: copied from graphdata.c
/* this function counts the number of components that the specified graph
 * has, and returns this number */
int calc_graph_components(Graph* graph)
{
  int order[graph->numvert];
  for (int i = 0; i < graph->numvert; i++)
    order[i] = 0;

  int j = 0;
  for (int i = 0; i < graph->numvert; i++) {
    if (order[i] == 0) {
      j++;
      component_dfs(graph, order, i, j);
    }
  }

  return j;
}

/// dhu: copied from graphdata.c
/* this function performs depth-first-search for detecting a cutpoint
 * returns CUTPNT_EXIST on first detection of a cutpoint, or returns
 * CUTPNT_NOTEXIST */
int cutpoint_dfs(Graph* graph, int vert, int back[], int dfsnumber[], int *dfnum)
{
  int loop;
  int nextvert;
  
  (*dfnum)++;
  dfsnumber[vert] = *dfnum;
  back[vert] = *dfnum;
  for (loop = 0; loop < graph->degree[vert]; loop++) {
    nextvert = graph->nbr[vert][loop];
    if (dfsnumber[nextvert] == 0) {
      if ( cutpoint_dfs(graph, nextvert, back, dfsnumber, dfnum) == CUTPNT_EXIST) {
        return CUTPNT_EXIST; 
      }
      if (back[nextvert] >= dfsnumber[vert]) {
        /* vertex 'vert' is an articulation point */
        return CUTPNT_EXIST;
      }
      else {
        if (back[nextvert] < back[vert]) 
          back[vert] = back[nextvert];
      }
    } 
    else {
      if (dfsnumber[nextvert] < back[vert]) 
	back[vert] = dfsnumber[nextvert];
    }
  }

  return CUTPNT_NOTEXIST;

}

/// dhu: copied from graphdata.c
/* this function checks to see if a cutpoint (articulation point) exists, 
 * and returns CUTPNT_EXIST if one is found, or CUTPNT_NOTEXIST if one does 
 * _not_ exist in the specified graph
 *
 * this code originally written by Culberson, from Baase
 * assumption:  vertex #1 has degree > 0 */
int check_graph_cutpoints(Graph* graph)
{
  int dfsnumber[graph->numvert];
  int back[graph->numvert];		/* record shallowest back edge point */
  int dfnum;			/* current depth first number */

  int loop;
  int nextvert;
  for (loop = 0; loop < graph->numvert; loop++)
    dfsnumber[loop] = 0;

  /* handle root case */
  back[0] = 1;
  dfsnumber[0] = 1;
  nextvert = graph->nbr[0][0];
  dfnum = 1;
  if (cutpoint_dfs(graph, nextvert, back, dfsnumber, &dfnum) == CUTPNT_EXIST) 
    return CUTPNT_EXIST;

  /* check if root vertex is a cutpoint */
  for (loop = 1; loop < graph->degree[0]; loop++) {
    if (dfsnumber[graph->nbr[0][loop]] == 0) {
      /* root vertex 0 is a cutpoint */
      return CUTPNT_EXIST;
    }
  }

  return CUTPNT_NOTEXIST;
}


//<MMattar> Created the whole function
// Return the bin in which the distance between beginning and end nodes is located
// assumes that the edge exists and that distances are normalized
int which_bin(Graph* graph, DBSeq* myDBSeq, int beginning, int end){

	// define the factor to find the most significant digit of the nodes
	int graph_factor = graph->numvert/myDBSeq->num_label;

	int a = beginning/graph_factor; // retrieve the most significant digit of the current node
	int b = end/graph_factor; // retrieve the most significant digit of the next node
	return myDBSeq->bin[a][b];
}
//</MMattar>


/* copied from backtrack.c
 * backtrack algorithm to find hamiltonian cycle, 
 *
 * function returns HC_FOUND if HC was found, HC_NOT_EXIST if no HC
 *   was found (final return of HC_NOT_EXIST means no HC exists)
 *
 * parameters:
 *   graph -> graph to solve
 *   pstart, pend, plenth, path, graphpath -> specify the hamiltonian path 
 *     being constructed
 *   nodecount -> # of 'nodes' (vertices) expanded
 *      = # of calls to the recursive routine
 *
 * This is the NON-PRUNING VERSION */
//<MMattar> Included parameter myDBSeq and included a conditional statement to check if binning is used
int calc_noprune_bt_alg(
  Graph* graph,
  int *pstart,
  int *pend,
  int *plength,
  vector<Path>& path,
  vector<GraphPath>& graphpath,
  int* nodecount,
  DBSeq* myDBSeq)
{
  int loop;
  int oldpend;
  int newvert;
 
  
  //<MMattar> Included a time-out functionality
  if ((time(NULL) - starttime) > 1)
  {
	  //cout << endl << endl << (time(NULL) - starttime) << endl << endl << endl;
	  cout << "Unable to find a valid Hamiltonian Circuit." << endl;
	  cout << "Elapsed time: " << (time(NULL) - starttime) << " seconds. Trying again..." << endl;
	  start(myDBSeq->num_label, myDBSeq->word_len, myDBSeq); // re-start the program and try again
  }
  //</MMattar>
  
  
  (*nodecount)++;

  /* if have hamiltonian path, then try to convert into a cycle */
  if (*plength == graph->numvert) 
  {
	  return(hc_path_to_cycle(graph, path, graphpath, pstart, pend, *plength) );
  }

  /* loop through neighbours of current endpoint, trying each in turn */

  //<MMattar>
  if (myDBSeq->num_bins != 0){
	  // figure out where in the guide function we are
	  int test_bin = 0;
	  // now we run num_bins loops, one for each bin
	  for (int b=0; b<myDBSeq->num_bins; b++){ // for each bin
		  for (loop = 0; loop < graph->degree[path[*pend].gvert]; loop++){
			  /* if neighbor is not in path, then add to path and recurse */
			  newvert = graph->nbr[path[*pend].gvert][loop];

			  // checks if newvert is in the preferred bin
			  test_bin = which_bin(graph, myDBSeq, path[*pend].gvert, newvert);
			  if((test_bin==(myDBSeq->pref_bin[*plength-1][b][0])) || (test_bin==(myDBSeq->pref_bin[*plength-1][b][1]))){
				  if (graphpath[newvert].pathpos == -1){

					  // saves important sequences for printing out if in debug mode
					  myDBSeq->bins_used[*plength-1] = which_bin(graph, myDBSeq, path[*pend].gvert, newvert);
					  myDBSeq->trans_seq[*plength-1] = myDBSeq->nm1[(path[*pend].gvert)/(graph->numvert/myDBSeq->num_label)][newvert/(graph->numvert/myDBSeq->num_label)];
					  myDBSeq->trans_seq2[*plength-1] = myDBSeq->nm2[(path[*pend].gvert)/(graph->numvert/myDBSeq->num_label)][newvert/(graph->numvert/myDBSeq->num_label)];
					  myDBSeq->trans_seq3[*plength-1] = myDBSeq->nm3[(path[*pend].gvert)/(graph->numvert/myDBSeq->num_label)][newvert/(graph->numvert/myDBSeq->num_label)];

					  /* add current neighbor to path and recurse */
					  oldpend = *pend;
					  add_vert_to_path(path, graphpath, pstart, pend, plength, newvert);
					  if (calc_noprune_bt_alg(graph, pstart, pend, plength, path, graphpath, nodecount, myDBSeq) == HC_FOUND){
						  return(HC_FOUND);
					  }
					  /* current try was bad, so backup (remove current neighbour from path) */
					  remove_endvert_from_path(path, graphpath, pstart, pend, plength, oldpend);
				  }
			  }
		  }  /* end of loop through neighbours of endpoint */
	  }
  }
  else{
	  for (loop = 0; loop < graph->degree[path[*pend].gvert]; loop++){
		  /* if neighbor is not in path, then add to path and recurse */
		  newvert = graph->nbr[path[*pend].gvert][loop];

		  if (graphpath[newvert].pathpos == -1){
			  /* add current neighbor to path and recurse */
			  oldpend = *pend;
			  add_vert_to_path(path, graphpath, pstart, pend, plength, newvert);

			  if (calc_noprune_bt_alg(graph, pstart, pend, plength, path, graphpath, nodecount, myDBSeq) == HC_FOUND){
				  return(HC_FOUND);
			  }

			  /* current try was bad, so backup (remove current neighbour from path) */
			  remove_endvert_from_path(path, graphpath, pstart, pend, plength, oldpend);
		  }
	  }  /* end of loop through neighbours of endpoint */
  }
  //</MMattar>

	return (HC_NOT_EXIST);

}  /* end of calc_noprune_bt_alg() */

/* Copied from hamcycle.c
 * this function trys to transform a given (hamiltonian) path to 
 * a (hamiltonian) cycle
 * note that the path does not have to be hamiltonian for the algorithm
 * to work.  it will just try to form a cycle from the given path
 *
 * returns HC_FOUND if successfull, HC_NOT_EXIST if unsuccessfull
 *
 * if HC_FOUND, the given path is modified so that the end vertex has an
 *   edge to the start vertex.
 *
 * note that in modifying the path, an edge will be 'removed' from the path
 * (between new end and predecessor of old end in new path)
 * this edge _cannot_ be forced, since both vertices must be at least
 * degree 3 in order for this transformation to work */
int hc_path_to_cycle(
  Graph* graph,
  vector<Path>& path,
  vector<GraphPath>& graphpath,
  int *pstart, 
  int *pend,
  int plength)
{
  int loop;
  int tempvert, curvert;
  int tempnum;
 
  /* check if have a simple hamiltonian path */
  if ( check_if_edge(graph, path[*pend].gvert, path[*pstart].gvert) 
       == EDGE_EXIST) 
  {
    return(HC_FOUND);
  }

  /* for each neighbour of end vertex, check to see if next vertex to it in
   * path has edge with start.  if so, then can construct a cycle
   */
  for (loop = 0; loop < graph->degree[path[*pend].gvert]; loop++) {
    tempvert = graph->nbr[path[*pend].gvert][loop];

    /* get next vertex in path */
    tempnum = graphpath[tempvert].pathpos;
    tempnum = path[tempnum].next;
    curvert = path[tempnum].gvert;

    /* check to see if this vertex has an edge with the start */
    if (check_if_edge(graph, curvert, path[*pstart].gvert) == EDGE_EXIST)
    {
      /* do reversal of path */
      hc_reverse_path(path, graphpath, pend, graphpath[tempvert].pathpos, 
			plength);

      return(HC_FOUND);
    }

  }  /* end of loop through neighbours of end-vertex on path */

  return(HC_NOT_EXIST);

}  /* end of hc_path_to_cycle() */
  

/* Copied from hamcycle.c
 * function to add a vertex to the current path
 *   vert = new vertex
 *   assumes vertices stored in path[] from 0 to *plength-1 */
void add_vert_to_path(
  vector<Path>& path,
  vector<GraphPath>& graphpath,
  int *pstart,
  int *pend,
  int *plength,
  int vert)
{

  path[*pend].next = *plength;
  path[*plength].gvert = vert;
  path[*plength].next = -1;
  graphpath[vert].pathpos = *plength;
  *pend = *plength;
  (*plength)++;

}  /* end of add_vert_to_path() */

/* Copied from hamcycle.c
 * function to remove the end vertex of the current path
 *   must specify second-from-the-end vertex (oldend)
 *   assumes vertices stored in path[] from 0 to *plength-1
 *
 * @@ this will break if the path order is messed with.  to fix, must
 * move vertex on the end (at location *plength-1) to the location of the
 * deleted endpoint. */
void remove_endvert_from_path(
  vector<Path>& path,
  vector<GraphPath>& graphpath,
  int *pstart,
  int *pend,
  int *plength,
  int oldend)
{
  graphpath[path[*pend].gvert].pathpos = -1;
  path[*pend].gvert = -1;
  path[*pend].next = -1;
  *pend = oldend;
  (*plength)--;
  path[*pend].next = -1;
}  /* end of remove_endvert_from_path() */

/* Copied from hamcycle.c
 * rotational transformation to reverse path in cycle.
 * ie: A-B-C-D, edge from D-A
 * new path: A-D-C-B  (B is new end of path)
 *
 * endpathv = path index of end-of-path, which is changed
 * revpathv = path index of start of reverse path (A)
 * plength = length of path
 *
 * included in hamcycle.c because it is used in the generic
 *   path_to_cycle() function */
void hc_reverse_path(
  vector<Path>& path,
  vector<GraphPath>& graphpath,
  int *endpathv,
  int revpathv,
  int plength)
{
  int i, j;
  int tempnum;
  int newend;

  newend = path[revpathv].next;
  path[revpathv].next = *endpathv;

  i = newend;
  j = path[i].next;

  do
  {    
    /* make j point to i, but save what j points to first */
    tempnum = path[j].next;
    path[j].next = i;

    /* advance i and j down path */
    i = j;
    j = tempnum;

  } while (j != -1);  /* stop when i = oldend, j points to nothing */

  /* set up new end of path */
  path[newend].next = -1;
  *endpathv = newend;
  graphpath[path[newend].gvert].ended = plength;   /* for posa's only */

}  /* end of hc_reverse_path() */
