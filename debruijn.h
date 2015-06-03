/************************************************************
 * Program that generates de Bruijn sequences. 
 * This program is based on Hamiltonian Cycle Program written by Basil Vandegriend 
 * For more details, please refer to:
 * http://webdocs.cs.ualberta.ca/~joe/Theses/vandegriend.html

 * Modified by Dongbo Hu
 * dongbo@mail.med.upenn.edu
 *
 * Modified by Marcelo Mattar (09/24/2010)
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

#ifndef DEBRUIJN_H
#define DEBRUIJN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <strings.h>
#include "getopt.h"
#include <time.h>  // added by dhu
#include <sstream>
#include <string>
#include <vector>

/* options data structure:  
 *   note that this is a global, so is accessible by the entire
 *   software package.
 *
 *   note also that this structure contains option data structures
 *     for the different modules, so that the options.h file must
 *     be included after the other module .h files.  */

//////////// backtrack.h /////////////////////////

/* degsortflag parameters */
#define DEGSORT_RAND 0
#define DEGSORT_MIN 1
#define DEGSORT_MAX 2

#define RUN_NORMAL 0
#define RUN_TIMELIMIT 1

/* defines for prunelevel in hc_do_pruning() */
#define HC_PRUNE_NONE  0x0
#define HC_PRUNE_BASIC 0x1
#define HC_PRUNE_CYC   0x2
#define HC_PRUNE_CONNECT 0x4
#define HC_PRUNE_CUTPOINT 0x8
#define HC_PRUNE_ALL (HC_PRUNE_BASIC | HC_PRUNE_CYC | HC_PRUNE_CONNECT | \
HC_PRUNE_CUTPOINT)

/* defines for the possible return values from the posa algorithm
 * hc_not_found - no cycle found, but one may exist
 * hc_not_exist - no cycle can exist on this graph */
#define HC_FOUND 0
#define HC_NOT_EXIST 1
#define HC_NOT_FOUND 2

/* defines for return values for hc_verify_solution */
#define HC_NOT_VERIFY 0
#define HC_VERIFY 1

/* values for selectflag parameter for select_initvert() */
#define INITVERT_RANDOM 0
#define INITVERT_MAXDEG 1
#define INITVERT_RANDEG 2
#define INITVERT_FIRST 3

/* return values for hc_check_timelimit() */
#define HC_QUIT 0
#define HC_CONTINUE 1

struct BackTrackOpt {
  BackTrackOpt() { init(); }
  void init();

  int initvertflag;
  int degsortflag;
  int pruneoptflag;
  int restart_increment;
  int max_nodes;
};

//////////////// heuristic.h /////////////////////////////

/* visitflag parameters */
#define VISIT_SMART 0
#define VISIT_RAND 1

/* completeflag parameters (for converting from line to cycle, to 
 * complete the solution) */
#define COMPLETE_NORM 0
#define COMPLETE_SMART 1

/* edgepruneflag parameters:  how many edges to remove initially */
#define EDGEPRUNE_ALL 0
#define EDGEPRUNE_BASIC 1
#define EDGEPRUNE_NONE 2

/* cycleextend parameters */
#define NOCYCLEEXTEND 0
#define CYCLEEXTEND 1
/* defines for return values for hc_verify_solution */
#define HC_NOT_VERIFY 0
#define HC_VERIFY 1

/* defines for prunelevel in hc_init_pruning() */
#define INIT_PRUNE_ALL 0
#define INIT_PRUNE_3CYC 1

struct HeuristicOpt {
  HeuristicOpt() { init(); }
  void init();

  int visitflag;
  int completeflag;
  int cycleextendflag;
};

//////// main.h //////////////////

/* math type defines */
#define EPSILON (0.0000001)
#define PI (3.14159)
#define LARGEVAL (10000000)

/* string type defines */
#define STRLEN (256)


/* return codes */
#define RET_OK 0
#define RET_ERROR 1

/* defines for print_info_summary() function */
#define PRINT_INFO_HEADER 	0
#define PRINT_INFO_ALL 		1

/////// graphdata.h ///////////

#define MAXVERT 1600
#define MAXDEGREE 50

/* graph file information */
#define GFILE_EXT ".graph"
#define GFILE_COMMENT "#"

/* defines for function return values */
#define CUTPNT_NOTEXIST 0
#define CUTPNT_EXIST 1

/* defines for 2 functions: rm_edge_graph, check_if_edge */
#define EDGE_NOTEXIST 	0
#define EDGE_REMOVE 	1
#define EDGE_EXIST	1

struct Graph {
  Graph() { init(); }
  void init();
  void setup(int num_label, int word_len);
  void printNeighbors();

  int numvert; 	/* total number of vertices */
  vector<int> degree;	/* degree (number of connections) of each vertex */

  vector<vector<int> > nbr;     /* neighbors of each vertex */

  int numedges;			/* total number of edges */
  float meandeg;		/* mean vertex degree */
  float stddevdeg;		/* deviation in vertex degree */
  int mindeg;			/* minimum vertex degree */
  int maxdeg;			/* maximum vertex degree */

  vector<int> deghistogram;	/* histogram of vertex degrees */

  int solve;			/* solvable flag for current problem */
				/* use HC_FOUND, HC_NOT_EXIST */
};


//<MMattar> Created the structure DBSeq
// contains information about the DeBruijn sequence
struct DBSeq {
  DBSeq() { init(); }
  void init();
  void setup(int num_label, int word_len, int num_bins, char* NeuralModel1_filename, char* NeuralModel2_filename, char* NeuralModel3_filename, char* GuideFunction_filename);

  // relevant variables for the carry-over design
  int num_label; // number of labels (k)
  int word_len; // word length; order of the sequence (n)
  int num_bins; // number of bins (B)
  int num_neuralmodels; // number of neural models used
  vector<double> guide_function; // guide function for selecting the bins
  vector<int> seq_bin; // sequence of bins, or, equivalently, discretization(binning) of the guide_function
  vector<double> trans_seq; // sequence of transitions obtained after running the algorithm
  vector<double> trans_seq2; // sequence of transitions obtained after running the algorithm
  vector<double> trans_seq3; // sequence of transitions obtained after running the algorithm
  vector<int> bins_used; // sequence of bins used when constructing the final sequence
  vector<vector<double> > nm1; // table of neural model
  vector<vector<double> > nm2; // table of neural model
  vector<vector<double> > nm3; // table of neural model
  vector<vector<double> > bin; // table of bins
  vector<vector<vector<int> > > pref_bin; // 3-dimensional matrix of preferred bins (k^n x B x 2)

};
//</MMattar>


///////// stats.h ////////////////

#define MAXTRIALTESTS 10
#define MAXGRAPHTESTS 10000

struct Stat {
  float ave;
  float stddev;
};

struct TrialStat {
  TrialStat() {
    result = HC_NOT_FOUND;
    time = 0;
    nodes = 0;
    edgeprune = 0;
    initprune = 0;
    retries = 0;
  }

  int result;  /* = HC_FOUND, HC_NOT_FOUND, HC_NOT_EXIST */
  float time;
  int nodes;
  int edgeprune;
  int initprune;
  int retries;
};

struct GraphStat {
  GraphStat() : trial(MAXTRIALTESTS) { }

  vector<TrialStat> trial;
  int graphham;      /* = HC_FOUND, HC_NOT_FOUND, HC_NOT_EXIST */
  int biconnected;   /* = 1 if biconnected, = 0 if not */
  int mindeg2;	     /* = 1 if min degree >= 2, = 0 if not */
};

struct ExpStat {
  ExpStat() : graph(MAXGRAPHTESTS) { }

  vector<GraphStat> graph; 
  /* node statistics */
  Stat nodes_tot;
  Stat nodes_ham;
  Stat nodes_noham;
  Stat nodes_nofound;

  Stat noderatio_tot;
  Stat noderatio_ham;
  Stat noderatio_noham;
  Stat noderatio_nofound;

  /* heuristic algorithm statistics */
  Stat algsuccess;
  Stat tmsuccess;
  Stat tmfail;
  Stat tmexpect;
  
  /* backtrack algorithm statistics */
  Stat tmham;
  Stat tmnoham;
  Stat tmtotal;
  Stat perham;
  Stat perbiconnect;
  Stat permindeg2;
};

//////////// hamcycle.h /////////////////////////


/* to be indexed by path position */
struct Path {
  int gvert;    /* graph vertex # */
  int next;     /* points to location in this array of next vertex in path */
};

/* to be indexed by graph vertex # */
struct GraphPath {
  int pathpos; /* position in path (-1 if not in path) */
  int ended;   /* length of path ended on this vertex (to avoid repeats) */
               /* for posa's only */
};

struct Edge {
  int v1;
  int v2;
};

struct EdgeStack {
  EdgeStack() {
    pointer = 0;
  }

  std::vector<Edge> stack;
  int pointer;
};

////////////  GraghGen.h ///////////////////////

#define CONNECT_NEAR 0
#define CONNECT_FAR 1

#define GRAPH_WRAP 0
#define GRAPH_NOWRAP 1

#define MAXBOARDSIZE 100

#define HAM_DONTCARE 0
#define HAM_ENSURE 1

#define MAXNUMADDPATHS 2

struct GraphGenOpt {  
  GraphGenOpt() { init(); }
  void init();

  /* common options */
  int nvertex;		/* number of vertices */
  int makeham;		/* make sure graph is hamilton (HAM_ENSURE or HAM_DONTCARE) */

  /* geometric graph options */
  float dist;		/* distance */
  int dim;		/* # of dimensions */
  int dflag;		/* distance flag = CONNECT_NEAR or CONNECT_FAR */
  int wrapflag;		/* wrap flag = GRAPH_WRAP or GRAPH_NOWRAP */
  int mindeg;		/* minimum degree that graph is forced to have */
  
  /* degreebound graph options */
  int degsize;		/* largest degree with non-zero percentage */
  float degpercent[MAXDEGREE];	/* array of percentages for each degree */

  /* knighttour parameters */
  int board1;
  int board2;
  int move1;
  int move2;

  /* crossroads parameters */
  int numsubgraphs;

  /* random graph parameters */
  float meandeg;	/* mean vertex degree */
  float degconst;	/* meandeg = degconst * (ln n + ln ln n) */

  /* add-cycle graph parameters */
  float numcycles;

  /* add-path graph parameters */
  float pathlengths[MAXNUMADDPATHS];

  /* ICCS graph parameters */
  int indsetsize;
};


///////// options.h ///////////////////

/* maximum length of option words */
#define OPTLEN	15

/* algorithms */
#define ALG_NOSOLVE 	0
#define ALG_NOPRUNE_BT	1
#define ALG_BACKTRACK	2
#define ALG_POSA_HEUR	3

#define NUM_ALG_OPT	4

#ifdef IN_OPTIONS_FILE
char opt_alg_str[NUM_ALG_OPT][OPTLEN] = {
                "NOSOLVE",
                "noprune_bt",
                "backtrack",
                "posa_heur" };
#else
extern char opt_alg_str[NUM_ALG_OPT][OPTLEN];
#endif

/* graph gen type */
#define GEN_NOGRAPH	0
#define GEN_GEOMETRIC 	1
#define GEN_DEGREEBOUND	2
#define GEN_KNIGHTTOUR	3
#define GEN_CROSSROADS	4
#define GEN_RANDOM	5
#define GEN_ADDCYCLE	6
#define GEN_ADDPATH	7
#define GEN_ICCS	8

#define NUM_GEN_OPT	9

#ifdef IN_OPTIONS_FILE
char opt_gen_str[NUM_GEN_OPT][OPTLEN] = {
                "NOSOLVE",
                "GEOMETRIC",
                "DEGREEBOUND",
                "KNIGHTTOUR",
		"CROSSROADS",
		"RANDOM",
		"ADDCYCLE",
		"ADDPATH",
		"ICCS"};
#else
extern char opt_gen_str[NUM_GEN_OPT][OPTLEN];
#endif

/* report flags - indicate which statistics to calculate/report */
#define REPORT_NONE	0
#define REPORT_GRAPH	0x01
#define REPORT_ALG	0x02
#define REPORT_SOLUTION	0x04
#define REPORT_OPTIONS	0x08
#define REPORT_SUMMARY	0x10

/* savegraph flag */
#define NOSAVEGRAPH 0
#define SAVEGRAPH 1

/* return codes from read_next_word() */
#define READ_OK		0x0
#define READ_STRLONG	0x1
#define READ_EOL	0x2
#define READ_EOF	0x4

/* flags to indicate type of word in option file */
#define WORD_NULL	0
#define WORD_EMPTY	1
#define WORD_ARG	2
#define WORD_PARM	3
#define WORD_COMMENT	4
#define WORD_OTHER	5

struct SeqOpt {
  SeqOpt() { init(); }
  void init();

  int algorithm;

  /* maximum amount of time to let algorithm run for */
  int alg_timelimit;

  HeuristicOpt heur_alg;
  BackTrackOpt bt_alg;
  GraphGenOpt graphgen;
  
  int graphgentype;  
  int rng_seed;      /* random number generator seed */
  int report_flags;  /* indicate which things to report */

  /* options for where to send the output */
  std::string log_str;
  std::string sol_str;
  std::string stat_str;
  std::string options_str;
  std::string summary_str;
};

// Some constants
struct SeqConst {
  static const int num_graph_tests = 1;
  static const int num_instance_tests = 1;
};

///////// my functions //////////////
void setOptions();
void usage();


//<MMattar> Included the parameter DBSeq* myDBSeq
void start(int num_label, int word_len, DBSeq* myDBSeq);
int test_hc_alg(int num_label, int word_len, Graph* graph, DBSeq* myDBSeq, TrialStat* trialstats);
int master_backtrack_alg(Graph* graph, DBSeq* myDBSeq, TrialStat* trialstats, int solution[]);
int calc_noprune_bt_alg(Graph* graph, int *pstart, int *pend, int *plength,
			vector<Path>& path, vector<GraphPath>& graphpath, int *nodecount, DBSeq* myDBSeq);
//</MMattar>


bool chkSolution(int num_label, int word_len, int solution[]);

void test_graph_properties(Graph* graph, GraphStat* graphstat);
void printSeq(int solution[], int n_vert, int factor);
void printSeq_v(int solution[], int n_vert, int factor, int word_len);

int select_initvertex(Graph* graph, int selectflag);
int hc_do_pruning(Graph* graph, int *prune, int prunelevel, EdgeStack* edgestack);
void push_edge_to_stack(int v1, int v2, EdgeStack* edges);
int extend_forced_path(int curvert, Graph *graph, int* used, int *prune, EdgeStack* edgestack);
int hc_verify_solution(Graph* graph, int solution[]);
int rm_edge_graph(Graph* graph,int x, int y);
int check_if_edge(Graph* graph, int x, int y);
void component_dfs(Graph* graph, int order[], int v, int c);
int calc_graph_components(Graph* graph);
int cutpoint_dfs(Graph* graph, int vert, int back[], int dfsnumber[], int *dfnum);
int check_graph_cutpoints(Graph* graph);

int hc_path_to_cycle(Graph* graph, vector<Path>& path, vector<GraphPath>& graphpath, 
		     int *pstart, int *pend, int plength);

void add_vert_to_path(vector<Path>& path, vector<GraphPath>& graphpath, 
		      int *pstart, int *pend, int *plength, int vert);

void remove_endvert_from_path(
  vector<Path>& path,
  vector<GraphPath>& graphpath,
  int *pstart,
  int *pend,
  int *plength,
  int oldend);

void hc_reverse_path(
  vector<Path>& path,
  vector<GraphPath>& graphpath,
  int *endpathv,
  int revpathv,
  int plength);

/* Function template that converts a non-string type into C++ string
 * It is mainly used to convert a numeric value (such as int32 to string). */
template<typename T>  
string num2str(T inputNum) 
{
  string foo;
  stringstream out;
  out << inputNum;
  foo = out.str();

  return foo;
}

#endif
