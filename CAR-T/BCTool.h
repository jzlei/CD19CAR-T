#ifndef BCTOOL_H
#define BCTOOL_H

//Definition of functions
#define _SQUARE(x) ((x)*(x))
#define _MAX(x,y) ((x)>=(y) ? (x):(y))
#define _MIN(x,y) ((x)<=(y) ? (x):(y))
#define _Heaviside(x) ((x)>(0) ? (1.0):(0.0))

//Definition of constants
#define ComLength 100
#define StrLength 1024
#define UNITMAX 2147483647  
#define MAXCELL 1000000
#define NUMREACT 0
#define NUMVAR 5
#define NUMRAND 2
#define PI 3.14159
#define POSITIVE 0
#define HSC 0
#define PREC 1
#define TD 2

struct IMD{
	char mdcrd[ComLength];	// trajectory file.
	char cellpar[ComLength]; // cell parameter file.
	double N;                  // Number of all cells (in 10^6).
    int N0;                 // Number of cells in the pool.
	int seed;				// Seed of the random numbers.
	double dt;              // Dynamical parameters.
	double T1;              // Dynamical parameters (period of one cell cycle).
	int ntpr;               // Output the result, if ntpr = 0, then no output.
	int ntpx;                // Output the path result for every ntpx steps. if ntpx = 0, then no result is outputed. 
};					    // MD parameters;

struct DCell{
    int type;   // type = 0: unchanged; type = 1: division; type = 2: cell death; type = 3: CD19 mutation; type = 4: differentiation
    int mutanttype1, mutanttype2;   // Cell type of daughter cells, 0 for wild type, 1 for mutant.
    int celltype1, celltype2;       // Cell tyoe of daughter cells, 0 for HSC, 1 for precursor cells, 2 for ternimal differentiated cells.
    double X1[NUMVAR];  // Protein number state of the daughter cells.
    double X2[NUMVAR];
    double q1[NUMVAR];     // Epigenetic state of the daughter cells.
    double q2[NUMVAR];
};

#endif /* defined (___main__) */


