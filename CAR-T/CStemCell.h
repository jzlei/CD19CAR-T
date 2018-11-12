#ifndef CSTEMCELL_H
#define CSTEMCELL_H

#include "Cell.h"

struct PARStemCell{
    double beta0, kappa0, m;
    double mu0, mu1, mu2;
    double alpha34, alpha19, alpha22, alphaLNK, alpha123;
    double gamma19, gamma22, gamma123;
    double p0, plost;  // p0: the probability of CD19 mutation; plost: the probability of CD19 lost
    double theta;
    double X0, n0;
    double X1, n1;
    double a123,a34,a22;
    double cd34b;
};   // System parameters


class CStemCell:public CCell{

private:

    double CD34, CD19, CD22, LNK, CD123;         // The gene expression levels _X[0]=CD34, _X[1]=CD19, _X[2]=CD22, _X[3]=LNK, _X[4]=CD123
    double _q[5];       // The epigenetic states. i=0 for CD34, i=1 for CD19, i=2 for CD22, i=3 for LNK, i=4 for CD123
    double _CARSignal;  // The strength of CAR-19 signal.
    
	PARStemCell _par;
 
    double GetnextEpi(int i, double y, double R);      // Get the epigenetic state of the daughter cell given the value y after DNA duplication.
    double fdeathrate(double R19, double R123);          // Get the death rate of cells.
    double fbeta(long N0, double R19, double R123);                      // The proliferation rate
    double fkappa();                      // The differentiation rate

	void ReadPar(char fpar[]);
    void SetDefaultPar();
    int GetCellType(double x);     // Get the type of a cell.
    
	// Begin Stochastic simulation
	void GetRandForces(double b[][NUMRAND], double x[],int k);
	void GetTrends(double a[], double x[]);
	// End Stochastic simulation	
	
	// Begin Gillespie algorithm
	void Propensities(double a[]);
	void Update(int mu);
	// End of Gillespie algorithm
    
    void OutPutParameters();
    
	
public:
    int _mutant;   // marker of CD19 mutation: 0:  WT, 1: CD19 mutant.
    int _celltype;  // marker for cell type, 0: for Stem cell; 1: for precursor cell; 2: terminal differentiated cell
    double _plost;  // The probability of CD19 lost.


	CStemCell();
	~CStemCell();
	
	bool Initialized(int k, char fpar[]);
    DCell CellFateDecision(long N0, double R19, double R123, double dt);
    
    void Pre_OneStep();
    bool MitosisQ;    

    
    
	friend class CSystem;
};

#endif
