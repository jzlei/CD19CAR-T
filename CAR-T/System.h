#ifndef CSYSTEM_H
#define CSYSTEM_H

class CSystem{
private:
    int _NumPoolCell;   // Number of cells in the simulation pool.
    double _NumCell;		// Number of total cells in the animal (in unit 10^6).
    int _MaxNumCell;    // Maximal cell number.
    double _R;          // CAR-T activity
    double _R19, _R123; // Activities of _R19 and _R123
    
    double _Prolif;       // Proliferation rate
    int _N0, _N2, _N1, _N3, _N4;
	CStemCell *_cells;	// Cells.

	void OutPutCells(double t);
	void OutPutSys(int step);
	void OutPutSys(FILE *fp);
    void RunOneStep(double t);
    double CARTActivity(double t);  // CAR-T Activity;
    
    DCell nextcell;
	
public:
	CSystem();
	~CSystem();
	
	CSystem(int N0);

	CStemCell& operator()(int indx); // Edit to return the type of your own class.

	bool Initialized();
	void Run();
    
    bool SystemUpdate(double t);
};

#endif
