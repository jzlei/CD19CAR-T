#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "string.h"
#include <iostream>
using namespace std;

#include "BCTool.h"
#include "Cell.h"
#include "CStemCell.h"
#include "System.h"
#include "Random.h"

extern struct IMD _MD;
extern double _R0, _tau1, _tau2, _ratio19, _rho;
extern CRandom Rand;

CSystem::CSystem()
{
	_NumCell = 0;
}

CSystem::CSystem(int N0)
{
	_NumCell = _MD.N;
    _Prolif=1.0;
    
    _NumPoolCell = N0;
    _MaxNumCell = MAXCELL;
    
    _cells = new CStemCell[_MaxNumCell];
}

CSystem::~CSystem()
{
     delete _cells;
}

bool CSystem::Initialized()
{
	int k;
	k=1;
	do
	{
		if((*this)(k).Initialized(k, _MD.cellpar))
		{
			k++;
		}
		else
		{
		}
	}while(k<=_MaxNumCell);

	return(true);
}

bool CSystem::SystemUpdate(double t)
{
    int k;
    int Ntemp;
    double *X34, *X19, *X22, *LNK, *X123;
    double *p, *q, *r, *w, *z;
    int *celltype;
    int *mutanttype;
    X34 = new double[2*_NumPoolCell+1];
    X19 = new double[2*_NumPoolCell+1];
    X22 = new double[2*_NumPoolCell+1];
    X123 = new double[2*_NumPoolCell+1];
    LNK = new double[2*_NumPoolCell+1];
    p = new double[2*_NumPoolCell+1];
    q = new double[2*_NumPoolCell+1];
    r = new double[2*_NumPoolCell+1];
    w = new double[2*_NumPoolCell+1];
    z = new double[2*_NumPoolCell+1];
    celltype = new int[2*_NumPoolCell+1];
    mutanttype = new int[2*_NumPoolCell+1];
    Ntemp=0;        // Number of cells after a cycle.
    _N0=0;          // Number of cells in the pool after cell fate decision.
    _N2=0;          // Number of death cells removed from TA cells
    _N1=0;          // Number of CD19 mutant cells
    _N3=0;          // Number of divided cells
    _N4=0;          // Number of terminal differentiated cells
    
    for (k=1; k<=_NumPoolCell; k++)
    {
        nextcell=(*this)(k).CellFateDecision(_NumCell, _R19, _R123, _MD.dt);
        switch(nextcell.type){
            case 4:                 // Cells with ternimal differentiation
                _N4++;
            case 0:                 // Cells remain unchanged
            case 3:                 // Cells with CD19 mutation
                Ntemp++;
                X34[Ntemp] = nextcell.X1[0];
                X19[Ntemp] = nextcell.X1[1];
                X22[Ntemp] = nextcell.X1[2];
                LNK[Ntemp] = nextcell.X1[3];
                X123[Ntemp] = nextcell.X1[4];
                p[Ntemp] = nextcell.q1[0];
                q[Ntemp] = nextcell.q1[1];
                r[Ntemp] = nextcell.q1[2];
                w[Ntemp] = nextcell.q1[3];
                z[Ntemp] = nextcell.q1[4];
                celltype[Ntemp] = nextcell.celltype1;
                mutanttype[Ntemp] = nextcell.mutanttype1;
                _N0++;
                _N1++;
                break;
            case 1:                 // Cell division
                Ntemp++;
                X34[Ntemp] = nextcell.X1[0];
                X19[Ntemp] = nextcell.X1[1];
                X22[Ntemp] = nextcell.X1[2];
                LNK[Ntemp] = nextcell.X1[3];
                X123[Ntemp] = nextcell.X1[4];
                p[Ntemp] = nextcell.q1[0];
                q[Ntemp] = nextcell.q1[1];
                r[Ntemp] = nextcell.q1[2];
                w[Ntemp] = nextcell.q1[3];
                z[Ntemp] = nextcell.q1[4];
                celltype[Ntemp] = nextcell.celltype1;
                mutanttype[Ntemp] = nextcell.mutanttype1;

                Ntemp++;
                X34[Ntemp] = nextcell.X2[0];
                X19[Ntemp] = nextcell.X2[1];
                X22[Ntemp] = nextcell.X2[2];
                LNK[Ntemp] = nextcell.X2[3];
                X123[Ntemp] = nextcell.X2[4];
                p[Ntemp] = nextcell.q2[0];
                q[Ntemp] = nextcell.q2[1];
                r[Ntemp] = nextcell.q2[2];
                w[Ntemp] = nextcell.q2[3];
                z[Ntemp] = nextcell.q2[4];
                celltype[Ntemp] = nextcell.celltype2;
                mutanttype[Ntemp] = nextcell.mutanttype2;

                _N0 = _N0+2;
                _N3 = _N3+1;
                
                break;
            case 2:                 // Cell death
                _N2++;
                break;
        }
        
    }
    
    _Prolif = 1.0*Ntemp/_NumPoolCell;
    _NumCell =_NumCell * _Prolif;
    
    if(Ntemp==0)
    {
        return(0);
    }
    
    int i;
    double p0;
    if(t > 30 && _Prolif > 1 && _NumCell>0.1)
    {
        _MaxNumCell = 10000;
    }

    p0 = 1.0*_MaxNumCell/Ntemp;
    k=0;
    for (i=1; i<=Ntemp; i++)
    {
        if(Rand()<p0 && k<_MaxNumCell)
        {
            k=k+1;
            (*this)(k)._X[0] = X34[k];
            (*this)(k)._X[1] = X19[k];
            (*this)(k)._X[2] = X22[k];
            (*this)(k)._X[3] = LNK[k];
            (*this)(k)._X[4] = X123[k];
            (*this)(k)._q[0] = p[k];
            (*this)(k)._q[1] = q[k];
            (*this)(k)._q[2] = r[k];
            (*this)(k)._q[3] = w[k];
            (*this)(k)._q[4] = z[k];
            if(Rand()<(*this)(k)._plost && (*this)(k)._q[1]>0.6)
            {
//                printf("%f\n",(*this)(k)._plost);
                (*this)(k)._q[1] = Rand(0,0.2);
                (*this)(k)._X[1] = (*this)(k)._q[1];
            }
            (*this)(k)._celltype = celltype[k];
            (*this)(k)._mutant = mutanttype[k];
            if((*this)(k)._mutant==1) _N1++;
        }
        if (k==_MaxNumCell)
        {
            break;
        }
    }
    _NumPoolCell = k;
    return(1);
}

void CSystem::Run()
{
    double t;
    int step;
    int k;
    FILE *fmdsys;
    char fn[ComLength];
    sprintf(fn,"%s.dat",_MD.mdcrd);
    if((fmdsys = fopen(fn,"w"))==NULL)
    {
        cout<<"Cannot open the file mdsys."<<endl;
        exit(0);
    }
    
    step=0;
    OutPutSys(step);

    for (t=0; t<=_MD.T1; t=t+_MD.dt)
    {
        step=step+1;
        _R = CARTActivity(t/24.0);
        _R19 = _ratio19 * _R;
        _R123 = (1-_ratio19) * _R;
        for (k=1; k<=_NumPoolCell; k++){
            (*this)(k).Pre_OneStep();
        }
        if(SystemUpdate(t)){
            fprintf(fmdsys,"%f %f %f %d %d %d %d %d %d %f\n",t,_Prolif, _NumCell, _NumPoolCell, _N0, _N1, _N2, _N3, _N4,_R);
            if((_MD.ntpx>0) && (step%_MD.ntpx==0))
            {
                OutPutSys(step);
            }
        }
        else{
            printf("Tumor cells cleaned at day %4.2f\n",t/24);
            break;
        }
    }
    fclose(fmdsys);
}

void CSystem::RunOneStep(double t)
{
}

double CSystem::CARTActivity(double t)
{
    double y;
    y=_R0 * (_rho + (1.0-_rho)*(1.0 + pow(t/_tau1,0.2))/(1.0 + pow(t/_tau2, 10.0)));
      
    return(y);
}

void CSystem::OutPutSys(int step)
{
    FILE *fp;
    int k;
    char fnc[StrLength];
    
    sprintf(fnc,"%s-%d.dat",_MD.mdcrd,step);
    if((fp = fopen(fnc,"w"))==NULL)
    {
        cout<<"Cannot open the file fmdc."<<endl;
        exit(0);
    }

	
	for(k = 1; k<= _NumPoolCell; k++)
	{
	// Edit to change you output items.

        fprintf(fp,"%d %1d %1d %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",k, (*this)(k)._celltype, (*this)(k)._mutant, (*this)(k)._X[0], (*this)(k)._X[1], (*this)(k)._X[2], (*this)(k)._X[3], (*this)(k)._X[4], (*this)(k)._q[0], (*this)(k)._q[1], (*this)(k)._q[2], (*this)(k)._q[3], (*this)(k)._q[4]);
	}
	 
    fclose(fp);
}

void CSystem::OutPutCells(double t)
{
	int k;
	for(k = 1; k<= _NumPoolCell; k++)
	{
        (*this)(k).OutPut(t);
	}
}


CStemCell& CSystem::operator()(int indx)
{
if(indx>0 && indx<=_MaxNumCell)
{
	return *(_cells+(indx-1));
}
else
{
	cout<<"Err<< CSystem () >>Dimensional error"<<endl;
	exit(0);
}
}
