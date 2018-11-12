#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "string.h"
#include <iostream>
using namespace std;

#include "Random.h"
#include "BCTool.h"
#include "CStemCell.h"

extern struct IMD _MD;
extern CRandom Rand;
extern double _ratio19;


CStemCell::CStemCell()
{
}

CStemCell::~CStemCell()
{
}

bool CStemCell::Initialized(int k, char fpar[])
{
    SetDefaultPar();
	ReadPar(fpar);
	_t = 0.0;
	_cellid=k;
	_NumVar=5;
	_NumRand=0;
	_X = new double[_NumVar];
    _X[0] = CD34;
    _X[1] = CD19;
    _X[2] = CD22;
    _X[3] = LNK;
    _q[0] = CD34;
    _q[1] = CD19;
    _q[2] = CD22;
    _q[3] = LNK;
    
    _celltype = GetCellType(_X[0]);
    _mutant = 0;
    return true;
}

int CStemCell::GetCellType(double x) // x: CD34
{
    int celltype;
    
    if (x > 0.6){
        celltype = HSC;        // HSC
    } else {
        celltype = PREC;      // Precursor
    }
    return(celltype);
}

void CStemCell::ReadPar(char fpar[])
{
	FILE *fp;
	char str[StrLength], *pst;
	if((fp = fopen(fpar,"r"))==NULL)
	{
		cout<<"Cannot open the cell parameter input file."<<endl;
		exit(0);
	}
	rewind(fp);
	while(!feof(fp))
	{
		fgets(str,StrLength,fp);
		if(str[0]=='#'){ continue;}

        if((pst=strstr(str,"beta0="))!=NULL)
	{
		_par.beta0=atof(pst+6);
	}
        if((pst=strstr(str,"cd34b="))!=NULL)
        {
            _par.cd34b=atof(pst+6);
        }
        if((pst=strstr(str,"alpha19="))!=NULL)
        {
            _par.alpha19=atof(pst+8);
        }
        if((pst=strstr(str,"kappa0="))!=NULL)
        {
            _par.kappa0=atof(pst+7);
        }
        if((pst=strstr(str,"mu0="))!=NULL)
        {
            _par.mu0=atof(pst+4);
        }
        if((pst=strstr(str,"mu1="))!=NULL)
        {
            _par.mu1=atof(pst+4);
        }
        if((pst=strstr(str,"mu2="))!=NULL)
        {
            _par.mu2=atof(pst+4);
        }
        if((pst=strstr(str,"m="))!=NULL)
        {
            _par.m=atof(pst+2);
        }
        if((pst=strstr(str,"p0="))!=NULL)
        {
            _par.p0=atof(pst+3);
        }
        if((pst=strstr(str,"plost="))!=NULL)
        {
            _par.plost=atof(pst+6);
        }
        if((pst=strstr(str,"X0="))!=NULL)
        {
            _par.X0=atof(pst+3);
        }
        if((pst=strstr(str,"n0="))!=NULL)
        {
            _par.n0=atof(pst+3);
        }
        if((pst=strstr(str,"X1="))!=NULL)
        {
            _par.X1=atof(pst+3);
        }
        if((pst=strstr(str,"n1="))!=NULL)
        {
            _par.n1=atof(pst+3);
        }
        if((pst=strstr(str,"theta="))!=NULL)
        {
            _par.theta=atof(pst+6);
        }
        if((pst=strstr(str,"gamma19v="))!=NULL)
        {
            _par.gamma19=atof(pst+9);
        }
        if((pst=strstr(str,"gamma22v="))!=NULL)
        {
            _par.gamma22=atof(pst+9);
        }
        if((pst=strstr(str,"gamma123v="))!=NULL)
        {
            _par.gamma123=atof(pst+10);
        }
        if((pst=strstr(str,"CD34="))!=NULL)
        {
            CD34=atof(pst+5);
        }
        if((pst=strstr(str,"CD19="))!=NULL)
        {
            CD19=atof(pst+5);
        }
        if((pst=strstr(str,"CD22="))!=NULL)
        {
            CD22=atof(pst+5);
        }
        if((pst=strstr(str,"LNK="))!=NULL)
        {
            LNK=atof(pst+4);
        }
        if((pst=strstr(str,"a123="))!=NULL)
        {
            _par.a123=atof(pst+5);
        }
        if((pst=strstr(str,"a34="))!=NULL)
        {
            _par.a34=atof(pst+4);
        }
        if((pst=strstr(str,"a22="))!=NULL)
        {
            _par.a22=atof(pst+4);
        }
	}
	fclose(fp);
    
    _plost = _par.plost;

//    OutPutParameters();
    
}

void CStemCell::OutPutParameters()
{
    printf("beta0=%f\n", _par.beta0);
    printf("mu0=%f, mu1=%f, m2=%f\n",_par.mu0, _par.mu1, _par.mu2);
    printf("alpha19=%f\n",_par.alpha19);
    printf("kappa0=%f\n",_par.kappa0);
    printf("p0=%f\n",_par.p0);
    printf("m=%f\n",_par.m);
    printf("gamma19=%f, gamma22=%f, gamma123=%f\n", _par.gamma19, _par.gamma22,_par.gamma123);
    printf("a34=%f, a123=%f, a22=%f\n", _par.a34, _par.a123, _par.a22);
    printf("X0=%f, X1=%f\n", _par.X0, _par.X1);
    printf("ratio19=%f\n", _ratio19);
    printf("alpha34=%f\n", _par.alpha34);
}
void CStemCell::SetDefaultPar()
{
    CD34 = Rand(0.01,0.1);
    CD19 = Rand(0.6,1.0);
    CD22 = Rand(0.01,0.1);
    LNK = Rand(0.6,1.0);
    CD123 = Rand(0.01,0.2);
    _par.plost=0;
    _par.p0=0;
    _par.a34=0;
    _par.a123=0;
    _par.a22=0;
    _par.mu2=0;
    _par.gamma123=0;
}

void CStemCell::GetRandForces(double b[][NUMRAND], double x[], int k)
{
}

void CStemCell::GetTrends(double a[], double x[])
{
}

DCell CStemCell::CellFateDecision(long N0, double R19, double R123, double dt)
{
    int i;
    double mu;
    double beta;
    double kappa;
    double p0;
    double _rand;
    
    _rand = Rand();
    _CARSignal = (1.0/((1 + pow(CD34/_par.X0, _par.n0))*(1 + pow(CD123/_par.X1, _par.n1)))) * R19 * (_par.gamma19 * CD19)/(1 + _par.gamma19 * CD19 + _par.gamma22 * CD22);
    mu = fdeathrate(R19, R123) * dt;
    beta = fbeta(N0,R19, R123) * dt;
    kappa = fkappa() * dt;
    p0 = _par.p0 * dt;
    
    _par.alpha22 = 1.54 + 0.25 * CD19 + _par.a22 * _CARSignal;
    _par.alpha123 = 1.54 + 0.45 * CD34 + _par.a123 * _CARSignal;
    _par.alphaLNK = 2.1 - 0.3 * _CARSignal;  
    _par.alpha34 = _par.cd34b + 0.16 * CD19 + _par.a34 * _CARSignal;
    
    DCell nextcell;
    
    nextcell.type = 0;  //Remain unchanged.
    for (i=0;i<_NumVar;i++)
    {
        nextcell.X1[i] = _X[i];
        nextcell.q1[i] = _q[i];
    }
    nextcell.celltype1 = _celltype;
    nextcell.mutanttype1 = _mutant;

    if (_rand < mu)  // Cell death
    {
        nextcell.type = 2;
    }
    else{
        if (_rand < mu+beta)   // Cell division
        {
            nextcell.type = 1;  // Cell division
            for (i=0;i<_NumVar;i++)
            {
                nextcell.q1[i] = GetnextEpi(i, _q[i], R19);
                nextcell.q2[i] = GetnextEpi(i, _q[i], R19);
            
                nextcell.X1[i] = 1.0*nextcell.q1[i];
                nextcell.X2[i] = 1.0*nextcell.q2[i];
            }
            nextcell.celltype1 = GetCellType(nextcell.X1[0]);
            nextcell.celltype2 = GetCellType(nextcell.X2[0]);
            nextcell.mutanttype1 = _mutant;
            nextcell.mutanttype2 = _mutant;
        }
        else{
            if (_celltype == 0 && _rand < mu+beta+p0)    // CD19 mutation
            {
                nextcell.type = 3;  //CD19 mutation
                for (i=0;i<_NumVar;i++)
                {
                    nextcell.X1[i] = _X[i];
                    nextcell.q1[i] = _q[i];
                }
                nextcell.mutanttype1 = 1;
                nextcell.celltype1 = _celltype;
            }
            else if (_celltype == PREC && _rand < kappa+mu+beta+p0)    // Terminal differentiation
            {
                nextcell.type = 4;      // Differentiation
                for (i=0; i<_NumVar; i++)
                {
                    nextcell.X1[i] = _X[i];
                    nextcell.q1[i] = _q[i];
                }
                nextcell.mutanttype1 = _mutant;
                nextcell.celltype1 = TD;
            }
        }
    }
    return(nextcell);
}

double CStemCell::GetnextEpi(int i, double u, double R)
{
    double y=0;
    double a, b;
    double z;
    switch (i) {
        case 0: // For CD34
            y = 0.08 + 1.06 * pow(_par.alpha34 * u, 2.2)/(1+pow(_par.alpha34 * u, 2.2));
            break;
        case 1: // For CD19
            y = 0.05 + 1.06 * pow(_par.alpha19 * u, 1.6)/(1+pow(_par.alpha19 * u, 1.6));
            break;
        case 2: // For CD22
            y = 0.04 + 0.96 * pow(_par.alpha22 * u, 1.6)/(1+pow(_par.alpha22 * u, 1.6));
            break;
        case 3: // For LNK
            y = 0.1 + 0.85 * pow(_par.alphaLNK * u, 2.0)/(1+pow(_par.alphaLNK * u, 2.0));
            break;
        case 4: // For CD123
            y = 0.05 + 1.02 * pow(_par.alpha123 * u, 1.8)/(1+pow(_par.alpha123 * u, 1.8));
            break;
    }
    a = (_par.m  - 1)*y;
    b = (_par.m - 1)*(1-y);
    z = Rand.BetaDistribution(a, b);
    return(z);
}

void CStemCell::Update(int mu)
{
}

void CStemCell::Propensities(double a[])
{
}

double CStemCell::fbeta(long N0, double R19, double R123)
{
    double beta;
    
    if(_celltype == TD){
        beta = 0;
    }
    else {
        beta = _par.beta0*(1.0/(1.0+(1.0*N0/_par.theta)))*((5.8*CD34+pow(2.2*CD34,6.0))/(1+pow(3.75*CD34,6.0)));   // Proliferation rate of cells
    }
    return(beta);
}

double CStemCell::fkappa()
{
    double kappa;
    kappa = _par.kappa0*(1.0/(1.0 + pow(4.0*CD34, 6.0)));
    return(kappa);
}

double CStemCell::fdeathrate(double R19, double R123)
{
    double mu;
    if(_mutant == 1 )    // CD19 mutant cells or HSC
    {
        mu = _par.mu0;
    } else {
        mu = _par.mu0 + _par.mu1* _CARSignal + _par.mu2 * (_par.gamma123 * CD123)/(1.0 + _par.gamma123 * CD123) * R123;
    }
    return(mu);
}

void CStemCell::Pre_OneStep()
{
    CD34=_X[0];
    CD19=_X[1];
    CD22=_X[2];
    LNK=_X[3];
    CD123=_X[4];
}



