#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include <cfloat>
#include <cstring>

#define debug(s)
#define debugmsg(s)

using namespace std;

/* as.cpp */
typedef struct {
	double var; // Var(f(X))
	double M;   // E(f(X))
	double *M11; // E(f(X1)f(X2))
	double *Mb; // E(f(X2))
	int B; // nombre de replic bootstrap
	bool bootstrapErrs; // bootstrap sur errs ?
	double *sigmaEst; // estimation de l'écart type de l'estimateur
} moments_t;

typedef struct {
	double *min;
	double *max;
} ival_t;

extern int AsTaille;
extern int vfTaille;
extern bool ErreurSortie;
extern bool TOTAL;

double Phi(double x);
double Probit(double p);
void GenEchantillon(double f(double*,int,double*), int p, int N, ival_t *iv, bool genErrs);
void mom_new(moments_t *mom, int p);
void mom_init(moments_t *mom, int p);
void mom_delete(moments_t *mom);
void sampler(double *x, int p, ival_t *iv);
void estiSobol2_regr(double f(double *,int, double *), 
							int p, int N, moments_t *mom, ival_t *iv,
							double *min, double *max, double *minsigma, double *maxsigma);
void estiSobol2_regr_replica(double f(double *,int, double *), 
							int p, int N, int* indEchantillons, ival_t *iv,
							double *minreplis, double *maxreplis);
void estiSobol2_regr_bc(double f(double*,int,double*), int p, int N, int B, ival_t *iv, double *min, double *max, double *valmin, double *valmax, double risque=0.05) ;
void estiSobol2_boot_bc(double f(double *,int,double *), int p, int N, int B, ival_t *iv,double *s, double *inf, double *sup);
void estiSobol2_boot(double f(double *,int,double *), int p, int N, double *s, moments_t *mom, ival_t *iv, double *min=NULL, double *max=NULL, double *minsigma=NULL, double *maxsigma=NULL);
void estiSobol2_asympt(double f(double *,int,double *), int p, int N, double *s, moments_t *mom, ival_t *iv);
void estiSobolT_asympt(double f(double *,int,double *), int p, int N, double *s, moments_t *mom, ival_t *iv);
void benchMultifidSobol(double f(double *,int, double *), double f2(double *,int,double *), int p, int N, double *s, moments_t *mom, ival_t *iv);
void multiFid_optim(double nf, double nc, double expoCout, double sigmac, double sigmaf, double P, double alpha);
void estiSobol2erreurs(int p, moments_t *mom, int N, double *errs, double *vf, int *indEchantillons, double *inf, double *sup);
void estiSobolRL(double f(double*,int,double*), int p, int N, double *s, moments_t *mom, ival_t *iv);
void estiSobol_bootErr(double f(double*,int,double*), int p, int N, double *s, moments_t *mom, ival_t *iv);
void estiSobol2_OPTbfgs_bc(double f(double*,int,double*), int p, int N, int B, ival_t *iv, double *min, double *max, double *valmin, double *valmax, double risque=.05);
void estiSobol_bfgs(int N, int j, int p, int *indEchantillons, double *valopti, double signe=1, bool dump=false, double *xsave=NULL) ;
void calcPPV_(int N, double *X, int j, int p,int *PPV, int *iPPV, double *distPPV);
void calcPPV(int N, int p);
double calcIntLap(int N,int j, int *PPV, int *iPPV, double *distPPV, double *Y, double *valSmooth=NULL);
double regularite(int N, int j, double *x, double *valSmooth);

/* benchas.cpp */
void benchTousAs(int idVariable, int p, int N, int B, ival_t *iv, double VV, double f(double *,int, double *), ostream& o);
double benchTousAs2(int idVariable, int p, int N, int B, int NReplic, ival_t *iv, double VV, double f(double *,int, double *), int* estiTest, ostream& o) ;


