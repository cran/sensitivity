// Author: Alexandre Janon <alexandre.janon at imag.fr>
//
#include <as.h>

extern "C" {
	void Rglue_set_as(int* p, int* N, double *As1, double *As2);
	void Rglue_fill_vf(double *vfsrc, double *errssrc);
	void Rglue_fill_as(int* p, int* N, double *AsXX1src, double *AsXX2src);
	void Rglue_sobol_regr_bc(int* Bptr, double *min, double *max, double *valmin, double *valmax, double *risque);
	void Rglue_next_input(double *curInput);
	void Rglue_reset_input(void);
	void Rglue_get_design(double *design);
	void Rglue_sobol_BFGS(int *Bptr, double *min, double *max, double *valmin, double *valmax,double *risque);
	void Rglue_set_lambda0_r(double *l0, double *r);
	void Rglue_sobol_BFGS_replica(double *valopti, int* jptr, int* signe, int *indEchantillons, double *xsave, double *regular);
	void Rglue_get_partial_design(int* jptr, double* XX);
	void Rglue_regularite(int* p, int* N, double* r, double* x, double* y, double* reg);
	void Rglue_getPN(int* p, int* N);
}

extern double *vf, *errs, *AsXX1, *AsXX2, *distPPV, *invsumkern;
extern int *PPV, *iPPV, szPPV;
extern int AsTaille, vfTaille;
int pGlobal, NGlobal, iGlobal;
ival_t ivGlobal;
extern double lambda0, rPPV;

void Rglue_getPN(int *p, int* N) {
	*p=pGlobal;
	*N=NGlobal;
}

void Rglue_fill_as(int *p, int *N, double *AsXX1src, double *AsXX2src) {
	pGlobal=*p;
	iGlobal=0;
	NGlobal=*N;
	if(AsXX1 && AsXX2) {
		delete[] AsXX1;
		delete[] AsXX2;
		delete[] PPV;
		delete[] iPPV;
		delete[] distPPV;
		delete[] invsumkern;
		PPV=NULL;
		iPPV=NULL;
		distPPV=NULL;
		invsumkern=NULL;
	}
	AsXX1=new double[pGlobal*NGlobal];
	AsXX2=new double[pGlobal*NGlobal];
	memcpy(AsXX1, AsXX1src, sizeof(double)*NGlobal*pGlobal);
	memcpy(AsXX2, AsXX2src, sizeof(double)*NGlobal*pGlobal);
	ivGlobal.min=new double[*p];
	ivGlobal.max=new double[*p];
	for(int i=0;i<*p;i++) {
		ivGlobal.min[i]=0;
		ivGlobal.max[i]=1;
	}
	calcPPV(*N,*p);

}

void Rglue_set_as(int* pptr, int* Nptr, double *X1, double *X2) {
	int p=*pptr;
	int N=*Nptr;
	pGlobal=p;
	iGlobal=0;
	NGlobal=N;
	if(AsXX1 && AsXX2) {
		delete[] AsXX1;
		delete[] AsXX2;
		delete[] PPV;
		delete[] iPPV;
		delete[] distPPV;
		delete[] invsumkern;
		PPV=iPPV=NULL;
		distPPV=invsumkern=NULL;
	}
	AsXX1=new double[p*N];
	AsXX2=new double[p*N];
	memcpy(AsXX1, X1, p*N*sizeof(double));
	memcpy(AsXX2, X2, p*N*sizeof(double));
	calcPPV(N,p);
}

void Rglue_next_input(double *curInput) {
	int j=iGlobal % (pGlobal+1); 
	int n=iGlobal/(pGlobal+1);
	if(j == 0) {
		memcpy(curInput, AsXX1+n*pGlobal, pGlobal*sizeof(double));
	} else {
		memcpy(curInput, AsXX2+n*pGlobal, pGlobal*sizeof(double));
		curInput[j-1]=*(AsXX1+n*pGlobal+(j-1));
	}
	iGlobal++;
}

void Rglue_get_design(double *design) {
	for(int i=0;i<NGlobal*(pGlobal+1);i++) 
		Rglue_next_input(design+i*pGlobal);
}

void Rglue_reset_input() {
	iGlobal=0;
}

void Rglue_fill_vf(double *vfsrc, double *errssrc) {
	int sz=(pGlobal+1)*NGlobal;
	if(vf && errs) {
		delete[] vf;
		delete[] errs;
	}
	vf=new double[sz];
	errs=new double[sz];
	AsTaille=vfTaille=pGlobal*NGlobal;
	memcpy(vf, vfsrc, sz*sizeof(double));
	memcpy(errs, errssrc, sz*sizeof(double));
	Rglue_reset_input();
}

void Rglue_sobol_regr_bc(int* Bptr, double *min, double *max, double *valmin, double *valmax, double *risque) {
	estiSobol2_regr_bc(NULL, pGlobal, NGlobal, *Bptr, &ivGlobal, min, max, valmin, valmax, *risque);
}

void Rglue_sobol_BFGS(int *Bptr, double *min, double *max, double *valmin, double *valmax, double *risque) {
	estiSobol2_OPTbfgs_bc(NULL, pGlobal, NGlobal, *Bptr, &ivGlobal, min, max, valmin, valmax, *risque);
}

void Rglue_get_partial_design(int* jptr, double* XX) {
	int p=pGlobal, N=NGlobal, j=*jptr;
		memcpy(XX, AsXX1, N*p*sizeof(double));
		memcpy(XX+N*p, AsXX2, N*p*sizeof(double));
		for(int i=0;i<N;i++)
			XX[N*p+i*p+j]=AsXX1[i*p+j];
}

void Rglue_set_lambda0_r(double *l0, double *r) {
	lambda0=*l0;
	rPPV=*r;
	debug(lambda0);
	debug(rPPV);
}

void Rglue_sobol_BFGS_replica(double *valopti, int* jptr, int* signe, int *indEchantillons, double *xsave, double *regular) {
	int *IE=(*indEchantillons==-1 ? NULL : indEchantillons);
	double *xsv=(*xsave==-1 ? NULL : xsave);
	estiSobol_bfgs(NGlobal, *jptr, pGlobal, IE, valopti, *signe, false, xsv);
	if(xsv)
		*regular=regularite(NGlobal, *jptr, xsv, NULL);
}

extern double rPPV;

void Rglue_regularite(int* p, int* N, double* r, double* x, double* y, double *reg) {
	double rppvold=rPPV;
	int* PPV=new int[szPPV];
	int* iPPV=new int[*N];
	double *distPPV=new double[szPPV];

	rPPV=*r;
	calcPPV_(*N, x, 0, *p, PPV, iPPV, distPPV);
	*reg=calcIntLap(*N, 0, PPV, iPPV, distPPV, y, NULL);
	delete[] PPV;
	delete[] iPPV;
	delete[] distPPV;

	rPPV=rppvold;
}

