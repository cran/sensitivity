// Author: Alexandre Janon <alexandre.janon at imag.fr>

#include <as.h>

#define TOL 1e-5
bool TOTAL=false;

void mom_new(moments_t *mom, int p) {
	mom->B=1;
	mom->M11=new double[p];
	mom->Mb=new double[p];
	mom->sigmaEst=new double[p];
	mom->bootstrapErrs=false;
}

void mom_init(moments_t *mom, int p) {
	mom->var=0;
	mom->M=0;
	for(int i=0;i<p;i++) {
		mom->M11[i]=mom->Mb[i]=0;
		mom->sigmaEst[i]=-1;
	}
}

void mom_add(double f0, moments_t *mom) {
	mom->M += f0;
	mom->var += f0*f0;
}

void mom_add_M1k(double f0, double fj, int j, moments_t *mom) {
	mom->M11[j] += f0*fj;
	mom->Mb[j] += fj;
}

void mom_final(moments_t *mom, int p, int N) {
	mom->M /= N;
	mom->var = mom->var/N-(mom->M)*(mom->M);

	for(int i=0;i<p;i++) {
		mom->M11[i] /= N;
		mom->Mb[i] /= N;
	}
}

void mom_delete(moments_t *mom) {
	delete[] mom->M11;
	delete[] mom->Mb;
	delete[] mom->sigmaEst;
}

void reinit(time_t& germ) {
	srand(germ);
}

void sampler(double *x, int p, ival_t *iv) {
	for(int i=0 ; i<p ; i++)
		x[i]=rand()/(double)RAND_MAX;

	if(iv) 
		for(int i=0;i<p;i++)
			x[i]=iv->min[i]+(iv->max[i]-iv->min[i])*x[i];
}

void reechelle(double *AS, int N, int p, ival_t *iv) {
	for(int i=0;i<N;i++)
		for(int k=0;k<p;k++) 
			AS[k+i*p]=iv->min[k]+(iv->max[k]-iv->min[k])*AS[k+i*p];
}

double *AsXX1=NULL;
double *AsXX2=NULL;
double *vf=NULL;
double *errs=NULL;
int AsTaille=0;
int vfTaille=0;
int *PPV=NULL;
int *iPPV;
double *distPPV;
double *invsumkern;
////

double calculeSumResidus(int N, int p, int j, double xm, double ym, double a) {
	double sumres=0;
	for(int i=0;i<N;i++) {
		double res=vf[i*(p+1)+j+1]-(ym+a*(vf[i*(p+1)]-xm));
		sumres+=res*res;
	}
	return(sumres);
}

void calculeSumResidusErr(int N, int p, int j, int *indEchantillons, double xm, double ym, double a, double& valinf, double& valsup) {
	valinf=valsup=0;
	for(int r=0;r<N;r++) {
		int i=indEchantillons[r];
		double res1=vf[i*(p+1)+j+1]-errs[i*(p+1)+j+1]-(ym+a*(vf[i*(p+1)]-errs[i*(p+1)]-xm));
		double res2=vf[i*(p+1)+j+1]-errs[i*(p+1)+j+1]-(ym+a*(vf[i*(p+1)]+errs[i*(p+1)]-xm));
		double res3=vf[i*(p+1)+j+1]+errs[i*(p+1)+j+1]-(ym+a*(vf[i*(p+1)]-errs[i*(p+1)]-xm));
		double res4=vf[i*(p+1)+j+1]+errs[i*(p+1)+j+1]-(ym+a*(vf[i*(p+1)]+errs[i*(p+1)]-xm));

		if( (res1>0 && res2>0 && res3>0 && res4>0) || (res1<0 && res2<0 && res3<0 && res4<0)) {
			res1=fabs(res1); res2=fabs(res2); res3=fabs(res3); res4=fabs(res4);
			double inf=(res1<res2 && res1<res3 && res1<res4 ? res1 :
						  (res2<res1 && res2<res3 && res2<res4 ? res2 :
						  (res3<res1 && res3<res2 && res3<res4 ? res3 : res4)));
			double sup=(res1>res2 && res1>res3 && res1>res4 ? res1 :
						  (res2>res1 && res2>res3 && res2>res4 ? res2 :
						  (res3>res1 && res3>res2 && res3>res4 ? res3 : res4)));
			valinf+=inf*inf;
			valsup+=sup*sup;
		} else {
			res1=fabs(res1); res2=fabs(res2); res3=fabs(res3); res4=fabs(res4);
			double sup=(res1>res2 && res1>res3 && res1>res4 ? res1 :
						  (res2>res1 && res2>res3 && res2>res4 ? res2 :
						  (res3>res1 && res3>res2 && res3>res4 ? res3 : res4)));
			valsup+=sup*sup;
		}
	}
}

void calculeABG(int N, int p, int j, int *indEchantillons, double xm, double ym, double& ainf, double& binf, double& ginf, double& asup, double& bsup, double& gsup) {
	double inf0, inf1, infdem;
	double sup0, sup1, supdem;
	calculeSumResidusErr(N,p,j,indEchantillons,xm,ym,0.25,inf0,sup0);
	calculeSumResidusErr(N,p,j,indEchantillons,xm,ym,0.5,infdem,supdem);
	calculeSumResidusErr(N,p,j,indEchantillons,xm,ym,1,inf1,sup1);

	ainf=(16./3.)*inf0 - 8*infdem + (8./3.)*inf1;
	binf=-8*inf0 + 10*infdem -  2*inf1;
	ginf=(8./3.)*inf0-2*infdem+(1./3.)*inf1;
	asup=(16./3.)*sup0 - 8*supdem + (8./3.)*sup1;
	bsup=-8*sup0 + 10*supdem -  2*sup1;
	gsup=(8./3.)*sup0-2*supdem+(1./3.)*sup1;
}

void estiSobol2_regr(double f(double *,int, double *), 
							int p, int N, moments_t *mom, ival_t *iv,
							double *min, double *max, double *minsigma, double *maxsigma) {
	int *indEchantillons=new int[N];
	int B=mom->B;

	memset(min, 0, p*sizeof(double));
	memset(max, 0, p*sizeof(double));
	if(minsigma && maxsigma) {
		memset(minsigma, 0, p*sizeof(double));
		memset(maxsigma, 0, p*sizeof(double));
	}

	for(int k=0;k<B;k++) {
		for(int r=0;r<N;r++) 
			indEchantillons[r]=(k==0 ? r : rand()*(N/(double)RAND_MAX));

		double xm=0;
		double errxm=0;
		for(int r=0;r<N;r++) {
			int i=indEchantillons[r];
			xm += vf[i*(p+1)];
			errxm += errs[i*(p+1)];
		}
		xm /= double(N);
		errxm /= double(N);

		for(int j=0;j<p;j++) {
			if(iv && fabs(iv->max[j]-iv->min[j])<DBL_EPSILON) continue;

			double ym=0;
			double errym=0;
			for(int r=0;r<N;r++) {
				int i=indEchantillons[r];
				ym += vf[i*(p+1)+j+1];
				errym += errs[i*(p+1)+j+1];
			}
			ym /= double(N);
			errym /= double(N);

			double ainf, asup, binf, bsup, ginf, gsup;

			double Smin=-1, Smax=-1;
			int NITMAX=1;
			int NIT=1, NIT2=1;
			for(double xxm=xm-errxm ;  xxm<=xm+errxm-TOL && NIT<NITMAX; xxm+=.01, NIT++) {
				for(double yym=ym-errym ; yym<=ym+errym-TOL && NIT2<NITMAX ; yym+=.01, NIT2++) {
					double delta;
					calculeABG(N,p,j,indEchantillons,xxm, yym, ainf, binf, ginf, asup, bsup, gsup);
					delta=2*sqrt( (asup-ainf)*(gsup-ginf) );
						double minn=-(binf+delta)/(2*ainf);
						double maxn=-(bsup-delta)/(2*asup);
						if(Smin<0 || minn<Smin) Smin=minn;
						if(Smax<0 || maxn>Smax) Smax=maxn;
				}
			}

			min[j] += Smin;
			max[j] += Smax;
			

			if(minsigma && maxsigma) {
				minsigma[j] += Smin*Smin;
				maxsigma[j] += Smax*Smax;
			}
		}
	}

	for(int j=0;j<p;j++) {
		if(iv && fabs(iv->max[j]-iv->min[j])<DBL_EPSILON) continue;
		min[j] /= double(B);
		max[j] /= double(B);
	 { debug(min[j]); debug(max[j]); }
	 	if(minsigma && maxsigma) {
			minsigma[j] = sqrt(minsigma[j]/double(B)-min[j]*min[j]);
			maxsigma[j] = sqrt(maxsigma[j]/double(B)-max[j]*max[j]);
		}
	}

	delete[] indEchantillons;
}

void estiSobol2_regr_replica(double f(double *,int, double *), 
							int p, int N, int* indEchantillons, ival_t *iv,
							double *minreplis, double *maxreplis) {
	double xm=0;
	double errxm=0;
	for(int r=0;r<N;r++) {
		int i=indEchantillons[r];
		xm += vf[i*(p+1)];
		errxm += errs[i*(p+1)];
	}
	xm /= double(N);
	errxm /= double(N);
	for(int j=0;j<p;j++) {
		double ym=0;
		double errym=0;
		for(int r=0;r<N;r++) {
			int i=indEchantillons[r];
			ym += vf[i*(p+1)+j+1];
			errym += errs[i*(p+1)+j+1];
		}
		ym /= double(N);
		errym /= double(N);

		double ainf, asup, binf, bsup, ginf, gsup;
		double Smin, Smax;

		calculeABG(N,p,j,indEchantillons,xm, ym, ainf, binf, ginf, asup, bsup, gsup);
		double delta=2*sqrt( (asup-ainf)*(gsup-ginf) );
		Smin=-(binf+delta)/(2*ainf);
		Smax=-(bsup-delta)/(2*asup);

		int NIT=2;
		for(double xxm=xm-errxm ; xxm<=xm+errxm-TOL ; xxm+=2*errxm/NIT) {
			for(double yym=ym-errym ; yym<=ym+errym-TOL ; yym+=2*errxm/NIT) {
				calculeABG(N,p,j,indEchantillons,xxm, yym, ainf, binf, ginf, asup, bsup, gsup);
				double delta=2*sqrt( (asup-ainf)*(gsup-ginf) );
				double minn=-(binf+delta)/(2*ainf);
				double maxn=-(bsup-delta)/(2*asup);
				if(Smin<0 || minn<Smin) Smin=minn;
				if(Smax<0 || maxn>Smax) Smax=maxn;
			}
		}

		minreplis[j]=Smin;
		maxreplis[j]=Smax;
	}
}

int compare(const void *a, const void *b) {
	return(*((double*)a)>*((double*)b));
}

double Phi(double x) {
	return(.5*erfc(-x/sqrt(2.)));
}

double Probit(double p) {
	if(fabs(p)<DBL_EPSILON)
		return(-10000);

	double x=0, xold=0;
	do {
		double e=Phi(x)-p;
		double u=e*sqrt(8*atan(1.))*exp(x*x/2.);
		xold=x;
		x=x-u/(1+x*u/2.);
	} while(fabs(x-xold)>1e-6);
	return(x);
}

void bc(double *replis, int B, double *min, double *max, double risque=.05) {
	double alpha=risque/2;
	double thetachap=replis[0];
	qsort(replis, B, sizeof(double), compare);
	int cdf=0;
	for(; cdf<B && replis[cdf]<=thetachap ; cdf++);
	double z0=Probit(double(cdf)/double(B));
	double zalpha=Probit(alpha);
	double zcalpha=Probit(1-alpha);
	double q1=Phi(2*z0+zalpha);
	double q2=Phi(2*z0+zcalpha);
	int A=int(double(B)*q1);
	int AA=int(double(B)*q2);
	A=(A<0 ? 0 : A);
	A=(A>=B ? B-1 : A);
	AA=(AA<0 ? 0 : AA);
	AA=(AA>=B ? B-1 : AA);
	*min=replis[A];
	*max=replis[AA];
}

void estiSobol2_regr_bc(double f(double*,int,double*), int p, int N, int B, ival_t *iv, double *min, double *max, double *valmin, double *valmax, double risque) {
	double *lminreplis=new double[B*p];
	double *lmaxreplis=new double[B*p];

	for(int k=0;k<B;k++) {
		unsigned int seedp=k;
		double *minreplis=new double[p];
		double *maxreplis=new double[p];
		int *indEchantillons=new int[N];
		for(int r=0;r<N;r++) 
			indEchantillons[r]=(k==0 ? r : rand()*(N/(double)RAND_MAX));

			estiSobol2_regr_replica(f,p,N,indEchantillons,iv,minreplis,maxreplis);
			for(int j=0;j<p;j++) {
				lminreplis[j*B+k]=minreplis[j];
				lmaxreplis[j*B+k]=maxreplis[j];
				if(k==0) {
					valmin[j]=minreplis[j];
					valmax[j]=maxreplis[j];
				}
			}
			delete[] minreplis;
			delete[] maxreplis;
			delete[] indEchantillons;
		}
	

	for(int j=0;j<p;j++) {
		double minn, maxn;
		bc(lminreplis+j*B, B, &minn, &maxn,risque);
		min[j]=minn;
		bc(lmaxreplis+j*B, B, &minn, &maxn,risque);
		max[j]=maxn;
	}

	delete[] lmaxreplis;
	delete[] lminreplis;
}

void estiSobol2_replica(double f(double *,int,double *), int p, int N, int* indEchantillons, ival_t *iv, double *replis) {
	moments_t mom;
	mom_new(&mom,p);
	mom_init(&mom, p);
	for(int r=0;r<N;r++) {
		int i=indEchantillons[r]; 
		mom_add(vf[i*(p+1)], &mom);
		for(int j=0 ; j<p ; j++) 
				mom_add_M1k(vf[i*(p+1)],vf[i*(p+1)+j+1], j, &mom);
	}
	mom_final(&mom, p, N);
	for(int j=0 ; j<p ; j++) {
			replis[j]=(mom.M11[j]-mom.M*mom.Mb[j])/mom.var;
	}
	mom_delete(&mom);
}

double *lreplis=NULL;

////

double lambda0=0;
double lambda=lambda0;
double rPPV=.1;
int szPPV=100*1048576;

void calcPPV_(int N, double *X, int plan, int p, int *PPV, int *iPPV, double *distPPV) {
	if(lambda0<=0) return;
	double dist;
	int nvois=0;
	for(int i=0;i<N;i++) {
		iPPV[plan*N+i]=(plan*N+i==0 ? 0 : iPPV[plan*N+i-1]);
		for(int j=0;j<N;j++) {
			if(j==i) continue;
			dist=0;
			for(int d=0;d<p;d++) {
				dist += (X[i*p+d]-X[j*p+d])*(X[i*p+d]-X[j*p+d]);
			}
			dist /= rPPV;
			if(dist<1) {
				double kern=(1-dist)*(1-dist);
				//cout << "dans le plan " << plan << ", " << j << " est voisin de " << i << endl;
				PPV[iPPV[plan*N+i]]=j; 
				distPPV[iPPV[plan*N+i]]=kern; 
				//debug(distPPV[iPPV[plan*N+i]]);
				iPPV[plan*N+i]++;  
				nvois++;
				if(iPPV[plan*N+i]>szPPV) { debugmsg("increase szPPV or decrease winsmooth"); debug(szPPV);}
			}
		}
	}
	//debug(double(nvois)/double(N));
	debug(double(nvois)/double(N*N));
}

double calcInvsumkern2(int N, int j, int *PPV, int *iPPV, double *distPPV, int i) {
	double sumkern=1;
	for(int k=(j*N+i==0 ? 0 : iPPV[j*N+i-1]);k<iPPV[j*N+i];k++) {
	//	cout << "dans le plan " << j << ", " << PPV[k] << " est voisin de " << i << "; poids=" << distPPV[k] << endl;
		sumkern += distPPV[k];
	}
	return(1./sumkern);
}

void calcInvsumkern(int N, int p, int *PPV, int *iPPV, double *distPPV, double *invsumkern) {
        for(int j=0;j<p;j++) {
                for(int i=0;i<N;i++) {
                        invsumkern[j*N+i]=1;
                        for(int k=(j*N+i==0 ? 0 : iPPV[j*N+i-1]);k<iPPV[j*N+i];k++) {
                        //	cout << "dans le plan " << j << ", " << PPV[k] << " est voisin de " << i << "; poids=" << distPPV[k] << endl;
                                invsumkern[j*N+i] += distPPV[k];
                        }
                        invsumkern[j*N+i]=1./invsumkern[j*N+i];
                }
        }
}


void calcPPV(int N, int p) {
	if(PPV || lambda0<=0) return;
	double *XX=new double[N*2*p]; // 2N points à p coords
	iPPV=new int[N*2*p]; // p tableaux iPPV de taille 2N
	szPPV=N*p*1000;
	PPV=new int[szPPV];
	distPPV=new double[szPPV];
   invsumkern=new double[N*2*p];
	for(int j=0;j<p;j++) {
		debug(j);
		if(j==0)
			memcpy(XX, AsXX1, N*p*sizeof(double));
		memcpy(XX+N*p, AsXX2, N*p*sizeof(double));
		for(int i=0;i<N;i++)
			XX[N*p+i*p+j]=AsXX1[i*p+j];
		calcPPV_(2*N, XX, j, p, PPV, iPPV, distPPV);
	}
   calcInvsumkern(2*N, p, PPV, iPPV, distPPV, invsumkern);
	delete[] XX;
}

void calcPPVindEch(int N, int p, int* indEchantillons) {
	double *XX=new double[N*2*p]; // 2N points à p coords

	//debugmsg("calcPPVindEch...");
	for(int j=0;j<p;j++) {
		if(j==0) { // copie du 1er plan extrait suivant indEchantillons dans XX
			for(int ii=0;ii<N;ii++) {
				int i=indEchantillons[ii];
				memcpy(XX+ii*p, AsXX1+i*p, p*sizeof(double));
			}
		}
		// copie du 2eme plan j extrait dans XX+N*p
		for(int ii=0;ii<N;ii++) {
			int i=indEchantillons[ii];
			memcpy(XX+N*p+ii*p, AsXX2+i*p, p*sizeof(double));
			// échange jème variable
			XX[N*p+ii*p+j]=AsXX1[i*p+j];
		}
		calcPPV_(2*N, XX, j, p, PPV, iPPV, distPPV);
	}
   calcInvsumkern(2*N, p, PPV, iPPV, distPPV, invsumkern);
	delete[] XX;
}

double calcIntLap(int N, int j, int *PPV, int *iPPV, double *distPPV, double *Y, double *valSmooth) {
	double r=0;
	double sumkern, smoothed;
	for(int i=0;i<N;i++) {
		smoothed=Y[i];
		sumkern=1;
		for(int k=(j*N+i==0 ? 0 : iPPV[j*N+i-1]);k<iPPV[j*N+i];k++) {
	//		cout << "dans le plan " << j << ", " << PPV[k] << " est voisin de " << i << "; poids=" << distPPV[k] << endl;
			smoothed += Y[PPV[k]]*distPPV[k];
			sumkern += distPPV[k];
		}
		smoothed/=sumkern;
		r+=(smoothed-Y[i])*(smoothed-Y[i]);
		if(valSmooth)
			valSmooth[i]=smoothed;
	}
	return(r/double(N));
}

double gradLap(int N, int j, int *PPV, int *iPPV, double *distPPV, int l, double *Y, double *valSmooth, double *invsumkern) {
	double grad=2*(Y[l]-valSmooth[l])*(1-invsumkern[j*N+l]);
	for(int k=(j*N+l==0 ? 0 : iPPV[j*N+l-1]);k<iPPV[j*N+l];k++) {
	//	cout << "dans le plan " << j << ", " << PPV[k] << " est voisin de " << l << "; poids=" << distPPV[k] << endl;
		grad -= 2*(Y[PPV[k]]-valSmooth[PPV[k]])*distPPV[k]*invsumkern[j*N+PPV[k]];
	}
	return(grad/N);
}

double regularite(int N, int j, double *x, double *valSmooth) {
	return(calcIntLap(2*N, j, PPV, iPPV, distPPV, x, valSmooth));
}

static double mu=1;

double calcSobol_bfgs(int N, double *x, int j, int p, double& M12, double& MC1, double& M1, double& M2, double *valSmooth) {
	M12=MC1=M1=M2=0;
	for(int i=0;i<N;i++) {
		M12 += x[i]*x[N+i];
		MC1 += x[i]*x[i];
		M1 += x[i];
		M2 += x[N+i];
	}
	M12 /= N; MC1 /= N;
	M1 /= N; M2 /= N;
	double penal=0;
	if(lambda0>0) penal=lambda*regularite(N, j, x, valSmooth);
//	debug(regularite(N,j,x,valSmooth));
	return((M12-M1*M2)/(MC1-M1*M1)+penal);
}

double gradSobol(int N, double *x, int j, int p, double& M12, double& MC1, double& M1, double& M2, int l, double *valSmooth) {
	double grad;
	if(l>=N) {
		grad=(x[l-N]/N-M1/N)/(MC1-M1*M1);
	} else {
		grad=((x[l+N]/N-M2/N)*(MC1-M1*M1)-2*(M12-M1*M2)*(x[l]/N-M1/N))/((MC1-M1*M1)*(MC1-M1*M1));
	}
	double gradregul=0;
	if(lambda0>0) 
		gradregul=gradLap(2*N,j,PPV,iPPV,distPPV,l,x,valSmooth,invsumkern);
	return(grad+lambda*gradregul);
}

extern "C" {
#include <R.h>
	void F77_NAME(setulb)(int* n, int *m, double *x, double *l, double *u, int *nbd, double *f, double *g, double *factr, double *pgtol, double *wa, int *iwa, char *task, int *iprint, char *csave, int *lsave, int *isave, double *dsave);
}

void estiSobol_bfgs(int N, int j, int p, int *indEchantillons, double *valopti, double signe, bool dump, double *xsave) {
	int dim=2*N;
	int *nbd=new int[dim];
	double *valSmooth=new double[dim];
	double M12, MC1, M1, M2;
	double *l=new double[dim];
	double *u=new double[dim];
	double *x;
	if(xsave)
		x=xsave;
	else
		x=new double[dim];
	double *grad=new double[dim];
	if(indEchantillons) {
		for(int i=0;i<N;i++) {
			l[i]=vf[indEchantillons[i]*(p+1)]-errs[indEchantillons[i]*(p+1)];
			u[i]=vf[indEchantillons[i]*(p+1)]+errs[indEchantillons[i]*(p+1)];
			nbd[i]=2;
		}
		for(int i=N;i<2*N;i++) {
			l[i]=vf[indEchantillons[i-N]*(p+1)+j+1]-errs[indEchantillons[i-N]*(p+1)+j+1];
			u[i]=vf[indEchantillons[i-N]*(p+1)+j+1]+errs[indEchantillons[i-N]*(p+1)+j+1];
			nbd[i]=2;
		}
	} else {
		for(int i=0;i<N;i++) {
			l[i]=vf[i*(p+1)]-errs[i*(p+1)];
			u[i]=vf[i*(p+1)]+errs[i*(p+1)];
			nbd[i]=2;
		}
		for(int i=N;i<2*N;i++) {
			l[i]=vf[(i-N)*(p+1)+j+1]-errs[(i-N)*(p+1)+j+1];
			u[i]=vf[(i-N)*(p+1)+j+1]+errs[(i-N)*(p+1)+j+1];
			nbd[i]=2;
		}
	}
	double factr=1e-3;
	double pgtol=1e-3;
	int m=20;
	int nmax=dim, mmax=m;
	double *wa=new double[(2*mmax+4)*nmax+12*mmax*mmax+12*mmax];
	int *iwa=new int[3*dim];
	char task[60];
	int iprint=-1;
	char csave[60];
	int lsave[4];
	double dsave[29];
	int isave[44];

	int Nrestart=1;

	for(int K=0;K<Nrestart;K++) {
		memset(task,' ',60*sizeof(char));
		task[0]='S';task[1]='T';task[2]='A';task[3]='R';task[4]='T';
	//	memset(csave,'X',60*sizeof(char));

		lambda=lambda0*signe;

		for(int jj=0;jj<dim;jj++) 
			x[jj]=l[jj]+(u[jj]-l[jj])*double(rand())/double(RAND_MAX);

		do {
		//	printf("%s\n", csave);
			F77_NAME(setulb)(&dim,&m,x,l,u,nbd,valopti,grad,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave);

			if(task[0]=='F' && task[1]=='G') {
				*valopti=signe*calcSobol_bfgs(N,x,j,p,M12,MC1,M1,M2,valSmooth);
				for(int jj=0;jj<dim;jj++) {
					grad[jj]=signe*gradSobol(N,x,j,p,M12,MC1,M1,M2,jj,valSmooth);
				}
			} else if(task[0]=='N' && task[1]=='E' && task[2]=='W') {
				continue;
			} else {
				//debug(task);
				break;
			}
			} while(1);

		lambda=0;
		*valopti=calcSobol_bfgs(N,x,j,p,M12,MC1,M1,M2,valSmooth);
//		debug(regularite(N,j,x,valSmooth));

		//debug(*valopti);
	}

	if(dump || (getenv("BFGSDUMP")&&j==atoi(getenv("BFGSDUMP"))) ) {
		ofstream d("/tmp/dump");
		ofstream dvf("/tmp/vfdump");
		for(int i=0;i<N;i++) {
			for(int jj=0;jj<p;jj++) {
				d << AsXX1[i*p+jj] << " ";
				dvf <<  AsXX1[i*p+jj] << " ";
			}
			d << x[i] << endl;
			dvf << vf[i*(p+1)] << endl;
		}
		for(int i=0;i<N;i++) {
			for(int jj=0;jj<p;jj++) {
				if(j==jj) {
					d << AsXX1[i*p+jj] << " ";
					dvf << AsXX1[i*p+jj] << " ";
				} else {
					d << AsXX2[i*p+jj] << " ";
					dvf << AsXX2[i*p+jj] << " ";
				}
			}
			d << x[i+N] << endl;
			dvf << vf[i*(p+1)+(j+1)] << endl;
		}
	}

	delete[] wa;
	delete[] iwa;
	delete[] u;
	delete[] l;
	delete[] nbd;
	if(!xsave)
		delete[] x;
	delete[] valSmooth;
	delete[] grad;
}

int varInt=-1;

void estiSobol2_OPTbfgs_bc(double f(double*,int,double*), int p, int N, int B, ival_t *iv, double *min, double *max, double *valmin, double *valmax, double risque) {

	// replications bootstrap
	int *indEchantillon=new int[N];
	double *repliMin=new double[p];
	double *repliMax=new double[p];
	double *lminreplis;
	double *lmaxreplis;
	ofstream ff;
	if(getenv("DUMPBOOT"))
		ff.open("/tmp/dumpboot");
	moments_t mom;
	mom.sigmaEst=NULL;
	lminreplis=new double[B*p];
	lmaxreplis=new double[B*p];
	for(int k=0 ; k<B;k++) {
	//	if(!(k%10)) 
	//		printf("\nComputing replication %d out of %d... ", k+1, B);
		for(int l=0;l<N;l++) indEchantillon[l]=(k==0 ? l : rand()*(N/(double)RAND_MAX)); 
		 debug(k);
		if(k!=0)
			calcPPVindEch(N,p,indEchantillon);

		for(int j=0;j<p;j++) {
			if(varInt<0 || (varInt>=0 && j==varInt)) {
				estiSobol_bfgs(N, j, p, indEchantillon, &repliMin[j]);
				estiSobol_bfgs(N, j, p, indEchantillon, &repliMax[j], -1);
			}
		}

		if(k==0) {
			for(int j=0;j<p;j++) {
				valmin[j] = repliMin[j];
				valmax[j] = repliMax[j];
				debug(valmin[j]);
				debug(valmax[j]);
			}
		}
		for(int j=0;j<p;j++) {
			lminreplis[k+j*B]=repliMin[j];
			lmaxreplis[k+j*B]=repliMax[j];
			if(getenv("DUMPBOOT") && atoi(getenv("DUMPBOOT"))==j)
				ff << repliMax[j] << endl;
		}
	}
	//printf("\n");

	if(mom.sigmaEst) {
		for(int j=0;j<p;j++) {
			if(iv && iv->max[j]-iv->min[j]<DBL_EPSILON) continue;
			double minn, maxn;
			min[j]=valmin[j]-2*mom.sigmaEst[j];
			max[j]=valmax[j]+2*mom.sigmaEst[j];
		}
		mom_delete(&mom);
	} else {
		for(int j=0;j<p;j++) {
			double minn, maxn;
			bc(lminreplis+j*B, B, &minn, &maxn, risque);
			min[j]=minn;
			bc(lmaxreplis+j*B, B, &minn, &maxn, risque);
			max[j]=maxn;
		}
	}

	delete[] indEchantillon;
	delete[] repliMin;
	delete[] repliMax;
	delete[] lminreplis;
	delete[] lmaxreplis;
}


