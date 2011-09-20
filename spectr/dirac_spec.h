#include<QCDenvNODEpt.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

using namespace std;


inline int eigv_cmp(const void*, const void*);


class eigsysws {

public:
    latt *z;
    int n;
    gsl_matrix_complex *matrix;
    gsl_matrix_complex *eigvec;
    gsl_vector *eigval;
    gsl_eigen_hermv_workspace *eigws;

    eigsysws(latt*, int);
    
    ~eigsysws();

    void braket_set(int, int, Cplx);

    inline void sort(){ gsl_eigen_genhermv_sort(eigval, eigvec, GSL_EIGEN_SORT_VAL_ASC); }

    Cplx braket_get(int, int);

    inline double eigval_get(int i) { return (gsl_vector_get(eigval, i)); }

    Cplx eigvec_get(int, int);
    
    void diagonalize();

    void ket(SpinColor_fld*, dgsp&, eigv*, int);

    void tosmart(SpinColor_fld*, SpinColor_fld*, dgsp&, eigv*);

    void tostupid(SpinColor_fld*, SpinColor_fld*, dgsp&, eigv*);

    void prout( gsl_matrix_complex*);

    void rebase(eigsysws&, dgsp&, int);
    
};


class diracws {
public:
    latt *Z;
    SpinColor_fld *q, *tq, *sq, *p, *sp;
    SpinColor_fld *inres, *tmpres;
    fftw_plan plan[3];

    diracws(latt*);

    ~diracws();

/********************** DIRAC OPERATOR IN POSITION SPACE **********************/
/* WARNING: you can safely use this method _only_ with _non_ overlapping (in  */
/* particular: non coinciding) memory areas for input and output spincolor    */
/* fields, since for each spincolor it reads the two neighborhood ones in the */
/* four directions.                                                           */
/******************************************************************************/
    SpinColor_fld *dirac_q(SpinColor_fld*, SpinColor_fld*, ptGluon_fld&, int, int, int);
    
/***************** Dirac Operator in Momentum Space *****************/

    SpinColor_fld *dirac_p(SpinColor_fld*, SpinColor_fld*, int, int, int);

/***************** Fast Fourier Transform *****************/

    void fft(int);

/************* DIRAC DAGGERED * DIRAC OPERATOR IN MOMENTUM SPACE **************/
    void dd(ptGluon_fld&, int);

    // in e out POSSONO essere sovrapposti
    void dd0inv(SpinColor_fld*, SpinColor_fld*, dgsp&, eigv*, int);
    
    // in e out POSSONO essere sovrapposti
    void dd1inv(SpinColor_fld*, SpinColor_fld*, dgsp&, eigv*, int, eigsysws*);

    // in e out POSSONO essere sovrapposti
    void res0(SpinColor_fld*, SpinColor_fld*, ptGluon_fld&, dgsp&, eigv*, int, int);

    // in e out POSSONO essere sovrapposti
    void res1(SpinColor_fld*, SpinColor_fld*, ptGluon_fld&, eigsysws*, dgsp&, eigv*, int, int);

};



// Per visualizzare i dati
void header(dgsp*, int, int);

void prline(dgsp*, eigv*, int, int, int);

void rawline(dgsp*, eigv*, int, int, int, FILE*);

void tail(dgsp*, eigv*, int, double*, double*, double*, int);



// Conteggio del numero di sottospazi di degenerazione
void conta_sottospazi(eigv*, int&, int, int, int&, int);

// Crea l struttura dati relativa al sottospazio di degenerazione
void degen_subspc(eigv*, dgsp*, int, int, int, int&, int);

void eigenproblem(eigsysws&, dgsp&, diracws&, eigv*, ptGluon_fld&, int);

void non_degen_solver(diracws&, eigv*, eigsysws&, dgsp&, ptGluon_fld&, int, int);
