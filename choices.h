/* Condizioni al contorno in lattice.h  */
/* E' l'unico altroparametro da settare */

#define dim 4
#define allocORD 12
extern int PTORD;


#define _SQRT_BETA_EXPANSION_

#define NF 4
#define NC 3
#define MBARE 0

#define WIL 0
#define TLSYM 1
#define IWA 2
#define DBW2 3
#define GAUGE_ACTION WIL //Alternatives: WIL, TLSYM, IWA, DBW2

#ifdef __PARALLEL_OMP__
#define ntt  1
#define ntx  1
#define nty  1
#define ntz  1
#define NTHR  (ntx*nty*ntz*ntt)

//#define NTHR 8
#else
#define NTHR 1
#endif

#define __ZERO_FLAG__ 1



