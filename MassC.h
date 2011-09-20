#include<iostream>

#if GAUGE_ACTION == WIL
#if NF == 0
double mcpt[6 + 1] = {4., 0., -2.6057, 0., -4.5203, 0., -13.11};
#elif NF == 1
double mcpt[6 + 1] = {4., 0., -2.6057, 0., -4.4064, 0., -12.487};
#elif NF == 2
//double mcpt[6 + 1] = {4., 0., 0., 0., 0., 0., 0.};
double mcpt[6 + 1] = {4., 0., -2.6057, 0., -4.2925, 0., -11.78};
#elif NF == 3
double mcpt[6 + 1] = {4., 0., -2.6057, 0., -4.1786, 0., -11.02};
#elif NF == 4
double mcpt[6 + 1] = {4., 0., -2.6057, 0., -4.0647, 0., -10.24};
#endif

#elif GAUGE_ACTION == TLSYM
#if NF == 0
double mcpt[6 + 1] = {4., 0., -2.0489, 0., -2.14  , 0., 0};
#elif NF == 1
double mcpt[6 + 1] = {4., 0., -2.0489, 0., -2.0618, 0., 0};
#elif NF == 2
double mcpt[6 + 1] = {4., 0., -2.0489, 0., -1.9836, 0., -3.99};
#elif NF == 3
double mcpt[6 + 1] = {4., 0., -2.0489, 0., -1.9053, 0., 0};
#elif NF == 4
double mcpt[6 + 1] = {4., 0., -2.0489, 0., -1.8271, 0., 0};
#endif

#elif GAUGE_ACTION == IWA
#if NF == 0
double mcpt[6 + 1] = {4., 0., -1.3209, 0., -0.3583, 0., 0};
#elif NF == 1
double mcpt[6 + 1] = {4., 0., -1.3209, 0., -0.3165, 0., 0};
#elif NF == 2
double mcpt[6 + 1] = {4., 0., -1.3209, 0., -0.2747, 0., 0};
#elif NF == 3
double mcpt[6 + 1] = {4., 0., -1.3209, 0., -0.2329, 0., 0};
#elif NF == 4
double mcpt[6 + 1] = {4., 0., -1.3209, 0., -0.1911, 0., 0};
#endif

#elif GAUGE_ACTION == DBW2
#if NF == 0
double mcpt[6 + 1] = {4., 0., -0.5832, 0., 0.1836, 0., 0};
#elif NF == 1
double mcpt[6 + 1] = {4., 0., -0.5832, 0., 0.1988, 0., 0};
#elif NF == 2
double mcpt[6 + 1] = {4., 0., -0.5832, 0., 0.214 , 0., 0};
#elif NF == 3
double mcpt[6 + 1] = {4., 0., -0.5832, 0., 0.2292, 0., 0};
#elif NF == 4
double mcpt[6 + 1] = {4., 0., -0.5832, 0., 0.2444, 0., 0};
#endif
#endif



