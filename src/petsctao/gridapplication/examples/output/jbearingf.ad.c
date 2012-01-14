/************************** DISCLAIMER ********************************/
/*                                                                    */
/*   This file was generated on 08/09/02 10:02:12 by the version of   */
/*   ADIC compiled on  08/07/00 18:06:31                              */
/*                                                                    */
/*   ADIC was prepared as an account of work sponsored by an          */
/*   agency of the United States Government and the University of     */
/*   Chicago.  NEITHER THE AUTHOR(S), THE UNITED STATES GOVERNMENT    */
/*   NOR ANY AGENCY THEREOF, NOR THE UNIVERSITY OF CHICAGO, INCLUDING */
/*   ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY, EXPRESS  */
/*   OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR */
/*   THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION OR  */
/*   PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE */
/*   PRIVATELY OWNED RIGHTS.                                          */
/*                                                                    */
/**********************************************************************/
#include "ad_deriv.h"
#include <math.h>
#include "adintrinsics.h"
typedef struct  {
InactiveDouble  ecc;
InactiveDouble  b;
InactiveDouble  *wq, *wl;
InactiveDouble  hx, hy;
InactiveDouble  area;
int  mx, my;
}
AppCtx;
int  ad_JBearLocalFunction(int  coor[2],DERIV_TYPE  x[4],DERIV_TYPE  *f,void   *ptr) {
AppCtx  *user = (AppCtx * )ptr;
InactiveDouble  hx, hy, area;
InactiveDouble  *wq, *wl;
DERIV_TYPE  avgWq, sqGrad, avgWV;
DERIV_TYPE  dvdx, dvdy, fl, fu;
int  i;
    double  ad_loc_0;
    double  ad_loc_1;
    double  ad_adj_0;
    double  ad_adj_1;
    double  ad_loc_2;
    double  ad_loc_3;
    double  ad_loc_4;
    double  ad_loc_5;
    double  ad_adj_2;
    hx = user->hx;
    hy = user->hy;
    area = user->area;
    wq = user->wq;
    wl = user->wl;
    i = coor[0];
    {
        ad_loc_0 = DERIV_val(x[0]) - DERIV_val(x[1]);
        ad_loc_1 = ad_loc_0 / hx;
        ad_adj_0 = 1.000000000000000e+00 / hx;
        ad_adj_1 =  -ad_adj_0;
        ad_grad_axpy_2(&(dvdx), ad_adj_0, &(x[0]), ad_adj_1, &(x[1]));
        DERIV_val(dvdx) = ad_loc_1;
    }
    {
        ad_loc_0 = DERIV_val(x[0]) - DERIV_val(x[2]);
        ad_loc_1 = ad_loc_0 / hy;
        ad_adj_0 = 1.000000000000000e+00 / hy;
        ad_adj_1 =  -ad_adj_0;
        ad_grad_axpy_2(&(dvdy), ad_adj_0, &(x[0]), ad_adj_1, &(x[2]));
        DERIV_val(dvdy) = ad_loc_1;
    }
    {
        ad_loc_0 = DERIV_val(dvdx) * DERIV_val(dvdx);
        ad_loc_1 = DERIV_val(dvdy) * DERIV_val(dvdy);
        ad_loc_2 = ad_loc_0 + ad_loc_1;
        ad_grad_axpy_4(&(sqGrad), DERIV_val(dvdx), &(dvdx), DERIV_val(dvdx), &(dvdx), DERIV_val(dvdy), &(dvdy), DERIV_val(dvdy), &(dvdy));
        DERIV_val(sqGrad) = ad_loc_2;
    }
    {
        ad_loc_0 = 2.0 * wq[i];
        ad_loc_1 = ad_loc_0 + wq[i + 1];
        ad_loc_2 = ad_loc_1 / 3.0;
        ad_grad_axpy_0(&(avgWq));
        DERIV_val(avgWq) = ad_loc_2;
    }
    {
        ad_loc_0 = wl[i] * DERIV_val(x[0]);
        ad_loc_1 = wl[i + 1] * DERIV_val(x[1]);
        ad_loc_2 = ad_loc_0 + ad_loc_1;
        ad_loc_3 = wl[i] * DERIV_val(x[2]);
        ad_loc_4 = ad_loc_2 + ad_loc_3;
        ad_loc_5 = ad_loc_4 / 3.0;
        ad_adj_0 = wl[i] * 3.333333333333333e-01;
        ad_adj_1 = wl[i + 1] * 3.333333333333333e-01;
        ad_adj_2 = wl[i] * 3.333333333333333e-01;
        ad_grad_axpy_3(&(avgWV), ad_adj_2, &(x[0]), ad_adj_1, &(x[1]), ad_adj_0, &(x[2]));
        DERIV_val(avgWV) = ad_loc_5;
    }
    {
        ad_loc_0 = 0.5 * DERIV_val(avgWq);
        ad_loc_1 = ad_loc_0 * DERIV_val(sqGrad);
        ad_loc_2 = ad_loc_1 - DERIV_val(avgWV);
        ad_adj_0 = 0.5 * DERIV_val(sqGrad);
        ad_grad_axpy_3(&(fl), ad_adj_0, &(avgWq), ad_loc_0, &(sqGrad), -1.000000000000000e+00, &(avgWV));
        DERIV_val(fl) = ad_loc_2;
    }
    {
        ad_loc_0 = DERIV_val(x[3]) - DERIV_val(x[2]);
        ad_loc_1 = ad_loc_0 / hx;
        ad_adj_0 = 1.000000000000000e+00 / hx;
        ad_adj_1 =  -ad_adj_0;
        ad_grad_axpy_2(&(dvdx), ad_adj_0, &(x[3]), ad_adj_1, &(x[2]));
        DERIV_val(dvdx) = ad_loc_1;
    }
    {
        ad_loc_0 = DERIV_val(x[3]) - DERIV_val(x[1]);
        ad_loc_1 = ad_loc_0 / hy;
        ad_adj_0 = 1.000000000000000e+00 / hy;
        ad_adj_1 =  -ad_adj_0;
        ad_grad_axpy_2(&(dvdy), ad_adj_0, &(x[3]), ad_adj_1, &(x[1]));
        DERIV_val(dvdy) = ad_loc_1;
    }
    {
        ad_loc_0 = DERIV_val(dvdx) * DERIV_val(dvdx);
        ad_loc_1 = DERIV_val(dvdy) * DERIV_val(dvdy);
        ad_loc_2 = ad_loc_0 + ad_loc_1;
        ad_grad_axpy_4(&(sqGrad), DERIV_val(dvdx), &(dvdx), DERIV_val(dvdx), &(dvdx), DERIV_val(dvdy), &(dvdy), DERIV_val(dvdy), &(dvdy));
        DERIV_val(sqGrad) = ad_loc_2;
    }
    {
        ad_loc_0 = 2.0 * wq[i + 1];
        ad_loc_1 = ad_loc_0 + wq[i];
        ad_loc_2 = ad_loc_1 / 3.0;
        ad_grad_axpy_0(&(avgWq));
        DERIV_val(avgWq) = ad_loc_2;
    }
    {
        ad_loc_0 = wl[i + 1] * DERIV_val(x[1]);
        ad_loc_1 = wl[i] * DERIV_val(x[2]);
        ad_loc_2 = ad_loc_0 + ad_loc_1;
        ad_loc_3 = wl[i + 1] * DERIV_val(x[3]);
        ad_loc_4 = ad_loc_2 + ad_loc_3;
        ad_loc_5 = ad_loc_4 / 3.0;
        ad_adj_0 = wl[i + 1] * 3.333333333333333e-01;
        ad_adj_1 = wl[i] * 3.333333333333333e-01;
        ad_adj_2 = wl[i + 1] * 3.333333333333333e-01;
        ad_grad_axpy_3(&(avgWV), ad_adj_2, &(x[1]), ad_adj_1, &(x[2]), ad_adj_0, &(x[3]));
        DERIV_val(avgWV) = ad_loc_5;
    }
    {
        ad_loc_0 = 0.5 * DERIV_val(avgWq);
        ad_loc_1 = ad_loc_0 * DERIV_val(sqGrad);
        ad_loc_2 = ad_loc_1 - DERIV_val(avgWV);
        ad_adj_0 = 0.5 * DERIV_val(sqGrad);
        ad_grad_axpy_3(&(fu), ad_adj_0, &(avgWq), ad_loc_0, &(sqGrad), -1.000000000000000e+00, &(avgWV));
        DERIV_val(fu) = ad_loc_2;
    }
    {
        ad_loc_0 = DERIV_val(fl) + DERIV_val(fu);
        ad_loc_1 = area * ad_loc_0;
        ad_grad_axpy_2(&(*f), area, &(fl), area, &(fu));
        DERIV_val(*f) = ad_loc_1;
    }
    return 0;
}
void   ad_AD_Init(int  arg0) {
    ad_AD_GradInit(arg0);

}
void   ad_AD_Final() {
    ad_AD_GradFinal();

}
