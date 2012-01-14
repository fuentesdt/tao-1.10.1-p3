/************************** DISCLAIMER ********************************/
/*                                                                    */
/*   This file was generated on 08/09/02 09:55:56 by the version of   */
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
InactiveDouble  param;
InactiveDouble  hx, hy, area;
int  mx, my;
}
AppCtx;
int  ad_CombLocalFunction(int  coor[2],DERIV_TYPE  x[4],DERIV_TYPE  *f,void   *ptr) {
AppCtx  *user = (AppCtx * )ptr;
InactiveDouble  lambda, hx, hy, area;
DERIV_TYPE  fquad, fexp, dvdx, dvdy;
DERIV_TYPE  ad_var_0, ad_var_1, ad_var_2, ad_var_3, ad_var_4, ad_var_5;
double  ad_adji_0;
    double  ad_loc_0;
    double  ad_loc_1;
    double  ad_adj_0;
    double  ad_adj_1;
    double  ad_loc_2;
    double  ad_loc_3;
    double  ad_loc_4;
    double  ad_adj_2;
    double  ad_adj_3;

        static int g_filenum = 0;
        if (g_filenum == 0) {
            adintr_ehsfid(&g_filenum, __FILE__, "ad_CombLocalFunction");
        }
            lambda = user->param;
    hx = user->hx;
    hy = user->hy;
    area = user->area;
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
        ad_grad_axpy_4(&(fquad), DERIV_val(dvdx), &(dvdx), DERIV_val(dvdx), &(dvdx), DERIV_val(dvdy), &(dvdy), DERIV_val(dvdy), &(dvdy));
        DERIV_val(fquad) = ad_loc_2;
    }
     DERIV_val(ad_var_0) = exp( DERIV_val(x[0]));  /*exp*/
      ad_adji_0 = DERIV_val(ad_var_0);
    {
        ad_grad_axpy_1(&(ad_var_0), ad_adji_0, &(x[0]));
    }
     DERIV_val(ad_var_1) = exp( DERIV_val(x[1]));  /*exp*/
      ad_adji_0 = DERIV_val(ad_var_1);
    {
        ad_grad_axpy_1(&(ad_var_1), ad_adji_0, &(x[1]));
    }
     DERIV_val(ad_var_2) = exp( DERIV_val(x[2]));  /*exp*/
      ad_adji_0 = DERIV_val(ad_var_2);
    {
        ad_grad_axpy_1(&(ad_var_2), ad_adji_0, &(x[2]));
    }
    {
        ad_loc_0 = DERIV_val(ad_var_0) + DERIV_val(ad_var_1);
        ad_loc_1 = ad_loc_0 + DERIV_val(ad_var_2);
        ad_grad_axpy_3(&(fexp), 1.000000000000000e+00, &(ad_var_0), 1.000000000000000e+00, &(ad_var_1), 1.000000000000000e+00, &(ad_var_2));
        DERIV_val(fexp) = ad_loc_1;
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
        ad_loc_3 = DERIV_val(fquad) + ad_loc_2;
        ad_grad_axpy_5(&(fquad), 1.000000000000000e+00, &(fquad), DERIV_val(dvdx), &(dvdx), DERIV_val(dvdx), &(dvdx), DERIV_val(dvdy), &(dvdy), DERIV_val(dvdy), &(dvdy));
        DERIV_val(fquad) = ad_loc_3;
    }
     DERIV_val(ad_var_3) = exp( DERIV_val(x[1]));  /*exp*/
      ad_adji_0 = DERIV_val(ad_var_3);
    {
        ad_grad_axpy_1(&(ad_var_3), ad_adji_0, &(x[1]));
    }
     DERIV_val(ad_var_4) = exp( DERIV_val(x[2]));  /*exp*/
      ad_adji_0 = DERIV_val(ad_var_4);
    {
        ad_grad_axpy_1(&(ad_var_4), ad_adji_0, &(x[2]));
    }
     DERIV_val(ad_var_5) = exp( DERIV_val(x[3]));  /*exp*/
      ad_adji_0 = DERIV_val(ad_var_5);
    {
        ad_grad_axpy_1(&(ad_var_5), ad_adji_0, &(x[3]));
    }
    {
        ad_loc_0 = DERIV_val(ad_var_3) + DERIV_val(ad_var_4);
        ad_loc_1 = ad_loc_0 + DERIV_val(ad_var_5);
        ad_loc_2 = DERIV_val(fexp) + ad_loc_1;
        ad_grad_axpy_4(&(fexp), 1.000000000000000e+00, &(fexp), 1.000000000000000e+00, &(ad_var_3), 1.000000000000000e+00, &(ad_var_4), 1.000000000000000e+00, &(ad_var_5));
        DERIV_val(fexp) = ad_loc_2;
    }
    {
        ad_loc_0 = 0.5 * DERIV_val(fquad);
        ad_loc_1 = lambda * DERIV_val(fexp);
        ad_loc_2 = ad_loc_1 / 3.0;
        ad_loc_3 = ad_loc_0 - ad_loc_2;
        ad_loc_4 = area * ad_loc_3;
        ad_adj_0 =  -area;
        ad_adj_1 = 3.333333333333333e-01 * ad_adj_0;
        ad_adj_2 = lambda * ad_adj_1;
        ad_adj_3 = 0.5 * area;
        ad_grad_axpy_2(&(*f), ad_adj_3, &(fquad), ad_adj_2, &(fexp));
        DERIV_val(*f) = ad_loc_4;
    }
    return 0;
}
void   ad_AD_Init(int  arg0) {
    ad_AD_GradInit(arg0);

}
void   ad_AD_Final() {
    ad_AD_GradFinal();

}
