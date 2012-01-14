/************************** DISCLAIMER ********************************/
/*                                                                    */
/*   This file was generated on 02/21/03 17:34:11 by the version of   */
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
#define ad_GRAD_MAX 16
typedef struct  {
int  mx, my;
InactiveDouble  hx, hy;
InactiveDouble  twopivornum;
InactiveDouble  tkappa;
}
AppCtx;
int  ad_GinzLandLocalFunction(int  coor[2],DERIV_TYPE  xx[16],DERIV_TYPE  *f,void   *ptr) {
AppCtx  *user = (AppCtx * )ptr;
InactiveDouble  hx = user->hx, hy = user->hy;
InactiveDouble  tkappa = user->tkappa;
InactiveDouble  twopivornum = user->twopivornum;
int  i, mx = user->mx, my = user->my, sqn = mx * my;
DERIV_TYPE  arg;
DERIV_TYPE  d1, d2;
DERIV_TYPE  x1, x2, xy;
DERIV_TYPE  delsq, fcond, fkin, ffield;
DERIV_TYPE  x[16];
int  ad_var_0;
DERIV_TYPE  ad_var_1, ad_var_2, ad_var_3, ad_var_4, ad_var_5, ad_var_6, ad_var_7, ad_var_8, ad_var_9, ad_var_10, ad_var_11, ad_var_12;
double  ad_adji_0;
    double  ad_loc_0;
    double  ad_loc_1;
    double  ad_loc_2;
    double  ad_adj_0;
    double  ad_adj_1;
    double  ad_loc_3;
    double  ad_adj_2;
    double  ad_adj_3;
    double  ad_loc_4;
    double  ad_adj_4;
    double  ad_loc_5;

        static int g_filenum = 0;
        if (g_filenum == 0) {
            adintr_ehsfid(&g_filenum, __FILE__, "ad_GinzLandLocalFunction");
        }
            for (i = 0;     i < 16;     )    {
        {
            ad_grad_axpy_copy(&(x[i]), &(xx[i]));
            DERIV_val(x[i]) = DERIV_val(xx[i]);
        }
        ad_var_0 = i++;
    }
    if (coor[0] == mx - 1)     {
        {
            ad_loc_0 = twopivornum * coor[1];
            ad_loc_1 = ad_loc_0 / my;
            ad_grad_axpy_0(&(arg));
            DERIV_val(arg) = ad_loc_1;
        }
     DERIV_val(ad_var_1) = cos( DERIV_val(arg)); /*cos*/
      ad_adji_0 = -sin( DERIV_val(arg));
        {
            ad_grad_axpy_1(&(ad_var_1), ad_adji_0, &(arg));
        }
     DERIV_val(ad_var_2) = sin( DERIV_val(arg)); /*sin*/
      ad_adji_0 = cos( DERIV_val(arg));
        {
            ad_grad_axpy_1(&(ad_var_2), ad_adji_0, &(arg));
        }
        {
            ad_loc_0 = DERIV_val(x[4]) * DERIV_val(ad_var_1);
            ad_loc_1 = DERIV_val(x[5]) * DERIV_val(ad_var_2);
            ad_loc_2 = ad_loc_0 - ad_loc_1;
            ad_adj_0 =  -DERIV_val(x[5]);
            ad_adj_1 =  -DERIV_val(ad_var_2);
            ad_grad_axpy_4(&(d1), DERIV_val(ad_var_1), &(x[4]), DERIV_val(x[4]), &(ad_var_1), ad_adj_1, &(x[5]), ad_adj_0, &(ad_var_2));
            DERIV_val(d1) = ad_loc_2;
        }
     DERIV_val(ad_var_3) = sin( DERIV_val(arg)); /*sin*/
      ad_adji_0 = cos( DERIV_val(arg));
        {
            ad_grad_axpy_1(&(ad_var_3), ad_adji_0, &(arg));
        }
     DERIV_val(ad_var_4) = cos( DERIV_val(arg)); /*cos*/
      ad_adji_0 = -sin( DERIV_val(arg));
        {
            ad_grad_axpy_1(&(ad_var_4), ad_adji_0, &(arg));
        }
        {
            ad_loc_0 = DERIV_val(x[4]) * DERIV_val(ad_var_3);
            ad_loc_1 = DERIV_val(x[5]) * DERIV_val(ad_var_4);
            ad_loc_2 = ad_loc_0 + ad_loc_1;
            ad_grad_axpy_4(&(d2), DERIV_val(ad_var_3), &(x[4]), DERIV_val(x[4]), &(ad_var_3), DERIV_val(ad_var_4), &(x[5]), DERIV_val(x[5]), &(ad_var_4));
            DERIV_val(d2) = ad_loc_2;
        }
        {
            ad_grad_axpy_copy(&(x[4]), &(d1));
            DERIV_val(x[4]) = DERIV_val(d1);
        }
        {
            ad_grad_axpy_copy(&(x[5]), &(d2));
            DERIV_val(x[5]) = DERIV_val(d2);
        }
        {
            ad_loc_0 = my * hy;
            ad_loc_1 = twopivornum / ad_loc_0;
            ad_loc_2 = DERIV_val(x[7]) + ad_loc_1;
            ad_grad_axpy_1(&(x[7]), 1.000000000000000e+00, &(x[7]));
            DERIV_val(x[7]) = ad_loc_2;
        }
    }
    {
        ad_loc_0 = DERIV_val(x[0]) * DERIV_val(x[0]);
        ad_loc_1 = DERIV_val(x[1]) * DERIV_val(x[1]);
        ad_loc_2 = ad_loc_0 + ad_loc_1;
        ad_grad_axpy_4(&(delsq), DERIV_val(x[0]), &(x[0]), DERIV_val(x[0]), &(x[0]), DERIV_val(x[1]), &(x[1]), DERIV_val(x[1]), &(x[1]));
        DERIV_val(delsq) = ad_loc_2;
    }
    {
        ad_loc_0 =  -DERIV_val(delsq);
        ad_loc_1 = DERIV_val(delsq) * DERIV_val(delsq);
        ad_loc_2 = ad_loc_1 / 2.0;
        ad_loc_3 = ad_loc_0 + ad_loc_2;
        ad_adj_0 = DERIV_val(delsq) * 5.000000000000000e-01;
        ad_adj_1 = DERIV_val(delsq) * 5.000000000000000e-01;
        ad_grad_axpy_3(&(fcond), -1.000000000000000e+00, &(delsq), ad_adj_1, &(delsq), ad_adj_0, &(delsq));
        DERIV_val(fcond) = ad_loc_3;
    }
    {
        ad_loc_0 = hx * DERIV_val(x[2]);
        ad_grad_axpy_1(&(d1), hx, &(x[2]));
        DERIV_val(d1) = ad_loc_0;
    }
     DERIV_val(ad_var_5) = cos( DERIV_val(d1)); /*cos*/
      ad_adji_0 = -sin( DERIV_val(d1));
    {
        ad_grad_axpy_1(&(ad_var_5), ad_adji_0, &(d1));
    }
     DERIV_val(ad_var_6) = sin( DERIV_val(d1)); /*sin*/
      ad_adji_0 = cos( DERIV_val(d1));
    {
        ad_grad_axpy_1(&(ad_var_6), ad_adji_0, &(d1));
    }
    {
        ad_loc_0 = DERIV_val(x[0]) * DERIV_val(ad_var_5);
        ad_loc_1 = DERIV_val(x[4]) - ad_loc_0;
        ad_loc_2 = DERIV_val(x[1]) * DERIV_val(ad_var_6);
        ad_loc_3 = ad_loc_1 + ad_loc_2;
        ad_adj_0 =  -DERIV_val(x[0]);
        ad_adj_1 =  -DERIV_val(ad_var_5);
        ad_grad_axpy_5(&(x1), 1.000000000000000e+00, &(x[4]), ad_adj_1, &(x[0]), ad_adj_0, &(ad_var_5), DERIV_val(ad_var_6), &(x[1]), DERIV_val(x[1]), &(ad_var_6));
        DERIV_val(x1) = ad_loc_3;
    }
     DERIV_val(ad_var_7) = cos( DERIV_val(d1)); /*cos*/
      ad_adji_0 = -sin( DERIV_val(d1));
    {
        ad_grad_axpy_1(&(ad_var_7), ad_adji_0, &(d1));
    }
     DERIV_val(ad_var_8) = sin( DERIV_val(d1)); /*sin*/
      ad_adji_0 = cos( DERIV_val(d1));
    {
        ad_grad_axpy_1(&(ad_var_8), ad_adji_0, &(d1));
    }
    {
        ad_loc_0 = DERIV_val(x[1]) * DERIV_val(ad_var_7);
        ad_loc_1 = DERIV_val(x[5]) - ad_loc_0;
        ad_loc_2 = DERIV_val(x[0]) * DERIV_val(ad_var_8);
        ad_loc_3 = ad_loc_1 - ad_loc_2;
        ad_adj_0 =  -DERIV_val(x[0]);
        ad_adj_1 =  -DERIV_val(ad_var_8);
        ad_adj_2 =  -DERIV_val(x[1]);
        ad_adj_3 =  -DERIV_val(ad_var_7);
        ad_grad_axpy_5(&(x2), 1.000000000000000e+00, &(x[5]), ad_adj_3, &(x[1]), ad_adj_2, &(ad_var_7), ad_adj_1, &(x[0]), ad_adj_0, &(ad_var_8));
        DERIV_val(x2) = ad_loc_3;
    }
    {
        ad_loc_0 = DERIV_val(x1) * DERIV_val(x1);
        ad_loc_1 = DERIV_val(x2) * DERIV_val(x2);
        ad_loc_2 = ad_loc_0 + ad_loc_1;
        ad_loc_3 = hx * hx;
        ad_loc_4 = ad_loc_2 / ad_loc_3;
        ad_adj_0 = 1.000000000000000e+00 / ad_loc_3;
        ad_adj_1 = DERIV_val(x2) * ad_adj_0;
        ad_adj_2 = DERIV_val(x2) * ad_adj_0;
        ad_adj_3 = DERIV_val(x1) * ad_adj_0;
        ad_adj_4 = DERIV_val(x1) * ad_adj_0;
        ad_grad_axpy_4(&(fkin), ad_adj_4, &(x1), ad_adj_3, &(x1), ad_adj_2, &(x2), ad_adj_1, &(x2));
        DERIV_val(fkin) = ad_loc_4;
    }
    {
        ad_loc_0 = hy * DERIV_val(x[3]);
        ad_grad_axpy_1(&(d2), hy, &(x[3]));
        DERIV_val(d2) = ad_loc_0;
    }
     DERIV_val(ad_var_9) = cos( DERIV_val(d2)); /*cos*/
      ad_adji_0 = -sin( DERIV_val(d2));
    {
        ad_grad_axpy_1(&(ad_var_9), ad_adji_0, &(d2));
    }
     DERIV_val(ad_var_10) = sin( DERIV_val(d2)); /*sin*/
      ad_adji_0 = cos( DERIV_val(d2));
    {
        ad_grad_axpy_1(&(ad_var_10), ad_adji_0, &(d2));
    }
    {
        ad_loc_0 = DERIV_val(x[0]) * DERIV_val(ad_var_9);
        ad_loc_1 = DERIV_val(x[8]) - ad_loc_0;
        ad_loc_2 = DERIV_val(x[1]) * DERIV_val(ad_var_10);
        ad_loc_3 = ad_loc_1 + ad_loc_2;
        ad_adj_0 =  -DERIV_val(x[0]);
        ad_adj_1 =  -DERIV_val(ad_var_9);
        ad_grad_axpy_5(&(x1), 1.000000000000000e+00, &(x[8]), ad_adj_1, &(x[0]), ad_adj_0, &(ad_var_9), DERIV_val(ad_var_10), &(x[1]), DERIV_val(x[1]), &(ad_var_10));
        DERIV_val(x1) = ad_loc_3;
    }
     DERIV_val(ad_var_11) = cos( DERIV_val(d2)); /*cos*/
      ad_adji_0 = -sin( DERIV_val(d2));
    {
        ad_grad_axpy_1(&(ad_var_11), ad_adji_0, &(d2));
    }
     DERIV_val(ad_var_12) = sin( DERIV_val(d2)); /*sin*/
      ad_adji_0 = cos( DERIV_val(d2));
    {
        ad_grad_axpy_1(&(ad_var_12), ad_adji_0, &(d2));
    }
    {
        ad_loc_0 = DERIV_val(x[1]) * DERIV_val(ad_var_11);
        ad_loc_1 = DERIV_val(x[9]) - ad_loc_0;
        ad_loc_2 = DERIV_val(x[0]) * DERIV_val(ad_var_12);
        ad_loc_3 = ad_loc_1 - ad_loc_2;
        ad_adj_0 =  -DERIV_val(x[0]);
        ad_adj_1 =  -DERIV_val(ad_var_12);
        ad_adj_2 =  -DERIV_val(x[1]);
        ad_adj_3 =  -DERIV_val(ad_var_11);
        ad_grad_axpy_5(&(x2), 1.000000000000000e+00, &(x[9]), ad_adj_3, &(x[1]), ad_adj_2, &(ad_var_11), ad_adj_1, &(x[0]), ad_adj_0, &(ad_var_12));
        DERIV_val(x2) = ad_loc_3;
    }
    {
        ad_loc_0 = DERIV_val(x1) * DERIV_val(x1);
        ad_loc_1 = DERIV_val(x2) * DERIV_val(x2);
        ad_loc_2 = ad_loc_0 + ad_loc_1;
        ad_loc_3 = hy * hy;
        ad_loc_4 = ad_loc_2 / ad_loc_3;
        ad_loc_5 = DERIV_val(fkin) + ad_loc_4;
        ad_adj_0 = 1.000000000000000e+00 / ad_loc_3;
        ad_adj_1 = DERIV_val(x2) * ad_adj_0;
        ad_adj_2 = DERIV_val(x2) * ad_adj_0;
        ad_adj_3 = DERIV_val(x1) * ad_adj_0;
        ad_adj_4 = DERIV_val(x1) * ad_adj_0;
        ad_grad_axpy_5(&(fkin), 1.000000000000000e+00, &(fkin), ad_adj_4, &(x1), ad_adj_3, &(x1), ad_adj_2, &(x2), ad_adj_1, &(x2));
        DERIV_val(fkin) = ad_loc_5;
    }
    {
        ad_loc_0 = DERIV_val(x[7]) - DERIV_val(x[3]);
        ad_loc_1 = ad_loc_0 / hx;
        ad_loc_2 = DERIV_val(x[10]) - DERIV_val(x[2]);
        ad_loc_3 = ad_loc_2 / hy;
        ad_loc_4 = ad_loc_1 - ad_loc_3;
        ad_adj_0 =  -(1.000000000000000e+00 / hy);
        ad_adj_1 =  -ad_adj_0;
        ad_adj_2 = 1.000000000000000e+00 / hx;
        ad_adj_3 =  -ad_adj_2;
        ad_grad_axpy_4(&(xy), ad_adj_2, &(x[7]), ad_adj_3, &(x[3]), ad_adj_0, &(x[10]), ad_adj_1, &(x[2]));
        DERIV_val(xy) = ad_loc_4;
    }
    {
        ad_loc_0 = DERIV_val(xy) * DERIV_val(xy);
        ad_loc_1 = tkappa * tkappa;
        ad_loc_2 = ad_loc_0 * ad_loc_1;
        ad_adj_0 = DERIV_val(xy) * ad_loc_1;
        ad_adj_1 = DERIV_val(xy) * ad_loc_1;
        ad_grad_axpy_2(&(ffield), ad_adj_1, &(xy), ad_adj_0, &(xy));
        DERIV_val(ffield) = ad_loc_2;
    }
    {
        ad_loc_0 = DERIV_val(fcond) + DERIV_val(fkin);
        ad_loc_1 = ad_loc_0 + DERIV_val(ffield);
        ad_loc_2 = ad_loc_1 / sqn;
        ad_adj_0 = 1.000000000000000e+00 / sqn;
        ad_grad_axpy_3(&(*f), ad_adj_0, &(fcond), ad_adj_0, &(fkin), ad_adj_0, &(ffield));
        DERIV_val(*f) = ad_loc_2;
    }
    return 0;
}
void   ad_AD_Init(int  arg0) {
    ad_AD_GradInit(arg0);

}
void   ad_AD_Final() {
    ad_AD_GradFinal();

}
