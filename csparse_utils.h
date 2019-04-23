//
// Created by John Lambert-Admin on 4/20/19.
//

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

# include <stdlib.h>
# include <limits.h>
# include <math.h>
# include <stdio.h>
# include <time.h>

# include "csparse/csparse.h"



typedef struct problem_struct
{
    cs *A ;
    cs *C ;
    int sym ;
    double *x ;
    double *b ;
    double *r ;
} problem ;



problem *free_problem ( problem *Prob ) ;


/* 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise */
static int is_sym ( cs *A )
{
//    return 1;
    int is_upper, is_lower, j, p, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i ;
    if (m != n) return (0) ;
    is_upper = 1 ;
    is_lower = 1 ;
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            if (Ai [p] > j) is_upper = 0 ;
            if (Ai [p] < j) is_lower = 0 ;
        }
    }
    return (is_upper ? 1 : (is_lower ? -1 : 0)) ;
}

/* true for off-diagonal entries */
static int dropdiag ( int i, int j, double aij, void *other )
{
    return (i != j);
}

/* C = A + triu(A,1)' */
static cs *make_sym ( cs *A )
{
    cs *AT, *C ;
    AT = cs_transpose (A, 1) ;		/* AT = A' */
    cs_fkeep (AT, &dropdiag, NULL) ;	/* drop diagonal entries from AT */
    C = cs_add (A, AT, 1, 1) ;		/* C = A+AT */
    cs_spfree (AT) ;
    return (C) ;
}

/* create a right-hand-side */
static void rhs (double *x, double *b, int m)
{
    int i ;
//    for (i = 0 ; i < m ; i++) b [i] = 1 + ((double) i) / m ;

    b[0] = -0.12887603;
    b[1] = 4.9677217;
    b[2] = -3.16320475;
    b[3] = 0.3012804;

    for (i = 0 ; i < m ; i++) x [i] = b [i] ;


}

/* infinity-norm of x */
static double norm (double *x, int n)
{
    int i ;
    double normx = 0 ;
    for (i = 0 ; i < n ; i++) normx = CS_MAX (normx, fabs (x [i])) ;
    return (normx) ;
}

/* compute residual, norm(A*x-b,inf) / (norm(A,1)*norm(x,inf) + norm(b,inf)) */
static void resid (int ok, cs *A, double *x, double *b, double *r)
{
    int i, m, n ;
    if (!ok) { printf ("    (failed)\n") ; return ; }
    m = A->m ; n = A->n ;
    for (i = 0 ; i < m ; i++) r [i] = -b [i] ;	    /* r = -b */
    cs_gaxpy (A, x, r) ;			    /* r = r + A*x  */
    printf ("resid: %8.2e\n",
            norm (r,m) / ((n == 0) ? 1 : (cs_norm (A) * norm (x,n) + norm (b,m)))) ;
}

static double tic (void)
{
    return (clock () / (double) CLOCKS_PER_SEC) ;
}

static double toc (double t)
{
    double s = tic ();
    return (CS_MAX (0, s-t));
}

static void print_order (int order)
{
    switch (order)
    {
        case -1: printf ("natural    ") ; break ;
        case  0: printf ("amd(A+A')  ") ; break ;
        case  1: printf ("amd(S'*S)  ") ; break ;
        case  2: printf ("amd(A'*A)  ") ; break ;
    }
}


/* free a problem */
problem *free_problem (problem *Prob)
{
    if (!Prob) return (NULL) ;
    cs_spfree (Prob->A) ;
    if (Prob->sym) cs_spfree (Prob->C) ;
    cs_free (Prob->b) ;
    cs_free (Prob->x) ;
    cs_free (Prob->r) ;
    return (problem *) (cs_free (Prob)) ;
}


#ifdef __cplusplus
}
#endif
