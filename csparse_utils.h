//
// Created by John Lambert-Admin on 4/20/19.
//

#pragma once

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




problem *get_problem (FILE *f, double tol) ;
int chol_ordering_demo (problem *Prob) ;
problem *free_problem (problem *Prob) ;


/* cs_demo2: read a matrix and solve a linear system */
void convert_eigen_to_csparse (void)
{
    problem *Prob = get_problem (stdin, 1e-14) ;
    chol_ordering_demo (Prob) ;
    free_problem (Prob) ;
}


/* cs_demo2: read a matrix and solve a linear system */
void convert_csparse_to_eigen (void)
{
//    problem *Prob = get_problem (stdin, 1e-14) ;
//    demo2 (Prob) ;
//    free_problem (Prob) ;
}





/* read a problem from a file */
problem *get_problem (FILE *f, double tol)
{
    cs *T, *A, *C ;
    int sym, m, n, mn, nz1, nz2 ;
    problem *Prob ;
    Prob = cs_calloc (1, sizeof (problem)) ;
    if (!Prob) return (NULL) ;
    T = cs_load (f) ;			/* load triplet matrix T from a file */
    Prob->A = A = cs_triplet (T) ;	/* A = compressed-column form of T */
    cs_spfree (T) ;			/* clear T */
    if (!cs_dupl (A)) return (free_problem (Prob)) ; /* sum up duplicates */
    Prob->sym = sym = is_sym (A) ;	/* determine if A is symmetric */
    m = A->m ; n = A->n ;
    mn = CS_MAX (m,n) ;
    nz1 = A->p [n] ;
    cs_dropzeros (A) ;			/* drop zero entries */
    nz2 = A->p [n] ;
    if (tol > 0) cs_droptol (A, tol) ;	/* drop tiny entries (just to test) */
    Prob->C = C = sym ? make_sym (A) : A ;  /* C = A + triu(A,1)', or C=A */
    if (!C) return (free_problem (Prob)) ;
    printf ("\n--- Matrix: %d-by-%d, nnz: %d (sym: %d: nnz %d), norm: %8.2e\n",
            m, n, A->p [n], sym, sym ? C->p [n] : 0, cs_norm (C)) ;
    if (nz1 != nz2) printf ("zero entries dropped: %d\n", nz1 - nz2) ;
    if (nz2 != A->p [n]) printf ("tiny entries dropped: %d\n", nz2 - A->p [n]) ;
    Prob->b = cs_malloc (mn, sizeof (double)) ;
    Prob->x = cs_malloc (mn, sizeof (double)) ;
    Prob->r = cs_malloc (mn, sizeof (double)) ;
    return ((!Prob->b || !Prob->x || !Prob->r) ? free_problem (Prob) : Prob) ;
}


/* solve a linear system using Cholesky with various orderings */
int chol_ordering_demo (problem *Prob)
{
    cs *A, *C ;
    double *b, *x, *r,  t, tol ;
    int k, m, n, ok, order, nb, ns, *R, *S, *rr, sprank ;
    csd *D ;
    if (!Prob) return (0) ;
    A = Prob->A ; C = Prob->C ; b = Prob->b ; x = Prob->x ; r = Prob->r ;
    m = A->m ; n = A->n ;
    tol = Prob->sym ? 0.001 : 1 ;		/* partial pivoting tolerance */
    D = cs_dmperm (C) ;				/* dmperm analysis */
    if (!D) return (0) ;
    nb = D->nb ; R = D->R ; S = D->S ; rr = D->rr ;
    sprank = rr [3] ;
    for (ns = 0, k = 0 ; k < nb ; k++)
    {
        ns += ((R [k+1] == R [k]+1) && (S [k+1] == S [k]+1)) ;
    }
    printf ("blocks: %d singletons: %d structural rank: %d\n", nb, ns, sprank) ;
    cs_dfree (D) ;

    if (!Prob->sym) return (1) ;
    for (order = -1 ; order <= 0 ; order++)	/* natural and amd(A+A') */
    {
        if (order == -1 && m > 1000) continue ;
        printf ("Chol ") ;
        print_order (order) ;
        rhs (x, b, m) ;				/* compute right-hand-side */
        t = tic () ;
        ok = cs_cholsol (C, x, order) ;		/* solve Ax=b with Cholesky */
        printf ("time: %8.2f ", toc (t)) ;
        resid (ok, C, x, b, r) ;		/* print residual */
    }
    return (1) ;
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

/* infinity-norm of x */
static double norm (double *x, int n)
{
    int i ;
    double normx = 0 ;
    for (i = 0 ; i < n ; i++) normx = CS_MAX (normx, fabs (x [i])) ;
    return (normx) ;
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
    return (cs_free (Prob)) ;
}
