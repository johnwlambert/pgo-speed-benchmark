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



/* CS_LOAD loads a triplet matrix from Eigen. */
cs *eigen_arr_to_csparse(double *A, size_t m, size_t n)
{
    int i_1d;
    double x ;
    cs *T ;
    T = cs_spalloc (0, 0, 1, 1, 1) ;

    size_t row = 0;
    size_t col = 0;

    for (i_1d=0; i_1d < m*n; i_1d++)
    {
        row = i_1d / n;
        col = i_1d % n;
        if (A[i_1d] != 0.0)
        {
            if (!cs_entry (T, row, col, A[i_1d] )) return (cs_spfree (T)) ;
        }
    }
    return (T) ;
}


problem *load_csparse_format(double *A_arr, double *b_arr, size_t num_rows, size_t num_cols)
{
    cs *T, *A, *C ;
    int sym, m, n, mn, nz1, nz2 ;
    problem *Prob ;
    Prob = (problem *) cs_calloc (1, sizeof (problem)) ;
    if (!Prob) return (NULL) ;
    T = eigen_arr_to_csparse(A_arr, num_rows, num_cols) ;			/* load triplet matrix T from Eigen "A_" */
    Prob->A = A = cs_triplet (T) ;	/* A = compressed-column form of T */
    cs_spfree (T) ;			/* clear T */
    if (!cs_dupl (A)) return (free_problem (Prob)) ; /* sum up duplicates */
    Prob->sym = sym = is_sym (A) ;	/* determine if A is symmetric */
    m = A->m ; n = A->n ;
    mn = CS_MAX (m,n) ;
    nz1 = A->p [n] ;
    cs_dropzeros (A) ;			/* drop zero entries */
    nz2 = A->p [n] ;
    Prob->C = C = sym ? make_sym (A) : A ;  /* C = A + triu(A,1)', or C=A */
    if (!C) return (free_problem (Prob)) ;
    printf ("\n--- Matrix: %d-by-%d, nnz: %d (sym: %d: nnz %d), norm: %8.2e\n",
            m, n, A->p [n], sym, sym ? C->p [n] : 0, cs_norm (C)) ;
    if (nz1 != nz2) printf ("zero entries dropped: %d\n", nz1 - nz2) ;
    if (nz2 != A->p [n]) printf ("tiny entries dropped: %d\n", nz2 - A->p [n]) ;
    Prob->b = (double *) cs_malloc (mn, sizeof (double)) ;

    for (size_t i = 0 ; i < num_rows ; i++) Prob->b[i] = b_arr [i] ;

    Prob->x =  (double *) cs_malloc (mn, sizeof (double)) ;
    Prob->r =  (double *) cs_malloc (mn, sizeof (double)) ;
    return ((!Prob->b || !Prob->x || !Prob->r) ? free_problem (Prob) : Prob) ;
}








/* solve a linear system using Cholesky, LU, and QR, with various orderings
 *
 * QR way too slow... QR   amd(A'*A)  time:     0.36 resid: 1.53e-16
 * */
void solve_with_cs_chol (problem *Prob) //, std::string solver_type)
{
    cs *A, *C ;
    double *b, *x, *r,  t, tol ;
    int k, m, n, ok, order, nb, ns, *R, *S, *rr, sprank ;
    csd *D ;
    if (!Prob) return ;
    A = Prob->A ; C = Prob->C ; b = Prob->b ; x = Prob->x ; r = Prob->r ;
    m = A->m ; n = A->n ;
    tol = Prob->sym ? 0.001 : 1 ;		/* partial pivoting tolerance */
    D = cs_dmperm (C) ;				/* dmperm analysis */
    if (!D) return ;
    nb = D->nb ; R = D->R ; S = D->S ; rr = D->rr ;
    sprank = rr [3] ;
    for (ns = 0, k = 0 ; k < nb ; k++)
    {
        ns += ((R [k+1] == R [k]+1) && (S [k+1] == S [k]+1)) ;
    }
    printf ("blocks: %d singletons: %d structural rank: %d\n", nb, ns, sprank) ;
    cs_dfree (D) ;


//    for (order = -1 ; order <= 2 ; order += 3)	/* natural and amd(A'*A) */
//    {
//        if (order == -1 && m > 1000) continue ;
//        printf ("QR   ") ;
//        print_order (order) ;
//        rhs (x, b, m) ;				/* compute right-hand-side */
//        t = tic () ;
//        ok = cs_qrsol (C, x, order) ;		/* min norm(Ax-b) with QR */
//        printf ("time: %8.2f ", toc (t)) ;
//
//        resid (ok, C, x, b, r) ;		/* print residual */
//    }


    if (m != n || sprank < n) return ;	/* return if rect. or singular*/
    for (order = -1 ; order <= 2 ; order++)	/* try all orderings */
    {
        if (order == -1 && m > 1000) continue ;
        printf ("LU   ") ;
        print_order (order) ;

        rhs (x, b, m) ;				/* compute right-hand-side */


        for(size_t idx=0; idx < 4; idx++){
            printf("b=%f\n", b[idx]);
        }

        t = tic () ;
        ok = cs_lusol (C, x, order, tol) ;	/* solve Ax=b with LU */

        printf("x\n");
        for(size_t idx=0; idx < 4; idx++){
            printf("x=%f\n", x[idx]);
        }

        printf ("time: %8.2f ", toc (t)) ;
        resid (ok, C, x, b, r) ;		/* print residual */
    }
    if (!Prob->sym) return ;
    for (order = -1 ; order <= 0 ; order++)	/* natural and amd(A+A') */
    {
        if (order == -1 && m > 1000) continue ;
        printf ("Chol ") ;
        print_order (order) ;
        rhs (x, b, m) ;				/* compute right-hand-side */
        t = tic () ;
        ok = cs_cholsol (C, x, order) ;		/* solve Ax=b with Cholesky */

        printf("x\n");
        for(size_t idx=0; idx < m; idx++){
            printf("%f\n", x[idx]);
        }

        printf ("time: %8.2f ", toc (t)) ;
        resid (ok, C, x, b, r) ;		/* print residual */
    }
    // CONVERT X BACK TO EIGEN NOW
}


void solve_c(double *A, double*b, double*x,  size_t num_rows, size_t num_cols)
{
    problem *Prob = load_csparse_format (A, b, num_rows, num_cols) ;
    solve_with_cs_chol (Prob) ;
    free_problem (Prob) ;
}



#ifdef __cplusplus
}
#endif
