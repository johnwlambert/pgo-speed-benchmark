//
// Created by John Lambert-Admin on 4/22/19.
//

#pragma once


#include <string>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

//#include "csparse/csparse.h"
//#include "csparse_utils.h"
#include "pose_graph.h" // to get MatrixXd, VectorXd

class SparseLinSolver{

    public:
        SparseLinSolver(MatrixXd & A, VectorXd & b, VectorXd & x) : A_(A), b_(b), x_(x)
        {}

        void solve_eigen_w_csparse(std::string solver_type)
        {
//            problem *Prob = load_csparse_format () ;
//            solve_with_cs_chol (Prob, solver_type) ;
//            free_problem (Prob) ;
        }

private:
    MatrixXd A_;
    VectorXd b_;
    VectorXd x_;

//    /* CS_LOAD loads a triplet matrix from Eigen. */
//    cs *eigen_to_csparse()
//    {
//        int i, j ;
//        double x ;
//        cs *T ;
//        T = cs_spalloc (0, 0, 1, 1, 1) ;
//        for (i = 0; i < A_.rows(); i++)
//        {
//            for (i = 0; i < A_.cols() ; i++)
//            {
//                if (A_(i,j) != 0.0)
//                {
//                    if (!cs_entry (T, i, j, A_(i,j) )) return (cs_spfree (T)) ;
//                }
//            }
//        }
//        return (T) ;
//    }
//
//
//    problem *load_csparse_format()
//    {
//        cs *T, *A, *C ;
//        int sym, m, n, mn, nz1, nz2 ;
//        problem *Prob ;
//        Prob = (problem *) cs_calloc (1, sizeof (problem)) ;
//        if (!Prob) return (NULL) ;
//        T = eigen_to_csparse() ;			/* load triplet matrix T from Eigen "A_" */
//        Prob->A = A = cs_triplet (T) ;	/* A = compressed-column form of T */
//        cs_spfree (T) ;			/* clear T */
//        if (!cs_dupl (A)) return (free_problem (Prob)) ; /* sum up duplicates */
//        Prob->sym = sym = is_sym (A) ;	/* determine if A is symmetric */
//        m = A->m ; n = A->n ;
//        mn = CS_MAX (m,n) ;
//        nz1 = A->p [n] ;
//        cs_dropzeros (A) ;			/* drop zero entries */
//        nz2 = A->p [n] ;
//        Prob->C = C = sym ? make_sym (A) : A ;  /* C = A + triu(A,1)', or C=A */
//        if (!C) return (free_problem (Prob)) ;
//        printf ("\n--- Matrix: %d-by-%d, nnz: %d (sym: %d: nnz %d), norm: %8.2e\n",
//                m, n, A->p [n], sym, sym ? C->p [n] : 0, cs_norm (C)) ;
//        if (nz1 != nz2) printf ("zero entries dropped: %d\n", nz1 - nz2) ;
//        if (nz2 != A->p [n]) printf ("tiny entries dropped: %d\n", nz2 - A->p [n]) ;
//        Prob->b = (double *) cs_malloc (mn, sizeof (double)) ;
//
//        Prob->x =  (double *) cs_malloc (mn, sizeof (double)) ;
//        Prob->r =  (double *) cs_malloc (mn, sizeof (double)) ;
//        return ((!Prob->b || !Prob->x || !Prob->r) ? free_problem (Prob) : Prob) ;
//    }
//
//
//    /* solve a linear system using Cholesky, LU, and QR, with various orderings */
//    void solve_with_cs_chol (problem *Prob, std::string solver_type)
//    {
//        cs *A, *C ;
//        double *b, *x, *r,  t, tol ;
//        int k, m, n, ok, order, nb, ns, *R, *S, *rr, sprank ;
//        csd *D ;
//        if (!Prob) return ;
//        A = Prob->A ; C = Prob->C ; b = Prob->b ; x = Prob->x ; r = Prob->r ;
//        m = A->m ; n = A->n ;
//        tol = Prob->sym ? 0.001 : 1 ;		/* partial pivoting tolerance */
//        D = cs_dmperm (C) ;				/* dmperm analysis */
//        if (!D) return ;
//        nb = D->nb ; R = D->R ; S = D->S ; rr = D->rr ;
//        sprank = rr [3] ;
//        for (ns = 0, k = 0 ; k < nb ; k++)
//        {
//            ns += ((R [k+1] == R [k]+1) && (S [k+1] == S [k]+1)) ;
//        }
//        printf ("blocks: %d singletons: %d structural rank: %d\n", nb, ns, sprank) ;
//        cs_dfree (D) ;
//        for (order = -1 ; order <= 2 ; order += 3)	/* natural and amd(A'*A) */
//        {
//            if (order == -1 && m > 1000) continue ;
//            printf ("QR   ") ;
//            print_order (order) ;
//            rhs (x, b, m) ;				/* compute right-hand-side */
//            t = tic () ;
//            ok = cs_qrsol (C, x, order) ;		/* min norm(Ax-b) with QR */
//            printf ("time: %8.2f ", toc (t)) ;
//
//            resid (ok, C, x, b, r) ;		/* print residual */
//        }
//        if (m != n || sprank < n) return ;	/* return if rect. or singular*/
//        for (order = -1 ; order <= 2 ; order++)	/* try all orderings */
//        {
//            if (order == -1 && m > 1000) continue ;
//            printf ("LU   ") ;
//            print_order (order) ;
//
//            rhs (x, b, m) ;				/* compute right-hand-side */
//
//
//            for(size_t idx=0; idx < 4; idx++){
//                printf("b=%f\n", b[idx]);
//            }
//
//            t = tic () ;
//            ok = cs_lusol (C, x, order, tol) ;	/* solve Ax=b with LU */
//
//            printf("x\n");
//            for(size_t idx=0; idx < 4; idx++){
//                printf("x=%f\n", x[idx]);
//            }
//
//            printf ("time: %8.2f ", toc (t)) ;
//            resid (ok, C, x, b, r) ;		/* print residual */
//        }
//        if (!Prob->sym) return ;
//        for (order = -1 ; order <= 0 ; order++)	/* natural and amd(A+A') */
//        {
//            if (order == -1 && m > 1000) continue ;
//            printf ("Chol ") ;
//            print_order (order) ;
//            rhs (x, b, m) ;				/* compute right-hand-side */
//            t = tic () ;
//            ok = cs_cholsol (C, x, order) ;		/* solve Ax=b with Cholesky */
//
//            printf("x\n");
//            for(size_t idx=0; idx < m; idx++){
//                printf("%f\n", x[idx]);
//            }
//
//            printf ("time: %8.2f ", toc (t)) ;
//            resid (ok, C, x, b, r) ;		/* print residual */
//        }
//        // CONVERT X BACK TO EIGEN NOW
//    }
};


