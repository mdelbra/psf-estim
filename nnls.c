/*----------------------------------------------------------------------------

 "Point Spread Function Estimation from a Random Target"

 Copyright 2010-2011 mauricio delbracio (mdelbra@gmail.com)

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

 ----------------------------------------------------------------------------*/

/**
 * @file nnls.c
 * @brief library code with numerical algorithms for solving least squares
 *        and non-negative least squares
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 * @date Nov 24, 2011
 */


/*This filed was updated in version 1.2*/
/*
 Double precision functions were added, there were some problems
 due to the float precision, now, all *d functions convert the
 input A and b matrix to double matrix and do all the computations with doubles
 finally the result is truncated into a float to be compatible with the other
 part of the program
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nnls.h"
#include <cblas.h>
#include <float.h>

/** Buffer Size */
#define BUFFER_SIZE 113337

/** If the absolute value is less than EPS_ZERO consider it is zero */
#define EPS_ZERO 1e-6


/*Wrapper functions to use LAPACK*/
/*
 *SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
 $                   WORK, LWORK, IWORK, INFO )
 *
 *  -- LAPACK driver routine (version 3.2.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     June 2010
 *
 *     .. Scalar Arguments ..
 INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
 DOUBLE PRECISION   RCOND
 *     ..
 *     .. Array Arguments ..
 INTEGER            IWORK( * )
 DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
 *     ..
 *
 *  Purpose
 *  =======
 *
 *  DGELSD computes the minimum-norm solution to a real linear least
 *  squares problem:
 *      minimize 2-norm(| b - A*x |)
 *  using the singular value decomposition (SVD) of A. A is an M-by-N
 *  matrix which may be rank-deficient.
 *
 *  Several right hand side vectors b and solution vectors x can be
 *  handled in a single call; they are stored as the columns of the
 *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
 *  matrix X.
 *
 *  The problem is solved in three steps:
 *  (1) Reduce the coefficient matrix A to bidiagonal form with
 *      Householder transformations, reducing the original problem
 *      into a "bidiagonal least squares problem" (BLS)
 *  (2) Solve the BLS using a divide and conquer approach.
 *  (3) Apply back all the Householder tranformations to solve
 *      the original least squares problem.
 *
 *  The effective rank of A is determined by treating as zero those
 *  singular values which are less than RCOND times the largest singular
 *  value.
 *
 *  The divide and conquer algorithm makes very mild assumptions about
 *  floating point arithmetic. It will work on machines with a guard
 *  digit in add/subtract, or on those binary machines without guard
 *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
 *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
 *  without guard digits, but we know of none.
 *
 *  Arguments
 *  =========
 *
 *  M       (input) INTEGER
 *          The number of rows of A. M >= 0.
 *
 *  N       (input) INTEGER
 *          The number of columns of A. N >= 0.
 *
 *  NRHS    (input) INTEGER
 *          The number of right hand sides, i.e., the number of columns
 *          of the matrices B and X. NRHS >= 0.
 *
 *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
 *          On entry, the M-by-N matrix A.
 *          On exit, A has been destroyed.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,M).
 *
 *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
 *          On entry, the M-by-NRHS right hand side matrix B.
 *          On exit, B is overwritten by the N-by-NRHS solution
 *          matrix X.  If m >= n and RANK = n, the residual
 *          sum-of-squares for the solution in the i-th column is given
 *          by the sum of squares of elements n+1:m in that column.
 *
 *  LDB     (input) INTEGER
 *          The leading dimension of the array B. LDB >= max(1,max(M,N)).
 *
 *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
 *          The singular values of A in decreasing order.
 *          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
 *
 *  RCOND   (input) DOUBLE PRECISION
 *          RCOND is used to determine the effective rank of A.
 *          Singular values S(i) <= RCOND*S(1) are treated as zero.
 *          If RCOND < 0, machine precision is used instead.
 *
 *  RANK    (output) INTEGER
 *          The effective rank of A, i.e., the number of singular values
 *          which are greater than RCOND*S(1).
 *
 *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *
 *  LWORK   (input) INTEGER
 *          The dimension of the array WORK. LWORK must be at least 1.
 *          The exact minimum amount of workspace needed depends on M,
 *          N and NRHS. As long as LWORK is at least
 *              12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,
 *          if M is greater than or equal to N or
 *              12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,
 *          if M is less than N, the code will execute correctly.
 *          SMLSIZ is returned by ILAENV and is equal to the maximum
 *          size of the subproblems at the bottom of the computation
 *          tree (usually about 25), and
 *             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
 *          For good performance, LWORK should generally be larger.
 *
 *          If LWORK = -1, then a workspace query is assumed; the routine
 *          only calculates the optimal size of the WORK array, returns
 *          this value as the first entry of the WORK array, and no error
 *          message related to LWORK is issued by XERBLA.
 *
 *  IWORK   (workspace) INTEGER array, dimension (MAX(1,LIWORK))
 *          LIWORK >= max(1, 3 * MINMN * NLVL + 11 * MINMN),
 *          where MINMN = MIN( M,N ).
 *          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK.
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value.
 *          > 0:  the algorithm for computing the SVD failed to converge;
 *                if INFO = i, i off-diagonal elements of an intermediate
 *                bidiagonal form did not converge to zero.
 *
 *  Further Details
 *  ===============
 *
 *  Based on contributions by
 *     Ming Gu and Ren-Cang Li, Computer Science Division, University of
 *       California at Berkeley, USA
 *     Osni Marques, LBNL/NERSC, USA
 *
 *  =====================================================================
 */
static long sgelsd(int m, int n, int nrhs,
                   float *a, int lda,  float *b, int ldb,
                   float *s, float rcond, int *rank,
                   float *work, int lwork, int *iwork)
{
    extern void sgelsd_(const int *m, const int *n, const int *nrhs,
                        float *a, const int *lda, float *b,
                        const int *ldb, float *s, const float *rcond,
                        int *rank, float *work, int *lwork, int *iwork,
                        int *info);
    int info;
    sgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank,
            work, &lwork, iwork, &info);
    return info;
}

static long dgelsd(int m, int n, int nrhs,
                   double *a, int lda,  double *b, int ldb,
                   double *s, double rcond, int *rank,
                   double *work, int lwork, int *iwork)
{
    extern void dgelsd_(const int *m, const int *n, const int *nrhs,
                double *a, const int *lda, double *b,
                const int *ldb, double *s, const double *rcond,
                int *rank, double *work, int *lwork, int *iwork,
            int *info);
    int info;
    dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank,
            work, &lwork, iwork, &info);
    return info;
}

/**
 * @brief Solve Least Squares problem x such that Ax = b.
 * @param A  - Array cointaining matrix 'A' elements (column major)
 * @param b  - Array of observed values 'b'
 * @param m  - number of rows of 'A'
 * @param n  - number of columns of 'A'
 * @return Array 'x' with the solution
 */
float *solve_ls(float *A, float *b, int m, int n)
{
    float *work;
    float *s;
    float *x;
    int i ;

    float lwork;
    int iwork[BUFFER_SIZE];
    int rank;

    /*Least Squares */
    s = (float *) malloc(n * sizeof(float));

    /*Do a query to know the optimum BufferSize */
    sgelsd(m, n, 1, A, m, b, m, s, EPS_ZERO, &rank, &lwork, -1, iwork);

    work = (float *) malloc(lwork * sizeof(float));
    sgelsd(m, n, 1, A, m, b, m, s, EPS_ZERO, &rank, work, (int) lwork,
                  iwork);

    free((void *) s);
    free((void *) work);

    x  = (float *) malloc(n * sizeof(float));
    for(i=0;i<n;i++)
        x[i] = b[i];

    return x;

}


/**
 * @brief Solve Least Squares problem x such that Ax = b and then
 *        thresholds the solution with 'th' such that x<=th is 0.
 * @param A  - Array cointaining matrix 'A' elements (column major)
 * @param b  - Array of observed values 'b'
 * @param m  - number of rows of 'A'
 * @param n  - number of columns of 'A'
 * @param th - final threshold.
 * @return Array 'x' with the solution.
 */
float *solve_ls_th(float *A, float *b, int m, int n, float th)
{

    float *x;
    int i;

    x = solve_ls(A, b, m, n);

    for (i = 0; i < n; i++)
        x[i] = (x[i] > th) ? x[i] : 0;

    return x;

}


/**
 * @brief Solve Non-Negative Least Squares problem x such that Ax = b, x>=0
 * @details (Reference) Portugal, Judice and Vicente, A comparison of block
 pivoting and interior point algorithms for linear least squares
 problems with nonnegative variables, Mathematics of Computation,
 63(1994), pp. 625-643
 * @param A  - Array cointaining matrix 'A' elements (column major)
 * @param b  - Array of observed values 'b'
 * @param m  - number of rows of 'A'
 * @param n  - number of columns of 'A'
 * @return Array 'x' with the solution
 */
float *solve_nnls(float *A, float *b, int m, int n)
{

    /*Input:
     A:      [mxn] matrix
     b:      [mx1] vector
     Output
     x:      solution
     */



    float *xk, *yk, *uk, *vk, *C, *d, *AtA, *Atb, *error;
    int k, i, j;
    float tol1 = 1e-4;
    float tol2 = 1e-1;

    float mk, theta1, theta2, theta;
    float *work;
    float *s;
    float lwork;
    int iwork[BUFFER_SIZE];
    int rank;

    float *Aaux;
    float *baux, *xini;

    char ready = 0;
    char theta1_empty;
    char theta2_empty;

    int niter = 30;

    xk = (float *) malloc(n * sizeof(float));
    yk = (float *) malloc(n * sizeof(float));
    uk = (float *) malloc(n * sizeof(float));
    vk = (float *) malloc(n * sizeof(float));

    C = (float *) malloc((m + n) * n * sizeof(float));

    d = (float *) malloc((m + n) * sizeof(float));
    AtA = (float *) malloc(n * n * sizeof(float));
    Atb = (float *) malloc(n * sizeof(float));
    error = (float *) malloc(n * sizeof(float));

    Aaux = (float *) malloc(m * n * sizeof(float));
    baux = (float *) malloc(m * sizeof(float));

    /*for the Least Squares */
    s = (float*)malloc(n*sizeof(float));

    /*AtA = A'*A;  C:= alpha*A*A' + beta*C, */
    cblas_sgemm(CblasColMajor,CblasTrans, CblasNoTrans, n , n ,
                m , 1.0 , A , m , A , m, 0.0, AtA, n);

    /*Atb = A'*b */
    cblas_sgemv(CblasColMajor, CblasTrans, m, n, 1.0, A, m,
                b, 1, 0.0, Atb, 1);


    /*Step 0 - Initialization */
    k = 0;
    ready = 0;

    /*I am going to initialize with the result of the ls squares */
    for (i = 0; i < m * n; i++)
        Aaux[i] = A[i];

    for (i = 0; i < m; i++)
        baux[i] = b[i];

    xini = solve_ls(Aaux, baux, m, n);

    for (i = 0; i < n; i++)
        xini[i] = xini[i] <= 0.00001 ? 0.00001 : xini[i];


    for (i = 0; i < n; i++)
    {
        xk[i] = xini[i];    /*old 1 */
        yk[i] = 1;
    }


    while (!ready)
    {
        k++;

        printf("    ->Iteration k= %d\n", k);
        /*Step 1 */
        mk = cblas_sdot(n, xk, 1, yk, 1) / (n * n);


        /*Generate the fixed submatrix  of C =  C = [A; Xk_sqrt_inv*Yk_sqrt];*/
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                C[i + j * (m + n)] = A[i + j * m];

        /*Update the bottom part of matrix C
         * with Xk_sqrt_inv*Yk_sqrt diagonal matrix */
        /* C = [A; Xk_sqrt_inv*Yk_sqrt]; */
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                if (i == j)
                    C[(i + m) + j * (m + n)] = sqrt(yk[i] / xk[i]);
                else
                    C[(i + m) + j * (m + n)] = 0;
            }
        /*d = [b-A*Xk*e; Xk_sqrt_inv*Yk_sqrt_inv*mk*e ]; */
        /*top part of dtop = b - Ax */
        /*initialize d = -b */
        for (i = 0; i < m; i++)
            d[i] = b[i];

        /*d = b - Ax; */
        cblas_sgemv(CblasColMajor, CblasNoTrans, m, n, -1.0, A, m,
                    xk, 1, 1.0, d, 1);


        /*fill in the bottom part of d
         *  Xk_sqrt_inv*Yk_sqrt_inv*mk*e */
        for (i = 0; i < n; i++)
            d[i + m] = mk / sqrt(yk[i] * xk[i]);


        /*Find uk : C*uk = d */

        /*If it is the first iteration I find the optimum
         * buffer size in order to speed up
         */

        if (k == 1)
        {
            /*Do a query to know the optimum BufferSize */
            sgelsd(m+n, n, 1, C, m+n, d, m+n, s, EPS_ZERO,
                          &rank, &lwork, -1, iwork);
            work = (float *) malloc(lwork * sizeof(float));
        }

        /*Find uk : C*uk = d */
        sgelsd(m+n, n, 1, C, m+n, d, m+n, s, EPS_ZERO, &rank, work,
                      (int)lwork, iwork);

        for (i = 0; i < n; i++)
            uk[i] = d[i];


        /*Calculate  vk = -Yk*e + Xk_inv*mk*e - Xk_inv*Yk*uk; */
        for (i = 0; i < n; i++)
            vk[i] = -yk[i] + mk / xk[i] - yk[i] * uk[i] / xk[i];


        /*Step2 */

        /*Finding theta1 -  T1 = min(-xk(uk < 0) ./ uk(uk < 0)); */
        theta1 = DBL_MAX;
        theta1_empty = 1;
        for (i = 0; i < n; i++)
            if ((uk[i] < 0) && (-xk[i] / uk[i] < theta1))
            {
                theta1 = -xk[i] / uk[i];
                theta1_empty = 0;
            }

        /*Finding theta2 -   T2 = min(-yk(vk < 0) ./ vk(vk < 0)); */
        theta2 = DBL_MAX;
        theta2_empty = 1;
        for (i = 0; i < n; i++)
            if ((vk[i] < 0) && (-yk[i] / vk[i] < theta2))
            {
                theta2 = -yk[i] / vk[i];
                theta2_empty = 0;
            }

        /*theta = 0.99995 * min(T1,T2) */
        theta = (theta1 < theta2) ? theta1 : theta2;
        theta = 0.99995 * theta;

        if (theta1_empty || theta2_empty)
        {
            theta = 0;
            ready = 1;
        }

        /* xk = xk + theta*uk;
         * yk = yk + theta*vk;
         */


        for (i = 0; i < n; i++)
        {
            xk[i] = xk[i] + theta * uk[i];
            yk[i] = yk[i] + theta * vk[i];
        }

        /*Step 3 */


        /*Computing error = AtAx - Atb -yk */
        /*initialize error = -Atb */
        for (i = 0; i < n; i++)
            error[i] = Atb[i];


        /*error = AtAx - Atb */
        cblas_ssymv(CblasColMajor, CblasUpper, n, 1.0, AtA, n, xk, 1,
                    -1.0, error, 1);

        /*error = AtAx - Atb - yk */
        for (i = 0; i < n; i++)
            error[i] = error[i] - yk[i];

        if ((cblas_sdot(n, xk, 1, yk, 1) < tol1) &&
            (cblas_snrm2(n, error, 1) < tol2))
            ready = 1;
        else if (k == niter) ready = 1;
    }

    /*Cleaning */
    free((void *) yk);
    free((void *) uk);
    free((void *) vk);
    free((void *) C);
    free((void *) d);
    free((void *) AtA);
    free((void *) Atb);
    free((void *) error);
    free((void *) Aaux);
    free((void *) baux);
    free((void *) work);
    free((void *) xini);
    free((void *) s);

    return xk;
    /*Non-negative Least Squares */
}


/* Added in Version 1.2 */
/* Single precision wasn't enough...*/

/**
 * @brief Solve Least Squares problem (double precision) x such that Ax = b.
 * @param A  - Array cointaining matrix 'A' elements (column major)
 * @param b  - Array of observed values 'b'
 * @param m  - number of rows of 'A'
 * @param n  - number of columns of 'A'
 * @return Array 'x' with the solution
 */
float *solve_lsd(float *Af, float *bf, int m, int n)
{
    double *work;
    double *s;
    float *x;
    int i ;

    double lwork;
    int iwork[BUFFER_SIZE];
    int rank, info;

    double *A, *b;

    A = (double *) malloc(m*n*sizeof(double));
    b = (double *) malloc(m*sizeof(double));


    for(i=0;i<m*n;i++)
        A[i] = (double) Af[i];

    for(i=0;i<m;i++)
        b[i] = (double) bf[i];


    /*Least Squares */
    s = (double *) malloc(n * sizeof(double));

    /*Do a query to know the optimum BufferSize */
    info =
    dgelsd(m, n, 1, A, m, b, m, s, EPS_ZERO, &rank, &lwork, -1, iwork);

    work = (double *) malloc(lwork * sizeof(double));
    info =
    dgelsd(m, n, 1, A, m, b, m, s, EPS_ZERO, &rank, work, (int) lwork,
           iwork);

    free((void *) s);
    free((void *) work);
    free((void *) A);
    free((void *) b);

    x  = (float *) malloc(n * sizeof(float));
    for(i=0;i<n;i++)
        x[i] = (float) b[i];

    return x;

}



/**
 * @brief Solve Least Squares problem (double precision) x such that Ax = b
 *        and then thresholds the solution with 'th' such that x<=th is 0.
 * @param A  - Array cointaining matrix 'A' elements (column major)
 * @param b  - Array of observed values 'b'
 * @param m  - number of rows of 'A'
 * @param n  - number of columns of 'A'
 * @param th - final threshold.
 * @return Array 'x' with the solution.
 */
float *solve_lsd_th(float *A, float *b, int m, int n, float th)
{

    float *x;
    int i;

    x = solve_lsd(A, b, m, n);

    for (i = 0; i < n; i++)
        x[i] = (x[i] > th) ? x[i] : 0;

    return x;

}


/**
 * @brief Solve Non-Negative Least Squares problem (double precision) x such
 * that Ax = b, x>=0.
 * @details (Reference) Portugal, Judice and Vicente, A comparison of block
 pivoting and interior point algorithms for linear least squares
 problems with nonnegative variables, Mathematics of Computation,
 63(1994), pp. 625-643
 * @param A  - Array cointaining matrix 'A' elements (column major)
 * @param b  - Array of observed values 'b'
 * @param m  - number of rows of 'A'
 * @param n  - number of columns of 'A'
 * @return Array 'x' with the solution
 */
float *solve_nnlsd(float *A, float *b, int m, int n)
{

    /*Input:
     A:      [mxn] matrix
     b:      [mx1] vector
     Output
     x:      solution
     */


    float *x;
    double *xk, *yk, *uk, *vk, *C, *d, *AtA, *Atb, *error;
    int k, i, j;
    double tol1 = 1e-4;
    double tol2 = 1e-1;

    double mk, theta1, theta2, theta;
    double *work;
    double *s;
    double lwork;
    int iwork[BUFFER_SIZE];
    int info;
    int rank;

    float *Aaux;
    float *baux, *xini;

    double *Ad;
    double *bd;

    char ready = 0;
    char theta1_empty;
    char theta2_empty;

    int niter = 30;

    xk = (double *) malloc(n * sizeof(double));
    yk = (double *) malloc(n * sizeof(double));
    uk = (double *) malloc(n * sizeof(double));
    vk = (double *) malloc(n * sizeof(double));

    C = (double *) malloc((m + n) * n * sizeof(double));

    d = (double *) malloc((m + n) * sizeof(double));
    AtA = (double *) malloc(n * n * sizeof(double));
    Atb = (double *) malloc(n * sizeof(double));
    error = (double *) malloc(n * sizeof(double));

    Aaux = (float *) malloc(m * n * sizeof(float));
    baux = (float *) malloc(m * sizeof(float));

    Ad = (double *) malloc(m * n * sizeof(double));
    bd = (double *) malloc(m * sizeof(double));



    for (i = 0; i < m * n; i++)
        Ad[i] = (double) A[i];

    for (i = 0; i < m; i++)
        bd[i] = (double) b[i];


    /*for the Least Squares */
    s = (double*)malloc(n*sizeof(double));

    /*AtA = A'*A;  C:= alpha*A*A' + beta*C, */
    cblas_dgemm(CblasColMajor,CblasTrans, CblasNoTrans, n , n ,
                m , 1.0 , Ad , m , Ad , m, 0.0, AtA, n);

    /*Atb = A'*b */
    cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, Ad, m,
                bd, 1, 0.0, Atb, 1);


    /*Step 0 - Initialization */
    k = 0;
    ready = 0;

    /*I am going to initialize with the result of the ls squares */
    for (i = 0; i < m * n; i++)
        Aaux[i] = A[i];

    for (i = 0; i < m; i++)
        baux[i] = b[i];

    xini = solve_lsd(Aaux, baux, m, n);

    for (i = 0; i < n; i++)
        xini[i] = xini[i] <= 0.00001 ? 0.00001 : xini[i];


    for (i = 0; i < n; i++)
    {
        xk[i] = (double) xini[i];    /*old 1 */
        yk[i] = 1;
    }


    while (!ready)
    {
        k++;

        printf("    ->Iteration k= %d\n", k);
        /*Step 1 */
        mk = cblas_ddot(n, xk, 1, yk, 1) / (n * n);


        /*Generate the fixed submatrix  of C =  C = [A; Xk_sqrt_inv*Yk_sqrt];*/
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                C[i + j * (m + n)] =  Ad[i + j * m];

        /*Update the bottom part of matrix C
         * with Xk_sqrt_inv*Yk_sqrt diagonal matrix */
        /* C = [A; Xk_sqrt_inv*Yk_sqrt]; */
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                if (i == j)
                    C[(i + m) + j * (m + n)] = sqrt(yk[i] / xk[i]);
                else
                    C[(i + m) + j * (m + n)] = 0;
            }
        /*d = [b-A*Xk*e; Xk_sqrt_inv*Yk_sqrt_inv*mk*e ]; */
        /*top part of dtop = b - Ax */
        /*initialize d = -b */
        for (i = 0; i < m; i++)
            d[i] = b[i];

        /*d = b - Ax; */
        cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, -1.0, Ad, m,
                    xk, 1, 1.0, d, 1);


        /*fill in the bottom part of d
         *  Xk_sqrt_inv*Yk_sqrt_inv*mk*e */
        for (i = 0; i < n; i++)
            d[i + m] = mk / sqrt(yk[i] * xk[i]);


        /*Find uk : C*uk = d */

        /*If it is the first iteration I find the optimum
         * buffer size in order to speed up
         */

        if (k == 1)
        {
            /*Do a query to know the optimum BufferSize */
            info = dgelsd(m+n, n, 1, C, m+n, d, m+n, s, EPS_ZERO,
                          &rank, &lwork, -1, iwork);
            work = (double *) malloc(lwork * sizeof(double));
        }

        /*Find uk : C*uk = d */
        info = dgelsd(m+n, n, 1, C, m+n, d, m+n, s, EPS_ZERO, &rank, work,
                      (int)lwork, iwork);

        for (i = 0; i < n; i++)
            uk[i] = d[i];


        /*Calculate  vk = -Yk*e + Xk_inv*mk*e - Xk_inv*Yk*uk; */
        for (i = 0; i < n; i++)
            vk[i] = -yk[i] + mk / xk[i] - yk[i] * uk[i] / xk[i];


        /*Step2 */

        /*Finding theta1 -  T1 = min(-xk(uk < 0) ./ uk(uk < 0)); */
        theta1 = DBL_MAX;
        theta1_empty = 1;
        for (i = 0; i < n; i++)
            if ((uk[i] < 0) && (-xk[i] / uk[i] < theta1))
            {
                theta1 = -xk[i] / uk[i];
                theta1_empty = 0;
            }

        /*Finding theta2 -   T2 = min(-yk(vk < 0) ./ vk(vk < 0)); */
        theta2 = DBL_MAX;
        theta2_empty = 1;
        for (i = 0; i < n; i++)
            if ((vk[i] < 0) && (-yk[i] / vk[i] < theta2))
            {
                theta2 = -yk[i] / vk[i];
                theta2_empty = 0;
            }


        /*theta = 0.99995 * min(T1,T2) */
        theta = (theta1 < theta2) ? theta1 : theta2;
        theta = 0.99995 * theta;

        if (theta1_empty || theta2_empty)
        {
            theta = 0;
            ready = 1;
        }

        /* xk = xk + theta*uk;
         * yk = yk + theta*vk;
         */


        for (i = 0; i < n; i++)
        {
            xk[i] = xk[i] + theta * uk[i];
            yk[i] = yk[i] + theta * vk[i];
        }

        /*Step 3 */


        /*Computing error = AtAx - Atb -yk */
        /*initialize error = -Atb */
        for (i = 0; i < n; i++)
            error[i] = Atb[i];


        /*error = AtAx - Atb */
        cblas_dsymv(CblasColMajor, CblasUpper, n, 1.0, AtA, n, xk, 1,
                    -1.0, error, 1);

        /*error = AtAx - Atb - yk */
        for (i = 0; i < n; i++)
            error[i] = error[i] - yk[i];

        if ((cblas_ddot(n, xk, 1, yk, 1) < tol1) &&
            (cblas_dnrm2(n, error, 1) < tol2))
            ready = 1;
        else if (k == niter) ready = 1;
    }


    x  = (float *) malloc(n * sizeof(float));
    for(i=0;i<n;i++)
        x[i] = (float) xk[i];


    /*Cleaning */
    free((void *) yk);
    free((void *) uk);
    free((void *) vk);
    free((void *) C);
    free((void *) d);
    free((void *) AtA);
    free((void *) Atb);
    free((void *) error);
    free((void *) Aaux);
    free((void *) baux);
    free((void *) Ad);
    free((void *) bd);
    free((void *) work);
    free((void *) xini);
    free((void *) s);
    free((void *) xk);

    return x;

    /*Non-negative Least Squares */
}
