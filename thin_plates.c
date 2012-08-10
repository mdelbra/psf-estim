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

/*Version 1.2 24 November 2011*/


/**
 * @file thin_plates.c
 * @brief library code to estimate/evaluate thin plates splines.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include "thin_plates.h"


/*---------------------------------------------------------------------------*/
/* LAPACK Wrapping functions */


/*
* subroutine DGEQRF    (    INTEGER     M,
* INTEGER     N,
* DOUBLE PRECISION,dimension( lda, * )     A,
* INTEGER     LDA,
* DOUBLE PRECISION,dimension( * )     TAU,
* DOUBLE PRECISION,dimension( * )     WORK,
* INTEGER     LWORK,
* INTEGER     INFO
* )
*
*
*
*
*  Purpose
*  =======
*
*  DGEQRF computes a QR factorization of a real M-by-N matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of min(m,n) elementary reflectors (see Further
*          Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is
*          the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*
* A two-dimensional Fortran array declared as
*
*    DOUBLE PRECISION A(LDA, N)
*
*  is a contiguous piece of LDA X N float-words of memory, stored in
*  column-major order: elements in a column are contiguous, and elements
*  within a row are separated by a stride of LDA float-words.
*/
static long sgeqrf(long m, long n, float *a, long lda, float *tau,
                   float *work, long lwork)
{
    extern void sgeqrf_(const long *m, const long *n, float *a,
                        const long *lda, float *tau, float *work,
                        const long *lwork, long *info);
    long info;
    sgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}



/*
* subroutine to re-build matrix Q from Lapack decomposition.
*SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
     INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
     DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORGQR generates an M-by-N real matrix Q with orthonormal columns,
*  which is defined as the first N columns of a product of K elementary
*  reflectors of order M
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by DGEQRF in the first k columns of its array
*          argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*/

static long sorgqr(long m, long n, long k, float *a, long lda, float *tau,
                   float *work, long lwork)
{
    extern void sorgqr_(const long *m, const long *n, const long *k,
                        float *a, const long *lda, float *tau, float *work,
                        const long *lwork, long *info);
    long info;
    sorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}



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

static long sgelsd(long m, long n, long nrhs,
                   float *a, long lda, float *b, long ldb,
                   float *s, float rcond, long *rank,
                   float *work, long lwork, long *iwork)
{
    extern void sgelsd_(const long *m, const long *n, const long *nrhs,
                        float *a, const long *lda, float *b,
                        const long *ldb, float *s, const float *rcond,
                        long *rank, float *work, long *lwork, long *iwork,
                        long *info);
    long info;
    sgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank,
            work, &lwork, iwork, &info);
    return info;
}


/** @brief Error/Exit print a message and exit.
 *  @param msg
 */
static void error(char *msg)
{
    fprintf(stderr, "ThinPlate Error: %s\n", msg);
    exit(EXIT_FAILURE);
}


/**
 * @brief Free memory used in ThinPlate 'tp'
 * @param tp
 */
void free_thinPlate(ThinPlate tp)
{
    if (tp == NULL || tp->xc == NULL || tp->yc == NULL
        || tp->coef_x == NULL || tp->coef_y == NULL)
        error("free_pointList: invalid pointList input.");

    free((void *) tp->xc);
    free((void *) tp->yc);
    free((void *) tp->coef_x);
    free((void *) tp->coef_y);
    free((void *) tp);
}


/**
 * @brief Create new ThinPlate for 'nc' points
 * @param nc - number of points
 * @return created ThinPlate
 */
static ThinPlate new_thinPlate(int nc)
{

    ThinPlate tp;

    if (nc < 3)
        error("new_pointList: 'nc' must be at least three.");

    tp = (ThinPlate) malloc(sizeof(struct thinPlateStruct));
    if (tp == NULL)
        error("not enough memory.");
    tp->nc = nc;
    tp->xc = (float *) malloc(tp->nc * sizeof(float));
    tp->yc = (float *) malloc(tp->nc * sizeof(float));
    tp->coef_x = (float *) malloc(tp->nc * sizeof(float));
    tp->coef_y = (float *) malloc(tp->nc * sizeof(float));
    if (tp->xc == NULL || tp->yc == NULL || tp->coef_x == NULL
        || tp->coef_y == NULL)
        error("not enough memory.");

    return tp;
}



/**
 * @brief Calculate a new ThinPlate between 'Pin' and 'Pout' points and
 *        regularization parameter 'lambda'
 * @param Pin - Array of input 2D points
 * @param Pout - Array of output 2D points
 * @param k   - number of points
 * @param lambda - Regularization parameter in [0,Inf)
                   0 - no regularization
                   Inf - affine transform between points
 * @return created ThinPlate
 */
ThinPlate calculate_thinPlate(float *Pin, float *Pout, int k, float lambda)
{

    /* Pin, Pout : 2k float arrays. x[0], x[1],...y[0],...
     * k : number of points
     * lambda : regularization parameter - 0 no reg.
     * */

    int i, j, info;
    float work[BUFFER_SIZE];
    float tau[3];
    float *A, *s, *R, *Phi, *Q1, *Q2, *X, *Y;
    float *aux, *aux2, *aux3, *d, *C;

    long iwork[BUFFER_SIZE];
    long rank = 0;


    float r;

    ThinPlate tp;

    /* X is of size k x k in order to house full matrix Q from QR decomp. */
    X = (float *) malloc(k * k * sizeof(float));
    Y = (float *) malloc(D * k * sizeof(float));
    A = (float *) malloc(D * k * sizeof(float));
    R = (float *) malloc(D * D * sizeof(float));
    Phi = (float *) malloc(k * k * sizeof(float));


    aux = (float *) malloc((k - D) * D * sizeof(float));
    aux2 = (float *) malloc(k * (k - D) * sizeof(float));

    /*Initizalize in 0 */
    aux3 = (float *) calloc((k - D) * (k - D), sizeof(float));

    d = (float *) malloc(D * D * sizeof(float));
    C = (float *) malloc(k * D * sizeof(float));

    s = (float *) malloc((k - D) * sizeof(float));


    /*
     * X data - input data
     * Y data - output data
     * k - number of points
     * D - dimension 2d, D = 2;
     */

    /*
     * Recall: LAPACK is column major.
     */

    for (i = 0; i < k; i++)
    {
        /* The first column is the x coord */
        X[i] = Pin[2 * i];
        Y[i] = Pout[2 * i];

        /* The second column is the y coord */
        X[i + k] = Pin[2 * i + 1];
        Y[i + k] = Pout[2 * i + 1];

        /* The third column is always one - homogeneaous coord */
        X[i + 2 * k] = 1.0;
        Y[i + 2 * k] = 1.0;
    }

    /* Fill in the Kernel matrix */
    for (i = 0; i < k; i++)
        for (j = 0; j < k; j++)
        {
            r = (Pin[2 * i] - Pin[2 * j]) * (Pin[2 * i] - Pin[2 * j])
                + (Pin[2 * i + 1] - Pin[2 * j + 1]) * (Pin[2 * i + 1]
                - Pin[2 * j + 1]);
            if (r > EPS_ZERO)
                Phi[i + j * k] = r * log(r);
            else
                Phi[i + j * k] = 0;
        }

    /* QR Factorization of matrix X, size(X) = k x 3  */

    /*Use a very big Buffer of size lwork = BUFFER_SIZE
     * it works but is not optimal*/
    info = sgeqrf(k, D, X, k, tau, work, BUFFER_SIZE);

    /* Get R matrix - is the top left submatrix saved in X */
    for (i = 0; i < D; i++)
        for (j = 0; j < D; j++)
        {
            if (j <= i)
                R[j + D * i] = X[j + k * i];
            else
                R[j + D * i] = 0;
        }


    /* Get Q matrix - Q matrix will be saved in X */

    /* Use a very big Buffer of size lwork = BUFFER_SIZE
     * it works but is not optimal.
     */
    info = sorgqr(k, k, D, X, k, tau, work, BUFFER_SIZE);


    /* Get Q1 and Q2 sub-matrices of Q
     * Q1 - first D columns, Q2 - from column D+1 to the end.
     */
    Q1 = X;
    Q2 = X + D * k;

    /*
     * c = Q2/(Q2'*Phi*Q2 + lambda*eye(k-D-1))*Q2'*Y;
     * d = R\(Q1'*(Y - Phi*c));
     */

    /*
       Q2 - k x (k-D)
       Q1 - k x D
       R  - D x D
       Y  - k x D
       X  - k x D
     */

    /*
     * c = Q2/(Q2'*Phi*Q2 + lambda*eye(k-D-1))*Q2'*Y;
     * d = R\(Q1'*(Y - Phi*c));
     */

    /* aux2 = PhiQ2 */
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                k, k - D, k, 1, Phi, k, Q2, k, 0.0, aux2, k);

    /*Initialize aux3 = lambda*I */
    for (i = 0; i < k - D; i++)
        aux3[i + (k - D) * i] = lambda;


    /* aux3 = Q2^t Phi Q2 + lambda I */
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                k - D, k - D, k, 1, Q2, k, aux2, k, 1.0, aux3, k - D);

    /* aux =  Q2^t Y   --  (k-D) x D */
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                k - D, D, k, 1, Q2, k, Y, k, 0.0, aux, k - D);

    /* aux =  inv(aux3)*aux   --  (k-D) x D */
    info =
    sgelsd(k - D, k - D, D, aux3, k - D, aux, k - D, s, EPS_ZERO,
           &rank, work, BUFFER_SIZE, iwork);
    if (info != 0)
        fprintf(stderr, "failure with error %d\n", info);

    /* C  = Q2 inv(Q2^tPhi Q2 + lambda eye(k-D-1)) Q2^t Y */
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                k, D, k - D, 1, Q2, k, aux, k - D, 0.0, C, k);

    /* d = R\(Q1'*(Y - Phi*c))   */

    /* Y = Y - Phi*c) */
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                k, D, k, -1.0, Phi, k, C, k, 1.0, Y, k);

    /*d =  Q1^t Y */
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                D, D, k, 1, Q1, k, Y, k, 0.0, d, D);
    /* d =  inv(R)*d   --  D x D */
    info = sgelsd(D, D, D, R, D, d, D, s, EPS_ZERO, &rank,
                  work, BUFFER_SIZE, iwork);
    if (info != 0)
        fprintf(stderr, "failure with error %d\n", info);

    /* Save the information in a ThinPlate Structure */
    tp = new_thinPlate(k);
    for (i = 0; i < k; i++)
    {
        tp->xc[i] = Pin[2 * i];
        tp->yc[i] = Pin[2 * i + 1];
        tp->coef_x[i] = C[i];
        tp->coef_y[i] = C[i + k];
        tp->lambda = lambda;
    }
    tp->affine[0] = d[0];
    tp->affine[1] = d[1];
    tp->affine[2] = d[2];
    tp->affine[3] = d[3];
    tp->affine[4] = d[4];
    tp->affine[5] = d[5];



    /* Cleaning the house... */
    free((void *) X);
    free((void *) Y);
    free((void *) A);
    free((void *) R);
    free((void *) Phi);
    free((void *) aux);
    free((void *) aux2);
    free((void *) aux3);
    free((void *) d);
    free((void *) C);
    free((void *) s);

    return tp;
}


/**
 * @brief Evaluate a ThinPlate in an array of 'np' points 'Pin' and
 *        return the result in the array 'Pout'
 * @param tp - input ThinPlate structure
 * @param Pin - Array of input 2D points
 * @param Pout - Array of output 2D points
 * @param np   - number of points
 * @return EXIT_SUCCESS
 */
int evaluate_thinPlate(ThinPlate tp, float *Pin, float *Pout, int np)
{
    int i, j;
    float xPo, yPo, xP, yP, r2;


    for (i = 0; i < np; i++)
    {
        /*Read the point */
        xP = Pin[2 * i];
        yP = Pin[2 * i + 1];

        /*Affine Part */
        xPo = tp->affine[0] * xP + tp->affine[1] * yP + tp->affine[2];
        yPo = tp->affine[3] * xP + tp->affine[4] * yP + tp->affine[5];

        /*Kernel part */
        for (j = 0; j < tp->nc; j++)
        {
            r2 = (xP - tp->xc[j]) * (xP - tp->xc[j])
                + (yP - tp->yc[j]) * (yP -tp->yc[j]);
            if (r2 > EPS_ZERO)
            {
                xPo += tp->coef_x[j] * r2 * log(r2);
                yPo += tp->coef_y[j] * r2 * log(r2);
            }
        }
        Pout[2 * i] = xPo;
        Pout[2 * i + 1] = yPo;
    }

    return EXIT_SUCCESS;
}
