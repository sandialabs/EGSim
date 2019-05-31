      SUBROUTINE LEV2PCHOL( UPLO, N, A, LDA, PIV, RANK, TOL, WORK,
     $                      INFO )
*
*     Modified to include pivoting for semidefinite matrices by
*     Craig Lucas, University of Manchester. January, 2004
*
*     Original LAPACK routine DPOTF2
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   TOL
      INTEGER            INFO, LDA, N, RANK
      CHARACTER          UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( 2*N )
      INTEGER            PIV( N )
*     ..
*
*  Purpose
*  =======
*
*  LEV2PCHOL computes the Cholesky factorization with complete
*  pivoting of a real symmetric positive semidefinite matrix A.
*
*  The factorization has the form
*     P' * A * P = U' * U ,  if UPLO = 'U',
*     P' * A * P = L  * L',  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular, and
*  P is stored as vector PIV.
*
*  This algorithm does not attempt to check that A is positive
*  semidefinite. This version of the algorithm calls level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n by n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization as above.
*
*  PIV     (output) INTEGER array, dimension (N)
*          PIV is such that the nonzero entries are P( PIV(K), K ) = 1.
*
*  RANK    (output) INTEGER
*          The rank of A given by the number of steps the algorithm 
*          completed.
*
*  TOL     (input) DOUBLE PRECISION
*          User defined tolerance. If TOL < 0, then N*U*MAX( A( K,K ) )
*          will be used. The algorithm terminates at the (K-1)st step
*          if the pivot <= TOL.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  WORK    DOUBLE PRECISION array, dimension (2*N)
*          Work space.
*
*  INFO    (output) INTEGER
*          < 0: if INFO = -K, the K-th argument had an illegal value
*          = 0  algorithm completed successfully.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AJJ, DSTOP, DTEMP, U
      INTEGER            ITEMP, J, P, PVT
      LOGICAL            UPPER
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      LOGICAL            LSAME
      EXTERNAL           DLAMCH, LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLAS_DMAX_VAL, DGEMV, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'LEV2PC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Initialize PIV
*
      DO 10 P = 1, N
         PIV( P ) = P
   10 CONTINUE
*
*     Get unit roundoff
*
      U = DLAMCH( 'E' )
*
*     Compute stopping value
*
      CALL BLAS_DMAX_VAL( N, A( 1, 1 ), LDA+1, PVT, DTEMP )
      AJJ = A( PVT, PVT )
      IF( AJJ.EQ.ZERO ) THEN
         RANK = 0
         GO TO 80
      END IF
*
*     Compute stopping value if not supplied
*
      IF( TOL.LT.ZERO ) THEN
         DSTOP = N*U*AJJ
      ELSE
         DSTOP = TOL
      END IF
*
*     Set first half of WORK to zero, holds dot products
*
      DO 20 P = 1, N
         WORK( P ) = 0
   20 CONTINUE
*
      IF( UPPER ) THEN
*
*        Compute the Cholesky factorization P' * A * P = U' * U
*
         DO 40 J = 1, N
*
*        Find pivot, test for exit, else swap rows and columns
*        Update dot products, compute possible pivots which are
*        stored in the second half of WORK
*
            DO 30 P = J, N
*
               IF( J.GT.1 ) THEN
                  WORK( P ) = WORK( P ) + A( J-1, P )**2
               END IF
               WORK( N+P ) = A( P, P ) - WORK( P )
*
   30       CONTINUE
*
            IF( J.GT.1 ) THEN
               CALL BLAS_DMAX_VAL( N-J+1, WORK( N+J ), 1, ITEMP, DTEMP )
               PVT = ITEMP + J - 1
               AJJ = WORK( N+PVT )
               IF( AJJ.LE.DSTOP ) THEN
                  A( J, J ) = AJJ
                  GO TO 70
               END IF
            END IF
*
            IF( J.NE.PVT ) THEN
*
*              Pivot OK, so can now swap pivot rows and columns
*
               A( PVT, PVT ) = A( J, J )
               CALL DSWAP( J-1, A( 1, J ), 1, A( 1, PVT ), 1 )
               IF( PVT.LT.N )
     $            CALL DSWAP( N-PVT, A( J, PVT+1 ), LDA,
     $                        A( PVT, PVT+1 ), LDA )
               CALL DSWAP( PVT-J-1, A( J, J+1 ), LDA, A( J+1, PVT ), 1 )
*
*              Swap dot products and PIV
*
               DTEMP = WORK( J )
               WORK( J ) = WORK( PVT )
               WORK( PVT ) = DTEMP
               ITEMP = PIV( PVT )
               PIV( PVT ) = PIV( J )
               PIV( J ) = ITEMP
            END IF
*
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of row J
*
            IF( J.LT.N ) THEN
               CALL DGEMV( 'Trans', J-1, N-J, -ONE, A( 1, J+1 ), LDA,
     $                     A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
*
   40    CONTINUE
*
      ELSE
*
*        Compute the Cholesky factorization P' * A * P = L * L'
*
         DO 60 J = 1, N
*
*        Find pivot, test for exit, else swap rows and columns
*        Update dot products, compute possible pivots which are
*        stored in the second half of WORK
*
            DO 50 P = J, N
*
               IF( J.GT.1 ) THEN
                  WORK( P ) = WORK( P ) + A( P, J-1 )**2
               END IF
               WORK( N+P ) = A( P, P ) - WORK( P )
*
   50       CONTINUE
*
            IF( J.GT.1 ) THEN
               CALL BLAS_DMAX_VAL( N-J+1, WORK( N+J ), 1, ITEMP, DTEMP )
               PVT = ITEMP + J - 1
               AJJ = WORK( N+PVT )
               IF( AJJ.LE.DSTOP ) THEN
                  A( J, J ) = AJJ
                  GO TO 70
               END IF
            END IF
*
            IF( J.NE.PVT ) THEN
*
*              Pivot OK, so can now swap pivot rows and columns
*
               A( PVT, PVT ) = A( J, J )
               CALL DSWAP( J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA )
               IF( PVT.LT.N )
     $            CALL DSWAP( N-PVT, A( PVT+1, J ), 1, A( PVT+1, PVT ),
     $                        1 )
               CALL DSWAP( PVT-J-1, A( J+1, J ), 1, A( PVT, J+1 ), LDA )
*
*              Swap dot products and PIV
*
               DTEMP = WORK( J )
               WORK( J ) = WORK( PVT )
               WORK( PVT ) = DTEMP
               ITEMP = PIV( PVT )
               PIV( PVT ) = PIV( J )
               PIV( J ) = ITEMP
            END IF
*
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of column J
*
            IF( J.LT.N ) THEN
               CALL DGEMV( 'No tran', N-J, J-1, -ONE, A( J+1, 1 ), LDA,
     $                     A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
*
   60    CONTINUE
*
      END IF
*
*     Ran to completion, A has full rank
*
      RANK = N
*
      GO TO 80
   70 CONTINUE
*
*     Rank is number of steps completed
*
      RANK = J - 1
*
   80 CONTINUE
      RETURN
*
*     End of LEV2PCHOL
*
      END

