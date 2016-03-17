module module_fftpack

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: rffti
    public :: rfftf
    public :: rfftb
    public :: ezffti
    public :: ezfftf
    public :: ezfftb
    public :: sinti
    public :: sint
    public :: costi
    public :: cost
    public :: sinqi
    public :: sinqf
    public :: sinqb
    public :: cosqi
    public :: cosqf
    public :: cosqb
    public :: cffti
    public :: cfftf
    public :: cfftb

contains
    !
    !     file fftpack.f
    !
    !
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !     *                                                               *
    !     *                  copyright (c) 2005 by UCAR                   *
    !     *                                                               *
    !     *       University Corporation for Atmospheric Research         *
    !     *                                                               *
    !     *                      all rights reserved                      *
    !     *                                                               *
    !     *                    FISHPACK90  version 1.1                    *
    !     *                                                               *
    !     *                 A Package of Fortran 77 and 90                *
    !     *                                                               *
    !     *                Subroutines and Example Programs               *
    !     *                                                               *
    !     *               for Modeling Geophysical Processes              *
    !     *                                                               *
    !     *                             by                                *
    !     *                                                               *
    !     *        John Adams, Paul Swarztrauber and Roland Sweet         *
    !     *                                                               *
    !     *                             of                                *
    !     *                                                               *
    !     *         the National Center for Atmospheric Research          *
    !     *                                                               *
    !     *                Boulder, Colorado  (80307)  U.S.A.             *
    !     *                                                               *
    !     *                   which is sponsored by                       *
    !     *                                                               *
    !     *              the National Science Foundation                  *
    !     *                                                               *
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !
    ! LATEST REVISION
    ! ---------------
    !     June 2004      (VERSION 5.0) FORTRAN 90 CHANGES
    !
    ! PURPOSE
    ! -------
    !     THIS PACKAGE CONSISTS OF PROGRAMS WHICH PERFORM FAST FOURIER
    !     TRANSFORMS FOR BOTH COMPLEX AND REAL PERIODIC SEQUENCES AND
    !     CERTAIN OTHER SYMMETRIC SEQUENCES THAT ARE LISTED BELOW.
    !
    ! USAGE
    ! -----
    !     1.   RFFTI     INITIALIZE  RFFTF AND RFFTB
    !     2.   RFFTF     FORWARD TRANSFORM OF A REAL PERIODIC SEQUENCE
    !     3.   RFFTB     BACKWARD TRANSFORM OF A REAL COEFFICIENT ARRAY
    !
    !     4.   EZFFTI    INITIALIZE EZFFTF AND EZFFTB
    !     5.   EZFFTF    A SIMPLIFIED REAL PERIODIC FORWARD TRANSFORM
    !     6.   EZFFTB    A SIMPLIFIED REAL PERIODIC BACKWARD TRANSFORM
    !
    !     7.   SINTI     INITIALIZE SINT
    !     8.   SINT      SINE TRANSFORM OF A REAL ODD SEQUENCE
    !
    !     9.   COSTI     INITIALIZE COST
    !     10.  COST      COSINE TRANSFORM OF A REAL EVEN SEQUENCE
    !
    !     11.  SINQI     INITIALIZE SINQF AND SINQB
    !     12.  SINQF     FORWARD SINE TRANSFORM WITH ODD WAVE NUMBERS
    !     13.  SINQB     UNNORMALIZED INVERSE OF SINQF
    !
    !     14.  COSQI     INITIALIZE COSQF AND COSQB
    !     15.  COSQF     FORWARD COSINE TRANSFORM WITH ODD WAVE NUMBERS
    !     16.  COSQB     UNNORMALIZED INVERSE OF COSQF
    !
    !     17.  CFFTI     INITIALIZE CFFTF AND CFFTB
    !     18.  CFFTF     FORWARD TRANSFORM OF A COMPLEX PERIODIC SEQUENCE
    !     19.  CFFTB     UNNORMALIZED INVERSE OF CFFTF
    !
    ! SPECIAL CONDITIONS
    ! ------------------
    !     BEFORE CALLING ROUTINES EZFFTB AND EZFFTF FOR THE FIRST TIME,
    !     OR BEFORE CALLING EZFFTB AND EZFFTF WITH A DIFFERENT LENGTH,
    !     USERS MUST INITIALIZE BY CALLING ROUTINE EZFFTI.
    !
    ! I/O
    ! ---
    !     NONE
    !
    ! PRECISION
    ! ---------
    !     NONE
    !
    ! REQUIRED LIBRARY FILES
    ! ----------------------
    !     NONE
    !
    ! LANGUAGE
    ! --------
    !     FORTRAN
    !
    ! HISTORY
    ! -------
    !     DEVELOPED AT NCAR IN BOULDER, COLORADO BY PAUL N. SWARZTRAUBER
    !     OF THE SCIENTIFIC COMPUTING DIVISION.  RELEASED ON NCAR'S PUBLIC
    !     SOFTWARE LIBRARIES IN JANUARY 1980.  MODIFIED MAY 29, 1985 TO
    !     INCREASE EFFICIENCY.
    !
    ! PORTABILITY
    ! -----------
    !     FORTRAN 77
    !
    ! **********************************************************************
    !
    !     SUBROUTINE RFFTI(N, wsave)
    !
    !     SUBROUTINE RFFTI INITIALIZES THE ARRAY wsave WHICH IS USED IN
    !     BOTH RFFTF AND RFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
    !     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
    !     STORED IN wsave.
    !
    !     INPUT PARAMETER
    !
    !     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.
    !
    !     OUTPUT PARAMETER
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
    !             THE SAME WORK ARRAY CAN BE USED FOR BOTH RFFTF AND RFFTB
    !             AS LONG AS N REMAINS UNCHANGED. DIFFERENT wsave ARRAYS
    !             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
    !             wsave MUST NOT BE CHANGED BETWEEN CALLS OF RFFTF OR RFFTB.
    !
    ! **********************************************************************
    !
    !     SUBROUTINE RFFTF(N, R, wsave)
    !
    !     SUBROUTINE RFFTF COMPUTES THE FOURIER COEFFICIENTS OF A REAL
    !     PERODIC SEQUENCE (FOURIER ANALYSIS). THE TRANSFORM IS DEFINED
    !     BELOW AT OUTPUT PARAMETER R.
    !
    !     INPUT PARAMETERS
    !
    !     N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
    !             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
    !             N MAY CHANGE SO LONG AS DIFFERENT WORK ARRAYS ARE PROVIDED
    !
    !     R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
    !             TO BE TRANSFORMED
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
    !             IN THE PROGRAM THAT CALLS RFFTF. THE wsave ARRAY MUST BE
    !             INITIALIZED BY CALLING SUBROUTINE RFFTI(N, wsave) AND A
    !             DIFFERENT wsave ARRAY MUST BE USED FOR EACH DIFFERENT
    !             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
    !             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
    !             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
    !             THE SAME wsave ARRAY CAN BE USED BY RFFTF AND RFFTB.
    !
    !
    !     OUTPUT PARAMETERS
    !
    !     R       R(1) = THE SUM FROM I=1 TO I=N OF R(I)
    !
    !             IF N IS EVEN SET L =N/2   , IF N IS ODD SET L = (N+1)/2
    !
    !               THEN FOR K = 2, ..., L
    !
    !                  R(2*K-2) = THE SUM FROM I = 1 TO I = N OF
    !
    !                       R(I)*COS((K-1)*(I-1)*2*PI/N)
    !
    !                  R(2*K-1) = THE SUM FROM I = 1 TO I = N OF
    !
    !                      -R(I)*SIN((K-1)*(I-1)*2*PI/N)
    !
    !             IF N IS EVEN
    !
    !                  R(N) = THE SUM FROM I = 1 TO I = N OF
    !
    !                       (-1)**(I-1)*R(I)
    !
    !      *****  NOTE
    !                  THIS TRANSFORM IS UNNORMALIZED SINCE A CALL OF RFFTF
    !                  FOLLOWED BY A CALL OF RFFTB WILL MULTIPLY THE INPUT
    !                  SEQUENCE BY N.
    !
    !     wsave   CONTAINS RESULTS WHICH MUST NOT BE DESTROYED BETWEEN
    !             CALLS OF RFFTF OR RFFTB.
    !
    !
    ! **********************************************************************
    !
    !     SUBROUTINE RFFTB(N, R, wsave)
    !
    !     SUBROUTINE RFFTB COMPUTES THE REAL PERODIC SEQUENCE FROM ITS
    !     FOURIER COEFFICIENTS (FOURIER SYNTHESIS). THE TRANSFORM IS DEFINED
    !     BELOW AT OUTPUT PARAMETER R.
    !
    !     INPUT PARAMETERS
    !
    !     N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
    !             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
    !             N MAY CHANGE SO LONG AS DIFFERENT WORK ARRAYS ARE PROVIDED
    !
    !     R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
    !             TO BE TRANSFORMED
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
    !             IN THE PROGRAM THAT CALLS RFFTB. THE wsave ARRAY MUST BE
    !             INITIALIZED BY CALLING SUBROUTINE RFFTI(N, wsave) AND A
    !             DIFFERENT wsave ARRAY MUST BE USED FOR EACH DIFFERENT
    !             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
    !             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
    !             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
    !             THE SAME wsave ARRAY CAN BE USED BY RFFTF AND RFFTB.
    !
    !
    !     OUTPUT PARAMETERS
    !
    !     R       FOR N EVEN AND FOR I = 1, ..., N
    !
    !                  R(I) = R(1)+(-1)**(I-1)*R(N)
    !
    !                       PLUS THE SUM FROM K=2 TO K=N/2 OF
    !
    !                        2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
    !
    !                       -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
    !
    !             FOR N ODD AND FOR I = 1, ..., N
    !
    !                  R(I) = R(1) PLUS THE SUM FROM K=2 TO K=(N+1)/2 OF
    !
    !                       2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
    !
    !                      -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
    !
    !      *****  NOTE
    !                  THIS TRANSFORM IS UNNORMALIZED SINCE A CALL OF RFFTF
    !                  FOLLOWED BY A CALL OF RFFTB WILL MULTIPLY THE INPUT
    !                  SEQUENCE BY N.
    !
    !     wsave   CONTAINS RESULTS WHICH MUST NOT BE DESTROYED BETWEEN
    !             CALLS OF RFFTB OR RFFTF.
    !
    !
    ! **********************************************************************
    !
    !     SUBROUTINE EZFFTI(N, wsave)
    !
    !     SUBROUTINE EZFFTI INITIALIZES THE ARRAY wsave WHICH IS USED IN
    !     BOTH EZFFTF AND EZFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
    !     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
    !     STORED IN wsave.
    !
    !     INPUT PARAMETER
    !
    !     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.
    !
    !     OUTPUT PARAMETER
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
    !             THE SAME WORK ARRAY CAN BE USED FOR BOTH EZFFTF AND EZFFTB
    !             AS LONG AS N REMAINS UNCHANGED. DIFFERENT wsave ARRAYS
    !             ARE REQUIRED FOR DIFFERENT VALUES OF N.
    !
    !
    ! **********************************************************************
    !
    !     SUBROUTINE EZFFTF(N, R, AZERO, A, B, wsave)
    !
    !     SUBROUTINE EZFFTF COMPUTES THE FOURIER COEFFICIENTS OF A REAL
    !     PERODIC SEQUENCE (FOURIER ANALYSIS). THE TRANSFORM IS DEFINED
    !     BELOW AT OUTPUT PARAMETERS AZERO, A AND B. EZFFTF IS A SIMPLIFIED
    !     BUT SLOWER VERSION OF RFFTF.
    !
    !     INPUT PARAMETERS
    !
    !     N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
    !             IS MUST EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
    !
    !     R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
    !             TO BE TRANSFORMED. R IS NOT DESTROYED.
    !
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
    !             IN THE PROGRAM THAT CALLS EZFFTF. THE wsave ARRAY MUST BE
    !             INITIALIZED BY CALLING SUBROUTINE EZFFTI(N, wsave) AND A
    !             DIFFERENT wsave ARRAY MUST BE USED FOR EACH DIFFERENT
    !             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
    !             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
    !             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
    !             THE SAME wsave ARRAY CAN BE USED BY EZFFTF AND EZFFTB.
    !
    !     OUTPUT PARAMETERS
    !
    !     AZERO   THE SUM FROM I=1 TO I=N OF R(I)/N
    !
    !     A, B     FOR N EVEN B(N/2)=0. AND A(N/2) IS THE SUM FROM I=1 TO
    !             I=N OF (-1)**(I-1)*R(I)/N
    !
    !             FOR N EVEN DEFINE KMAX=N/2-1
    !             FOR N ODD  DEFINE KMAX=(N-1)/2
    !
    !             THEN FOR  K=1, ..., KMAX
    !
    !                  A(K) EQUALS THE SUM FROM I=1 TO I=N OF
    !
    !                       2./N*R(I)*COS(K*(I-1)*2*PI/N)
    !
    !                  B(K) EQUALS THE SUM FROM I=1 TO I=N OF
    !
    !                       2./N*R(I)*SIN(K*(I-1)*2*PI/N)
    !
    !
    ! **********************************************************************
    !
    !     SUBROUTINE EZFFTB(N, R, AZERO, A, B, wsave)
    !
    !     SUBROUTINE EZFFTB COMPUTES A REAL PERODIC SEQUENCE FROM ITS
    !     FOURIER COEFFICIENTS (FOURIER SYNTHESIS). THE TRANSFORM IS
    !     DEFINED BELOW AT OUTPUT PARAMETER R. EZFFTB IS A SIMPLIFIED
    !     BUT SLOWER VERSION OF RFFTB.
    !
    !     INPUT PARAMETERS
    !
    !     N       THE LENGTH OF THE OUTPUT ARRAY R.  THE METHOD IS MOST
    !             EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
    !
    !     AZERO   THE CONSTANT FOURIER COEFFICIENT
    !
    !     A, B     ARRAYS WHICH CONTAIN THE REMAINING FOURIER COEFFICIENTS
    !             THESE ARRAYS ARE NOT DESTROYED.
    !
    !             THE LENGTH OF THESE ARRAYS DEPENDS ON WHETHER N IS EVEN OR
    !             ODD.
    !
    !             IF N IS EVEN N/2    LOCATIONS ARE REQUIRED
    !             IF N IS ODD (N-1)/2 LOCATIONS ARE REQUIRED
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
    !             IN THE PROGRAM THAT CALLS EZFFTB. THE wsave ARRAY MUST BE
    !             INITIALIZED BY CALLING SUBROUTINE EZFFTI(N, wsave) AND A
    !             DIFFERENT wsave ARRAY MUST BE USED FOR EACH DIFFERENT
    !             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
    !             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
    !             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
    !             THE SAME wsave ARRAY CAN BE USED BY EZFFTF AND EZFFTB.
    !
    !
    !     OUTPUT PARAMETERS
    !
    !     R       IF N IS EVEN DEFINE KMAX=N/2
    !             IF N IS ODD  DEFINE KMAX=(N-1)/2
    !
    !             THEN FOR I=1, ..., N
    !
    !                  R(I)=AZERO PLUS THE SUM FROM K=1 TO K=KMAX OF
    !
    !                  A(K)*COS(K*(I-1)*2*PI/N)+B(K)*SIN(K*(I-1)*2*PI/N)
    !
    !     ********************* COMPLEX NOTATION **************************
    !
    !             FOR J=1, ..., N
    !
    !             R(J) EQUALS THE SUM FROM K=-KMAX TO K=KMAX OF
    !
    !                  C(K)*EXP(I*K*(J-1)*2*PI/N)
    !
    !             WHERE
    !
    !                  C(K) = .5*CMPLX(A(K), -B(K))   FOR K=1, ..., KMAX
    !
    !                  C(-K) = CONJG(C(K))
    !
    !                  C(0) = AZERO
    !
    !                       AND I=SQRT(-1)
    !
    !     *************** AMPLITUDE - PHASE NOTATION ***********************
    !
    !             FOR I=1, ..., N
    !
    !             R(I) EQUALS AZERO PLUS THE SUM FROM K=1 TO K=KMAX OF
    !
    !                  ALPHA(K)*COS(K*(I-1)*2*PI/N+BETA(K))
    !
    !             WHERE
    !
    !                  ALPHA(K) = SQRT(A(K)*A(K)+B(K)*B(K))
    !
    !                  COS(BETA(K))=A(K)/ALPHA(K)
    !
    !                  SIN(BETA(K))=-B(K)/ALPHA(K)
    !
    ! **********************************************************************
    !
    !     SUBROUTINE SINTI(N, wsave)
    !
    !     SUBROUTINE SINTI INITIALIZES THE ARRAY wsave WHICH IS USED IN
    !     SUBROUTINE SINT. THE PRIME FACTORIZATION OF N TOGETHER WITH
    !     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
    !     STORED IN wsave.
    !
    !     INPUT PARAMETER
    !
    !     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
    !             IS MOST EFFICIENT WHEN N+1 IS A PRODUCT OF SMALL PRIMES.
    !
    !     OUTPUT PARAMETER
    !
    !     wsave   A WORK ARRAY WITH AT LEAST INT(2.5*N+15) LOCATIONS.
    !             DIFFERENT wsave ARRAYS ARE REQUIRED FOR DIFFERENT VALUES
    !             OF N. THE CONTENTS OF wsave MUST NOT BE CHANGED BETWEEN
    !             CALLS OF SINT.
    !
    ! **********************************************************************
    !
    !     SUBROUTINE SINT(N, X, wsave)
    !
    !     SUBROUTINE SINT COMPUTES THE DISCRETE FOURIER SINE TRANSFORM
    !     OF AN ODD SEQUENCE X(I). THE TRANSFORM IS DEFINED BELOW AT
    !     OUTPUT PARAMETER X.
    !
    !     SINT IS THE UNNORMALIZED INVERSE OF ITSELF SINCE A CALL OF SINT
    !     FOLLOWED BY ANOTHER CALL OF SINT WILL MULTIPLY THE INPUT SEQUENCE
    !     X BY 2*(N+1).
    !
    !     THE ARRAY wsave WHICH IS USED BY SUBROUTINE SINT MUST BE
    !     INITIALIZED BY CALLING SUBROUTINE SINTI(N, wsave).
    !
    !     INPUT PARAMETERS
    !
    !     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
    !             IS MOST EFFICIENT WHEN N+1 IS THE PRODUCT OF SMALL PRIMES.
    !
    !     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
    !
    !
    !     wsave   A WORK ARRAY WITH DIMENSION AT LEAST INT(2.5*N+15)
    !             IN THE PROGRAM THAT CALLS SINT. THE wsave ARRAY MUST BE
    !             INITIALIZED BY CALLING SUBROUTINE SINTI(N, wsave) AND A
    !             DIFFERENT wsave ARRAY MUST BE USED FOR EACH DIFFERENT
    !             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
    !             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
    !             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
    !
    !     OUTPUT PARAMETERS
    !
    !     X       FOR I=1, ..., N
    !
    !                  X(I)= THE SUM FROM K=1 TO K=N
    !
    !                       2*X(K)*SIN(K*I*PI/(N+1))
    !
    !                  A CALL OF SINT FOLLOWED BY ANOTHER CALL OF
    !                  SINT WILL MULTIPLY THE SEQUENCE X BY 2*(N+1).
    !                  HENCE SINT IS THE UNNORMALIZED INVERSE
    !                  OF ITSELF.
    !
    !     wsave   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
    !             DESTROYED BETWEEN CALLS OF SINT.
    !
    ! **********************************************************************
    !
    !     SUBROUTINE COSTI(N, wsave)
    !
    !     SUBROUTINE COSTI INITIALIZES THE ARRAY wsave WHICH IS USED IN
    !     SUBROUTINE COST. THE PRIME FACTORIZATION OF N TOGETHER WITH
    !     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
    !     STORED IN wsave.
    !
    !     INPUT PARAMETER
    !
    !     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
    !             IS MOST EFFICIENT WHEN N-1 IS A PRODUCT OF SMALL PRIMES.
    !
    !     OUTPUT PARAMETER
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
    !             DIFFERENT wsave ARRAYS ARE REQUIRED FOR DIFFERENT VALUES
    !             OF N. THE CONTENTS OF wsave MUST NOT BE CHANGED BETWEEN
    !             CALLS OF COST.
    !
    ! **********************************************************************
    !
    !     SUBROUTINE COST(N, X, wsave)
    !
    !     SUBROUTINE COST COMPUTES THE DISCRETE FOURIER COSINE TRANSFORM
    !     OF AN EVEN SEQUENCE X(I). THE TRANSFORM IS DEFINED BELOW AT OUTPUT
    !     PARAMETER X.
    !
    !     COST IS THE UNNORMALIZED INVERSE OF ITSELF SINCE A CALL OF COST
    !     FOLLOWED BY ANOTHER CALL OF COST WILL MULTIPLY THE INPUT SEQUENCE
    !     X BY 2*(N-1). THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER X
    !
    !     THE ARRAY wsave WHICH IS USED BY SUBROUTINE COST MUST BE
    !     INITIALIZED BY CALLING SUBROUTINE COSTI(N, wsave).
    !
    !     INPUT PARAMETERS
    !
    !     N       THE LENGTH OF THE SEQUENCE X. N MUST BE GREATER THAN 1.
    !             THE METHOD IS MOST EFFICIENT WHEN N-1 IS A PRODUCT OF
    !             SMALL PRIMES.
    !
    !     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15
    !             IN THE PROGRAM THAT CALLS COST. THE wsave ARRAY MUST BE
    !             INITIALIZED BY CALLING SUBROUTINE COSTI(N, wsave) AND A
    !             DIFFERENT wsave ARRAY MUST BE USED FOR EACH DIFFERENT
    !             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
    !             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
    !             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
    !
    !     OUTPUT PARAMETERS
    !
    !     X       FOR I=1, ..., N
    !
    !                 X(I) = X(1)+(-1)**(I-1)*X(N)
    !
    !                  + THE SUM FROM K=2 TO K=N-1
    !
    !                      2*X(K)*COS((K-1)*(I-1)*PI/(N-1))
    !
    !                  A CALL OF COST FOLLOWED BY ANOTHER CALL OF
    !                  COST WILL MULTIPLY THE SEQUENCE X BY 2*(N-1)
    !                  HENCE COST IS THE UNNORMALIZED INVERSE
    !                  OF ITSELF.
    !
    !     wsave   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
    !             DESTROYED BETWEEN CALLS OF COST.
    !
    ! **********************************************************************
    !
    !     SUBROUTINE SINQI(N, wsave)
    !
    !     SUBROUTINE SINQI INITIALIZES THE ARRAY wsave WHICH IS USED IN
    !     BOTH SINQF AND SINQB. THE PRIME FACTORIZATION OF N TOGETHER WITH
    !     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
    !     STORED IN wsave.
    !
    !     INPUT PARAMETER
    !
    !     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED. THE METHOD
    !             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
    !
    !     OUTPUT PARAMETER
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
    !             THE SAME WORK ARRAY CAN BE USED FOR BOTH SINQF AND SINQB
    !             AS LONG AS N REMAINS UNCHANGED. DIFFERENT wsave ARRAYS
    !             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
    !             wsave MUST NOT BE CHANGED BETWEEN CALLS OF SINQF OR SINQB.
    !
    ! **********************************************************************
    !
    !     SUBROUTINE SINQF(N, X, wsave)
    !
    !     SUBROUTINE SINQF COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
    !     WAVE DATA. THAT IS , SINQF COMPUTES THE COEFFICIENTS IN A SINE
    !     SERIES REPRESENTATION WITH ONLY ODD WAVE NUMBERS. THE TRANSFORM
    !     IS DEFINED BELOW AT OUTPUT PARAMETER X.
    !
    !     SINQB IS THE UNNORMALIZED INVERSE OF SINQF SINCE A CALL OF SINQF
    !     FOLLOWED BY A CALL OF SINQB WILL MULTIPLY THE INPUT SEQUENCE X
    !     BY 4*N.
    !
    !     THE ARRAY wsave WHICH IS USED BY SUBROUTINE SINQF MUST BE
    !     INITIALIZED BY CALLING SUBROUTINE SINQI(N, wsave).
    !
    !
    !     INPUT PARAMETERS
    !
    !     N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
    !             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
    !
    !     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
    !             IN THE PROGRAM THAT CALLS SINQF. THE wsave ARRAY MUST BE
    !             INITIALIZED BY CALLING SUBROUTINE SINQI(N, wsave) AND A
    !             DIFFERENT wsave ARRAY MUST BE USED FOR EACH DIFFERENT
    !             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
    !             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
    !             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
    !
    !     OUTPUT PARAMETERS
    !
    !     X       FOR I=1, ..., N
    !
    !                  X(I) = (-1)**(I-1)*X(N)
    !
    !                     + THE SUM FROM K=1 TO K=N-1 OF
    !
    !                     2*X(K)*SIN((2*I-1)*K*PI/(2*N))
    !
    !                  A CALL OF SINQF FOLLOWED BY A CALL OF
    !                  SINQB WILL MULTIPLY THE SEQUENCE X BY 4*N.
    !                  THEREFORE SINQB IS THE UNNORMALIZED INVERSE
    !                  OF SINQF.
    !
    !     wsave   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
    !             BE DESTROYED BETWEEN CALLS OF SINQF OR SINQB.
    !
    ! **********************************************************************
    !
    !     SUBROUTINE SINQB(N, X, wsave)
    !
    !     SUBROUTINE SINQB COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
    !     WAVE DATA. THAT IS , SINQB COMPUTES A SEQUENCE FROM ITS
    !     REPRESENTATION IN TERMS OF A SINE SERIES WITH ODD WAVE NUMBERS.
    !     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER X.
    !
    !     SINQF IS THE UNNORMALIZED INVERSE OF SINQB SINCE A CALL OF SINQB
    !     FOLLOWED BY A CALL OF SINQF WILL MULTIPLY THE INPUT SEQUENCE X
    !     BY 4*N.
    !
    !     THE ARRAY wsave WHICH IS USED BY SUBROUTINE SINQB MUST BE
    !     INITIALIZED BY CALLING SUBROUTINE SINQI(N, wsave).
    !
    !
    !     INPUT PARAMETERS
    !
    !     N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
    !             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
    !
    !     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
    !             IN THE PROGRAM THAT CALLS SINQB. THE wsave ARRAY MUST BE
    !             INITIALIZED BY CALLING SUBROUTINE SINQI(N, wsave) AND A
    !             DIFFERENT wsave ARRAY MUST BE USED FOR EACH DIFFERENT
    !             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
    !             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
    !             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
    !
    !     OUTPUT PARAMETERS
    !
    !     X       FOR I=1, ..., N
    !
    !                  X(I)= THE SUM FROM K=1 TO K=N OF
    !
    !                    4*X(K)*SIN((2K-1)*I*PI/(2*N))
    !
    !                  A CALL OF SINQB FOLLOWED BY A CALL OF
    !                  SINQF WILL MULTIPLY THE SEQUENCE X BY 4*N.
    !                  THEREFORE SINQF IS THE UNNORMALIZED INVERSE
    !                  OF SINQB.
    !
    !     wsave   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
    !             BE DESTROYED BETWEEN CALLS OF SINQB OR SINQF.
    !
    ! **********************************************************************
    !
    !     SUBROUTINE COSQI(N, wsave)
    !
    !     SUBROUTINE COSQI INITIALIZES THE ARRAY wsave WHICH IS USED IN
    !     BOTH COSQF AND COSQB. THE PRIME FACTORIZATION OF N TOGETHER WITH
    !     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
    !     STORED IN wsave.
    !
    !     INPUT PARAMETER
    !
    !     N       THE LENGTH OF THE ARRAY TO BE TRANSFORMED.  THE METHOD
    !             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
    !
    !     OUTPUT PARAMETER
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
    !             THE SAME WORK ARRAY CAN BE USED FOR BOTH COSQF AND COSQB
    !             AS LONG AS N REMAINS UNCHANGED. DIFFERENT wsave ARRAYS
    !             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
    !             wsave MUST NOT BE CHANGED BETWEEN CALLS OF COSQF OR COSQB.
    !
    ! **********************************************************************
    !
    !     SUBROUTINE COSQF(N, X, wsave)
    !
    !     SUBROUTINE COSQF COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
    !     WAVE DATA. THAT IS , COSQF COMPUTES THE COEFFICIENTS IN A COSINE
    !     SERIES REPRESENTATION WITH ONLY ODD WAVE NUMBERS. THE TRANSFORM
    !     IS DEFINED BELOW AT OUTPUT PARAMETER X
    !
    !     COSQF IS THE UNNORMALIZED INVERSE OF COSQB SINCE A CALL OF COSQF
    !     FOLLOWED BY A CALL OF COSQB WILL MULTIPLY THE INPUT SEQUENCE X
    !     BY 4*N.
    !
    !     THE ARRAY wsave WHICH IS USED BY SUBROUTINE COSQF MUST BE
    !     INITIALIZED BY CALLING SUBROUTINE COSQI(N, wsave).
    !
    !
    !     INPUT PARAMETERS
    !
    !     N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
    !             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
    !
    !     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15
    !             IN THE PROGRAM THAT CALLS COSQF. THE wsave ARRAY MUST BE
    !             INITIALIZED BY CALLING SUBROUTINE COSQI(N, wsave) AND A
    !             DIFFERENT wsave ARRAY MUST BE USED FOR EACH DIFFERENT
    !             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
    !             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
    !             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
    !
    !     OUTPUT PARAMETERS
    !
    !     X       FOR I=1, ..., N
    !
    !                  X(I) = X(1) PLUS THE SUM FROM K=2 TO K=N OF
    !
    !                     2*X(K)*COS((2*I-1)*(K-1)*PI/(2*N))
    !
    !                  A CALL OF COSQF FOLLOWED BY A CALL OF
    !                  COSQB WILL MULTIPLY THE SEQUENCE X BY 4*N.
    !                  THEREFORE COSQB IS THE UNNORMALIZED INVERSE
    !                  OF COSQF.
    !
    !     wsave   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
    !             BE DESTROYED BETWEEN CALLS OF COSQF OR COSQB.
    !
    ! **********************************************************************
    !
    !     SUBROUTINE COSQB(N, X, wsave)
    !
    !     SUBROUTINE COSQB COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
    !     WAVE DATA. THAT IS , COSQB COMPUTES A SEQUENCE FROM ITS
    !     REPRESENTATION IN TERMS OF A COSINE SERIES WITH ODD WAVE NUMBERS.
    !     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER X.
    !
    !     COSQB IS THE UNNORMALIZED INVERSE OF COSQF SINCE A CALL OF COSQB
    !     FOLLOWED BY A CALL OF COSQF WILL MULTIPLY THE INPUT SEQUENCE X
    !     BY 4*N.
    !
    !     THE ARRAY wsave WHICH IS USED BY SUBROUTINE COSQB MUST BE
    !     INITIALIZED BY CALLING SUBROUTINE COSQI(N, wsave).
    !
    !
    !     INPUT PARAMETERS
    !
    !     N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
    !             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
    !
    !     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
    !
    !     wsave   A WORK ARRAY THAT MUST BE DIMENSIONED AT LEAST 3*N+15
    !             IN THE PROGRAM THAT CALLS COSQB. THE wsave ARRAY MUST BE
    !             INITIALIZED BY CALLING SUBROUTINE COSQI(N, wsave) AND A
    !             DIFFERENT wsave ARRAY MUST BE USED FOR EACH DIFFERENT
    !             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
    !             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
    !             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
    !
    !     OUTPUT PARAMETERS
    !
    !     X       FOR I=1, ..., N
    !
    !                  X(I)= THE SUM FROM K=1 TO K=N OF
    !
    !                    4*X(K)*COS((2*K-1)*(I-1)*PI/(2*N))
    !
    !                  A CALL OF COSQB FOLLOWED BY A CALL OF
    !                  COSQF WILL MULTIPLY THE SEQUENCE X BY 4*N.
    !                  THEREFORE COSQF IS THE UNNORMALIZED INVERSE
    !                  OF COSQB.
    !
    !     wsave   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
    !             BE DESTROYED BETWEEN CALLS OF COSQB OR COSQF.
    !
    ! **********************************************************************
    !
    !     SUBROUTINE CFFTI(N, wsave)
    !
    !     SUBROUTINE CFFTI INITIALIZES THE ARRAY wsave WHICH IS USED IN
    !     BOTH CFFTF AND CFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
    !     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
    !     STORED IN wsave.
    !
    !     INPUT PARAMETER
    !
    !     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED
    !
    !     OUTPUT PARAMETER
    !
    !     wsave   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4*N+15
    !             THE SAME WORK ARRAY CAN BE USED FOR BOTH CFFTF AND CFFTB
    !             AS LONG AS N REMAINS UNCHANGED. DIFFERENT wsave ARRAYS
    !             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
    !             wsave MUST NOT BE CHANGED BETWEEN CALLS OF CFFTF OR CFFTB.
    !
    ! **********************************************************************
    !
    !     SUBROUTINE CFFTF(N, C, wsave)
    !
    !     SUBROUTINE CFFTF COMPUTES THE FORWARD COMPLEX DISCRETE FOURIER
    !     TRANSFORM (THE FOURIER ANALYSIS). EQUIVALENTLY , CFFTF COMPUTES
    !     THE FOURIER COEFFICIENTS OF A COMPLEX PERIODIC SEQUENCE.
    !     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER C.
    !
    !     THE TRANSFORM IS NOT NORMALIZED. TO OBTAIN A NORMALIZED TRANSFORM
    !     THE OUTPUT MUST BE DIVIDED BY N. OTHERWISE A CALL OF CFFTF
    !     FOLLOWED BY A CALL OF CFFTB WILL MULTIPLY THE SEQUENCE BY N.
    !
    !     THE ARRAY wsave WHICH IS USED BY SUBROUTINE CFFTF MUST BE
    !     INITIALIZED BY CALLING SUBROUTINE CFFTI(N, wsave).
    !
    !     INPUT PARAMETERS
    !
    !
    !     N      THE LENGTH OF THE COMPLEX SEQUENCE C. THE METHOD IS
    !            MORE EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES. N
    !
    !     C      A COMPLEX ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
    !
    !     wsave   A REAL WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4N+15
    !             IN THE PROGRAM THAT CALLS CFFTF. THE wsave ARRAY MUST BE
    !             INITIALIZED BY CALLING SUBROUTINE CFFTI(N, wsave) AND A
    !             DIFFERENT wsave ARRAY MUST BE USED FOR EACH DIFFERENT
    !             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
    !             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
    !             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
    !             THE SAME wsave ARRAY CAN BE USED BY CFFTF AND CFFTB.
    !
    !     OUTPUT PARAMETERS
    !
    !     C      FOR J=1, ..., N
    !
    !                C(J)=THE SUM FROM K=1, ..., N OF
    !
    !                      C(K)*EXP(-I*(J-1)*(K-1)*2*PI/N)
    !
    !                            WHERE I=SQRT(-1)
    !
    !     wsave   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
    !             DESTROYED BETWEEN CALLS OF SUBROUTINE CFFTF OR CFFTB
    !
    ! **********************************************************************
    !
    !     SUBROUTINE CFFTB(N, C, wsave)
    !
    !     SUBROUTINE CFFTB COMPUTES THE BACKWARD COMPLEX DISCRETE FOURIER
    !     TRANSFORM (THE FOURIER SYNTHESIS). EQUIVALENTLY , CFFTB COMPUTES
    !     A COMPLEX PERIODIC SEQUENCE FROM ITS FOURIER COEFFICIENTS.
    !     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER C.
    !
    !     A CALL OF CFFTF FOLLOWED BY A CALL OF CFFTB WILL MULTIPLY THE
    !     SEQUENCE BY N.
    !
    !     THE ARRAY wsave WHICH IS USED BY SUBROUTINE CFFTB MUST BE
    !     INITIALIZED BY CALLING SUBROUTINE CFFTI(N, wsave).
    !
    !     INPUT PARAMETERS
    !
    !
    !     N      THE LENGTH OF THE COMPLEX SEQUENCE C. THE METHOD IS
    !            MORE EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
    !
    !     C      A COMPLEX ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
    !
    !     wsave   A REAL WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4N+15
    !             IN THE PROGRAM THAT CALLS CFFTB. THE wsave ARRAY MUST BE
    !             INITIALIZED BY CALLING SUBROUTINE CFFTI(N, wsave) AND A
    !             DIFFERENT wsave ARRAY MUST BE USED FOR EACH DIFFERENT
    !             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
    !             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
    !             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
    !             THE SAME wsave ARRAY CAN BE USED BY CFFTF AND CFFTB.
    !
    !     OUTPUT PARAMETERS
    !
    !     C      FOR J=1, ..., N
    !
    !                C(J)=THE SUM FROM K=1, ..., N OF
    !
    !                      C(K)*EXP(I*(J-1)*(K-1)*2*PI/N)
    !
    !                            WHERE I=SQRT(-1)
    !
    !     wsave   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
    !             DESTROYED BETWEEN CALLS OF SUBROUTINE CFFTF OR CFFTB
    ! **********************************************************************
    subroutine EZFFTF(n, r, azero, a, b, wsave)
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer  :: n
        real , intent (out) :: azero
        real , intent (in) :: r(*)
        real , intent (out) :: a(*)
        real , intent (out) :: b(*)
        real  :: wsave(*)
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: ns2, ns2m
        real :: cf, cfm
        !-----------------------------------------------
        !
        if (n - 2 <= 0) then
            if (n - 2 /= 0) then
                azero = R(1)
                return
            end if
            azero = 0.5*(R(1)+R(2))
            a(1) = 0.5*(R(1)-R(2))
            return
        end if
        wsave(:n) = R(:n)
        call RFFTF (n, wsave, wsave(n+1))
        cf = 2./real(n)
        cfm = -cf
        azero = 0.5*cf*wsave(1)
        ns2 = (n + 1)/2
        ns2m = ns2 - 1
        a(:ns2m) = cf*wsave(2:ns2m*2:2)
        b(:ns2m) = cfm*wsave(3:ns2m*2+1:2)
        if (mod(n, 2) == 1) return
        a(ns2) = 0.5*cf*wsave(n)
        b(ns2) = 0.
        return
    end subroutine EZFFTF


    subroutine EZFFTB(n, r, azero, a, b, wsave)
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer  :: n
        real , intent (in) :: azero
        real  :: r(*)
        real , intent (in) :: a(*)
        real , intent (in) :: b(*)
        real  :: wsave(*)
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: ns2
        !-----------------------------------------------
        !
        if (n - 2 <= 0) then
            if (n - 2 /= 0) then
                r(1) = azero
                return
            end if
            r(1) = azero + A(1)
            r(2) = azero - A(1)
            return
        end if
        ns2 = (n - 1)/2
        r(2:ns2*2:2) = 0.5*A(:ns2)
        r(3:ns2*2+1:2) = -0.5*B(:ns2)
        r(1) = azero
        if (mod(n, 2) == 0) r(n) = A(ns2+1)
        call RFFTB (n, r, wsave(n+1))
        return
    end subroutine EZFFTB


    subroutine EZFFTI(n, wsave)
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer  :: n
        real  :: wsave(*)
        !-----------------------------------------------
        !
        if (n == 1) return
        call EZFFT1 (n, wsave(2*n+1), int(wsave(3*n+1:3*n+1), kind=ip) )

    end subroutine EZFFTI


    subroutine EZFFT1(n, wa, ifac)
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer , intent (in) :: n
        integer :: ifac(*)
        !integer , intent (in out) :: ifac(*)
        real , intent (in out) :: wa(*)
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer , dimension(4) :: ntryh
        integer:: nl, nf, j, ntry, nq, nr, i, is, nfm1, l1, k1, ip_rename, l2, ido, ipm, ii
        real :: tpi, argh, arg1, ch1, sh1, dch1, dsh1, ch1h
        !-----------------------------------------------
        data NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 4, 2, 3, 5/
        tpi = 8.0*ATAN(1.0)
        nl = n
        nf = 0
        j = 0
101 continue
    j = j + 1
    if (j - 4 <= 0) then
        ntry = NTRYH(j)
    else
        ntry = ntry + 2
    end if
104 continue
    nq = nl/ntry
    nr = nl - ntry*nq
    if (nr /= 0) go to 101
    nf = nf + 1
    ifac(nf+2) = ntry
    nl = nq
    if (ntry == 2) then
        if (nf /= 1) then
            ifac(nf+2:4:(-1)) = IFAC(nf+1:3:(-1))
            ifac(3) = 2
        end if
    end if
    if (nl /= 1) go to 104
    ifac(1) = n
    ifac(2) = nf
    argh = tpi/real(n)
    is = 0
    nfm1 = nf - 1
    l1 = 1
    if (nfm1 == 0) return
    do k1 = 1, nfm1
        ip_rename = IFAC(k1+2)
        l2 = l1*ip_rename
        ido = n/l2
        ipm = ip_rename - 1
        arg1 = real(l1)*argh
        ch1 = 1.
        sh1 = 0.
        dch1 = COS(arg1)
        dsh1 = SIN(arg1)
        do j = 1, ipm
            ch1h = dch1*ch1 - dsh1*sh1
            sh1 = dch1*sh1 + dsh1*ch1
            ch1 = ch1h
            i = is + 2
            wa(i-1) = ch1
            wa(i) = sh1
            if (ido >= 5) then
                do ii = 5, ido, 2
                    i = i + 2
                    wa(i-1) = ch1*WA(i-3) - sh1*WA(i-2)
                    wa(i) = ch1*WA(i-2) + sh1*WA(i-3)
                end do
            end if
            is = is + ido
        end do
        l1 = l2
    end do

end subroutine EZFFT1

subroutine COSTI(n, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: n
    real  :: wsave(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: nm1, np1, ns2, k, kc
    real :: pi, dt, fk
    !-----------------------------------------------
    !
    pi = acos( -1.0 )
    if (n <= 3) return
    nm1 = n - 1
    np1 = n + 1
    ns2 = n/2
    dt = pi/real(nm1)
    fk = 0.
    do k = 2, ns2
        kc = np1 - k
        fk = fk + 1.
        wsave(k) = 2.*SIN(fk*dt)
        wsave(kc) = 2.*COS(fk*dt)
    end do
    call RFFTI (nm1, wsave(n+1))
    return
end subroutine COSTI


subroutine COST(n, x, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: n
    real  :: x(*)
    real  :: wsave(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: nm1, np1, ns2, k, kc, modn, i
    real :: x1h, x1p3, tx2, c1, t1, t2, xim2, xi
    !-----------------------------------------------
    !
    nm1 = n - 1
    np1 = n + 1
    ns2 = n/2
    if (n - 2 >= 0) then
        if (n - 2 <= 0) then
            x1h = X(1) + X(2)
            x(2) = X(1) - X(2)
            x(1) = x1h
            return
        end if
        if (n <= 3) then
            x1p3 = X(1) + X(3)
            tx2 = X(2) + X(2)
            x(2) = X(1) - X(3)
            x(1) = x1p3 + tx2
            x(3) = x1p3 - tx2
            return
        end if
        c1 = X(1) - X(n)
        x(1) = X(1) + X(n)
        do k = 2, ns2
            kc = np1 - k
            t1 = X(k) + X(kc)
            t2 = X(k) - X(kc)
            c1 = c1 + wsave(kc)*t2
            t2 = wsave(k)*t2
            x(k) = t1 - t2
            x(kc) = t1 + t2
        end do
        modn = mod(n, 2)
        if (modn /= 0) x(ns2+1) = X(ns2+1) + X(ns2+1)
        call RFFTF (nm1, x, wsave(n+1))
        xim2 = X(2)
        x(2) = c1
        do i = 4, n, 2
            xi = X(i)
            x(i) = X(i-2) - X(i-1)
            x(i-1) = xim2
            xim2 = xi
        end do
        if (modn /= 0) x(n) = xim2
    end if
    return
end subroutine COST

subroutine SINTI(n, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: n
    real  :: wsave(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: ns2, np1, k
    real :: pi, dt
    !-----------------------------------------------
    !
    pi = acos( -1.0 )
    if (n <= 1) return
    ns2 = n/2
    np1 = n + 1
    dt = pi/real(np1)
    do k = 1, ns2
        wsave(k) = 2.*SIN(k*dt)
    end do
    call RFFTI (np1, wsave(ns2+1))
    return
end subroutine SINTI

subroutine SINT(n, x, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: x(*)
    real  :: wsave(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: np1, iw1, iw2, iw3
    !-----------------------------------------------
    !
    np1 = n + 1
    iw1 = n/2 + 1
    iw2 = iw1 + np1
    iw3 = iw2 + np1
    call SINT1 (n, x, wsave, wsave(iw1), wsave(iw2), int(wsave(iw3:iw3),kind=ip) )
    return
end subroutine SINT

subroutine SINT1(n, war, was, xh, x, ifac)

    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer, intent (in) :: n
    integer             :: ifac(*)
    real (wp)                :: war(*)
    real (wp),    intent (in) :: was(*)
    real (wp)                :: xh(*)
    real (wp)                :: x(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer         :: i, np1, ns2, k, kc, modn
    real (wp)            :: xhold, t1, t2
    !-----------------------------------------------

    xh(:n) = WAR(:n)
    war(:n) = X(:n)
    if (n - 2 <= 0) then
        if (n - 2 /= 0) then
            xh(1) = XH(1) + XH(1)
            go to 106
        end if
        xhold = sqrt(3.0)*(XH(1)+XH(2))
        xh(2) = sqrt(3.0)*(XH(1)-XH(2))
        xh(1) = xhold
        go to 106
    end if
    np1 = n + 1
    ns2 = n/2
    x(1) = 0.
    do k = 1, ns2
        kc = np1 - k
        t1 = XH(k) - XH(kc)
        t2 = WAS(k)*(XH(k)+XH(kc))
        x(k+1) = t1 + t2
        x(kc+1) = t2 - t1
    end do
    modn = mod(n, 2)
    if (modn /= 0) x(ns2+2) = 4.*XH(ns2+1)
    call RFFTF1 (np1, x, xh, war, ifac)
    xh(1) = 0.5*X(1)
    do i = 3, n, 2
        xh(i-1) = -X(i)
        xh(i) = XH(i-2) + X(i-1)
    end do
    if (modn == 0) xh(n) = -X(n+1)
106 continue
    x(:n) = WAR(:n)
    war(:n) = XH(:n)
    return
end subroutine SINT1

subroutine COSQI(n, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: wsave(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k
    real :: pih, dt, fk
    !-----------------------------------------------
    !
    pih = 2.0*ATAN(1.0)
    dt = pih/real(n)
    fk = 0.
    do k = 1, n
        fk = fk + 1.
        wsave(k) = COS(fk*dt)
    end do
    call RFFTI (n, wsave(n+1))
    return
end subroutine COSQI

subroutine COSQF(n, x, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: x(*)
    real  :: wsave(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    real ::  tsqx
        !-----------------------------------------------
    !
    if (n - 2 >= 0) then
        if (n - 2 > 0) go to 103
        tsqx = sqrt(2.0)*X(2)
        x(2) = X(1) - tsqx
        x(1) = X(1) + tsqx
    end if
    return
103 continue
    call COSQF1 (n, x, wsave, wsave(n+1))
    return
end subroutine COSQF

subroutine COSQF1(n, x, w, xh)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: x(*)
    real , intent (in) :: w(*)
    real  :: xh(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: ns2, np2, k, kc, modn, i
    real :: xim1
    !-----------------------------------------------
    ns2 = (n + 1)/2
    np2 = n + 2
    do k = 2, ns2
        kc = np2 - k
        xh(k) = X(k) + X(kc)
        xh(kc) = X(k) - X(kc)
    end do
    modn = mod(n, 2)
    if (modn == 0) xh(ns2+1) = X(ns2+1) + X(ns2+1)
    do k = 2, ns2
        kc = np2 - k
        x(k) = W(k-1)*XH(kc) + W(kc-1)*XH(k)
        x(kc) = W(k-1)*XH(k) - W(kc-1)*XH(kc)
    end do
    if (modn == 0) x(ns2+1) = W(ns2)*XH(ns2+1)
    call RFFTF (n, x, xh)
    do i = 3, n, 2
        xim1 = X(i-1) - X(i)
        x(i) = X(i-1) + X(i)
        x(i-1) = xim1
    end do
    return
end subroutine COSQF1

subroutine COSQB(n, x, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: x(*)
    real  :: wsave(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    real :: x1
    !-----------------------------------------------
    !
    if (n - 2 <= 0) then
        if (n - 2 /= 0) then
            x(1) = 4.*X(1)
            return
        end if
        x1 = 4.*(X(1)+X(2))
        x(2) = 2.0*sqrt(2.0)*(X(1)-X(2))
        x(1) = x1
        return
    end if
    call COSQB1 (n, x, wsave, wsave(n+1))
    return
end subroutine COSQB

subroutine COSQB1(n, x, w, xh)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: x(*)
    real , intent (in) :: w(*)
    real  :: xh(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: ns2, np2, i, modn, k, kc
    real :: xim1
    !-----------------------------------------------
    ns2 = (n + 1)/2
    np2 = n + 2
    do i = 3, n, 2
        xim1 = X(i-1) + X(i)
        x(i) = X(i) - X(i-1)
        x(i-1) = xim1
    end do
    x(1) = X(1) + X(1)
    modn = mod(n, 2)
    if (modn == 0) x(n) = X(n) + X(n)
    call RFFTB (n, x, xh)
    do k = 2, ns2
        kc = np2 - k
        xh(k) = W(k-1)*X(kc) + W(kc-1)*X(k)
        xh(kc) = W(k-1)*X(k) - W(kc-1)*X(kc)
    end do
    if (modn == 0) x(ns2+1) = W(ns2)*(X(ns2+1)+X(ns2+1))
    do k = 2, ns2
        kc = np2 - k
        x(k) = XH(k) + XH(kc)
        x(kc) = XH(k) - XH(kc)
    end do
    x(1) = X(1) + X(1)
    return
end subroutine COSQB1


subroutine SINQI(n, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: wsave(*)
    !-----------------------------------------------
    !
    call COSQI (n, wsave)
    return
end subroutine SINQI


subroutine SINQF(n, x, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: x(*)
    real  :: wsave(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: ns2, k, kc
    real :: xhold
    !-----------------------------------------------
    !
    if (n == 1) return
    ns2 = n/2
    do k = 1, ns2
        kc = n - k
        xhold = X(k)
        x(k) = X(kc+1)
        x(kc+1) = xhold
    end do
    call COSQF (n, x, wsave)
    x(2:n:2) = -X(2:n:2)
    return
end subroutine SINQF

subroutine SINQB(n, x, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: x(*)
    real  :: wsave(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: ns2, k, kc
    real :: xhold
    !-----------------------------------------------
    !
    if (n <= 1) then
        x(1) = 4.*X(1)
        return
    end if
    ns2 = n/2
    x(2:n:2) = -X(2:n:2)
    call COSQB (n, x, wsave)
    do k = 1, ns2
        kc = n - k
        xhold = X(k)
        x(k) = X(kc+1)
        x(kc+1) = xhold
    end do
    return
end subroutine SINQB


subroutine CFFTI(n, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: wsave(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: iw1, iw2
    !-----------------------------------------------
    !
    if (n == 1) return
    iw1 = n + n + 1
    iw2 = iw1 + n + n
    call CFFTI1 (n, wsave(iw1), int(wsave(iw2:iw2), kind=ip) )
    return
end subroutine CFFTI


subroutine CFFTI1(n, wa, ifac)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: n
    !integer , intent (in out) :: ifac(*)
    integer :: ifac(*)
    real , intent (in out) :: wa(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer , dimension(4) :: ntryh
    integer :: nl, nf, j, ntry, nq, nr, i, l1, k1, ip_rename, ld, l2, ido &
        , idot, ipm, i1, ii
    real :: tpi, argh, fi, argld, arg
    !-----------------------------------------------
    data NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 3, 4, 2, 5/
    nl = n
    nf = 0
    j = 0
101 continue
    j = j + 1
    if (j - 4 <= 0) then
        ntry = NTRYH(j)
    else
        ntry = ntry + 2
    end if
104 continue
    nq = nl/ntry
    nr = nl - ntry*nq
    if (nr /= 0) go to 101
    nf = nf + 1
    ifac(nf+2) = ntry
    nl = nq
    if (ntry == 2) then
        if (nf /= 1) then
            ifac(nf+2:4:(-1)) = IFAC(nf+1:3:(-1))
            ifac(3) = 2
        end if
    end if
    if (nl /= 1) go to 104
    ifac(1) = n
    ifac(2) = nf
    tpi = 8.0*ATAN(1.0)
    argh = tpi/real(n)
    i = 2
    l1 = 1
    do k1 = 1, nf
        ip_rename = IFAC(k1+2)
        ld = 0
        l2 = l1*ip_rename
        ido = n/l2
        idot = ido + ido + 2
        ipm = ip_rename - 1
        do j = 1, ipm
            i1 = i
            wa(i-1) = 1.
            wa(i) = 0.
            ld = ld + l1
            fi = 0.
            argld = real(ld)*argh
            do ii = 4, idot, 2
                i = i + 2
                fi = fi + 1.
                arg = fi*argld
                wa(i-1) = COS(arg)
                wa(i) = SIN(arg)
            end do
            if (ip_rename <= 5) cycle
            wa(i1-1) = WA(i-1)
            wa(i1) = WA(i)
        end do
        l1 = l2
    end do
    return
end subroutine CFFTI1


subroutine CFFTB(n, c, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: c(*)
    real  :: wsave(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: iw1, iw2
    !-----------------------------------------------
    !
    if (n == 1) return
    iw1 = n + n + 1
    iw2 = iw1 + n + n
    call CFFTB1 (n, c, wsave, wsave(iw1), int(wsave(iw2:iw2), kind=ip) )
    return
end subroutine CFFTB

subroutine CFFTB1(n, c, ch, wa, ifac)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: n
    integer , intent (in) :: ifac(*)
    real  :: c(*)
    real  :: ch(*)
    real  :: wa(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer:: nf, na, l1, iw, k1, ip_rename, l2, ido, idot, idl1, ix2, ix3, ix4, nac, n2
    !-----------------------------------------------
    nf = IFAC(2)
    na = 0
    l1 = 1
    iw = 1
    do k1 = 1, nf
        ip_rename = IFAC(k1+2)
        l2 = ip_rename*l1
        ido = n/l2
        idot = ido + ido
        idl1 = idot*l1
        if (ip_rename == 4) then
            ix2 = iw + idot
            ix3 = ix2 + idot
            if (na == 0) then
                call PASSB4 (idot, l1, c, ch, WA(iw), WA(ix2), WA(ix3))
            else
                call PASSB4 (idot, l1, ch, c, WA(iw), WA(ix2), WA(ix3))
            end if
            na = 1 - na
        else
            if (ip_rename == 2) then
                if (na == 0) then
                    call PASSB2 (idot, l1, c, ch, WA(iw))
                else
                    call PASSB2 (idot, l1, ch, c, WA(iw))
                end if
                na = 1 - na
            else
                if (ip_rename == 3) then
                    ix2 = iw + idot
                    if (na == 0) then
                        call PASSB3 (idot, l1, c, ch, WA(iw), WA(ix2))
                    else
                        call PASSB3 (idot, l1, ch, c, WA(iw), WA(ix2))
                    end if
                    na = 1 - na
                else
                    if (ip_rename == 5) then
                        ix2 = iw + idot
                        ix3 = ix2 + idot
                        ix4 = ix3 + idot
                        if (na == 0) then
                            call PASSB5 (idot, l1, c, ch, WA(iw), WA(ix2), &
                                WA(ix3), WA(ix4))
                        else
                            call PASSB5 (idot, l1, ch, c, WA(iw), WA(ix2), &
                                WA(ix3), WA(ix4))
                        end if
                        na = 1 - na
                    else
                        if (na == 0) then
                            call PASSB (nac, idot, ip_rename, l1, idl1, c, c, c, ch &
                                , ch, WA(iw))
                        else
                            call PASSB (nac, idot, ip_rename, l1, idl1, ch, ch, ch &
                                , c, c, WA(iw))
                        end if
                        if (nac /= 0) na = 1 - na
                    end if
                end if
            end if
        end if
        l1 = l2
        iw = iw + (ip_rename - 1)*idot
    end do
    if (na == 0) return
    n2 = n + n
    c(:n2) = CH(:n2)
    return
end subroutine CFFTB1


subroutine PASSB2(ido, l1, cc, ch, wa1)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, 2, l1)
    real , intent (out) :: ch(ido, l1, 2)
    real , intent (in) :: wa1(1)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, i
    real :: tr2, ti2
    !-----------------------------------------------
    if (ido <= 2) then
        ch(1, :, 1) = CC(1, 1, :) + CC(1, 2, :)
        ch(1, :, 2) = CC(1, 1, :) - CC(1, 2, :)
        ch(2, :, 1) = CC(2, 1, :) + CC(2, 2, :)
        ch(2, :, 2) = CC(2, 1, :) - CC(2, 2, :)
        return
    end if
    do k = 1, l1
        do i = 2, ido, 2
            ch(i-1, k, 1) = CC(i-1, 1, k) + CC(i-1, 2, k)
            tr2 = CC(i-1, 1, k) - CC(i-1, 2, k)
            ch(i, k, 1) = CC(i, 1, k) + CC(i, 2, k)
            ti2 = CC(i, 1, k) - CC(i, 2, k)
            ch(i, k, 2) = WA1(i-1)*ti2 + WA1(i)*tr2
            ch(i-1, k, 2) = WA1(i-1)*tr2 - WA1(i)*ti2
        end do
    end do
    return
end subroutine PASSB2


subroutine PASSB3(ido, l1, cc, ch, wa1, wa2)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, 3, l1)
    real , intent (out) :: ch(ido, l1, 3)
    real , intent (in) :: wa1(*)
    real , intent (in) :: wa2(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, i
    real::taur, taui, tr2, cr2, ti2, ci2, cr3, ci3, dr2, dr3, di2, di3
        !-----------------------------------------------
    parameter ( TAUI = 0.866025403784439 )
    parameter ( TAUR = -.5 )
    if (ido == 2) then
        do k = 1, l1
            tr2 = CC(1, 2, k) + CC(1, 3, k)
            cr2 = CC(1, 1, k) + taur*tr2
            ch(1, k, 1) = CC(1, 1, k) + tr2
            ti2 = CC(2, 2, k) + CC(2, 3, k)
            ci2 = CC(2, 1, k) + taur*ti2
            ch(2, k, 1) = CC(2, 1, k) + ti2
            cr3 = taui*(CC(1, 2, k)-CC(1, 3, k))
            ci3 = taui*(CC(2, 2, k)-CC(2, 3, k))
            ch(1, k, 2) = cr2 - ci3
            ch(1, k, 3) = cr2 + ci3
            ch(2, k, 2) = ci2 + cr3
            ch(2, k, 3) = ci2 - cr3
        end do
        return
    end if
    do k = 1, l1
        do i = 2, ido, 2
            tr2 = CC(i-1, 2, k) + CC(i-1, 3, k)
            cr2 = CC(i-1, 1, k) + taur*tr2
            ch(i-1, k, 1) = CC(i-1, 1, k) + tr2
            ti2 = CC(i, 2, k) + CC(i, 3, k)
            ci2 = CC(i, 1, k) + taur*ti2
            ch(i, k, 1) = CC(i, 1, k) + ti2
            cr3 = taui*(CC(i-1, 2, k)-CC(i-1, 3, k))
            ci3 = taui*(CC(i, 2, k)-CC(i, 3, k))
            dr2 = cr2 - ci3
            dr3 = cr2 + ci3
            di2 = ci2 + cr3
            di3 = ci2 - cr3
            ch(i, k, 2) = WA1(i-1)*di2 + WA1(i)*dr2
            ch(i-1, k, 2) = WA1(i-1)*dr2 - WA1(i)*di2
            ch(i, k, 3) = WA2(i-1)*di3 + WA2(i)*dr3
            ch(i-1, k, 3) = WA2(i-1)*dr3 - WA2(i)*di3
        end do
    end do
    return
end subroutine PASSB3


subroutine PASSB4(ido, l1, cc, ch, wa1, wa2, wa3)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, 4, l1)
    real , intent (out) :: ch(ido, l1, 4)
    real , intent (in) :: wa1(*)
    real , intent (in) :: wa2(*)
    real , intent (in) :: wa3(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, i
    real::ti1, ti2, tr4, ti3, tr1, tr2, ti4, tr3, cr3, ci3, cr2, cr4, ci2, ci4
    !-----------------------------------------------
    if (ido == 2) then
        do k = 1, l1
            ti1 = CC(2, 1, k) - CC(2, 3, k)
            ti2 = CC(2, 1, k) + CC(2, 3, k)
            tr4 = CC(2, 4, k) - CC(2, 2, k)
            ti3 = CC(2, 2, k) + CC(2, 4, k)
            tr1 = CC(1, 1, k) - CC(1, 3, k)
            tr2 = CC(1, 1, k) + CC(1, 3, k)
            ti4 = CC(1, 2, k) - CC(1, 4, k)
            tr3 = CC(1, 2, k) + CC(1, 4, k)
            ch(1, k, 1) = tr2 + tr3
            ch(1, k, 3) = tr2 - tr3
            ch(2, k, 1) = ti2 + ti3
            ch(2, k, 3) = ti2 - ti3
            ch(1, k, 2) = tr1 + tr4
            ch(1, k, 4) = tr1 - tr4
            ch(2, k, 2) = ti1 + ti4
            ch(2, k, 4) = ti1 - ti4
        end do
        return
    end if
    do k = 1, l1
        do i = 2, ido, 2
            ti1 = CC(i, 1, k) - CC(i, 3, k)
            ti2 = CC(i, 1, k) + CC(i, 3, k)
            ti3 = CC(i, 2, k) + CC(i, 4, k)
            tr4 = CC(i, 4, k) - CC(i, 2, k)
            tr1 = CC(i-1, 1, k) - CC(i-1, 3, k)
            tr2 = CC(i-1, 1, k) + CC(i-1, 3, k)
            ti4 = CC(i-1, 2, k) - CC(i-1, 4, k)
            tr3 = CC(i-1, 2, k) + CC(i-1, 4, k)
            ch(i-1, k, 1) = tr2 + tr3
            cr3 = tr2 - tr3
            ch(i, k, 1) = ti2 + ti3
            ci3 = ti2 - ti3
            cr2 = tr1 + tr4
            cr4 = tr1 - tr4
            ci2 = ti1 + ti4
            ci4 = ti1 - ti4
            ch(i-1, k, 2) = WA1(i-1)*cr2 - WA1(i)*ci2
            ch(i, k, 2) = WA1(i-1)*ci2 + WA1(i)*cr2
            ch(i-1, k, 3) = WA2(i-1)*cr3 - WA2(i)*ci3
            ch(i, k, 3) = WA2(i-1)*ci3 + WA2(i)*cr3
            ch(i-1, k, 4) = WA3(i-1)*cr4 - WA3(i)*ci4
            ch(i, k, 4) = WA3(i-1)*ci4 + WA3(i)*cr4
        end do
    end do
    return
end subroutine PASSB4


subroutine PASSB5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, 5, l1)
    real , intent (out) :: ch(ido, l1, 5)
    real , intent (in) :: wa1(*)
    real , intent (in) :: wa2(*)
    real , intent (in) :: wa3(*)
    real , intent (in) :: wa4(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, i
    real :: tr11, ti11, tr12, ti12, ti5, ti2, ti4, ti3, tr5, tr2, tr4 &
        , tr3, cr2, ci2, cr3, ci3, cr5, ci5, cr4, ci4, dr3, dr4, di3, &
        di4, dr5, dr2, di5, di2
        !-----------------------------------------------
    parameter ( TI12 = 0.587785252292473 )
    parameter ( TR12 = -.809016994374947 )
    parameter ( TI11 = 0.951056516295154 )
    parameter ( TR11 = 0.309016994374947 )
    if (ido == 2) then
        do k = 1, l1
            ti5 = CC(2, 2, k) - CC(2, 5, k)
            ti2 = CC(2, 2, k) + CC(2, 5, k)
            ti4 = CC(2, 3, k) - CC(2, 4, k)
            ti3 = CC(2, 3, k) + CC(2, 4, k)
            tr5 = CC(1, 2, k) - CC(1, 5, k)
            tr2 = CC(1, 2, k) + CC(1, 5, k)
            tr4 = CC(1, 3, k) - CC(1, 4, k)
            tr3 = CC(1, 3, k) + CC(1, 4, k)
            ch(1, k, 1) = CC(1, 1, k) + tr2 + tr3
            ch(2, k, 1) = CC(2, 1, k) + ti2 + ti3
            cr2 = CC(1, 1, k) + tr11*tr2 + tr12*tr3
            ci2 = CC(2, 1, k) + tr11*ti2 + tr12*ti3
            cr3 = CC(1, 1, k) + tr12*tr2 + tr11*tr3
            ci3 = CC(2, 1, k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            ch(1, k, 2) = cr2 - ci5
            ch(1, k, 5) = cr2 + ci5
            ch(2, k, 2) = ci2 + cr5
            ch(2, k, 3) = ci3 + cr4
            ch(1, k, 3) = cr3 - ci4
            ch(1, k, 4) = cr3 + ci4
            ch(2, k, 4) = ci3 - cr4
            ch(2, k, 5) = ci2 - cr5
        end do
        return
    end if
    do k = 1, l1
        do i = 2, ido, 2
            ti5 = CC(i, 2, k) - CC(i, 5, k)
            ti2 = CC(i, 2, k) + CC(i, 5, k)
            ti4 = CC(i, 3, k) - CC(i, 4, k)
            ti3 = CC(i, 3, k) + CC(i, 4, k)
            tr5 = CC(i-1, 2, k) - CC(i-1, 5, k)
            tr2 = CC(i-1, 2, k) + CC(i-1, 5, k)
            tr4 = CC(i-1, 3, k) - CC(i-1, 4, k)
            tr3 = CC(i-1, 3, k) + CC(i-1, 4, k)
            ch(i-1, k, 1) = CC(i-1, 1, k) + tr2 + tr3
            ch(i, k, 1) = CC(i, 1, k) + ti2 + ti3
            cr2 = CC(i-1, 1, k) + tr11*tr2 + tr12*tr3
            ci2 = CC(i, 1, k) + tr11*ti2 + tr12*ti3
            cr3 = CC(i-1, 1, k) + tr12*tr2 + tr11*tr3
            ci3 = CC(i, 1, k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            dr3 = cr3 - ci4
            dr4 = cr3 + ci4
            di3 = ci3 + cr4
            di4 = ci3 - cr4
            dr5 = cr2 + ci5
            dr2 = cr2 - ci5
            di5 = ci2 - cr5
            di2 = ci2 + cr5
            ch(i-1, k, 2) = WA1(i-1)*dr2 - WA1(i)*di2
            ch(i, k, 2) = WA1(i-1)*di2 + WA1(i)*dr2
            ch(i-1, k, 3) = WA2(i-1)*dr3 - WA2(i)*di3
            ch(i, k, 3) = WA2(i-1)*di3 + WA2(i)*dr3
            ch(i-1, k, 4) = WA3(i-1)*dr4 - WA3(i)*di4
            ch(i, k, 4) = WA3(i-1)*di4 + WA3(i)*dr4
            ch(i-1, k, 5) = WA4(i-1)*dr5 - WA4(i)*di5
            ch(i, k, 5) = WA4(i-1)*di5 + WA4(i)*dr5
        end do
    end do

end subroutine PASSB5


subroutine PASSB(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (out) :: nac
    integer , intent (in) :: ido
    integer , intent (in) :: ip
    integer , intent (in) :: l1
    integer , intent (in) :: idl1
    real , intent (in) :: cc(ido, ip, l1)
    real , intent (out) :: c1(ido, l1, ip)
    real , intent (in out) :: c2(idl1, ip)
    real , intent (in out) :: ch(ido, l1, ip)
    real , intent (in out) :: ch2(idl1, ip)
    real , intent (in) :: wa(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: idot, nt, ipp2, ipph, idp, j, jc, k, i, idl, inc, l, lc, idlj, idij, idj
    real :: war, wai
    !-----------------------------------------------
    idot = ido/2
    nt = ip*idl1
    ipp2 = ip + 2
    ipph = (ip + 1)/2
    idp = ip*ido
    !
    if (ido >= l1) then
        do j = 2, ipph
            jc = ipp2 - j
            ch(:, :, j) = CC(:, j, :) + CC(:, jc, :)
            ch(:, :, jc) = CC(:, j, :) - CC(:, jc, :)
        end do
        ch(:, :, 1) = CC(:, 1, :)
    else
        do j = 2, ipph
            jc = ipp2 - j
            ch(:, :, j) = CC(:, j, :) + CC(:, jc, :)
            ch(:, :, jc) = CC(:, j, :) - CC(:, jc, :)
        end do
        ch(:, :, 1) = CC(:, 1, :)
    end if
    idl = 2 - ido
    inc = 0
    do l = 2, ipph
        lc = ipp2 - l
        idl = idl + ido
        c2(:, l) = CH2(:, 1) + WA(idl-1)*CH2(:, 2)
        c2(:, lc) = WA(idl)*CH2(:, ip)
        idlj = idl
        inc = inc + ido
        do j = 3, ipph
            jc = ipp2 - j
            idlj = idlj + inc
            if (idlj > idp) idlj = idlj - idp
            war = WA(idlj-1)
            wai = WA(idlj)
            c2(:, l) = C2(:, l) + war*CH2(:, j)
            c2(:, lc) = C2(:, lc) + wai*CH2(:, jc)
        end do
    end do
    do j = 2, ipph
        ch2(:, 1) = CH2(:, 1) + CH2(:, j)
    end do
    do j = 2, ipph
        jc = ipp2 - j
        ch2(:idl1-1:2, j) = C2(:idl1-1:2, j) - C2(2:idl1:2, jc)
        ch2(:idl1-1:2, jc) = C2(:idl1-1:2, j) + C2(2:idl1:2, jc)
        ch2(2:idl1:2, j) = C2(2:idl1:2, j) + C2(:idl1-1:2, jc)
        ch2(2:idl1:2, jc) = C2(2:idl1:2, j) - C2(:idl1-1:2, jc)
    end do
    nac = 1
    if (ido == 2) return
    c1(1, :, 2:ip) = CH(1, :, 2:ip)
    c1(2, :, 2:ip) = CH(2, :, 2:ip)
    if (idot <= l1) then
        idij = 0
        do j = 2, ip
            idij = idij + 2
            do i = 4, ido, 2
                idij = idij + 2
                c1(i-1, :, j) = WA(idij-1)*CH(i-1, :, j) - WA(idij)*CH(i, :, j)
                c1(i, :, j) = WA(idij-1)*CH(i, :, j) + WA(idij)*CH(i-1, :, j)
            end do
        end do
        return
    end if
    idj = 2 - ido
    do j = 2, ip
        idj = idj + ido
        do k = 1, l1
            idij = idj
            c1(3:ido-1:2, k, j) = WA(idij+1:ido-3+idij:2)*CH(3:ido-1:2, k, j &
                ) - WA(idij+2:ido-2+idij:2)*CH(4:ido:2, k, j)
            c1(4:ido:2, k, j) = WA(idij+1:ido-3+idij:2)*CH(4:ido:2, k, j) + &
                WA(idij+2:ido-2+idij:2)*CH(3:ido-1:2, k, j)
        end do
    end do

end subroutine PASSB

subroutine CFFTF(n, c, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: c(*)
    real  :: wsave(*)
    !-----------------------------------------------

    if (n == 1) return

    associate( iw1 => n + n + 1 )
        associate( iw2 => iw1 + n + n)
            call CFFTF1 (n, c, wsave, wsave(iw1), int(wsave(iw2:iw2), kind=ip) )
        end associate
    end associate

end subroutine CFFTF

subroutine CFFTF1(n, c, ch, wa, ifac)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: n
    integer , intent (in) :: ifac(*)
    real  :: c(*)
    real  :: ch(*)
    real  :: wa(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer:: nf, na, l1, iw, k1, ip_rename, l2, ido, idot, idl1, ix2, ix3, ix4, nac, n2
    !-----------------------------------------------
    nf = IFAC(2)
    na = 0
    l1 = 1
    iw = 1
    do k1 = 1, nf
        ip_rename = IFAC(k1+2)
        l2 = ip_rename*l1
        ido = n/l2
        idot = ido + ido
        idl1 = idot*l1
        if (ip_rename == 4) then
            ix2 = iw + idot
            ix3 = ix2 + idot
            if (na == 0) then
                call PASSF4 (idot, l1, c, ch, WA(iw), WA(ix2), WA(ix3))
            else
                call PASSF4 (idot, l1, ch, c, WA(iw), WA(ix2), WA(ix3))
            end if
            na = 1 - na
        else
            if (ip_rename == 2) then
                if (na == 0) then
                    call PASSF2 (idot, l1, c, ch, WA(iw))
                else
                    call PASSF2 (idot, l1, ch, c, WA(iw))
                end if
                na = 1 - na
            else
                if (ip_rename == 3) then
                    ix2 = iw + idot
                    if (na == 0) then
                        call PASSF3 (idot, l1, c, ch, WA(iw), WA(ix2))
                    else
                        call PASSF3 (idot, l1, ch, c, WA(iw), WA(ix2))
                    end if
                    na = 1 - na
                else
                    if (ip_rename == 5) then
                        ix2 = iw + idot
                        ix3 = ix2 + idot
                        ix4 = ix3 + idot
                        if (na == 0) then
                            call PASSF5 (idot, l1, c, ch, WA(iw), WA(ix2), &
                                WA(ix3), WA(ix4))
                        else
                            call PASSF5 (idot, l1, ch, c, WA(iw), WA(ix2), &
                                WA(ix3), WA(ix4))
                        end if
                        na = 1 - na
                    else
                        if (na == 0) then
                            call PASSF (nac, idot, ip_rename, l1, idl1, c, c, c, ch &
                                , ch, WA(iw))
                        else
                            call PASSF (nac, idot, ip_rename, l1, idl1, ch, ch, ch &
                                , c, c, WA(iw))
                        end if
                        if (nac /= 0) na = 1 - na
                    end if
                end if
            end if
        end if
        l1 = l2
        iw = iw + (ip_rename - 1)*idot
    end do
    if (na == 0) return
    n2 = n + n
    c(:n2) = CH(:n2)

end subroutine CFFTF1

subroutine PASSF2(ido, l1, cc, ch, wa1)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, 2, l1)
    real , intent (out) :: ch(ido, l1, 2)
    real , intent (in) :: wa1(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, i
    real :: tr2, ti2
    !-----------------------------------------------
    if (ido <= 2) then
        ch(1, :, 1) = CC(1, 1, :) + CC(1, 2, :)
        ch(1, :, 2) = CC(1, 1, :) - CC(1, 2, :)
        ch(2, :, 1) = CC(2, 1, :) + CC(2, 2, :)
        ch(2, :, 2) = CC(2, 1, :) - CC(2, 2, :)
        return
    end if
    do k = 1, l1
        do i = 2, ido, 2
            ch(i-1, k, 1) = CC(i-1, 1, k) + CC(i-1, 2, k)
            tr2 = CC(i-1, 1, k) - CC(i-1, 2, k)
            ch(i, k, 1) = CC(i, 1, k) + CC(i, 2, k)
            ti2 = CC(i, 1, k) - CC(i, 2, k)
            ch(i, k, 2) = WA1(i-1)*ti2 - WA1(i)*tr2
            ch(i-1, k, 2) = WA1(i-1)*tr2 + WA1(i)*ti2
        end do
    end do

end subroutine PASSF2

subroutine PASSF3(ido, l1, cc, ch, wa1, wa2)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, 3, l1)
    real , intent (out) :: ch(ido, l1, 3)
    real , intent (in) :: wa1(*)
    real , intent (in) :: wa2(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, i
    real (wp), parameter :: taui = -sqrt( 3.0_wp )/2
    real (wp), parameter :: taur = -0.5_wp
    real (wp)            :: tr2, cr2, ti2, ci2, cr3, ci3, dr2, dr3, di2, di3
    !-----------------------------------------------

    if (ido == 2) then
        do k = 1, l1
            tr2 = CC(1, 2, k) + CC(1, 3, k)
            cr2 = CC(1, 1, k) + taur*tr2
            ch(1, k, 1) = CC(1, 1, k) + tr2
            ti2 = CC(2, 2, k) + CC(2, 3, k)
            ci2 = CC(2, 1, k) + taur*ti2
            ch(2, k, 1) = CC(2, 1, k) + ti2
            cr3 = taui*(CC(1, 2, k)-CC(1, 3, k))
            ci3 = taui*(CC(2, 2, k)-CC(2, 3, k))
            ch(1, k, 2) = cr2 - ci3
            ch(1, k, 3) = cr2 + ci3
            ch(2, k, 2) = ci2 + cr3
            ch(2, k, 3) = ci2 - cr3
        end do
        return
    end if
    do k = 1, l1
        do i = 2, ido, 2
            tr2 = CC(i-1, 2, k) + CC(i-1, 3, k)
            cr2 = CC(i-1, 1, k) + taur*tr2
            ch(i-1, k, 1) = CC(i-1, 1, k) + tr2
            ti2 = CC(i, 2, k) + CC(i, 3, k)
            ci2 = CC(i, 1, k) + taur*ti2
            ch(i, k, 1) = CC(i, 1, k) + ti2
            cr3 = taui*(CC(i-1, 2, k)-CC(i-1, 3, k))
            ci3 = taui*(CC(i, 2, k)-CC(i, 3, k))
            dr2 = cr2 - ci3
            dr3 = cr2 + ci3
            di2 = ci2 + cr3
            di3 = ci2 - cr3
            ch(i, k, 2) = WA1(i-1)*di2 - WA1(i)*dr2
            ch(i-1, k, 2) = WA1(i-1)*dr2 + WA1(i)*di2
            ch(i, k, 3) = WA2(i-1)*di3 - WA2(i)*dr3
            ch(i-1, k, 3) = WA2(i-1)*dr3 + WA2(i)*di3
        end do
    end do
    return
end subroutine PASSF3

subroutine PASSF4(ido, l1, cc, ch, wa1, wa2, wa3)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, 4, l1)
    real , intent (out) :: ch(ido, l1, 4)
    real , intent (in) :: wa1(*)
    real , intent (in) :: wa2(*)
    real , intent (in) :: wa3(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, i
    real::ti1, ti2, tr4, ti3, tr1, tr2, ti4, tr3, cr3, ci3, cr2, cr4, ci2, ci4
    !-----------------------------------------------
    if (ido == 2) then
        do k = 1, l1
            ti1 = CC(2, 1, k) - CC(2, 3, k)
            ti2 = CC(2, 1, k) + CC(2, 3, k)
            tr4 = CC(2, 2, k) - CC(2, 4, k)
            ti3 = CC(2, 2, k) + CC(2, 4, k)
            tr1 = CC(1, 1, k) - CC(1, 3, k)
            tr2 = CC(1, 1, k) + CC(1, 3, k)
            ti4 = CC(1, 4, k) - CC(1, 2, k)
            tr3 = CC(1, 2, k) + CC(1, 4, k)
            ch(1, k, 1) = tr2 + tr3
            ch(1, k, 3) = tr2 - tr3
            ch(2, k, 1) = ti2 + ti3
            ch(2, k, 3) = ti2 - ti3
            ch(1, k, 2) = tr1 + tr4
            ch(1, k, 4) = tr1 - tr4
            ch(2, k, 2) = ti1 + ti4
            ch(2, k, 4) = ti1 - ti4
        end do
        return
    end if
    do k = 1, l1
        do i = 2, ido, 2
            ti1 = CC(i, 1, k) - CC(i, 3, k)
            ti2 = CC(i, 1, k) + CC(i, 3, k)
            ti3 = CC(i, 2, k) + CC(i, 4, k)
            tr4 = CC(i, 2, k) - CC(i, 4, k)
            tr1 = CC(i-1, 1, k) - CC(i-1, 3, k)
            tr2 = CC(i-1, 1, k) + CC(i-1, 3, k)
            ti4 = CC(i-1, 4, k) - CC(i-1, 2, k)
            tr3 = CC(i-1, 2, k) + CC(i-1, 4, k)
            ch(i-1, k, 1) = tr2 + tr3
            cr3 = tr2 - tr3
            ch(i, k, 1) = ti2 + ti3
            ci3 = ti2 - ti3
            cr2 = tr1 + tr4
            cr4 = tr1 - tr4
            ci2 = ti1 + ti4
            ci4 = ti1 - ti4
            ch(i-1, k, 2) = WA1(i-1)*cr2 + WA1(i)*ci2
            ch(i, k, 2) = WA1(i-1)*ci2 - WA1(i)*cr2
            ch(i-1, k, 3) = WA2(i-1)*cr3 + WA2(i)*ci3
            ch(i, k, 3) = WA2(i-1)*ci3 - WA2(i)*cr3
            ch(i-1, k, 4) = WA3(i-1)*cr4 + WA3(i)*ci4
            ch(i, k, 4) = WA3(i-1)*ci4 - WA3(i)*cr4
        end do
    end do

end subroutine PASSF4

subroutine PASSF5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, 5, l1)
    real , intent (out) :: ch(ido, l1, 5)
    real , intent (in) :: wa1(*)
    real , intent (in) :: wa2(*)
    real , intent (in) :: wa3(*)
    real , intent (in) :: wa4(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, i
    real :: tr11, ti11, tr12, ti12, ti5, ti2, ti4, ti3, tr5, tr2, tr4 &
        , tr3, cr2, ci2, cr3, ci3, cr5, ci5, cr4, ci4, dr3, dr4, di3, &
        di4, dr5, dr2, di5, di2
        !-----------------------------------------------
    parameter ( TI12 = - 0.587785252292473 )
    parameter ( TR12 = -.809016994374947 )
    parameter ( TI11 = -.951056516295154 )
    parameter ( TR11 = 0.309016994374947 )
    if (ido == 2) then
        do k = 1, l1
            ti5 = CC(2, 2, k) - CC(2, 5, k)
            ti2 = CC(2, 2, k) + CC(2, 5, k)
            ti4 = CC(2, 3, k) - CC(2, 4, k)
            ti3 = CC(2, 3, k) + CC(2, 4, k)
            tr5 = CC(1, 2, k) - CC(1, 5, k)
            tr2 = CC(1, 2, k) + CC(1, 5, k)
            tr4 = CC(1, 3, k) - CC(1, 4, k)
            tr3 = CC(1, 3, k) + CC(1, 4, k)
            ch(1, k, 1) = CC(1, 1, k) + tr2 + tr3
            ch(2, k, 1) = CC(2, 1, k) + ti2 + ti3
            cr2 = CC(1, 1, k) + tr11*tr2 + tr12*tr3
            ci2 = CC(2, 1, k) + tr11*ti2 + tr12*ti3
            cr3 = CC(1, 1, k) + tr12*tr2 + tr11*tr3
            ci3 = CC(2, 1, k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            ch(1, k, 2) = cr2 - ci5
            ch(1, k, 5) = cr2 + ci5
            ch(2, k, 2) = ci2 + cr5
            ch(2, k, 3) = ci3 + cr4
            ch(1, k, 3) = cr3 - ci4
            ch(1, k, 4) = cr3 + ci4
            ch(2, k, 4) = ci3 - cr4
            ch(2, k, 5) = ci2 - cr5
        end do
        return
    end if
    do k = 1, l1
        do i = 2, ido, 2
            ti5 = CC(i, 2, k) - CC(i, 5, k)
            ti2 = CC(i, 2, k) + CC(i, 5, k)
            ti4 = CC(i, 3, k) - CC(i, 4, k)
            ti3 = CC(i, 3, k) + CC(i, 4, k)
            tr5 = CC(i-1, 2, k) - CC(i-1, 5, k)
            tr2 = CC(i-1, 2, k) + CC(i-1, 5, k)
            tr4 = CC(i-1, 3, k) - CC(i-1, 4, k)
            tr3 = CC(i-1, 3, k) + CC(i-1, 4, k)
            ch(i-1, k, 1) = CC(i-1, 1, k) + tr2 + tr3
            ch(i, k, 1) = CC(i, 1, k) + ti2 + ti3
            cr2 = CC(i-1, 1, k) + tr11*tr2 + tr12*tr3
            ci2 = CC(i, 1, k) + tr11*ti2 + tr12*ti3
            cr3 = CC(i-1, 1, k) + tr12*tr2 + tr11*tr3
            ci3 = CC(i, 1, k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            dr3 = cr3 - ci4
            dr4 = cr3 + ci4
            di3 = ci3 + cr4
            di4 = ci3 - cr4
            dr5 = cr2 + ci5
            dr2 = cr2 - ci5
            di5 = ci2 - cr5
            di2 = ci2 + cr5
            ch(i-1, k, 2) = WA1(i-1)*dr2 + WA1(i)*di2
            ch(i, k, 2) = WA1(i-1)*di2 - WA1(i)*dr2
            ch(i-1, k, 3) = WA2(i-1)*dr3 + WA2(i)*di3
            ch(i, k, 3) = WA2(i-1)*di3 - WA2(i)*dr3
            ch(i-1, k, 4) = WA3(i-1)*dr4 + WA3(i)*di4
            ch(i, k, 4) = WA3(i-1)*di4 - WA3(i)*dr4
            ch(i-1, k, 5) = WA4(i-1)*dr5 + WA4(i)*di5
            ch(i, k, 5) = WA4(i-1)*di5 - WA4(i)*dr5
        end do
    end do

end subroutine PASSF5

subroutine PASSF(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (out) :: nac
    integer , intent (in) :: ido
    integer , intent (in) :: ip
    integer , intent (in) :: l1
    integer , intent (in) :: idl1
    real , intent (in) :: cc(ido, ip, l1)
    real , intent (out) :: c1(ido, l1, ip)
    real , intent (in out) :: c2(idl1, ip)
    real , intent (in out) :: ch(ido, l1, ip)
    real , intent (in out) :: ch2(idl1, ip)
    real , intent (in) :: wa(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: idot, nt, ipp2, ipph, idp, j, jc, k, i, idl, inc, l, lc, idlj, idij, idj
    real :: war, wai
    !-----------------------------------------------
    idot = ido/2
    nt = ip*idl1
    ipp2 = ip + 2
    ipph = (ip + 1)/2
    idp = ip*ido
    !
    if (ido >= l1) then
        do j = 2, ipph
            jc = ipp2 - j
            ch(:, :, j) = CC(:, j, :) + CC(:, jc, :)
            ch(:, :, jc) = CC(:, j, :) - CC(:, jc, :)
        end do
        ch(:, :, 1) = CC(:, 1, :)
    else
        do j = 2, ipph
            jc = ipp2 - j
            ch(:, :, j) = CC(:, j, :) + CC(:, jc, :)
            ch(:, :, jc) = CC(:, j, :) - CC(:, jc, :)
        end do
        ch(:, :, 1) = CC(:, 1, :)
    end if
    idl = 2 - ido
    inc = 0
    do l = 2, ipph
        lc = ipp2 - l
        idl = idl + ido
        c2(:, l) = CH2(:, 1) + WA(idl-1)*CH2(:, 2)
        c2(:, lc) = -WA(idl)*CH2(:, ip)
        idlj = idl
        inc = inc + ido
        do j = 3, ipph
            jc = ipp2 - j
            idlj = idlj + inc
            if (idlj > idp) idlj = idlj - idp
            war = WA(idlj-1)
            wai = WA(idlj)
            c2(:, l) = C2(:, l) + war*CH2(:, j)
            c2(:, lc) = C2(:, lc) - wai*CH2(:, jc)
        end do
    end do
    do j = 2, ipph
        ch2(:, 1) = CH2(:, 1) + CH2(:, j)
    end do
    do j = 2, ipph
        jc = ipp2 - j
        ch2(:idl1-1:2, j) = C2(:idl1-1:2, j) - C2(2:idl1:2, jc)
        ch2(:idl1-1:2, jc) = C2(:idl1-1:2, j) + C2(2:idl1:2, jc)
        ch2(2:idl1:2, j) = C2(2:idl1:2, j) + C2(:idl1-1:2, jc)
        ch2(2:idl1:2, jc) = C2(2:idl1:2, j) - C2(:idl1-1:2, jc)
    end do
    nac = 1
    if (ido == 2) return
    nac = 0
    c2(:, 1) = CH2(:, 1)
    c1(1, :, 2:ip) = CH(1, :, 2:ip)
    c1(2, :, 2:ip) = CH(2, :, 2:ip)
    if (idot <= l1) then
        idij = 0
        do j = 2, ip
            idij = idij + 2
            do i = 4, ido, 2
                idij = idij + 2
                c1(i-1, :, j) = WA(idij-1)*CH(i-1, :, j) + WA(idij)*CH(i, :, j)
                c1(i, :, j) = WA(idij-1)*CH(i, :, j) - WA(idij)*CH(i-1, :, j)
            end do
        end do
        return
    end if
    idj = 2 - ido
    do j = 2, ip
        idj = idj + ido
        do k = 1, l1
            idij = idj
            c1(3:ido-1:2, k, j) = WA(idij+1:ido-3+idij:2)*CH(3:ido-1:2, k, j &
                ) + WA(idij+2:ido-2+idij:2)*CH(4:ido:2, k, j)
            c1(4:ido:2, k, j) = WA(idij+1:ido-3+idij:2)*CH(4:ido:2, k, j) - &
                WA(idij+2:ido-2+idij:2)*CH(3:ido-1:2, k, j)
        end do
    end do

end subroutine PASSF

subroutine RFFTI(n, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: wsave(*)
    !-----------------------------------------------
    !
    if (n == 1) return

    call RFFTI1 (n, wsave(n+1), int(wsave(2*n+1:2*n+1), kind=ip) )

end subroutine RFFTI

subroutine RFFTI1(n, wa, ifac)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: n
    !integer , intent (in out) :: ifac(*)
    integer :: ifac(*)
    real , intent (out) :: wa(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer , dimension(4) :: ntryh
    integer :: nl, nf, j, ntry, nq, nr, i, is, nfm1, l1, k1, ip_rename, &
        ld, l2, ido, ipm, ii
    real :: tpi, argh, argld, fi, arg
    !-----------------------------------------------
    data NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 4, 2, 3, 5/
    nl = n
    nf = 0
    j = 0
101 continue
    j = j + 1
    if (j - 4 <= 0) then
        ntry = NTRYH(j)
    else
        ntry = ntry + 2
    end if
104 continue
    nq = nl/ntry
    nr = nl - ntry*nq
    if (nr /= 0) go to 101
    nf = nf + 1
    ifac(nf+2) = ntry
    nl = nq
    if (ntry == 2) then
        if (nf /= 1) then
            ifac(nf+2:4:(-1)) = IFAC(nf+1:3:(-1))
            ifac(3) = 2
        end if
    end if
    if (nl /= 1) go to 104
    ifac(1) = n
    ifac(2) = nf
    tpi = 8.0*ATAN(1.0)
    argh = tpi/real(n)
    is = 0
    nfm1 = nf - 1
    l1 = 1
    if (nfm1 == 0) return
    do k1 = 1, nfm1
        ip_rename = IFAC(k1+2)
        ld = 0
        l2 = l1*ip_rename
        ido = n/l2
        ipm = ip_rename - 1
        do j = 1, ipm
            ld = ld + l1
            i = is
            argld = real(ld)*argh
            fi = 0.
            do ii = 3, ido, 2
                i = i + 2
                fi = fi + 1.
                arg = fi*argld
                wa(i-1) = COS(arg)
                wa(i) = SIN(arg)
            end do
            is = is + ido
        end do
        l1 = l2
    end do

end subroutine RFFTI1

subroutine RFFTB(n, r, wsave)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real  :: r(*)
    real  :: wsave(*)
    !-----------------------------------------------
    !
    if (n == 1) return

    call RFFTB1 (n, r, wsave, wsave(n+1), int(wsave(2*n+1:2*n+1), kind=ip) )

end subroutine RFFTB

subroutine RFFTB1(n, c, ch, wa, ifac)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: n
    integer , intent (in) :: ifac(*)
    real  :: c(*)
    real  :: ch(*)
    real  :: wa(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: nf, na, l1, iw, k1, ip_rename, l2, ido, idl1, ix2, ix3, ix4
    !-----------------------------------------------
    nf = IFAC(2)
    na = 0
    l1 = 1
    iw = 1
    do k1 = 1, nf
        ip_rename = IFAC(k1+2)
        l2 = ip_rename*l1
        ido = n/l2
        idl1 = ido*l1
        if (ip_rename == 4) then
            ix2 = iw + ido
            ix3 = ix2 + ido
            if (na == 0) then
                call RADB4 (ido, l1, c, ch, WA(iw), WA(ix2), WA(ix3))
            else
                call RADB4 (ido, l1, ch, c, WA(iw), WA(ix2), WA(ix3))
            end if
            na = 1 - na
        else
            if (ip_rename == 2) then
                if (na == 0) then
                    call RADB2 (ido, l1, c, ch, WA(iw))
                else
                    call RADB2 (ido, l1, ch, c, WA(iw))
                end if
                na = 1 - na
            else
                if (ip_rename == 3) then
                    ix2 = iw + ido
                    if (na == 0) then
                        call RADB3 (ido, l1, c, ch, WA(iw), WA(ix2))
                    else
                        call RADB3 (ido, l1, ch, c, WA(iw), WA(ix2))
                    end if
                    na = 1 - na
                else
                    if (ip_rename == 5) then
                        ix2 = iw + ido
                        ix3 = ix2 + ido
                        ix4 = ix3 + ido
                        if (na == 0) then
                            call RADB5 (ido, l1, c, ch, WA(iw), WA(ix2), WA( &
                                ix3), WA(ix4))
                        else
                            call RADB5 (ido, l1, ch, c, WA(iw), WA(ix2), WA( &
                                ix3), WA(ix4))
                        end if
                        na = 1 - na
                    else
                        if (na == 0) then
                            call RADBG(ido, ip_rename, l1, idl1, c, c, c, ch, ch, WA(iw))
                        else
                            call RADBG(ido, ip_rename, l1, idl1, ch, ch, ch, c, c, WA(iw))
                        end if
                        if (ido == 1) na = 1 - na
                    end if
                end if
            end if
        end if
        l1 = l2
        iw = iw + (ip_rename - 1)*ido
    end do
    if (na == 0) return
    c(:n) = CH(:n)

end subroutine RFFTB1

subroutine RADB2(ido, l1, cc, ch, wa1)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, 2, l1)
    real , intent (out) :: ch(ido, l1, 2)
    real , intent (in) :: wa1(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, idp2, i, ic
    real :: tr2, ti2
    !-----------------------------------------------
    ch(1, :, 1) = CC(1, 1, :) + CC(ido, 2, :)
    ch(1, :, 2) = CC(1, 1, :) - CC(ido, 2, :)
    if (ido - 2 >= 0) then
        if (ido - 2 /= 0) then
            idp2 = ido + 2
            do k = 1, l1
                do i = 3, ido, 2
                    ic = idp2 - i
                    ch(i-1, k, 1) = CC(i-1, 1, k) + CC(ic-1, 2, k)
                    tr2 = CC(i-1, 1, k) - CC(ic-1, 2, k)
                    ch(i, k, 1) = CC(i, 1, k) - CC(ic, 2, k)
                    ti2 = CC(i, 1, k) + CC(ic, 2, k)
                    ch(i-1, k, 2) = WA1(i-2)*tr2 - WA1(i-1)*ti2
                    ch(i, k, 2) = WA1(i-2)*ti2 + WA1(i-1)*tr2
                end do
            end do
            if (mod(ido, 2) == 1) return
        end if
        ch(ido, :, 1) = CC(ido, 1, :) + CC(ido, 1, :)
        ch(ido, :, 2) = -(CC(1, 2, :)+CC(1, 2, :))
    end if

end subroutine RADB2

subroutine RADB3(ido, l1, cc, ch, wa1, wa2)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, 3, l1)
    real , intent (out) :: ch(ido, l1, 3)
    real , intent (in) :: wa1(*)
    real , intent (in) :: wa2(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, idp2, i, ic
    real::taur, taui, tr2, cr2, ci3, ti2, ci2, cr3, dr2, dr3, di2, di3
        !-----------------------------------------------
    parameter ( TAUI = 0.866025403784439 )
    parameter ( TAUR = -.5 )
    do k = 1, l1
        tr2 = CC(ido, 2, k) + CC(ido, 2, k)
        cr2 = CC(1, 1, k) + taur*tr2
        ch(1, k, 1) = CC(1, 1, k) + tr2
        ci3 = taui*(CC(1, 3, k)+CC(1, 3, k))
        ch(1, k, 2) = cr2 - ci3
        ch(1, k, 3) = cr2 + ci3
    end do
    if (ido == 1) return
    idp2 = ido + 2
    do k = 1, l1
        do i = 3, ido, 2
            ic = idp2 - i
            tr2 = CC(i-1, 3, k) + CC(ic-1, 2, k)
            cr2 = CC(i-1, 1, k) + taur*tr2
            ch(i-1, k, 1) = CC(i-1, 1, k) + tr2
            ti2 = CC(i, 3, k) - CC(ic, 2, k)
            ci2 = CC(i, 1, k) + taur*ti2
            ch(i, k, 1) = CC(i, 1, k) + ti2
            cr3 = taui*(CC(i-1, 3, k)-CC(ic-1, 2, k))
            ci3 = taui*(CC(i, 3, k)+CC(ic, 2, k))
            dr2 = cr2 - ci3
            dr3 = cr2 + ci3
            di2 = ci2 + cr3
            di3 = ci2 - cr3
            ch(i-1, k, 2) = WA1(i-2)*dr2 - WA1(i-1)*di2
            ch(i, k, 2) = WA1(i-2)*di2 + WA1(i-1)*dr2
            ch(i-1, k, 3) = WA2(i-2)*dr3 - WA2(i-1)*di3
            ch(i, k, 3) = WA2(i-2)*di3 + WA2(i-1)*dr3
        end do
    end do

end subroutine RADB3

subroutine RADB4(ido, l1, cc, ch, wa1, wa2, wa3)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, 4, l1)
    real , intent (out) :: ch(ido, l1, 4)
    real , intent (in) :: wa1(*)
    real , intent (in) :: wa2(*)
    real , intent (in) :: wa3(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, idp2, i, ic
    real :: tr1, tr2, tr3, tr4, ti1, ti2, ti3, ti4, cr3, ci3, &
        cr2, cr4, ci2, ci4
        !-----------------------------------------------
    do k = 1, l1
        tr1 = CC(1, 1, k) - CC(ido, 4, k)
        tr2 = CC(1, 1, k) + CC(ido, 4, k)
        tr3 = CC(ido, 2, k) + CC(ido, 2, k)
        tr4 = CC(1, 3, k) + CC(1, 3, k)
        ch(1, k, 1) = tr2 + tr3
        ch(1, k, 2) = tr1 - tr4
        ch(1, k, 3) = tr2 - tr3
        ch(1, k, 4) = tr1 + tr4
    end do
    if (ido - 2 >= 0) then
        if (ido - 2 /= 0) then
            idp2 = ido + 2
            do k = 1, l1
                do i = 3, ido, 2
                    ic = idp2 - i
                    ti1 = CC(i, 1, k) + CC(ic, 4, k)
                    ti2 = CC(i, 1, k) - CC(ic, 4, k)
                    ti3 = CC(i, 3, k) - CC(ic, 2, k)
                    tr4 = CC(i, 3, k) + CC(ic, 2, k)
                    tr1 = CC(i-1, 1, k) - CC(ic-1, 4, k)
                    tr2 = CC(i-1, 1, k) + CC(ic-1, 4, k)
                    ti4 = CC(i-1, 3, k) - CC(ic-1, 2, k)
                    tr3 = CC(i-1, 3, k) + CC(ic-1, 2, k)
                    ch(i-1, k, 1) = tr2 + tr3
                    cr3 = tr2 - tr3
                    ch(i, k, 1) = ti2 + ti3
                    ci3 = ti2 - ti3
                    cr2 = tr1 - tr4
                    cr4 = tr1 + tr4
                    ci2 = ti1 + ti4
                    ci4 = ti1 - ti4
                    ch(i-1, k, 2) = WA1(i-2)*cr2 - WA1(i-1)*ci2
                    ch(i, k, 2) = WA1(i-2)*ci2 + WA1(i-1)*cr2
                    ch(i-1, k, 3) = WA2(i-2)*cr3 - WA2(i-1)*ci3
                    ch(i, k, 3) = WA2(i-2)*ci3 + WA2(i-1)*cr3
                    ch(i-1, k, 4) = WA3(i-2)*cr4 - WA3(i-1)*ci4
                    ch(i, k, 4) = WA3(i-2)*ci4 + WA3(i-1)*cr4
                end do
            end do
            if (mod(ido, 2) == 1) return
        end if
        do k = 1, l1
            ti1 = CC(1, 2, k) + CC(1, 4, k)
            ti2 = CC(1, 4, k) - CC(1, 2, k)
            tr1 = CC(ido, 1, k) - CC(ido, 3, k)
            tr2 = CC(ido, 1, k) + CC(ido, 3, k)
            ch(ido, k, 1) = tr2 + tr2
            ch(ido, k, 2) = sqrt(2.0)*(tr1 - ti1)
            ch(ido, k, 3) = ti2 + ti2
            ch(ido, k, 4) = -sqrt(2.0)*(tr1 + ti1)
        end do
    end if
    return
end subroutine RADB4

subroutine RADB5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, 5, l1)
    real , intent (out) :: ch(ido, l1, 5)
    real , intent (in) :: wa1(*)
    real , intent (in) :: wa2(*)
    real , intent (in) :: wa3(*)
    real , intent (in) :: wa4(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, idp2, i, ic
    real :: tr11, ti11, tr12, ti12, ti5, ti4, tr2, tr3, cr2, cr3, ci5 &
        , ci4, ti2, ti3, tr5, tr4, ci2, ci3, cr5, cr4, dr3, dr4, di3, &
        di4, dr5, dr2, di5, di2
        !-----------------------------------------------
    parameter ( TI12 = 0.587785252292473 )
    parameter ( TR12 = -.809016994374947 )
    parameter ( TI11 = 0.951056516295154 )
    parameter ( TR11 = 0.309016994374947 )
    do k = 1, l1
        ti5 = CC(1, 3, k) + CC(1, 3, k)
        ti4 = CC(1, 5, k) + CC(1, 5, k)
        tr2 = CC(ido, 2, k) + CC(ido, 2, k)
        tr3 = CC(ido, 4, k) + CC(ido, 4, k)
        ch(1, k, 1) = CC(1, 1, k) + tr2 + tr3
        cr2 = CC(1, 1, k) + tr11*tr2 + tr12*tr3
        cr3 = CC(1, 1, k) + tr12*tr2 + tr11*tr3
        ci5 = ti11*ti5 + ti12*ti4
        ci4 = ti12*ti5 - ti11*ti4
        ch(1, k, 2) = cr2 - ci5
        ch(1, k, 3) = cr3 - ci4
        ch(1, k, 4) = cr3 + ci4
        ch(1, k, 5) = cr2 + ci5
    end do
    if (ido == 1) return
    idp2 = ido + 2
    do k = 1, l1
        do i = 3, ido, 2
            ic = idp2 - i
            ti5 = CC(i, 3, k) + CC(ic, 2, k)
            ti2 = CC(i, 3, k) - CC(ic, 2, k)
            ti4 = CC(i, 5, k) + CC(ic, 4, k)
            ti3 = CC(i, 5, k) - CC(ic, 4, k)
            tr5 = CC(i-1, 3, k) - CC(ic-1, 2, k)
            tr2 = CC(i-1, 3, k) + CC(ic-1, 2, k)
            tr4 = CC(i-1, 5, k) - CC(ic-1, 4, k)
            tr3 = CC(i-1, 5, k) + CC(ic-1, 4, k)
            ch(i-1, k, 1) = CC(i-1, 1, k) + tr2 + tr3
            ch(i, k, 1) = CC(i, 1, k) + ti2 + ti3
            cr2 = CC(i-1, 1, k) + tr11*tr2 + tr12*tr3
            ci2 = CC(i, 1, k) + tr11*ti2 + tr12*ti3
            cr3 = CC(i-1, 1, k) + tr12*tr2 + tr11*tr3
            ci3 = CC(i, 1, k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            dr3 = cr3 - ci4
            dr4 = cr3 + ci4
            di3 = ci3 + cr4
            di4 = ci3 - cr4
            dr5 = cr2 + ci5
            dr2 = cr2 - ci5
            di5 = ci2 - cr5
            di2 = ci2 + cr5
            ch(i-1, k, 2) = WA1(i-2)*dr2 - WA1(i-1)*di2
            ch(i, k, 2) = WA1(i-2)*di2 + WA1(i-1)*dr2
            ch(i-1, k, 3) = WA2(i-2)*dr3 - WA2(i-1)*di3
            ch(i, k, 3) = WA2(i-2)*di3 + WA2(i-1)*dr3
            ch(i-1, k, 4) = WA3(i-2)*dr4 - WA3(i-1)*di4
            ch(i, k, 4) = WA3(i-2)*di4 + WA3(i-1)*dr4
            ch(i-1, k, 5) = WA4(i-2)*dr5 - WA4(i-1)*di5
            ch(i, k, 5) = WA4(i-2)*di5 + WA4(i-1)*dr5
        end do
    end do

end subroutine RADB5

subroutine RADBG(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: ip
    integer , intent (in) :: l1
    integer , intent (in) :: idl1
    real , intent (in) :: cc(ido, ip, l1)
    real , intent (in out) :: c1(ido, l1, ip)
    real , intent (in out) :: c2(idl1, ip)
    real , intent (in out) :: ch(ido, l1, ip)
    real , intent (in out) :: ch2(idl1, ip)
    real , intent (in) :: wa(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer:: idp2, nbd, ipp2, ipph, k, i, j, jc, j2, l, lc, is, idij
    real:: tpi, arg, dcp, dsp, ar1, ai1, ar1h, dc2, ds2, ar2, ai2, ar2h
    !-----------------------------------------------
    tpi = 8.0*ATAN(1.0)
    arg = tpi/real(ip)
    dcp = COS(arg)
    dsp = SIN(arg)
    idp2 = ido + 2
    nbd = (ido - 1)/2
    ipp2 = ip + 2
    ipph = (ip + 1)/2
    if (ido >= l1) then
        ch(:, :, 1) = CC(:, 1, :)
    else
        ch(:, :, 1) = CC(:, 1, :)
    end if
    do j = 2, ipph
        jc = ipp2 - j
        j2 = j + j
        ch(1, :, j) = CC(ido, j2-2, :) + CC(ido, j2-2, :)
        ch(1, :, jc) = CC(1, j2-1, :) + CC(1, j2-1, :)
    end do
    if (ido /= 1) then
        if (nbd >= l1) then
            do j = 2, ipph
                jc = ipp2 - j
                ch(2:ido-1:2, :, j) = CC(2:ido-1:2, 2*j-1, :) + CC(idp2-4: &
                    idp2-1-ido:(-2), 2*j-2, :)
                ch(2:ido-1:2, :, jc) = CC(2:ido-1:2, 2*j-1, :) - CC(idp2-4: &
                    idp2-1-ido:(-2), 2*j-2, :)
                ch(3:ido:2, :, j) = CC(3:ido:2, 2*j-1, :) - CC(idp2-3:idp2- &
                    ido:(-2), 2*j-2, :)
                ch(3:ido:2, :, jc) = CC(3:ido:2, 2*j-1, :) + CC(idp2-3:idp2- &
                    ido:(-2), 2*j-2, :)
            end do
        else
            do j = 2, ipph
                jc = ipp2 - j
                ch(2:ido-1:2, :, j) = CC(2:ido-1:2, 2*j-1, :) + CC(idp2-4: &
                    idp2-1-ido:(-2), 2*j-2, :)
                ch(2:ido-1:2, :, jc) = CC(2:ido-1:2, 2*j-1, :) - CC(idp2-4: &
                    idp2-1-ido:(-2), 2*j-2, :)
                ch(3:ido:2, :, j) = CC(3:ido:2, 2*j-1, :) - CC(idp2-3:idp2- &
                    ido:(-2), 2*j-2, :)
                ch(3:ido:2, :, jc) = CC(3:ido:2, 2*j-1, :) + CC(idp2-3:idp2- &
                    ido:(-2), 2*j-2, :)
            end do
        end if
    end if
    ar1 = 1.
    ai1 = 0.
    do l = 2, ipph
        lc = ipp2 - l
        ar1h = dcp*ar1 - dsp*ai1
        ai1 = dcp*ai1 + dsp*ar1
        ar1 = ar1h
        c2(:, l) = CH2(:, 1) + ar1*CH2(:, 2)
        c2(:, lc) = ai1*CH2(:, ip)
        dc2 = ar1
        ds2 = ai1
        ar2 = ar1
        ai2 = ai1
        do j = 3, ipph
            jc = ipp2 - j
            ar2h = dc2*ar2 - ds2*ai2
            ai2 = dc2*ai2 + ds2*ar2
            ar2 = ar2h
            c2(:, l) = C2(:, l) + ar2*CH2(:, j)
            c2(:, lc) = C2(:, lc) + ai2*CH2(:, jc)
        end do
    end do
    do j = 2, ipph
        ch2(:, 1) = CH2(:, 1) + CH2(:, j)
    end do
    do j = 2, ipph
        jc = ipp2 - j
        ch(1, :, j) = C1(1, :, j) - C1(1, :, jc)
        ch(1, :, jc) = C1(1, :, j) + C1(1, :, jc)
    end do
    if (ido /= 1) then
        if (nbd >= l1) then
            do j = 2, ipph
                jc = ipp2 - j
                ch(2:ido-1:2, :, j) = C1(2:ido-1:2, :, j) - C1(3:ido:2, :, jc)
                ch(2:ido-1:2, :, jc) = C1(2:ido-1:2, :, j) + C1(3:ido:2, :, jc)
                ch(3:ido:2, :, j) = C1(3:ido:2, :, j) + C1(2:ido-1:2, :, jc)
                ch(3:ido:2, :, jc) = C1(3:ido:2, :, j) - C1(2:ido-1:2, :, jc)
            end do
        else
            do j = 2, ipph
                jc = ipp2 - j
                ch(2:ido-1:2, :, j) = C1(2:ido-1:2, :, j) - C1(3:ido:2, :, jc)
                ch(2:ido-1:2, :, jc) = C1(2:ido-1:2, :, j) + C1(3:ido:2, :, jc)
                ch(3:ido:2, :, j) = C1(3:ido:2, :, j) + C1(2:ido-1:2, :, jc)
                ch(3:ido:2, :, jc) = C1(3:ido:2, :, j) - C1(2:ido-1:2, :, jc)
            end do
        end if
    end if
    if (ido == 1) return
    c2(:, 1) = CH2(:, 1)
    c1(1, :, 2:ip) = CH(1, :, 2:ip)
    if (nbd <= l1) then
        is = -ido
        do j = 2, ip
            is = is + ido
            idij = is
            do i = 3, ido, 2
                idij = idij + 2
                c1(i-1, :, j) = WA(idij-1)*CH(i-1, :, j) - WA(idij)*CH(i, :, j)
                c1(i, :, j) = WA(idij-1)*CH(i, :, j) + WA(idij)*CH(i-1, :, j)
            end do
        end do
    else
        is = -ido
        do j = 2, ip
            is = is + ido
            do k = 1, l1
                idij = is
                c1(2:ido-1:2, k, j) = WA(idij+1:ido-2+idij:2)*CH(2:ido-1:2, &
                    k, j) - WA(idij+2:ido-1+idij:2)*CH(3:ido:2, k, j)
                c1(3:ido:2, k, j) = WA(idij+1:ido-2+idij:2)*CH(3:ido:2, k, j) &
                    + WA(idij+2:ido-1+idij:2)*CH(2:ido-1:2, k, j)
            end do
        end do
    end if

end subroutine RADBG

subroutine RFFTF(n, r, wsave)

    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: n
    real (wp)     :: r(*)
    real (wp)     :: wsave(*)
    !-----------------------------------------------

    if (n == 1) return

    call RFFTF1 (n, r, wsave, wsave(n+1), int(wsave(2*n+1:2*n+1), kind=ip) )

end subroutine RFFTF

subroutine RFFTF1(n, c, ch, wa, ifac)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: n
    integer , intent (in) :: ifac(*)
    real  :: c(*)
    real  :: ch(*)
    real  :: wa(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer:: nf, na, l2, iw, k1, kh, ip_rename, l1, ido, idl1, ix2, ix3, ix4
    !-----------------------------------------------
    nf = IFAC(2)
    na = 1
    l2 = n
    iw = n
    do k1 = 1, nf
        kh = nf - k1
        ip_rename = IFAC(kh+3)
        l1 = l2/ip_rename
        ido = n/l2
        idl1 = ido*l1
        iw = iw - (ip_rename - 1)*ido
        na = 1 - na
        if (ip_rename == 4) then
            ix2 = iw + ido
            ix3 = ix2 + ido
            if (na == 0) then
                call RADF4 (ido, l1, c, ch, WA(iw), WA(ix2), WA(ix3))
                go to 110
            end if
            call RADF4 (ido, l1, ch, c, WA(iw), WA(ix2), WA(ix3))
            go to 110
        end if
        if (ip_rename == 2) then
            if (na == 0) then
                call RADF2 (ido, l1, c, ch, WA(iw))
                go to 110
            end if
            call RADF2 (ido, l1, ch, c, WA(iw))
            go to 110
        end if
104 continue
    if (ip_rename == 3) then
        ix2 = iw + ido
        if (na == 0) then
            call RADF3 (ido, l1, c, ch, WA(iw), WA(ix2))
            go to 110
        end if
        call RADF3 (ido, l1, ch, c, WA(iw), WA(ix2))
        go to 110
    end if
106 continue
    if (ip_rename == 5) then
        ix2 = iw + ido
        ix3 = ix2 + ido
        ix4 = ix3 + ido
        if (na == 0) then
            call RADF5(ido, l1, c, ch, WA(iw), WA(ix2), WA(ix3), WA(ix4))
            go to 110
        end if
        call RADF5(ido, l1, ch, c, WA(iw), WA(ix2), WA(ix3), WA(ix4))
        go to 110
    end if
    if (ido == 1) na = 1 - na
    if (na == 0) then
        call RADFG (ido, ip_rename, l1, idl1, c, c, c, ch, ch, WA(iw))
        na = 1
    else
        call RADFG (ido, ip_rename, l1, idl1, ch, ch, ch, c, c, WA(iw))
        na = 0
    end if
110 continue
    l2 = l1
end do
if (na == 1) return
c(:n) = CH(:n)

end subroutine RFFTF1

subroutine RADF2(ido, l1, cc, ch, wa1)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, l1, 2)
    real , intent (out) :: ch(ido, 2, l1)
    real , intent (in) :: wa1(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, idp2, i, ic
    real :: tr2, ti2
    !-----------------------------------------------
    ch(1, 1, :) = CC(1, :, 1) + CC(1, :, 2)
    ch(ido, 2, :) = CC(1, :, 1) - CC(1, :, 2)
    if (ido - 2 >= 0) then
        if (ido - 2 /= 0) then
            idp2 = ido + 2
            do k = 1, l1
                do i = 3, ido, 2
                    ic = idp2 - i
                    tr2 = WA1(i-2)*CC(i-1, k, 2) + WA1(i-1)*CC(i, k, 2)
                    ti2 = WA1(i-2)*CC(i, k, 2) - WA1(i-1)*CC(i-1, k, 2)
                    ch(i, 1, k) = CC(i, k, 1) + ti2
                    ch(ic, 2, k) = ti2 - CC(i, k, 1)
                    ch(i-1, 1, k) = CC(i-1, k, 1) + tr2
                    ch(ic-1, 2, k) = CC(i-1, k, 1) - tr2
                end do
            end do
            if (mod(ido, 2) == 1) return
        end if
        ch(1, 2, :) = -CC(ido, :, 2)
        ch(ido, 1, :) = CC(ido, :, 1)
    end if

end subroutine RADF2

subroutine RADF3(ido, l1, cc, ch, wa1, wa2)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, l1, 3)
    real , intent (out) :: ch(ido, 3, l1)
    real , intent (in) :: wa1(*)
    real , intent (in) :: wa2(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, idp2, i, ic
    real::taur, taui, cr2, dr2, di2, dr3, di3, ci2, tr2, ti2, tr3, ti3
        !-----------------------------------------------
    parameter ( TAUI = 0.866025403784439 )
    parameter ( TAUR = -.5 )
    do k = 1, l1
        cr2 = CC(1, k, 2) + CC(1, k, 3)
        ch(1, 1, k) = CC(1, k, 1) + cr2
        ch(1, 3, k) = taui*(CC(1, k, 3)-CC(1, k, 2))
        ch(ido, 2, k) = CC(1, k, 1) + taur*cr2
    end do
    if (ido == 1) return
    idp2 = ido + 2
    do k = 1, l1
        do i = 3, ido, 2
            ic = idp2 - i
            dr2 = WA1(i-2)*CC(i-1, k, 2) + WA1(i-1)*CC(i, k, 2)
            di2 = WA1(i-2)*CC(i, k, 2) - WA1(i-1)*CC(i-1, k, 2)
            dr3 = WA2(i-2)*CC(i-1, k, 3) + WA2(i-1)*CC(i, k, 3)
            di3 = WA2(i-2)*CC(i, k, 3) - WA2(i-1)*CC(i-1, k, 3)
            cr2 = dr2 + dr3
            ci2 = di2 + di3
            ch(i-1, 1, k) = CC(i-1, k, 1) + cr2
            ch(i, 1, k) = CC(i, k, 1) + ci2
            tr2 = CC(i-1, k, 1) + taur*cr2
            ti2 = CC(i, k, 1) + taur*ci2
            tr3 = taui*(di2 - di3)
            ti3 = taui*(dr3 - dr2)
            ch(i-1, 3, k) = tr2 + tr3
            ch(ic-1, 2, k) = tr2 - tr3
            ch(i, 3, k) = ti2 + ti3
            ch(ic, 2, k) = ti3 - ti2
        end do
    end do

end subroutine RADF3

subroutine RADF4(ido, l1, cc, ch, wa1, wa2, wa3)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, l1, 4)
    real , intent (out) :: ch(ido, 4, l1)
    real , intent (in) :: wa1(*)
    real , intent (in) :: wa2(*)
    real , intent (in) :: wa3(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, idp2, i, ic
    real :: hsqt2, tr1, tr2, cr2, ci2, cr3, ci3, cr4, ci4, tr4, ti1, &
        ti4, ti2, ti3, tr3
        !-----------------------------------------------
    parameter ( HSQT2 = 0.7071067811865475 )
    do k = 1, l1
        tr1 = CC(1, k, 2) + CC(1, k, 4)
        tr2 = CC(1, k, 1) + CC(1, k, 3)
        ch(1, 1, k) = tr1 + tr2
        ch(ido, 4, k) = tr2 - tr1
        ch(ido, 2, k) = CC(1, k, 1) - CC(1, k, 3)
        ch(1, 3, k) = CC(1, k, 4) - CC(1, k, 2)
    end do
    if (ido - 2 >= 0) then
        if (ido - 2 /= 0) then
            idp2 = ido + 2
            do k = 1, l1
                do i = 3, ido, 2
                    ic = idp2 - i
                    cr2 = WA1(i-2)*CC(i-1, k, 2) + WA1(i-1)*CC(i, k, 2)
                    ci2 = WA1(i-2)*CC(i, k, 2) - WA1(i-1)*CC(i-1, k, 2)
                    cr3 = WA2(i-2)*CC(i-1, k, 3) + WA2(i-1)*CC(i, k, 3)
                    ci3 = WA2(i-2)*CC(i, k, 3) - WA2(i-1)*CC(i-1, k, 3)
                    cr4 = WA3(i-2)*CC(i-1, k, 4) + WA3(i-1)*CC(i, k, 4)
                    ci4 = WA3(i-2)*CC(i, k, 4) - WA3(i-1)*CC(i-1, k, 4)
                    tr1 = cr2 + cr4
                    tr4 = cr4 - cr2
                    ti1 = ci2 + ci4
                    ti4 = ci2 - ci4
                    ti2 = CC(i, k, 1) + ci3
                    ti3 = CC(i, k, 1) - ci3
                    tr2 = CC(i-1, k, 1) + cr3
                    tr3 = CC(i-1, k, 1) - cr3
                    ch(i-1, 1, k) = tr1 + tr2
                    ch(ic-1, 4, k) = tr2 - tr1
                    ch(i, 1, k) = ti1 + ti2
                    ch(ic, 4, k) = ti1 - ti2
                    ch(i-1, 3, k) = ti4 + tr3
                    ch(ic-1, 2, k) = tr3 - ti4
                    ch(i, 3, k) = tr4 + ti3
                    ch(ic, 2, k) = tr4 - ti3
                end do
            end do
            if (mod(ido, 2) == 1) return
        end if
        do k = 1, l1
            ti1 = -hsqt2*(CC(ido, k, 2)+CC(ido, k, 4))
            tr1 = hsqt2*(CC(ido, k, 2)-CC(ido, k, 4))
            ch(ido, 1, k) = tr1 + CC(ido, k, 1)
            ch(ido, 3, k) = CC(ido, k, 1) - tr1
            ch(1, 2, k) = ti1 - CC(ido, k, 3)
            ch(1, 4, k) = ti1 + CC(ido, k, 3)
        end do
    end if

end subroutine RADF4

subroutine RADF5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: l1
    real , intent (in) :: cc(ido, l1, 5)
    real , intent (out) :: ch(ido, 5, l1)
    real , intent (in) :: wa1(*)
    real , intent (in) :: wa2(*)
    real , intent (in) :: wa3(*)
    real , intent (in) :: wa4(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: k, idp2, i, ic
    real :: tr11, ti11, tr12, ti12, cr2, ci5, cr3, ci4, dr2, di2, dr3 &
        , di3, dr4, di4, dr5, di5, cr5, ci2, cr4, ci3, tr2, ti2, tr3, &
        ti3, tr5, ti5, tr4, ti4
        !-----------------------------------------------
    parameter ( TI12 = 0.587785252292473 )
    parameter ( TR12 = -.809016994374947 )
    parameter ( TI11 = 0.951056516295154 )
    parameter ( TR11 = 0.309016994374947 )
    do k = 1, l1
        cr2 = CC(1, k, 5) + CC(1, k, 2)
        ci5 = CC(1, k, 5) - CC(1, k, 2)
        cr3 = CC(1, k, 4) + CC(1, k, 3)
        ci4 = CC(1, k, 4) - CC(1, k, 3)
        ch(1, 1, k) = CC(1, k, 1) + cr2 + cr3
        ch(ido, 2, k) = CC(1, k, 1) + tr11*cr2 + tr12*cr3
        ch(1, 3, k) = ti11*ci5 + ti12*ci4
        ch(ido, 4, k) = CC(1, k, 1) + tr12*cr2 + tr11*cr3
        ch(1, 5, k) = ti12*ci5 - ti11*ci4
    end do
    if (ido == 1) return
    idp2 = ido + 2
    do k = 1, l1
        do i = 3, ido, 2
            ic = idp2 - i
            dr2 = WA1(i-2)*CC(i-1, k, 2) + WA1(i-1)*CC(i, k, 2)
            di2 = WA1(i-2)*CC(i, k, 2) - WA1(i-1)*CC(i-1, k, 2)
            dr3 = WA2(i-2)*CC(i-1, k, 3) + WA2(i-1)*CC(i, k, 3)
            di3 = WA2(i-2)*CC(i, k, 3) - WA2(i-1)*CC(i-1, k, 3)
            dr4 = WA3(i-2)*CC(i-1, k, 4) + WA3(i-1)*CC(i, k, 4)
            di4 = WA3(i-2)*CC(i, k, 4) - WA3(i-1)*CC(i-1, k, 4)
            dr5 = WA4(i-2)*CC(i-1, k, 5) + WA4(i-1)*CC(i, k, 5)
            di5 = WA4(i-2)*CC(i, k, 5) - WA4(i-1)*CC(i-1, k, 5)
            cr2 = dr2 + dr5
            ci5 = dr5 - dr2
            cr5 = di2 - di5
            ci2 = di2 + di5
            cr3 = dr3 + dr4
            ci4 = dr4 - dr3
            cr4 = di3 - di4
            ci3 = di3 + di4
            ch(i-1, 1, k) = CC(i-1, k, 1) + cr2 + cr3
            ch(i, 1, k) = CC(i, k, 1) + ci2 + ci3
            tr2 = CC(i-1, k, 1) + tr11*cr2 + tr12*cr3
            ti2 = CC(i, k, 1) + tr11*ci2 + tr12*ci3
            tr3 = CC(i-1, k, 1) + tr12*cr2 + tr11*cr3
            ti3 = CC(i, k, 1) + tr12*ci2 + tr11*ci3
            tr5 = ti11*cr5 + ti12*cr4
            ti5 = ti11*ci5 + ti12*ci4
            tr4 = ti12*cr5 - ti11*cr4
            ti4 = ti12*ci5 - ti11*ci4
            ch(i-1, 3, k) = tr2 + tr5
            ch(ic-1, 2, k) = tr2 - tr5
            ch(i, 3, k) = ti2 + ti5
            ch(ic, 2, k) = ti5 - ti2
            ch(i-1, 5, k) = tr3 + tr4
            ch(ic-1, 4, k) = tr3 - tr4
            ch(i, 5, k) = ti3 + ti4
            ch(ic, 4, k) = ti4 - ti3
        end do
    end do

end subroutine RADF5

subroutine RADFG(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: ido
    integer , intent (in) :: ip
    integer , intent (in) :: l1
    integer , intent (in) :: idl1
    real , intent (out) :: cc(ido, ip, l1)
    real , intent (in out) :: c1(ido, l1, ip)
    real , intent (in out) :: c2(idl1, ip)
    real , intent (in out) :: ch(ido, l1, ip)
    real , intent (in out) :: ch2(idl1, ip)
    real , intent (in) :: wa(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer:: ipph, ipp2, idp2, nbd, j, k, is, idij, i, jc, l, lc
    real:: tpi, arg, dcp, dsp, ar1, ai1, ar1h, dc2, ds2, ar2, ai2, ar2h
    !-----------------------------------------------
    tpi = 8.0*ATAN(1.0)
    arg = tpi/real(ip)
    dcp = COS(arg)
    dsp = SIN(arg)
    ipph = (ip + 1)/2
    ipp2 = ip + 2
    idp2 = ido + 2
    nbd = (ido - 1)/2
    if (ido /= 1) then
        ch2(:, 1) = C2(:, 1)
        ch(1, :, 2:ip) = C1(1, :, 2:ip)
        if (nbd <= l1) then
            is = -ido
            do j = 2, ip
                is = is + ido
                idij = is
                do i = 3, ido, 2
                    idij = idij + 2
                    ch(i-1, :, j)=WA(idij-1)*C1(i-1, :, j)+WA(idij)*C1(i, :, j)
                    ch(i, :, j)=WA(idij-1)*C1(i, :, j)-WA(idij)*C1(i-1, :, j)
                end do
            end do
        else
            is = -ido
            do j = 2, ip
                is = is + ido
                do k = 1, l1
                    idij = is
                    ch(2:ido-1:2, k, j) = WA(idij+1:ido-2+idij:2)*C1(2:ido-1 &
                        :2, k, j) + WA(idij+2:ido-1+idij:2)*C1(3:ido:2, k, j)
                    ch(3:ido:2, k, j) = WA(idij+1:ido-2+idij:2)*C1(3:ido:2, k &
                        , j) - WA(idij+2:ido-1+idij:2)*C1(2:ido-1:2, k, j)
                end do
            end do
        end if
        if (nbd >= l1) then
            do j = 2, ipph
                jc = ipp2 - j
                c1(2:ido-1:2, :, j)=CH(2:ido-1:2, :, j)+CH(2:ido-1:2, :, jc)
                c1(2:ido-1:2, :, jc) = CH(3:ido:2, :, j) - CH(3:ido:2, :, jc)
                c1(3:ido:2, :, j) = CH(3:ido:2, :, j) + CH(3:ido:2, :, jc)
                c1(3:ido:2, :, jc) = CH(2:ido-1:2, :, jc) - CH(2:ido-1:2, :, j)
            end do
            go to 121
        end if
        do j = 2, ipph
            jc = ipp2 - j
            c1(2:ido-1:2, :, j) = CH(2:ido-1:2, :, j) + CH(2:ido-1:2, :, jc)
            c1(2:ido-1:2, :, jc) = CH(3:ido:2, :, j) - CH(3:ido:2, :, jc)
            c1(3:ido:2, :, j) = CH(3:ido:2, :, j) + CH(3:ido:2, :, jc)
            c1(3:ido:2, :, jc) = CH(2:ido-1:2, :, jc) - CH(2:ido-1:2, :, j)
        end do
        go to 121
    end if
    c2(:, 1) = CH2(:, 1)
121 continue
    do j = 2, ipph
        jc = ipp2 - j
        c1(1, :, j) = CH(1, :, j) + CH(1, :, jc)
        c1(1, :, jc) = CH(1, :, jc) - CH(1, :, j)
    end do
    !
    ar1 = 1.
    ai1 = 0.
    do l = 2, ipph
        lc = ipp2 - l
        ar1h = dcp*ar1 - dsp*ai1
        ai1 = dcp*ai1 + dsp*ar1
        ar1 = ar1h
        ch2(:, l) = C2(:, 1) + ar1*C2(:, 2)
        ch2(:, lc) = ai1*C2(:, ip)
        dc2 = ar1
        ds2 = ai1
        ar2 = ar1
        ai2 = ai1
        do j = 3, ipph
            jc = ipp2 - j
            ar2h = dc2*ar2 - ds2*ai2
            ai2 = dc2*ai2 + ds2*ar2
            ar2 = ar2h
            ch2(:, l) = CH2(:, l) + ar2*C2(:, j)
            ch2(:, lc) = CH2(:, lc) + ai2*C2(:, jc)
        end do
    end do
    do j = 2, ipph
        ch2(:, 1) = CH2(:, 1) + C2(:, j)
    end do
    !
    if (ido >= l1) then
        cc(:, 1, :) = CH(:, :, 1)
    else
        cc(:, 1, :) = CH(:, :, 1)
    end if
    cc(ido, 2:(ipph-1)*2:2, :) = TRANSPOSE(CH(1, :, 2:ipph))
    cc(1, 3:ipph*2-1:2, :) = TRANSPOSE(CH(1, :, ipp2-2:ipp2-ipph:(-1)))
    if (ido == 1) return
    if (nbd >= l1) then
        cc(2:ido-1:2, 3:ipph*2-1:2, :) = reshape(SOURCE = CH(2:ido-1:2, :, &
            2:ipph)+CH(2:ido-1:2, :, ipp2-2:ipp2-ipph:(-1)), SHAPE = (/(ido &
            -1)/2, ipph-1, l1/), ORDER = (/1, 3, 2/))
        cc(idp2-4:idp2-1-ido:(-2), 2:(ipph-1)*2:2, :) = reshape(SOURCE = &
            CH(2:ido-1:2, :, 2:ipph)-CH(2:ido-1:2, :, ipp2-2:ipp2-ipph:(-1)) &
            , SHAPE = (/(ido-1)/2, ipph-1, l1/), ORDER = (/1, 3, 2/))
        cc(3:ido:2, 3:ipph*2-1:2, :) = reshape(SOURCE = CH(3:ido:2, :, 2: &
            ipph)+CH(3:ido:2, :, ipp2-2:ipp2-ipph:(-1)), SHAPE = (/(ido-1)/ &
            2, ipph-1, l1/), ORDER = (/1, 3, 2/))
        cc(idp2-3:idp2-ido:(-2), 2:(ipph-1)*2:2, :) = reshape(SOURCE = CH &
            (3:ido:2, :, ipp2-2:ipp2-ipph:(-1))-CH(3:ido:2, :, 2:ipph), SHAPE &
            = (/(ido-1)/2, ipph-1, l1/), ORDER = (/1, 3, 2/))
        return
    end if

    cc(2:ido-1:2, 3:ipph*2-1:2, :) = reshape(SOURCE = CH(2:ido-1:2, :, 2: &
        ipph)+CH(2:ido-1:2, :, ipp2-2:ipp2-ipph:(-1)), SHAPE = (/(ido-1)/2 &
        , ipph-1, l1/), ORDER = (/1, 3, 2/))
    cc(idp2-4:idp2-1-ido:(-2), 2:(ipph-1)*2:2, :) = reshape(SOURCE = CH( &
        2:ido-1:2, :, 2:ipph)-CH(2:ido-1:2, :, ipp2-2:ipp2-ipph:(-1)), SHAPE &
        = (/(ido-1)/2, ipph-1, l1/), ORDER = (/1, 3, 2/))
    cc(3:ido:2, 3:ipph*2-1:2, :) = reshape(SOURCE = CH(3:ido:2, :, 2:ipph) &
        +CH(3:ido:2, :, ipp2-2:ipp2-ipph:(-1)), SHAPE = (/(ido-1)/2, ipph-1 &
        , l1/), ORDER = (/1, 3, 2/))
    cc(idp2-3:idp2-ido:(-2), 2:(ipph-1)*2:2, :) = reshape(SOURCE = CH(3: &
        ido:2, :, ipp2-2:ipp2-ipph:(-1))-CH(3:ido:2, :, 2:ipph), SHAPE = (/( &
        ido-1)/2, ipph-1, l1/), ORDER = (/1, 3, 2/))

end subroutine RADFG

end module module_fftpack
!     this function is define in the file comf.f
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June 2004 2004    fortran 90 updates
!-----------------------------------------------------------------------
!     END
