C  modified version in which subroutine names such as DSCAL are replaced by
C  DSCAL_ to avoid multiple definition errors when linking with AUTO
C
      DOUBLE PRECISION FUNCTION A02ABF(XXR,XXI)
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS THE ABSOLUTE VALUE OF A COMPLEX NUMBER VIA ROUTINE
C     NAME
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 XXI, XXR
C     .. Local Scalars ..
      DOUBLE PRECISION                 H, ONE, XI, XR, ZERO
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, SQRT
C     .. Data statements ..
      DATA                             ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
C
      XR = ABS(XXR)
      XI = ABS(XXI)
      IF (XI.LE.XR) GO TO 20
      H = XR
      XR = XI
      XI = H
   20 IF (XI.NE.ZERO) GO TO 40
      A02ABF = XR
      RETURN
   40 H = XR*SQRT(ONE+(XI/XR)**2)
      A02ABF = H
      RETURN
      END
      SUBROUTINE A02ACF(XXR,XXI,YYR,YYI,ZR,ZI)
C     MARK 2A RELEASE.  NAG COPYRIGHT 1973
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     DIVIDES ONE COMPLEX NUMBER BY A SECOND
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XXI, XXR, YYI, YYR, ZI, ZR
C     .. Local Scalars ..
      DOUBLE PRECISION  A, H, ONE
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              ONE/1.0D0/
C     .. Executable Statements ..
C
      IF (ABS(YYR).LE.ABS(YYI)) GO TO 20
      H = YYI/YYR
      A = ONE/(H*YYI+YYR)
      ZR = (XXR+H*XXI)*A
      ZI = (XXI-H*XXR)*A
      RETURN
   20 H = YYR/YYI
      A = ONE/(H*YYR+YYI)
      ZR = (H*XXR+XXI)*A
      ZI = (H*XXI-XXR)*A
      RETURN
      END
      SUBROUTINE F01AKF(N,K,L,A,IA,INTGER)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 8 REVISED. IER-248 (JUN 1980).
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. EXTENDED BLAS (JUNE 1986)
C
C     DIRHES
C     AUGUST 1ST, 1971 .
C     GIVEN THE UNSYMMETRIC MATRIX, A, STORED IN THE ARRAY A(N,N),
C     THIS SUBROUTINE REDUCES THE SUB-MATRIX OF ORDER L - K + 1,
C     WHICH STARTS AT THE ELEMENT A(K,K) AND FINISHES AT THE
C     ELEMENT A(L,L), TO HESSENBERG FORM, H, BY THE DIRECT
C     METHOD(AN = NH). THE MATRIX H IS OVERWRITTEN ON A WITH
C     DETAILS OF THE TRANSFORMATIONS (N) STORED IN THE REMAINING
C     TRIANGLE UNDER H AND IN ELEMENTS K TO L OF THE ARRAY
C     INTGER(N).
C     1ST AUGUST 1971
C
C     .. Scalar Arguments ..
      INTEGER           IA, K, L, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N)
      INTEGER           INTGER(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  X, Y
      INTEGER           I, J, K1, M
C     .. External Subroutines ..
      EXTERNAL          DGEMV_, DTRSV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      K1 = K + 1
      IF (K1.GT.L) RETURN
      DO 140 J = K1, N
         M = J
         X = 0.0D0
         IF (J.GT.L) GO TO 120
         DO 20 I = J, L
            IF (ABS(A(I,J-1)).LE.ABS(X)) GO TO 20
            X = A(I,J-1)
            M = I
   20    CONTINUE
         INTGER(J) = M
         IF (M.EQ.J) GO TO 80
C        INTERCHANGE ROWS AND COLUMNS OF A.
         DO 40 I = K, N
            Y = A(M,I)
            A(M,I) = A(J,I)
            A(J,I) = Y
   40    CONTINUE
         DO 60 I = 1, L
            Y = A(I,M)
            A(I,M) = A(I,J)
            A(I,J) = Y
   60    CONTINUE
   80    IF (X.NE.0.0D0 .AND. J.LT.L) THEN
            DO 100 I = J + 1, L
               A(I,J-1) = A(I,J-1)/X
  100       CONTINUE
            CALL DGEMV_('N',L,L-J,1.0D0,A(1,J+1),IA,A(J+1,J-1),1,1.0D0,
     *                 A(1,J),1)
         END IF
  120    CALL DTRSV('L','N','U',J-K,A(K+1,K),IA,A(K+1,J),1)
      IF (J.LT.L) CALL DGEMV_('N',L-J,J-K,-1.0D0,A(J+1,K),IA,A(K+1,J),
     *                      1,1.0D0,A(J+1,J),1)
  140 CONTINUE
      RETURN
      END
      SUBROUTINE F01APF(N,LOW,IUPP,INTGER,H,IH,V,IV)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     DIRTRANS
C     FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS IN THE ARRAY
C     V(N,N) FROM THE INFORMATION LEFT BY SUBROUTINE F01AKF
C     BELOW THE UPPER HESSENBERG MATRIX, H, IN THE ARRAY H(N,N)
C     AND IN THE INTEGER ARRAY INTGER(N).
C     1ST AUGUST 1971
C
C     .. Scalar Arguments ..
      INTEGER           IH, IUPP, IV, LOW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  H(IH,N), V(IV,N)
      INTEGER           INTGER(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  X
      INTEGER           I, I1, II, J, LOW1, M
C     .. Executable Statements ..
      DO 40 I = 1, N
         DO 20 J = 1, N
            V(I,J) = 0.0D0
   20    CONTINUE
         V(I,I) = 1.0D0
   40 CONTINUE
      LOW1 = LOW + 1
      IF (LOW1.GT.IUPP) RETURN
      DO 120 II = LOW1, IUPP
         I = LOW1 + IUPP - II
         I1 = I - 1
         IF (LOW1.GT.I1) GO TO 80
         DO 60 J = LOW1, I1
            V(I,J) = H(I,J-1)
   60    CONTINUE
   80    M = INTGER(I)
         IF (M.EQ.I) GO TO 120
         DO 100 J = LOW1, IUPP
            X = V(M,J)
            V(M,J) = V(I,J)
            V(I,J) = X
  100    CONTINUE
  120 CONTINUE
      RETURN
      END
      SUBROUTINE F01ATF(N,IB,A,IA,LOW,LHI,D)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     BALANCE
C     REDUCE THE NORM OF A(N,N) BY EXACT DIAGONAL SIMILARITY
C     TRANSFORMATIONS STORED IN D(N).
C     DECEMBER 1ST.,1971
C
C     .. Scalar Arguments ..
      INTEGER           IA, IB, LHI, LOW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), D(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  B2, C, F, G, R, S
      INTEGER           I, J, JJ, K, L
      LOGICAL           NOCONV
C     .. External Subroutines ..
      EXTERNAL          F01ATZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
      B2 = IB*IB
      L = 1
      K = N
   20 IF (K.LT.1) GO TO 100
C     SEARCH FOR ROWS ISOLATING AN EIGENVALUE AND PUSH THEM DOWN
      J = K + 1
      DO 60 JJ = 1, K
         J = J - 1
         R = 0.0D0
         DO 40 I = 1, K
            IF (I.EQ.J) GO TO 40
            R = R + ABS(A(J,I))
   40    CONTINUE
         IF (R.EQ.0.0D0) GO TO 80
   60 CONTINUE
      GO TO 100
   80 CALL F01ATZ(K,A,IA,D,K,L,N,J)
      K = K - 1
      GO TO 20
C     SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE AND PUSH THEM
C     LEFT.
  100 IF (L.GT.K) GO TO 180
      DO 140 J = L, K
         C = 0.0D0
         DO 120 I = L, K
            IF (I.EQ.J) GO TO 120
            C = C + ABS(A(I,J))
  120    CONTINUE
         IF (C.EQ.0.0D0) GO TO 160
  140 CONTINUE
      GO TO 180
  160 CALL F01ATZ(L,A,IA,D,K,L,N,J)
      L = L + 1
      GO TO 100
C     NOW BALANCE THE SUBMATRIX IN ROWS L THROUGH K.
  180 LOW = L
      LHI = K
      IF (L.GT.K) GO TO 220
      DO 200 I = L, K
         D(I) = 1.0D0
  200 CONTINUE
  220 NOCONV = .FALSE.
      IF (L.GT.K) GO TO 420
      DO 400 I = L, K
         C = 0.0D0
         R = 0.0D0
         DO 240 J = L, K
            IF (J.EQ.I) GO TO 240
            C = C + ABS(A(J,I))
            R = R + ABS(A(I,J))
  240    CONTINUE
         G = R/DBLE(IB)
         F = 1.0D0
         S = C + R
  260    IF (C.GE.G) GO TO 280
         F = F*DBLE(IB)
         C = C*B2
         GO TO 260
  280    G = R*DBLE(IB)
  300    IF (C.LT.G) GO TO 320
         F = F/DBLE(IB)
         C = C/B2
         GO TO 300
  320    IF (((C+R)/F).GE.(0.95D0*S)) GO TO 400
         G = 1.0D0/F
         D(I) = D(I)*F
         NOCONV = .TRUE.
         IF (L.GT.N) GO TO 360
         DO 340 J = L, N
            A(I,J) = A(I,J)*G
  340    CONTINUE
  360    DO 380 J = 1, K
            A(J,I) = A(J,I)*F
  380    CONTINUE
  400 CONTINUE
  420 IF (NOCONV) GO TO 220
      RETURN
      END
      SUBROUTINE F01ATZ(M,A,IA,D,K,L,N,J)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C
C     AUXILIARY ROUTINE CALLED BY F01ATF.
C     INTERCHANGES ELEMENTS 1 TO K OF COLUMNS J AND M,
C     AND ELEMENTS L TO N OF ROWS J AND M.
C     .. Scalar Arguments ..
      INTEGER           IA, J, K, L, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), D(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F
      INTEGER           I
C     .. Executable Statements ..
      D(M) = J
      IF (J.EQ.M) GO TO 60
      DO 20 I = 1, K
         F = A(I,J)
         A(I,J) = A(I,M)
         A(I,M) = F
   20 CONTINUE
      IF (L.GT.N) GO TO 60
      DO 40 I = L, N
         F = A(J,I)
         A(J,I) = A(M,I)
         A(M,I) = F
   40 CONTINUE
   60 RETURN
      END
      SUBROUTINE F01AUF(N,LOW,LHI,M,D,Z,IZ)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     BALBAK
C     BACKWARD TRANSFORMATION OF A SET OF RIGHT-HAND EIGENVECTORS
C     OF A BALANCED MATRIX INTO THE EIGENVECTORS OF THE ORIGINAL
C     MATRIX FROM WHICH THE BALANCED MATRIX WAS DERIVED BY A CALL
C     OF SUBROUTINE F01ATF.
C     DECEMBER 1ST.,1971
C
C     .. Scalar Arguments ..
      INTEGER           IZ, LHI, LOW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N), Z(IZ,M)
C     .. Local Scalars ..
      DOUBLE PRECISION  S
      INTEGER           I, II, J, K, LHI1, LOW1
C     .. Executable Statements ..
      IF (LOW.GT.LHI) GO TO 60
      DO 40 I = LOW, LHI
         S = D(I)
C        LEFT-HAND EIGENVECTORS ARE BACK TRANSFORMED IF THE
C        FOREGOING STATEMENT IS REPLACED BY S=1/D(I)
         DO 20 J = 1, M
            Z(I,J) = Z(I,J)*S
   20    CONTINUE
   40 CONTINUE
   60 I = LOW
      LOW1 = LOW - 1
      IF (LOW1.LT.1) GO TO 120
      DO 100 II = 1, LOW1
         I = I - 1
         K = D(I)
         IF (K.EQ.I) GO TO 100
         DO 80 J = 1, M
            S = Z(I,J)
            Z(I,J) = Z(K,J)
            Z(K,J) = S
   80    CONTINUE
  100 CONTINUE
  120 LHI1 = LHI + 1
      IF (LHI1.GT.N) RETURN
      DO 160 I = LHI1, N
         K = D(I)
         IF (K.EQ.I) GO TO 160
         DO 140 J = 1, M
            S = Z(I,J)
            Z(I,J) = Z(K,J)
            Z(K,J) = S
  140    CONTINUE
  160 CONTINUE
      RETURN
      END
      SUBROUTINE F01LZF(N,A,NRA,C,NRC,WANTB,B,WANTQ,WANTY,Y,NRY,LY,
     *                  WANTZ,Z,NRZ,NCZ,D,E,WORK1,WORK2,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (BIDIAG)
C
C     F01LZF RETURNS ALL OR PART OF THE FACTORIZATION OF THE
C     N*N UPPER TRIANGULAR MATRIX A GIVEN BY
C
C     A = Q*C*(P**T) ,
C
C     WHERE Q AND P ARE N*N ORTHOGONAL MATRICES AND C IS AN
C     N*N UPPER BIDIAGONAL MATRIX.
C
C     IF WANTB IS .TRUE. THEN B RETURNS (Q**T)*B.
C     IF WANTY IS .TRUE. THEN Y RETURNS Y*Q.
C     IF WANTZ IS .TRUE. THEN Z RETURNS (P**T)*Z.
C
C     INPUT PARAMETERS.
C
C     N     - ORDER OF THE MATRIX A.
C
C     A     - THE N*N UPPER TRIANGULAR MATRIX TO BE FACTORIZED. THE
C             STRICTLY LOWER TRIANGULAR PART OF A IS NOT REFERENCED.
C
C     NRA   - ROW DIMENSION OF A AS DECLARED IN THE CALLING PROGRAM.
C             NRA MUST BE AT LEAST N.
C
C     NRC   - ROW DIMENSION OF C AS DECLARED IN THE CALLING PROGRAM.
C             NRC MUST BE AT LEAST N.
C
C     WANTB - MUST BE .TRUE. IF (Q**T)*B IS REQUIRED.
C             IF WANTB IS .FALSE. THEN B IS NOT REFERENCED.
C
C     B     - AN N ELEMENT REAL VECTOR.
C
C     WANTQ - MUST BE .TRUE. IF DETAILS OF Q ARE TO BE
C             STORED BELOW THE BIDIAGONAL PART OF C.
C             IF WANTQ IS .FALSE. THEN THE LOWER TRIANGULAR
C             PART OF C IS NOT REFERENCED.
C
C     WANTY - MUST BE .TRUE. IF Y*Q IS REQUIRED.
C             IF WANTY IS .FALSE. THEN Y IS NOT REFERENCED.
C
C     Y     - AN LY*N REAL MATRIX.
C
C     NRY   - IF WANTY IS .TRUE. THEN NRY MUST BE THE ROW
C             DIMENSION OF Y AS DECLARED IN THE CALLING
C             PROGRAM AND MUST BE AT LEAST LY.
C
C     LY    - IF WANTY IS .TRUE. THEN LY MUST BE THE NUMBER
C             OF ROWS OF Y AND MUST BE AT LEAST 1.
C
C     WANTZ - MUST BE .TRUE. IF (P**T)*Z IS REQUIRED.
C             IF WANTZ IS .FALSE. THEN Z IS NOT REFERENCED.
C
C     Z     - AN N*NCZ REAL MATRIX.
C
C     NRZ   - IF WANTZ IS .TRUE. THEN NRZ MUST BE THE ROW
C             DIMENSION OF Z AS DECLARED IN THE CALLING
C             PROGRAM AND MUST BE AT LEAST N.
C
C     NCZ   - IF WANTZ IS .TRUE. THEN NCZ MUST BE THE
C             NUMBER OF COLUMNS OF Z AND MUST BE AT LEAST
C             1.
C
C     IFAIL - THE USUAL FAILURE PARAMETER. IF IN DOUBT SET
C             IFAIL TO ZERO BEFORE CALLING THIS ROUTINE.
C
C     OUTPUT PARAMETERS.
C
C     C     - N*N MATRIX CONTAINING THE UPPER BIDIAGONAL MATRIX B.
C             DETAILS OF P ARE STORED ABOVE THE BIDIAGONAL
C             PART OF C. UNLESS WANTQ IS .TRUE. THE
C             STRICTLY LOWER TRIANGULAR PART OF C IS NOT
C             REFERENCED.
C             THE ROUTINE MAY BE CALLED WITH C=A.
C
C     B     - IF WANTB IS .TRUE. THEN B WILL RETURN THE N ELEMENT
C             VECTOR (Q**T)*B.
C
C     Y     - IF WANTY IS .TRUE. THEN Y WILL RETURN THE
C             LY*N MATRIX Y*Q.
C
C     Z     - IF WANTZ IS .TRUE. THEN Z WILL RETURN THE N*NCZ MATRIX
C             (P**T)*Z.
C
C     D     - N ELEMENT VECTOR CONTAINING THE DIAGONAL ELEMENTS OF C
C             SUCH THAT D(I)=C(I,I), I=1,2,...,N.
C
C     E     - N ELEMENT VECTOR CONTAINING THE
C             SUPER-DIAGONAL ELEMENTS OF C SUCH THAT
C             E(I)=C(I-1,I), I=2,3,...,N. E(1) IS NOT
C             REFERENCED.
C
C     IFAIL - ON NORMAL RETURN IFAIL WILL BE ZERO.
C             IF AN INPUT PARAMETER IS INCORRECTLY SUPPLIED
C             THEN IFAIL IS SET TO UNITY. NO OTHER FAILURE
C             IS POSSIBLE.
C
C     WORKSPACE PARAMETERS.
C
C     WORK1
C     WORK2 - N ELEMENT REAL VECTORS.
C             IF WANTZ IS .FALSE. THEN WORK1 AND WORK2 ARE NOT
C             REFERENCED.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01LZF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LY, N, NCZ, NRA, NRC, NRY, NRZ
      LOGICAL           WANTB, WANTQ, WANTY, WANTZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NRA,N), B(N), C(NRC,N), D(N), E(N), WORK1(N),
     *                  WORK2(N), Y(NRY,N), Z(NRZ,NCZ)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG, CS, EPS, RSQTPS, SMALL, SN, SQTEPS, T, W, X
      INTEGER           I, IERR, J, JJ, JP1, K, KP1, KP2, NM2
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  F01LZZ, X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          F01LZZ, X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01LZW, F01LZX, F01LZY
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      IERR = IFAIL
      IF (IERR.EQ.0) IFAIL = 1
C
      IF (NRA.LT.N .OR. NRC.LT.N .OR. N.LT.1) GO TO 220
      IF (WANTY .AND. (NRY.LT.LY .OR. LY.LT.1)) GO TO 220
      IF (WANTZ .AND. (NRZ.LT.N .OR. NCZ.LT.1)) GO TO 220
C
      SMALL = X02AMF()
      BIG = 1.0D0/SMALL
      EPS = X02AJF()
      SQTEPS = SQRT(EPS)
      RSQTPS = 1.0D0/SQTEPS
C
      D(1) = A(1,1)
C
      DO 40 J = 1, N
         DO 20 I = 1, J
            C(I,J) = A(I,J)
   20    CONTINUE
   40 CONTINUE
C
      IFAIL = 0
      IF (N.EQ.1) RETURN
      IF (N.EQ.2) GO TO 200
C
C     START MAIN LOOP. K(TH) STEP PUTS ZEROS INTO K(TH) ROW OF C.
C
      NM2 = N - 2
      DO 180 K = 1, NM2
         KP1 = K + 1
C
C        SET UP PLANE ROTATION P(J,J+1) TO ANNIHILATE C(K,J+1).
C        THIS ROTATION INTRODUCES AN UNWANTED ELEMENT IN C(J+1,J)
C        WHICH IS STORED IN X.
C        J GOES N-1,N-2,...,K+1.
C
         J = N
         DO 100 JJ = K, NM2
            JP1 = J
            J = J - 1
            W = C(K,JP1)
C
            T = F01LZZ(C(K,J),W,SMALL,BIG)
C
            C(K,JP1) = T
            X = 0.0D0
C
            CALL F01LZW(T,CS,SN,SQTEPS,RSQTPS,BIG)
C
            IF ( .NOT. WANTZ) GO TO 60
            WORK1(J) = CS
            WORK2(J) = SN
C
   60       IF (T.EQ.0.0D0) GO TO 80
            C(K,J) = CS*C(K,J) + SN*W
C
C           NOW APPLY THE TRANSFORMATION P(J,J+1).
C
            CALL F01LZY(J-K,CS,SN,C(KP1,J),C(KP1,JP1))
C
            X = SN*C(JP1,JP1)
            C(JP1,JP1) = CS*C(JP1,JP1)
C
C           NOW SET UP PLANE ROTATION Q(J,J+1)**T TO ANNIHILATE
C           X=C(J+1,J).
C
   80       T = F01LZZ(C(J,J),X,SMALL,BIG)
C
            IF (WANTQ) C(JP1,K) = T
C
            CALL F01LZW(T,D(J),E(J),SQTEPS,RSQTPS,BIG)
C
            C(J,J) = D(J)*C(J,J) + E(J)*X
C
            IF (WANTY) CALL F01LZY(LY,D(J),E(J),Y(1,J),Y(1,JP1))
C
  100    CONTINUE
C
C        NOW APPLY THE TRANSFORMATIONS Q(J,J+1)**T AND FORM
C        (P(J,J+1)**T)*Z, J=N-1,N-2,...,K+1 COLUMN BY COLUMN
C
         KP2 = KP1 + 1
         DO 120 J = KP2, N
C
            CALL F01LZX(J-K,D(K),E(K),C(KP1,J))
C
  120    CONTINUE
C
         IF (WANTB) CALL F01LZX(N-K,D(K),E(K),B(KP1))
C
         IF ( .NOT. WANTZ) GO TO 160
         DO 140 J = 1, NCZ
C
            CALL F01LZX(N-K,WORK1(K),WORK2(K),Z(KP1,J))
C
  140    CONTINUE
C
  160    D(KP1) = C(KP1,KP1)
         E(KP1) = C(K,KP1)
C
  180 CONTINUE
C
  200 D(N) = C(N,N)
      E(N) = C(N-1,N)
      RETURN
C
  220 IFAIL = P01ABF(IERR,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F01LZW(T,C,S,SQTEPS,RSQTPS,BIG)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (COSSIN)
C
C     F01LZW RETURNS THE VALUES
C
C     C = COS(THETA)   AND   S = SIN(THETA)
C
C     FOR A GIVEN VALUE OF
C
C     T = TAN(THETA) .
C
C     C IS ALWAYS NON-NEGATIVE AND S HAS THE SAME SIGN AS T.
C
C     SQTEPS, RSQTPS AND BIG MUST BE SUCH THAT
C
C     SQTEPS = SQRT(X02AJF) , RSQTPS = 1.0/SQTEPS AND BIG =
C     1.0/X02AMF ,
C
C     WHERE X02AJF AND X02AMF ARE THE NUMBERS RETURNED FROM
C     ROUTINES X02AJF AND X02AMF RESPECTIVELY.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIG, C, RSQTPS, S, SQTEPS, T
C     .. Local Scalars ..
      DOUBLE PRECISION  ABST, TT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, SQRT
C     .. Executable Statements ..
      IF (T.NE.0.0D0) GO TO 20
      C = 1.0D0
      S = 0.0D0
      RETURN
C
   20 ABST = ABS(T)
      IF (ABST.LT.SQTEPS) GO TO 60
      IF (ABST.GT.RSQTPS) GO TO 80
C
      TT = ABST*ABST
      IF (ABST.GT.1.0D0) GO TO 40
C
      TT = 0.25D0*TT
      C = 0.5D0/SQRT(0.25D0+TT)
      S = C*T
      RETURN
C
   40 TT = 0.25D0/TT
      S = 0.5D0/SQRT(0.25D0+TT)
      C = S/ABST
      S = SIGN(S,T)
      RETURN
C
   60 C = 1.0D0
      S = T
      RETURN
C
   80 C = 0.0D0
      IF (ABST.LT.BIG) C = 1.0D0/ABST
      S = SIGN(1.0D0,T)
      RETURN
      END
      SUBROUTINE F01LZX(N,C,S,X)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (PLROT6)
C
C     F01LZX RETURNS THE N ELEMENT VECTOR
C
C     Y = R(1,2)*R(2,3)*...*R(N-1,N)*X ,
C
C     WHERE X IS AN N ELEMENT VECTOR AND R(J-1,J) IS A PLANE
C     ROTATION FOR THE (J-1,J)-PLANE.
C
C     Y IS OVERWRITTEN ON X.
C
C     THE N ELEMENT VECTORS C AND S MUST BE SUCH THAT THE
C     NON-IDENTITY PART OF R(J-1,J) IS GIVEN BY
C
C     R(J-1,J) = (  C(J)  S(J) ) .
C                ( -S(J)  C(J) )
C
C     C(1) AND S(1) ARE NOT REFERENCED.
C
C
C     N MUST BE AT LEAST 1. IF N=1 THEN AN IMMEDIATE RETURN TO
C     THE CALLING PROGRAM IS MADE.
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(N), S(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  W
      INTEGER           I, II, IM1
C     .. Executable Statements ..
      IF (N.EQ.1) RETURN
C
      I = N
      DO 20 II = 2, N
         IM1 = I - 1
         W = X(IM1)
         X(IM1) = C(I)*W + S(I)*X(I)
         X(I) = C(I)*X(I) - S(I)*W
         I = IM1
   20 CONTINUE
C
      RETURN
      END
      SUBROUTINE F01LZY(N,C,S,X,Y)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (PLROT8)
C
C     F01LZY FORMS THE N*2 MATRIX
C
C     Z = ( X  Y )*( C  -S ) ,
C                  ( S   C )
C
C     WHERE X AND Y ARE N ELEMENT VECTORS, C=COS(THETA) AND
C     S=SIN(THETA).
C
C     THE FIRST COLUMN OF Z IS OVERWRITTEN ON X AND THE SECOND
C     COLUMN OF Z IS OVERWRITTEN ON Y.
C
C
C     N MUST BE AT LEAST 1.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C, S
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  W
      INTEGER           I
C     .. Executable Statements ..
      DO 20 I = 1, N
         W = X(I)
         X(I) = C*W + S*Y(I)
         Y(I) = C*Y(I) - S*W
   20 CONTINUE
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION F01LZZ(A,B,SMALL,BIG)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (TANGNT)
C
C     F01LZZ RETURNS THE VALUE
C
C     F01LZZ = B/A .
C
C     SMALL AND BIG MUST BE SUCH THAT
C
C     SMALL = X02AMF     AND     BIG = 1.0/SMALL ,
C
C     WHERE X02AMF IS THE SMALL NUMBER RETURNED FROM ROUTINE
C     X02AMF.
C
C     IF B/A IS LESS THAN SMALL THEN F01LZZ IS RETURNED AS
C     ZERO AND IF B/A IS GREATER THAN BIG THEN F01LZZ IS
C     RETURNED AS SIGN(BIG,B).
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B, BIG, SMALL
C     .. Local Scalars ..
      DOUBLE PRECISION                 ABSA, ABSB, X
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, SIGN
C     .. Executable Statements ..
      F01LZZ = 0.0D0
      IF (B.EQ.0.0D0) RETURN
C
      ABSA = ABS(A)
      ABSB = ABS(B)
      X = 0.0D0
      IF (ABSA.GE.1.0D0) X = ABSA*SMALL
C
      IF (ABSB.LT.X) RETURN
C
      X = 0.0D0
      IF (ABSB.GE.1.0D0) X = ABSB*SMALL
C
      IF (ABSA.LE.X) GO TO 20
C
      F01LZZ = B/A
      RETURN
C
   20 F01LZZ = SIGN(BIG,B)
      RETURN
      END
      SUBROUTINE F01QCF(M,N,A,LDA,ZETA,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F01QCF  finds  the  QR factorization  of the real  m by n,  m .ge. n,
C  matrix A,  so that  A is reduced to upper triangular form by means of
C  orthogonal transformations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )   when   m.gt.n,
C           ( 0 )
C
C     A = Q*R       when   m = n,
C
C  where  Q  is an  m by m orthogonal matrix and  R  is an  n by n upper
C  triangular matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ), which is used  to introduce zeros into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',
C
C     u( k ) = ( zeta( k ) ),
C              (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k ) and z( k )  are chosen to annhilate the elements  below the
C  triangular part of  A.
C
C  The vector  u( k ) is returned in the kth element of  ZETA and in the
C  kth column of A, such that zeta( k ) is in ZETA( k ) and the elements
C  of  z( k ) are in  a( k + 1, k ), ..., a( m, k ).  The elements of  R
C  are returned in the upper triangular part of  A.
C
C  Q is given by
C
C     Q = ( Q( n )*Q( n - 1 )*...*Q( 1 ) )'.
C
C  3. Parameters
C     ==========
C
C  M      - INTEGER.
C
C           On entry, M must specify the number of rows of  A. M must be
C           at least  n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the number of columns of  A. N must
C           be  at  least zero. When  N = 0  then an immediate return is
C           effected.
C
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On exit, the  N by N upper triangular part of A will contain
C           the upper triangular matrix R and the  M by N strictly lower
C           triangular  part   of   A   will  contain  details   of  the
C           factorization as described above.
C
C  LDA    - INTEGER.
C
C           On entry, LDA  must  specify  the  leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  ZETA   - REAL             array of DIMENSION at least ( n ).
C
C           On exit,  ZETA( k )  contains the scalar  zeta( k )  for the
C           kth  transformation.  If  T( k ) = I  then  ZETA( k ) = 0.0,
C           otherwise  ZETA( k )  contains  zeta( k ) as described above
C           and  zeta( k ) is always in the range  ( 1.0, sqrt( 2.0 ) ).
C
C  IFAIL  - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1 to specify noisy soft failure or noisy hard failure  or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On successful  exit  IFAIL  will be  zero,  otherwise  IFAIL
C           will  be set to  -1  indicating that an  input parameter has
C           been  incorrectly  set. See  the  next section  for  further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  IFAIL = -1
C
C     One or more of the following conditions holds:
C
C        M   .lt. N
C        N   .lt. 0
C        LDA .lt. M
C
C  If  on  entry,  IFAIL  was  either  -1 or 0  then  further diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C        B := Q*B   and   B := Q'*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  NAG Library routine  F01QDF. The  operation  B := Q*B can be obtained
C  by the call:
C
C     IFAIL = 0
C     CALL F01QDF( 'No transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, IFAIL )
C
C  and  B := Q'*B  can be obtained by the call:
C
C     IFAIL = 0
C     CALL F01QDF( 'Transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, IFAIL )
C
C  In  both  cases  WORK  must be a  k  element array  that  is used  as
C  workspace. If  B  is a one-dimensional array (single column) then the
C  parameter  LDB  can be replaced by  M. See routine F01QDF for further
C  details.
C
C  The first k columns of the orthogonal matrix Q can either be obtained
C  by setting  B to the first k columns of the unit matrix and using the
C  first of the above two calls,  or by calling the  NAG Library routine
C  F01QEF, which overwrites the k columns of Q on the first k columns of
C  the array A.  Q is obtained by the call:
C
C     CALL F01QEF( 'Separate', M, N, K, A, LDA, ZETA, WORK, IFAIL )
C
C  As above WORK must be a k element array.  If K is larger than N, then
C  A must have been declared to have at least K columns.
C
C  Operations involving the matrix  R  can readily  be performed by  the
C  Level 2 BLAS  routines  DTRSV  and DTRMV  (see Chapter F06), but note
C  that no test for  near singularity  of  R  is incorporated in DTRSV .
C  If  R  is singular,  or nearly singular then the  NAG Library routine
C  F02WUF  can be  used to  determine  the  singular value decomposition
C  of  R.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 21-December-1985.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01QCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ZETA(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           IERR, K, LA
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMV_, DGER_, F06FRF, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IERR = 0
      IF (M.LT.N) CALL P01ABY(M,'M',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.M) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Perform the factorization.
C
      IF (N.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      LA = LDA
      DO 20 K = 1, MIN(M-1,N)
C
C        Use a  Householder reflection  to  zero the  kth column  of  A.
C        First set up the reflection.
C
         CALL F06FRF(M-K,A(K,K),A(K+1,K),1,ZERO,ZETA(K))
         IF ((ZETA(K).GT.ZERO) .AND. (K.LT.N)) THEN
            IF ((K+1).EQ.N) LA = M - K + 1
C
C           Temporarily  store  beta and  put  zeta( k )  in  a( k, k ).
C
            TEMP = A(K,K)
            A(K,K) = ZETA(K)
C
C           We now perform the operation  A := Q( k )*A.
C
C           Let  B  denote  the bottom  ( m - k + 1 ) by ( n - k )  part
C           of  A.
C
C           First form   work = B'*u.  ( work  is stored in the elements
C           ZETA( k + 1 ), ..., ZETA( n ). )
C
            CALL DGEMV_('Transpose',M-K+1,N-K,ONE,A(K,K+1),LA,A(K,K),1,
     *                 ZERO,ZETA(K+1),1)
C
C           Now form  B := B - u*work'.
C
            CALL DGER_(M-K+1,N-K,-ONE,A(K,K),1,ZETA(K+1),1,A(K,K+1),LA)
C
C           Restore beta.
C
            A(K,K) = TEMP
         END IF
   20 CONTINUE
C
C     Set the final  ZETA  when  m.eq.n.
C
      IF (M.EQ.N) ZETA(N) = ZERO
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01QCF. ( SGEQR )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
      SUBROUTINE F01QDF(TRANS,WHERET,M,N,A,LDA,ZETA,NCOLB,B,LDB,WORK,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F01QDF performs one of the transformations
C
C     B := Q*B   or   B := Q'*B,
C
C  where  B is an  m by ncolb real matrix and  Q is an m by m orthogonal
C  matrix, given as the product of  Householder transformation matrices.
C
C  This  routine  is  intended  for  use  following  NAG Fortran Library
C  routine  F01QCF.
C
C  2. Description
C     ===========
C
C  Q is assumed to be given by
C
C     Q = ( Q( n )*Q( n - 1 )*...*Q( 1 ) )',
C
C  Q( k ) being given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )'
C
C     u( k ) = ( zeta( k ) ),
C              (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C
C  z( k )  must  be  supplied  in  the  kth  column  of  A  in  elements
C  a( k + 1, k ), ..., a( m, k )  and  zeta( k ) must be supplied either
C  in  a( k, k ) or in zeta( k ),  depending upon the parameter  WHERET.
C
C  To obtain Q explicitly B may be set to I and premultiplied by Q. This
C  is more efficient than obtaining Q'.
C
C  3. Parameters
C     ==========
C
C  TRANS  - CHARACTER*1.
C
C           On entry, TRANS  specifies the operation to be performed  as
C           follows.
C
C           TRANS = 'N' or 'n'  ( No transpose )
C
C              Perform the operation  B := Q*B.
C
C           TRANS = 'T' or 't' or 'C' or 'c'  ( Transpose )
C
C              Perform the operation  B := Q'*B.
C
C           Unchanged on exit.
C
C  WHERET - CHARACTER*1.
C
C           On entry,  WHERET  specifies where the elements of  zeta are
C           to be found as follows.
C
C           WHERET = 'I' or 'i'   ( In A )
C
C              The elements of zeta are in A.
C
C           WHERET = 'S' or 's'   ( Separate )
C
C              The elements of zeta are separate from A, in ZETA.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at least n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be  at least zero. When  N = 0  then an immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  stricly lower  triangular
C           part of the array  A  must contain details of the matrix  Q.
C           In  addition, when  WHERET = 'I' or 'i'  then  the  diagonal
C           elements of A must contain the elements of zeta as described
C           under the argument  ZETA  below.
C
C           When  WHERET = 'S' or 's'  then the diagonal elements of the
C           array  A  are referenced, since they are used temporarily to
C           store the  zeta( k ), but they contain their original values
C           on return.
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must specify  the leading dimension  of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least m.
C
C           Unchanged on exit.
C
C  ZETA   - REAL             array of  DIMENSION  at least  ( n ),  when
C           WHERET = 'S' or 's'.
C
C           Before entry with  WHERET = 'S' or 's', the array  ZETA must
C           contain  the  elements  of  zeta.  If  ZETA( k ) = 0.0  then
C           T( k )  is assumed  to be  I otherwise  ZETA( k ) is assumed
C           to contain zeta( k ).
C
C           When WHERET = 'I' or 'i', the array  ZETA is not referenced.
C
C           Unchanged on exit.
C
C  NCOLB  - INTEGER.
C
C           On  entry, NCOLB  must specify  the number of columns of  B.
C           NCOLB  must  be  at  least  zero.  When  NCOLB = 0  then  an
C           immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - REAL             array of DIMENSION ( LDB, ncolb ).
C
C           Before entry, the leading  M by NCOLB  part of  the array  B
C           must  contain  the matrix to be  transformed.
C
C           On  exit,  B  is  overwritten  by  the  transformed  matrix.
C
C  LDB    - INTEGER.
C
C           On  entry, LDB  must specify  the  leading dimension of  the
C           array  B as declared in the calling (sub) program. LDB  must
C           be at least m.
C
C           Unchanged on exit.
C
C  WORK   - REAL             array of DIMENSION at least ( ncolb ).
C
C           Used as internal workspace.
C
C  IFAIL  - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1 to specify noisy soft failure or noisy hard failure  or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On  successful exit  IFAIL  will be  zero,  otherwise  IFAIL
C           will  be set to   -1  indicating that an input parameter has
C           been  incorrectly  set. See  the  next  section  for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  IFAIL = -1
C
C     One or more of the following conditions holds:
C
C        TRANS  .ne. 'N' or 'n' or 'T' or 't' or 'C' or 'c'
C        WHERET .ne. 'I' or 'i' or 'S' or 's'
C        M      .lt. N
C        N      .lt. 0
C        LDA    .lt. M
C        NCOLB  .lt. 0
C        LDB    .lt. M
C
C  If  on  entry,  IFAIL  was either  -1 or 0  then  further  diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 13-November-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01QDF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDB, M, N, NCOLB
      CHARACTER*1       TRANS, WHERET
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), WORK(*), ZETA(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP, ZETAK
      INTEGER           IERR, K, KK, LB
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMV_, DGER_, P01ABW, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IERR = 0
      IF ((TRANS.NE.'N') .AND. (TRANS.NE.'n') .AND. (TRANS.NE.'T')
     *    .AND. (TRANS.NE.'t') .AND. (TRANS.NE.'C') .AND. (TRANS.NE.'c')
     *    ) CALL P01ABW(TRANS,'TRANS',IFAIL,IERR,SRNAME)
      IF ((WHERET.NE.'I') .AND. (WHERET.NE.'i') .AND. (WHERET.NE.'S')
     *     .AND. (WHERET.NE.'s')) CALL P01ABW(WHERET,'WHERET',IFAIL,
     *    IERR,SRNAME)
      IF (M.LT.N) CALL P01ABY(M,'M',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.M) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (NCOLB.LT.0) CALL P01ABY(NCOLB,'NCOLB',IFAIL,IERR,SRNAME)
      IF (LDB.LT.M) CALL P01ABY(LDB,'LDB',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Perform the transformation.
C
      IF (MIN(N,NCOLB).EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      LB = LDB
      DO 20 KK = 1, N
         IF ((TRANS.EQ.'T') .OR. (TRANS.EQ.'t') .OR. (TRANS.EQ.'C')
     *       .OR. (TRANS.EQ.'c')) THEN
C
C           Q'*B = Q( n )*...*Q( 2 )*Q( 1 )*B,
C
            K = KK
         ELSE
C
C           Q*B  = Q( 1 )'*Q( 2 )'*...*Q( n )'*B,
C
            K = N + 1 - KK
         END IF
         IF ((WHERET.EQ.'S') .OR. (WHERET.EQ.'s')) THEN
            ZETAK = ZETA(K)
         ELSE
            ZETAK = A(K,K)
         END IF
         IF (ZETAK.GT.ZERO) THEN
            TEMP = A(K,K)
            A(K,K) = ZETAK
            IF (NCOLB.EQ.1) LB = M - K + 1
C
C           Let C denote the bottom ( m - k + 1 ) by ncolb part of B.
C
C           First form  work = C'*u.
C
            CALL DGEMV_('Transpose',M-K+1,NCOLB,ONE,B(K,1),LB,A(K,K),1,
     *                 ZERO,WORK,1)
C
C           Now form  C := C - u*work'.
C
            CALL DGER_(M-K+1,NCOLB,-ONE,A(K,K),1,WORK,1,B(K,1),LB)
C
C           Restore the diagonal element of A.
C
            A(K,K) = TEMP
         END IF
   20 CONTINUE
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01QDF. ( SGEAPQ )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
      SUBROUTINE F02AGF(A,IA,N,RR,RI,VR,IVR,VI,IVI,INTGER,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 14A REVISED. IER-685 (DEC 1989).
C
C     EIGENVALUES AND EIGENVECTORS OF REAL UNSYMMETRIC MATRIX
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AGF')
C     .. Scalar Arguments ..
      INTEGER           IA, IFAIL, IVI, IVR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), RI(N), RR(N), VI(IVI,N), VR(IVR,N)
      INTEGER           INTGER(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, D, MACHEP, MAX, SUM, TERM
      INTEGER           I, IB, ISAVE, J, K, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF, X02BHF
      EXTERNAL          X02AJF, P01ABF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          F01AKF, F01APF, F01ATF, F01AUF, F02AQF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      MACHEP = X02AJF()
      IB = X02BHF()
      CALL F01ATF(N,IB,A,IA,K,L,RR)
      CALL F01AKF(N,K,L,A,IA,INTGER)
      CALL F01APF(N,K,L,INTGER,A,IA,VR,IVR)
      CALL F01AUF(N,K,L,N,RR,VR,IVR)
      CALL F02AQF(N,1,N,MACHEP,A,IA,VR,IVR,RR,RI,INTGER,IFAIL)
      IF (IFAIL.EQ.0) GO TO 20
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
   20 DO 140 I = 1, N
         IF (RI(I).EQ.0.0D0) GO TO 60
         IF (RI(I).GT.0.0D0) GO TO 100
         DO 40 J = 1, N
            VR(J,I) = VR(J,I-1)
            VI(J,I) = -VI(J,I-1)
   40    CONTINUE
         GO TO 140
   60    DO 80 J = 1, N
            VI(J,I) = 0.0D0
   80    CONTINUE
         GO TO 140
  100    DO 120 J = 1, N
            VI(J,I) = VR(J,I+1)
  120    CONTINUE
  140 CONTINUE
      DO 280 I = 1, N
         SUM = 0.0D0
         MAX = 0.0D0
         DO 180 J = 1, N
            IF (ABS(VR(J,I)).LE.MAX) GO TO 160
            MAX = ABS(VR(J,I))
  160       IF (ABS(VI(J,I)).LE.MAX) GO TO 180
            MAX = ABS(VI(J,I))
  180    CONTINUE
         DO 200 J = 1, N
            VR(J,I) = VR(J,I)/MAX
            VI(J,I) = VI(J,I)/MAX
  200    CONTINUE
         MAX = 0.0D0
         DO 240 J = 1, N
            TERM = VR(J,I)**2 + VI(J,I)**2
            SUM = SUM + TERM
            IF (TERM.LE.MAX) GO TO 220
            MAX = TERM
            C = VR(J,I)
            D = -VI(J,I)
  220       CONTINUE
  240    CONTINUE
         SUM = SUM*(C**2+D**2)
         SUM = SQRT(SUM)
         DO 260 J = 1, N
            TERM = VR(J,I)
            VR(J,I) = (VR(J,I)*C-VI(J,I)*D)/SUM
            VI(J,I) = (D*TERM+C*VI(J,I))/SUM
  260    CONTINUE
  280 CONTINUE
      RETURN
      END
      SUBROUTINE F02AQF(N,LOW,UPP,MACHEP,H,IH,VECS,IVECS,WR,WI,CNT,
     *                  IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C
C     HQR2
C     FINDS THE EIGENVALUES AND EIGENVECTORS OF A REAL MATRIX
C     WHICH HAS BEEN REDUCED TO UPPER HESSENBERG FORM IN THE ARRAY
C     H(N,N) WITH THE ACCUMULATED TRANSFORMATIONS STORED IN
C     THE ARRAY VECS(N,N). THE REAL AND IMAGINARY PARTS OF THE
C     EIGENVALUES ARE FORMED IN THE ARRAYS WR, WI(N) AND THE
C     EIGENVECTORS ARE FORMED IN THE ARRAY VECS(N,N) WHERE
C     ONLY ONE COMPLEX VECTOR, CORRESPONDING TO THE ROOT WITH
C     POSITIVE IMAGINARY PART, IS FORMED FOR A COMPLEX PAIR. LOW
C     AND UPP ARE TWO INTEGERS PRODUCED IN BALANCING WHERE
C     EIGENVALUES ARE ISOLATED IN POSITIONS 1 TO LOW-1 AND UPP+1
C     TO N. IF BALANCING IS NOT USED LOW=1, UPP=N. MACHEP IS THE
C     RELATIVE MACHINE PRECISION. THE SUBROUTINE WILL FAIL IF
C     ALL EIGENVALUES TAKE MORE THAN 30*N ITERATIONS.
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AQF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  MACHEP
      INTEGER           IFAIL, IH, IVECS, LOW, N, UPP
C     .. Array Arguments ..
      DOUBLE PRECISION  H(IH,N), VECS(IVECS,N), WI(N), WR(N)
      INTEGER           CNT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  NORM, P, Q, R, RA, S, SA, T, U, VI, VR, W, X, Y,
     *                  Z
      INTEGER           EN, EN2, I, I1, II, ISAVE, ITN, ITS, J, JJ, K,
     *                  KK, L, LL, LOW1, M, M2, M3, MM, NA, NA1, NHS,
     *                  UPP1
      LOGICAL           NOTLAS
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF, X02AJF, X02AKF
      INTEGER           P01ABF
      EXTERNAL          A02ABF, X02AJF, X02AKF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          A02ACF, DGEMV_
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN, DBLE, SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
C     COMPUTE MATRIX NORM
      NORM = 0.0D0
      K = 1
      DO 40 I = 1, N
         DO 20 J = K, N
            NORM = NORM + ABS(H(I,J))
   20    CONTINUE
         K = I
   40 CONTINUE
      NHS = N*(N+1)/2 + N - 1
C     ISOLATED ROOTS
      IF (LOW.LE.1) GO TO 80
      J = LOW - 1
      DO 60 I = 1, J
         WR(I) = H(I,I)
         WI(I) = 0.0D0
         CNT(I) = 0
   60 CONTINUE
   80 IF (UPP.GE.N) GO TO 120
      J = UPP + 1
      DO 100 I = J, N
         WR(I) = H(I,I)
         WI(I) = 0.0D0
         CNT(I) = 0
  100 CONTINUE
  120 EN = UPP
      T = 0.0D0
      ITN = 30*N
  140 IF (EN.LT.LOW) GO TO 880
      ITS = 0
      NA = EN - 1
C     LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
  160 IF (LOW+1.GT.EN) GO TO 200
      LOW1 = LOW + 1
      DO 180 LL = LOW1, EN
         L = EN + LOW1 - LL
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S.LT.X02AKF()/X02AJF()) S = NORM/DBLE(NHS)
         IF (ABS(H(L,L-1)).LE.MACHEP*S) GO TO 220
  180 CONTINUE
  200 L = LOW
  220 X = H(EN,EN)
      IF (L.EQ.EN) GO TO 740
      Y = H(NA,NA)
      W = H(EN,NA)*H(NA,EN)
      IF (L.EQ.NA) GO TO 760
      IF (ITN.LE.0) GO TO 1500
C     FORM SHIFT
      IF ((ITS.NE.10) .AND. (ITS.NE.20)) GO TO 280
      T = T + X
      IF (LOW.GT.EN) GO TO 260
      DO 240 I = LOW, EN
         H(I,I) = H(I,I) - X
  240 CONTINUE
  260 S = ABS(H(EN,NA)) + ABS(H(NA,EN-2))
      X = 0.75D0*S
      Y = X
      W = -0.4375D0*S**2
  280 ITS = ITS + 1
      ITN = ITN - 1
C     LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS
      IF (L.GT.EN-2) GO TO 320
      EN2 = EN - 2
      DO 300 MM = L, EN2
         M = L + EN2 - MM
         Z = H(M,M)
         R = X - Z
         S = Y - Z
         P = (R*S-W)/H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - Z - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P/S
         Q = Q/S
         R = R/S
         IF (M.EQ.L) GO TO 320
         IF ((ABS(H(M,M-1))*(ABS(Q)+ABS(R))).LE.(MACHEP*ABS(P)
     *       *(ABS(H(M-1,M-1))+ABS(Z)+ABS(H(M+1,M+1))))) GO TO 320
  300 CONTINUE
  320 M2 = M + 2
      IF (M2.GT.EN) GO TO 360
      DO 340 I = M2, EN
         H(I,I-2) = 0.0D0
  340 CONTINUE
  360 M3 = M + 3
      IF (M3.GT.EN) GO TO 400
      DO 380 I = M3, EN
         H(I,I-3) = 0.0D0
  380 CONTINUE
  400 IF (M.GT.NA) GO TO 720
      DO 700 K = M, NA
         NOTLAS = K .NE. NA
         IF (K.EQ.M) GO TO 420
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X.EQ.0.0D0) GO TO 700
         P = P/X
         Q = Q/X
         R = R/X
  420    S = SQRT(P**2+Q**2+R**2)
         IF (P.LT.0.0D0) S = -S
         IF (K.NE.M) GO TO 440
         IF (L.NE.M) H(K,K-1) = -H(K,K-1)
         GO TO 460
  440    H(K,K-1) = -S*X
  460    P = P + S
         X = P/S
         Y = Q/S
         Z = R/S
         Q = Q/P
         R = R/P
C        ROW MODIFICATION
         IF (NOTLAS) GO TO 500
         DO 480 J = K, N
            P = H(K,J) + Q*H(K+1,J)
            H(K+1,J) = H(K+1,J) - P*Y
            H(K,J) = H(K,J) - P*X
  480    CONTINUE
         GO TO 540
  500    DO 520 J = K, N
            P = H(K,J) + Q*H(K+1,J) + R*H(K+2,J)
            H(K+2,J) = H(K+2,J) - P*Z
            H(K+1,J) = H(K+1,J) - P*Y
            H(K,J) = H(K,J) - P*X
  520    CONTINUE
  540    J = EN
         IF (K+3.LT.EN) J = K + 3
C        COLUMN MODIFICATION
         IF (NOTLAS) GO TO 580
         DO 560 I = 1, J
            P = X*H(I,K) + Y*H(I,K+1)
            H(I,K+1) = H(I,K+1) - P*Q
            H(I,K) = H(I,K) - P
  560    CONTINUE
         GO TO 620
  580    DO 600 I = 1, J
            P = X*H(I,K) + Y*H(I,K+1) + Z*H(I,K+2)
            H(I,K+2) = H(I,K+2) - P*R
            H(I,K+1) = H(I,K+1) - P*Q
            H(I,K) = H(I,K) - P
  600    CONTINUE
C        ACCUMULATE TRANSFORMATIONS
  620    IF (LOW.GT.UPP) GO TO 700
         IF (NOTLAS) GO TO 660
         DO 640 I = LOW, UPP
            P = X*VECS(I,K) + Y*VECS(I,K+1)
            VECS(I,K+1) = VECS(I,K+1) - P*Q
            VECS(I,K) = VECS(I,K) - P
  640    CONTINUE
         GO TO 700
  660    DO 680 I = LOW, UPP
            P = X*VECS(I,K) + Y*VECS(I,K+1) + Z*VECS(I,K+2)
            VECS(I,K+2) = VECS(I,K+2) - P*R
            VECS(I,K+1) = VECS(I,K+1) - P*Q
            VECS(I,K) = VECS(I,K) - P
  680    CONTINUE
  700 CONTINUE
  720 GO TO 160
C     ONE ROOT FOUND
  740 WR(EN) = X + T
      H(EN,EN) = WR(EN)
      WI(EN) = 0.0D0
      CNT(EN) = ITS
      EN = NA
      GO TO 140
C     TWO ROOTS FOUND
  760 P = (Y-X)/2.0D0
      Q = P**2 + W
      Z = SQRT(ABS(Q))
      X = X + T
      H(EN,EN) = X
      H(NA,NA) = Y + T
      CNT(EN) = -ITS
      CNT(NA) = ITS
      IF (Q.LT.0.0D0) GO TO 840
C     REAL PAIR
      IF (P.LT.0.0D0) Z = P - Z
      IF (P.GT.0.0D0) Z = P + Z
      WR(NA) = X + Z
      WR(EN) = WR(NA)
      IF (Z.NE.0.0D0) WR(EN) = X - W/Z
      WI(NA) = 0.0D0
      WI(EN) = 0.0D0
      X = H(EN,NA)
      R = A02ABF(X,Z)
      P = X/R
      Q = Z/R
C     ROW MODIFICATION
      DO 780 J = NA, N
         Z = H(NA,J)
         H(NA,J) = Q*Z + P*H(EN,J)
         H(EN,J) = Q*H(EN,J) - P*Z
  780 CONTINUE
C     COLUMN MODIFICATION
      DO 800 I = 1, EN
         Z = H(I,NA)
         H(I,NA) = Q*Z + P*H(I,EN)
         H(I,EN) = Q*H(I,EN) - P*Z
  800 CONTINUE
C     ACCUMULATE TRANSFORMATIONS
      DO 820 I = LOW, UPP
         Z = VECS(I,NA)
         VECS(I,NA) = Q*Z + P*VECS(I,EN)
         VECS(I,EN) = Q*VECS(I,EN) - P*Z
  820 CONTINUE
      GO TO 860
C     COMPLEX PAIR
  840 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = Z
      WI(EN) = -Z
  860 EN = EN - 2
      GO TO 140
C     ALL ROOTS FOUND NOW BACKSUBSTITUTE
  880 IF (NORM.EQ.0.0D0) GO TO 1480
      NORM = NORM*MACHEP
C     BACKSUBSTITUTION
      DO 1340 KK = 1, N
         EN = N + 1 - KK
         P = WR(EN)
         Q = WI(EN)
         NA = EN - 1
         IF (Q.NE.0.0D0) GO TO 1120
C        REAL VECTOR
         H(EN,EN) = 1.0D0
         IF (NA.LT.1) GO TO 1340
         DO 1100 II = 1, NA
            I = NA + 1 - II
            I1 = I - 1
            W = H(I,I) - P
            R = H(I,EN)
            IF (WI(I).GE.0.0D0) GO TO 900
            Z = W
            S = R
            GO TO 1100
  900       IF (WI(I).GT.0.0D0) GO TO 1020
C           MODIFICATION TO STOP OVERFLOW
            IF (W.NE.0.0D0) GO TO 940
            IF (ABS(R).LT.10.0D0*NORM) GO TO 960
            R = -R
            DO 920 J = 1, EN
               H(J,EN) = H(J,EN)*NORM
  920       CONTINUE
            GO TO 980
  940       R = -R/W
            GO TO 980
  960       R = -R/NORM
  980       H(I,EN) = R
            IF (I1.EQ.0) GO TO 1100
            DO 1000 J = 1, I1
               H(J,EN) = H(J,EN) + H(J,I)*R
 1000       CONTINUE
            GO TO 1100
C           SOLVE REAL EQUATIONS
 1020       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I)-P)**2 + WI(I)**2
            T = (X*S-Z*R)/Q
            H(I,EN) = T
            IF (ABS(X).GT.ABS(Z)) GO TO 1040
            R = (-S-Y*T)/Z
            GO TO 1060
 1040       R = (-R-W*T)/X
 1060       H(I+1,EN) = R
            IF (I1.EQ.0) GO TO 1100
            DO 1080 J = 1, I1
               H(J,EN) = (H(J,EN)+H(J,I+1)*R) + H(J,I)*T
 1080       CONTINUE
 1100    CONTINUE
C        END REAL VECTOR
         GO TO 1340
 1120    IF (Q.GT.0.0D0) GO TO 1340
C        COMPLEX VECTOR ASSOCIATED WITH LAMBDA=P-I*Q
         IF (ABS(H(EN,NA)).LE.ABS(H(NA,EN))) GO TO 1140
         R = Q/H(EN,NA)
         S = -(H(EN,EN)-P)/H(EN,NA)
         GO TO 1160
 1140    CALL A02ACF(0.0D0,-H(NA,EN),H(NA,NA)-P,Q,R,S)
 1160    H(EN,NA) = 0.0D0
         H(EN,EN) = 1.0D0
         H(NA,NA) = R
         H(NA,EN) = S
         IF (NA.LT.2) GO TO 1340
         NA1 = NA - 1
         DO 1180 J = 1, NA1
            H(J,EN) = H(J,EN) + H(J,NA)*S
            H(J,NA) = H(J,NA)*R
 1180    CONTINUE
         DO 1320 II = 1, NA1
            I = 1 + NA1 - II
            I1 = I - 1
            W = H(I,I) - P
            RA = H(I,NA)
            SA = H(I,EN)
            IF (WI(I).GE.0.0D0) GO TO 1200
            Z = W
            R = RA
            S = SA
            GO TO 1320
 1200       IF (WI(I).EQ.0.0D0) GO TO 1280
C           SOLVE COMPLEX EQUATIONS
            X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I)-P)**2 + WI(I)**2 - Q**2
            VI = (WR(I)-P)*2.0D0*Q
            IF ((VR.EQ.0.0D0) .AND. (VI.EQ.0.0D0))
     *          VR = MACHEP*NORM*(ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(Z))
            CALL A02ACF(X*R-Z*RA+Q*SA,X*S-Z*SA-Q*RA,VR,VI,T,U)
            IF (ABS(X).LE.ABS(Z)+ABS(Q)) GO TO 1220
            R = (-RA-W*T+Q*U)/X
            S = (-SA-W*U-Q*T)/X
            GO TO 1240
 1220       CALL A02ACF(-R-Y*T,-S-Y*U,Z,Q,R,S)
 1240       H(I,NA) = T
            H(I,EN) = U
            H(I+1,NA) = R
            H(I+1,EN) = S
            IF (I1.EQ.0) GO TO 1320
            DO 1260 J = 1, I1
               H(J,NA) = (H(J,NA)+H(J,I+1)*R) + H(J,I)*T
               H(J,EN) = (H(J,EN)+H(J,I+1)*S) + H(J,I)*U
 1260       CONTINUE
            GO TO 1320
 1280       CALL A02ACF(-RA,-SA,W,Q,R,S)
            H(I,NA) = R
            H(I,EN) = S
            IF (I1.EQ.0) GO TO 1320
            DO 1300 J = 1, I1
               H(J,NA) = H(J,NA) + H(J,I)*R
               H(J,EN) = H(J,EN) + H(J,I)*S
 1300       CONTINUE
 1320    CONTINUE
C        END COMPLEX VECTOR
 1340 CONTINUE
C     END BACKSUBSTITUTION
C     VECTORS OF ISOLATED ROOTS
      LOW1 = LOW - 1
      UPP1 = UPP + 1
      DO 1420 J = 1, N
         M = MIN(J,LOW1)
         IF (M.LT.1) GO TO 1380
         DO 1360 I = 1, M
            VECS(I,J) = H(I,J)
 1360    CONTINUE
 1380    IF (UPP1.GT.J) GO TO 1420
         DO 1400 I = UPP1, J
            VECS(I,J) = H(I,J)
 1400    CONTINUE
 1420 CONTINUE
C     MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C     VECTORS OF ORIGINAL FULL MATRIX
      DO 1460 JJ = LOW, N
         J = LOW + N - JJ
         M = MIN(J,UPP)
         DO 1440 I = LOW, UPP
            VECS(I,J) = VECS(I,M)*H(M,J)
 1440    CONTINUE
         M = M - 1
         IF (M+1.GE.LOW) CALL DGEMV_('N',UPP-LOW+1,M-LOW+1,1.0D0,
     *                              VECS(LOW,LOW),IVECS,H(LOW,J),1,
     *                              1.0D0,VECS(LOW,J),1)
 1460 CONTINUE
 1480 IFAIL = 0
      RETURN
 1500 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F02SZF(N,D,E,SV,WANTB,B,WANTY,Y,NRY,LY,WANTZ,Z,NRZ,NCZ,
     *                  WORK1,WORK2,WORK3,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-328 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. IER-518 (AUG 1986).
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (SVDBID)
C
C     F02SZF RETURNS PART OR ALL OF THE SINGULAR VALUE
C     DECOMPOSITION OF THE N*N UPPER BIDIAGONAL MATRIX A. THAT
C     IS, A IS FACTORIZED AS
C
C     A = Q*DIAG(SV)*(P**T) ,
C
C     WHERE Q AND P ARE N*N ORTHOGONAL MATRICES AND DIAG(SV)
C     IS AN N*N DIAGONAL MATRIX WITH NON-NEGATIVE DIAGONAL
C     ELEMENTS SV(1),SV(2),..., SV(N), THESE BEING THE
C     SINGULAR VALUES OF A.
C
C     IF WANTB IS .TRUE. THEN B RETURNS (Q**T)*B.
C     IF WANTY IS .TRUE. THEN Y RETURNS Y*Q.
C     IF WANTZ IS .TRUE. THEN Z RETURNS (P**T)*Z.
C
C     INPUT PARAMETERS.
C
C     N     - THE ORDER OF THE MATRIX. MUST BE AT LEAST 1.
C
C     D     - N ELEMENT VECTOR SUCH THAT D(I)=A(I,I), I=1,2,...,N.
C             D IS UNALTERED UNLESS ROUTINE IS CALLED WITH SV=D.
C
C     E     - N ELEMENT VECTOR SUCH THAT E(I)=A(I-1,I), I=2,3,...,N.
C             E(1) IS NOT REFERENCED.
C             E IS UNALTERED UNLESS ROUTINE IS CALLED WITH WORK1=E.
C
C     WANTB - MUST BE .TRUE. IF (Q**T)*B IS REQUIRED.
C             IF WANTB IS .FALSE. THEN B IS NOT REFERENCED.
C
C     B     - AN N ELEMENT REAL VECTOR.
C
C     WANTY - MUST BE .TRUE. IF Y*Q IS REQUIRED.
C             IF WANTY IS .FALSE. THEN Y IS NOT REFERENCED.
C
C     Y     - AN LY*N REAL MATRIX.
C
C     NRY   - IF WANTY IS .TRUE. THEN NRY MUST BE THE ROW
C             DIMENSION OF Y AS DECLARED IN THE CALLING
C             PROGRAM AND MUST BE AT LEAST LY.
C
C     LY    - IF WANTY IS .TRUE. THEN LY MUST BE THE NUMBER
C             OF ROWS OF Y AND MUST BE AT LEAST 1.
C
C     WANTZ - MUST BE .TRUE. IF (P**T)*Z IS REQUIRED.
C             IF WANTZ IS .FALSE. THEN Z IS NOT REFERENCED.
C
C     Z     - AN N*NCZ REAL MATRIX.
C
C     NRZ   - IF WANTZ IS .TRUE. THEN NRZ MUST BE THE ROW
C             DIMENSION OF Z AS DECLARED IN THE CALLING
C             PROGRAM AND MUST BE AT LEAST N.
C
C     NCZ   - IF WANTZ IS .TRUE. THEN NCZ MUST BE THE
C             NUMBER OF COLUMNS OF Z AND MUST BE AT LEAST
C             1.
C
C     IFAIL - THE USUAL FAILURE PARAMETER. IF IN DOUBT SET
C             IFAIL TO ZERO BEFORE CALLING F02SZF.
C
C     OUTPUT PARAMETERS.
C
C     SV    - N ELEMENT VECTOR CONTAINING THE SINGULAR
C             VALUES OF A. THEY ARE ORDERED SO THAT
C             SV(1).GE.SV(2).GE. ... .GE.SV(N). THE ROUTINE
C             MAY BE CALLED WITH SV=D.
C
C     B     - IF WANTB IS .TRUE. THEN B WILL RETURN THE N
C             ELEMENT VECTOR (Q**T)*B.
C
C     Y     - IF WANTY IS .TRUE. THEN Y WILL RETURN THE
C             LY*N MATRIX Y*Q.
C
C     Z     - IF WANTZ IS .TRUE. THEN Z WILL RETURN THE N*NCZ MATRIX
C             (P**T)*Z.
C
C     IFAIL - ON NORMAL RETURN IFAIL WILL BE ZERO.
C             IN THE UNLIKELY EVENT THAT THE QR-ALGORITHM
C             FAILS TO FIND THE SINGULAR VALUES IN 50*N
C             ITERATIONS THEN IFAIL WILL BE 2 OR MORE AND
C             SUCH THAT SV(1),SV(2),..,SV(IFAIL-1) MAY NOT
C             HAVE BEEN FOUND. SEE WORK1 BELOW. THIS
C             FAILURE IS NOT LIKELY TO OCCUR.
C             IF AN INPUT PARAMETER IS INCORRECTLY SUPPLIED
C             THEN IFAIL IS SET TO UNITY.
C
C     WORKSPACE PARAMETERS.
C
C     WORK1 - AN N ELEMENT VECTOR. IF E IS NOT REQUIRED ON
C             RETURN THEN THE ROUTINE MAY BE CALLED WITH
C             WORK1=E. WORK1(1) RETURNS THE TOTAL NUMBER OF
C             ITERATIONS TAKEN BY THE  QR-ALGORITHM. IF
C             IFAIL IS POSITIVE ON RETURN THEN THE MATRIX A
C             IS GIVEN  BY A=Q*C*(P**T) , WHERE C IS THE
C             UPPER BIDIAGONAL MATRIX WITH SV AS ITS
C             DIAGONAL AND WORK1 AS ITS SUPER-DIAGONAL.
C
C     WORK2
C     WORK3 - N ELEMENT VECTORS. IF WANTZ IS .FALSE. THEN WORK2 AND
C             WORK3 ARE NOT REFERENCED.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02SZF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LY, N, NCZ, NRY, NRZ
      LOGICAL           WANTB, WANTY, WANTZ
C     .. Array Arguments ..
      DOUBLE PRECISION  B(N), D(N), E(N), SV(N), WORK1(N), WORK2(N),
     *                  WORK3(N), Y(NRY,N), Z(NRZ,NCZ)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANORM, BIG, C, DK, DKM1, DL, EK, EKM1, EPS, F,
     *                  G, RSQTPS, S, SHUFT, SMALL, SQTEPS, SVI, T, X
      INTEGER           I, IERR, ITER, J, JJ, K, KK, L, LL, LM1, LP1,
     *                  MAXIT
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  F01LZZ, X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          F01LZZ, X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01LZW, F01LZY, F02SZZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
      IERR = IFAIL
      IF (IERR.EQ.0) IFAIL = 1
C
      IF (N.LT.1) GO TO 500
      IF (WANTY .AND. (NRY.LT.LY .OR. LY.LT.1)) GO TO 500
      IF (WANTZ .AND. (NRZ.LT.N .OR. NCZ.LT.1)) GO TO 500
C
      SMALL = X02AMF()
      BIG = 1.0D0/SMALL
      EPS = X02AJF()
      SQTEPS = SQRT(EPS)
      RSQTPS = 1.0D0/SQTEPS
C
      ITER = 0
      K = N
      SV(1) = D(1)
      ANORM = ABS(D(1))
      IF (N.EQ.1) GO TO 280
C
      DO 20 I = 2, N
         SV(I) = D(I)
         WORK1(I) = E(I)
         ANORM = MAX(ANORM,ABS(D(I)),ABS(E(I)))
   20 CONTINUE
C
      MAXIT = 50*N
      EPS = EPS*ANORM
C
C     MAXIT IS THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C     EPS WILL BE USED TO TEST FOR NEGLIGIBLE ELEMENTS.
C     START MAIN LOOP. ONE SINGULAR VALUE IS FOUND FOR EACH
C     VALUE OF K. K GOES IN OPPOSITE DIRECTION TO KK.
C
      DO 260 KK = 2, N
C
C        NOW TEST FOR SPLITTING. L GOES IN OPPOSITE DIRECTION TO LL.
C
   40    L = K
         DO 60 LL = 2, K
            IF (ABS(WORK1(L)).LE.EPS) GO TO 240
            L = L - 1
            IF (ABS(SV(L)).LT.EPS) GO TO 180
   60    CONTINUE
C
   80    IF (ITER.EQ.MAXIT) GO TO 280
C
C        MAXIT QR-STEPS WITHOUT CONVERGENCE. FAILURE.
C
         ITER = ITER + 1
C
C        NOW DETERMINE SHIFT.
C
         LP1 = L + 1
         DL = SV(L)
         DKM1 = SV(K-1)
         DK = SV(K)
         EKM1 = 0.0D0
         IF (K.NE.2) EKM1 = WORK1(K-1)
         EK = WORK1(K)
         F = (DKM1-DK)*(DKM1+DK) + (EKM1-EK)*(EKM1+EK)
         F = F/(2.0D0*EK*DKM1)
         G = ABS(F)
         IF (G.LE.RSQTPS) G = SQRT(1.0D0+F**2)
         IF (F.LT.0.0D0) G = -G
C
         SHUFT = EK*(EK-DKM1/(F+G))
         F = (DL-DK)*(DL+DK) - SHUFT
         X = DL*WORK1(LP1)
C
C        NOW PERFORM THE QR-STEP AND CHASE ZEROS.
C
         DO 140 I = LP1, K
C
            T = F01LZZ(F,X,SMALL,BIG)
C
            CALL F01LZW(T,C,S,SQTEPS,RSQTPS,BIG)
C
            IF (I.GT.LP1) WORK1(I-1) = C*F + S*X
            F = C*SV(I-1) + S*WORK1(I)
            WORK1(I) = C*WORK1(I) - S*SV(I-1)
            X = S*SV(I)
            SVI = C*SV(I)
C
            IF ( .NOT. WANTZ) GO TO 100
            WORK2(I) = C
            WORK3(I) = S
C
  100       T = F01LZZ(F,X,SMALL,BIG)
C
            CALL F01LZW(T,C,S,SQTEPS,RSQTPS,BIG)
C
            IF (WANTY) CALL F01LZY(LY,C,S,Y(1,I-1),Y(1,I))
C
            IF ( .NOT. WANTB) GO TO 120
            T = B(I)
            B(I) = C*T - S*B(I-1)
            B(I-1) = C*B(I-1) + S*T
C
  120       SV(I-1) = C*F + S*X
            F = C*WORK1(I) + S*SVI
            SV(I) = C*SVI - S*WORK1(I)
C
            IF (I.EQ.K) GO TO 140
            X = S*WORK1(I+1)
            WORK1(I+1) = C*WORK1(I+1)
C
  140    CONTINUE
C
         WORK1(K) = F
         IF ( .NOT. WANTZ) GO TO 40
         DO 160 J = 1, NCZ
C
            CALL F02SZZ(K-L+1,WORK2(L),WORK3(L),Z(L,J))
C
  160    CONTINUE
         GO TO 40
C
C        COME TO NEXT PIECE IF SV(L-1) IS NEGLIGIBLE. FORCE A SPLIT.
C
  180    LM1 = L
         L = L + 1
         X = WORK1(L)
         WORK1(L) = 0.0D0
         DO 220 I = L, K
C
            T = F01LZZ(SV(I),X,SMALL,BIG)
C
            CALL F01LZW(T,C,S,SQTEPS,RSQTPS,BIG)
C
            IF (WANTY) CALL F01LZY(LY,C,-S,Y(1,LM1),Y(1,I))
C
            IF ( .NOT. WANTB) GO TO 200
            T = B(I)
            B(I) = C*T + S*B(LM1)
            B(LM1) = C*B(LM1) - S*T
C
  200       SV(I) = C*SV(I) + S*X
            IF (I.EQ.K) GO TO 220
            X = -S*WORK1(I+1)
            WORK1(I+1) = C*WORK1(I+1)
C
  220    CONTINUE
C
C        IF WE COME HERE WITH L=K THEN A SINGULAR VALUE HAS BEEN
C        FOUND.
C
  240    IF (L.LT.K) GO TO 80
C
         K = K - 1
  260 CONTINUE
C
  280 IFAIL = K - 1
      WORK1(1) = ITER
C
C     NOW MAKE SINGULAR VALUES NON-NEGATIVE.
C     K WILL BE 1 UNLESS FAILURE HAS OCCURED.
C
      DO 320 J = K, N
         IF (SV(J).GE.0.0D0) GO TO 320
C
         SV(J) = -SV(J)
C
         IF (WANTB) B(J) = -B(J)
         IF ( .NOT. WANTY) GO TO 320
         DO 300 I = 1, LY
            Y(I,J) = -Y(I,J)
  300    CONTINUE
C
  320 CONTINUE
C
C     NOW SORT THE SINGULAR VALUES INTO DESCENDING ORDER.
C
      IF (WANTZ) JJ = 0
      DO 400 J = K, N
         S = 0.0D0
         L = J
C
         DO 340 I = J, N
            IF (SV(I).LE.S) GO TO 340
            S = SV(I)
            L = I
  340    CONTINUE
C
         IF (S.EQ.0.0D0) GO TO 420
         IF (WANTZ) WORK2(J) = L
         IF (L.EQ.J) GO TO 400
         IF (WANTZ) JJ = J
C
         SV(L) = SV(J)
         SV(J) = S
         IF ( .NOT. WANTY) GO TO 380
C
         DO 360 I = 1, LY
            T = Y(I,J)
            Y(I,J) = Y(I,L)
            Y(I,L) = T
  360    CONTINUE
C
  380    IF ( .NOT. WANTB) GO TO 400
         T = B(J)
         B(J) = B(L)
         B(L) = T
C
  400 CONTINUE
C
  420 IF ( .NOT. WANTZ) GO TO 480
      IF (JJ.EQ.0) GO TO 480
      DO 460 I = 1, NCZ
         DO 440 J = K, JJ
            L = WORK2(J)
            IF (J.EQ.L) GO TO 440
            T = Z(J,I)
            Z(J,I) = Z(L,I)
            Z(L,I) = T
  440    CONTINUE
  460 CONTINUE
C
  480 IF (IFAIL.EQ.0) RETURN
C
      IFAIL = IFAIL + 1
  500 IFAIL = P01ABF(IERR,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F02SZZ(N,C,S,X)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (PLRT10)
C
C     F02SZZ RETURNS THE N ELEMENT VECTOR
C
C     Y = R(N-1,N)*R(N-2,N-1)*...*R(1,2)*X ,
C
C     WHERE X IS AN N ELEMENT VECTOR AND R(J-1,J) IS A PLANE
C     ROTATION FOR THE (J-1,J)-PLANE.
C
C     Y IS OVERWRITTEN ON X.
C
C     THE N ELEMENT VECTORS C AND S MUST BE SUCH THAT THE
C     NON-IDENTITY PART OF R(J-1,J) IS GIVEN BY
C
C     R(J-1,J) = (  C(J)  S(J) ) .
C                ( -S(J)  C(J) )
C
C     C(1) AND S(1) ARE NOT REFERENCED.
C
C
C     N MUST BE AT LEAST 1. IF N=1 THEN AN IMMEDIATE RETURN TO
C     THE CALLING PROGRAM IS MADE.
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(N), S(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  W
      INTEGER           I
C     .. Executable Statements ..
      IF (N.EQ.1) RETURN
C
      DO 20 I = 2, N
         W = X(I-1)
         X(I-1) = C(I)*W + S(I)*X(I)
         X(I) = C(I)*X(I) - S(I)*W
   20 CONTINUE
C
      RETURN
      END
      SUBROUTINE F02WAF(M,N,A,NRA,WANTB,B,SV,WORK,LWORK,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 15 REVISED. IER-912 (APR 1991).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (SVDGN1)
C     Modified by Sven to replace calls to F01QAF and F02WAZ.
C     7-Feb-1991.
C
C     F02WAF RETURNS PART OF THE SINGULAR VALUE DECOMPOSITION
C     OF THE M*N (M.GE.N) MATRIX A GIVEN BY
C
C     A = Q*(D)*(P**T) ,
C           (0)
C
C     WHERE Q AND P ARE ORTHOGONAL MATRICES AND D IS AN N*N
C     DIAGONAL MATRIX WITH NON-NEGATIVE DIAGONAL ELEMENTS,
C     THESE BEING THE SINGULAR VALUES OF A.
C
C     P**T AND THE DIAGONAL ELEMENTS OF D ARE RETURNED.
C     IF WANTB IS .TRUE. THEN (Q**T)*B IS ALSO RETURNED.
C
C     INPUT PARAMETERS.
C
C     M     - NUMBER OF ROWS OF A. M MUST BE AT LEAST N.
C
C     N     - NUMBER OF COLUMNS OF A. N MUST BE AT LEAST
C             UNITY AND MUST NOT BE LARGER THAN THAN M.
C
C     A     - THE M*N MATRIX TO BE FACTORIZED.
C
C     NRA   - ROW DIMENSION OF A AS DECLARED IN THE CALLING PROGRAM.
C             NRA MUST BE AT LEAST M.
C
C     WANTB - MUST BE .TRUE. IF (Q**T)*B IS REQUIRED.
C             IF WANTB IS .FALSE. THEN B IS NOT REFERENCED.
C
C     B     - AN M ELEMENT VECTOR.
C
C     IFAIL - THE USUAL FAILURE PARAMETER. IF IN DOUBT SET
C             IFAIL TO ZERO BEFORE CALLING F02WAF.
C
C     OUTPUT PARAMETERS.
C
C     A     - THE TOP N*N PART OF A WILL CONTAIN THE N*N ORTHOGONAL
C             MATRIX P**T.
C             THE REMAINING (M-N)*N PART OF A IS USED FOR INTERNAL
C             WORKSPACE.
C
C     B     - IF WANTB IS .TRUE. THEN B IS OVERWRITTEN BY
C             THE M ELEMENT VECTOR (Q**T)*B.
C
C     SV    - N ELEMENT VECTOR CONTAINING THE SINGULAR
C             VALUES OF A. THEY ARE ORDERED SO THAT
C             SV(1).GE.SV(2).GE. ... .GE.SV(N).GE.0.
C
C     IFAIL - ON NORMAL RETURN IFAIL WILL BE ZERO.
C             IN THE UNLIKELY EVENT THAT THE QR-ALGORITHM
C             FAILS TO FIND THE SINGULAR VALUES IN 50*N
C             ITERATIONS THEN IFAIL IS SET TO 2.
C             IF AN INPUT PARAMETER IS INCORRECTLY SUPPLIED
C             THEN IFAIL IS SET TO UNITY.
C
C     WORKSPACE PARAMETERS.
C
C     WORK  - A 3*N ELEMENT VECTOR.
C             WORK(1) RETURNS THE TOTAL NUMBER OF ITERATIONS TAKEN
C             BY THE QR-ALGORITHM.
C
C     LWORK - THE LENGTH OF THE VECTOR WORK. LWORK MUST BE
C             AT LEAST 3*N.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02WAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LWORK, M, N, NRA
      LOGICAL           WANTB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NRA,N), B(M), SV(N), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           IERR, NP1, NPNP1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01LZF, F01QCF, F02SZF, F02WAY, F01QDF
C     .. Executable Statements ..
      IERR = IFAIL
      IF (IERR.EQ.0) IFAIL = 1
C
      IF (NRA.LT.M .OR. M.LT.N .OR. LWORK.LT.3*N .OR. N.LT.1)
     *    GO TO 20
C
      NP1 = N + 1
      NPNP1 = N + NP1
C
C     Call to F01QAF replaced by a call to F01QCF.
C     CALL F01QAF(M,N,A,NRA,A,NRA,WORK,IFAIL)
      CALL F01QCF(M,N,A,NRA,WORK,IFAIL)
C
C     Call to F02WAZ replaced by a call to F01QDF.
C     IF (WANTB) CALL F02WAZ(M,N,A,NRA,WORK,B,B)
      IF (WANTB) CALL F01QDF('Transpose','Separate',M,N,A,NRA,WORK,1,B,
     *                       M,WORK(N+1),IFAIL)
C
      CALL F01LZF(N,A,NRA,A,NRA,WANTB,B,.FALSE.,.FALSE.,WORK,1,1,
     *            .FALSE.,WORK,1,1,SV,WORK,WORK,WORK,IFAIL)
C
      CALL F02WAY(N,A,NRA,A,NRA)
C
      IFAIL = 1
      CALL F02SZF(N,SV,WORK,SV,WANTB,B,.FALSE.,WORK,1,1,.TRUE.,A,NRA,N,
     *            WORK,WORK(NP1),WORK(NPNP1),IFAIL)
C
      IF (IFAIL.EQ.0) RETURN
C
      IFAIL = 2
   20 IFAIL = P01ABF(IERR,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F02WAY(N,C,NRC,PT,NRPT)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. IER-519 (AUG 1986).
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (BIGVPT)
C
C     F02WAY RETURNS THE N*N ORTHOGONAL MATRIX P**T FOR THE
C     FACTORIZATION OF ROUTINE F01LZF.
C
C     DETAILS OF P MUST BE SUPPLIED IN THE N*N MATRIX C AS
C     RETURNED FROM ROUTINE F01LZF.
C
C     NRC AND NRPT MUST BE THE ROW DIMENSIONS OF C AND PT
C     RESPECTIVELY AS DECLARED IN THE CALLING PROGRAM AND MUST
C     EACH BE AT LEAST N.
C
C     THE ROUTINE MAY BE CALLED WITH PT=C.
C
C     .. Scalar Arguments ..
      INTEGER           N, NRC, NRPT
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NRC,N), PT(NRPT,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG, CS, RSQTPS, SN, SQTEPS, T
      INTEGER           I, J, K, KK, KM1, KP1
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          F01LZW, F01LZY
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      BIG = 1.0D0/X02AMF()
      SQTEPS = SQRT(X02AJF())
      RSQTPS = 1.0D0/SQTEPS
      DO 15 J = 3, N
         DO 5 I = 1, J - 2
            PT(I,J) = C(I,J)
    5    CONTINUE
   15 CONTINUE
C
      PT(N,N) = 1.0D0
      IF (N.EQ.1) RETURN
C
      PT(N-1,N) = 0.0D0
      PT(N,N-1) = 0.0D0
      PT(N-1,N-1) = 1.0D0
      IF (N.EQ.2) RETURN
C
      K = N
      DO 60 KK = 3, N
         KP1 = K
         K = K - 1
         KM1 = K - 1
         PT(KM1,K) = 0.0D0
C
         DO 20 J = KP1, N
            T = PT(KM1,J)
            PT(KM1,J) = 0.0D0
            IF (T.EQ.0.0D0) GO TO 20
C
            CALL F01LZW(-T,CS,SN,SQTEPS,RSQTPS,BIG)
C
            CALL F01LZY(N-KM1,CS,SN,PT(K,J-1),PT(K,J))
C
   20    CONTINUE
C
         PT(KM1,KM1) = 1.0D0
         DO 40 I = K, N
            PT(I,KM1) = 0.0D0
   40    CONTINUE
C
   60 CONTINUE
C
      RETURN
      END
      INTEGER FUNCTION F02WDY(N,SV,TOL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (IRKSVD)
C
C     F02WDY RETURNS THE RANK OF AN M*K MATRIX A FOLLOWING A
C     SINGULAR VALUE DECOMPOSITION OF A.
C
C     THE N=MIN(M,K) SINGULAR VALUES OF A MUST BE IN
C     DESCENDING ORDER IN THE N ELEMENT VECTOR SV. THEN F02WDY
C     RETURNS THE LARGEST INTEGER SUCH THAT
C
C     SV(F02WDY) .GT. TOL*SV(1) .
C
C     IF SV(1)=0 THEN F02WDY IS RETURNED AS ZERO.
C
C     IF TOL.LT.EPS OR TOL.GE.1 THEN THE VALUE EPS IS USED IN
C     PLACE OF TOL, WHERE EPS IS THE SMALLEST REAL FOR WHICH
C     1.0+EPS.GT.1.0 ON THE MACHINE. FOR MOST PROBLEMS THIS IS
C     UNREASONABLY SMALL AND TOL SHOULD BE CHOSEN TO
C     APPROXIMATE THE RELATIVE ERRORS IN THE ELEMENTS OF A.
C
C     IF INSTEAD SINGULAR VALUES BELOW SOME VALUE DELTA ARE TO
C     BE REGARDED AS ZERO THEN SUPPLY TOL AS DELTA/SV(1).
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION        TOL
      INTEGER                 N
C     .. Array Arguments ..
      DOUBLE PRECISION        SV(N)
C     .. Local Scalars ..
      DOUBLE PRECISION        DELTA, TL
      INTEGER                 I, IR
C     .. External Functions ..
      DOUBLE PRECISION        X02AJF
      EXTERNAL                X02AJF
C     .. Executable Statements ..
      TL = TOL
      DELTA = X02AJF()
      IF (TL.LT.DELTA .OR. TL.GE.1.0D0) TL = DELTA
C
      IR = 0
      DELTA = TL*SV(1)
C
      DO 20 I = 1, N
         IF (SV(I).LE.DELTA) GO TO 40
         IR = I
   20 CONTINUE
C
   40 F02WDY = IR
      RETURN
      END
      SUBROUTINE F04AEF(A,IA,B,IB,N,M,C,IC,WKSPCE,AA,IAA,BB,IBB,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Accurate solution of a set of real linear equations
C     with multiple right hand sides.
C     1st August 1971
C
C     Rewritten to call F07ADG, a modified version of LAPACK routine
C     SGETRF/F07ADF; new IFAIL exit inserted for illegal input
C     parameters; error messages inserted. February 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04AEF')
C     .. Scalar Arguments ..
      INTEGER           IA, IAA, IB, IBB, IC, IFAIL, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), AA(IAA,*), B(IB,*), BB(IBB,*), C(IC,*),
     *                  WKSPCE(*)
C     .. Local Scalars ..
      INTEGER           I, IERR, IFAIL1, INFO, J, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F04AHF, F07ADG
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      IERR = 0
      NREC = 0
      IF (N.LT.0) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (M.LT.0) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99998) M
      ELSE IF (IA.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99997) IA, N
      ELSE IF (IB.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99996) IB, N
      ELSE IF (IC.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99995) IC, N
      ELSE IF (IAA.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99994) IAA, N
      ELSE IF (IBB.LT.MAX(1,N)) THEN
         IERR = 3
         NREC = 1
         WRITE (P01REC,FMT=99993) IBB, N
      ELSE
C
C        Copy A to working array AA
C
         DO 40 I = 1, N
            DO 20 J = 1, N
               AA(I,J) = A(I,J)
   20       CONTINUE
   40    CONTINUE
C
         CALL F07ADG(N,N,AA,IAA,WKSPCE,INFO)
C
         IF (INFO.EQ.0) THEN
            IF (N.GT.0 .AND. M.GT.0) THEN
C
               IFAIL1 = 1
               CALL F04AHF(N,M,A,IA,AA,IAA,WKSPCE,B,IB,X02AJF(),C,IC,BB,
     *                     IBB,I,IFAIL1)
C
               IF (IFAIL1.NE.0) THEN
                  IERR = 2
                  NREC = 1
                  WRITE (P01REC,FMT=99991)
               END IF
C
            END IF
C
         ELSE
            IERR = 1
            NREC = 1
            WRITE (P01REC,FMT=99992)
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.0: N =',I16)
99998 FORMAT (1X,'** On entry, M.lt.0: M =',I16)
99997 FORMAT (1X,'** On entry, IA.lt.max(1,N): IA =',I16,', N =',I16)
99996 FORMAT (1X,'** On entry, IB.lt.max(1,N): IB =',I16,', N =',I16)
99995 FORMAT (1X,'** On entry, IC.lt.max(1,N): IC =',I16,', N =',I16)
99994 FORMAT (1X,'** On entry, IAA.lt.max(1,N): IAA =',I16,', N =',I16)
99993 FORMAT (1X,'** On entry, IBB.lt.max(1,N): IBB =',I16,', N =',I16)
99992 FORMAT (1X,'** Matrix A is approximately singular.')
99991 FORMAT (1X,'** Matrix A is too ill-conditioned.')
      END
      SUBROUTINE F04AHF(N,IR,A,IA,AA,IAA,P,B,IB,EPS,X,IX,BB,IBB,L,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     UNSYMACCSOLVE
C     SOLVES AX=B WHERE A IS AN N*N UNSYMMETRIC MATRIX AND B IS AN
C     N*IR MATRIX OF RIGHT HAND SIDES, USING THE SUBROUTINE F04AJF.
C     THE SUBROUTINE MUST BY PRECEDED BY F03AFF IN WHICH L AND U
C     ARE PRODUCED IN AA(I,J) AND THE INTERCHANGES IN P(I). THE
C     RESIDUALS BB=B-AX ARE CALCULATED AND AD=BB IS SOLVED, OVER-
C     WRITING D ON BB. THE REFINEMENT IS REPEATED, AS LONG AS THE
C     MAXIMUM CORRECTION AT ANY STAGE IS LESS THAN HALF THAT AT THE
C     PREVIOUS STAGE, UNTIL THE MAXIMUM CORRECTION IS LESS THAN 2
C     EPS TIMES THE MAXIMUM X. SETS IFAIL = 1 IF THE SOLUTION FAILS
C     TO IMPROVE, ELSE IFAIL = 0. L IS THE NUMBER OF ITERATIONS.
C     ADDITIONAL PRECISION INNERPRODUCTS ARE ABSOLUTELY NECESSARY.
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04AHF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           IA, IAA, IB, IBB, IFAIL, IR, IX, L, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), AA(IAA,N), B(IB,IR), BB(IBB,IR), P(N),
     *                  X(IX,IR)
C     .. Local Scalars ..
      DOUBLE PRECISION  BBMAX, D0, D1, D11, D2, XMAX
      INTEGER           I, ID2, IFAIL1, ISAVE, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F04AJF, X03AAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL1 = 0
      DO 40 J = 1, IR
         DO 20 I = 1, N
            X(I,J) = 0.0D0
            BB(I,J) = B(I,J)
   20    CONTINUE
   40 CONTINUE
      L = 0
      D0 = 0.0D0
   60 CALL F04AJF(N,IR,AA,IAA,P,BB,IBB)
      L = L + 1
      ID2 = 0
      D1 = 0.0D0
      DO 100 J = 1, IR
         DO 80 I = 1, N
            X(I,J) = X(I,J) + BB(I,J)
   80    CONTINUE
  100 CONTINUE
      DO 140 J = 1, IR
         XMAX = 0.0D0
         BBMAX = 0.0D0
         DO 120 I = 1, N
            IF (ABS(X(I,J)).GT.XMAX) XMAX = ABS(X(I,J))
            IF (ABS(BB(I,J)).GT.BBMAX) BBMAX = ABS(BB(I,J))
            CALL X03AAF(A(I,1),N*IA-I+1,X(1,J),(IR-J+1)
     *                  *IX,N,IA,1,-B(I,J),0.0D0,D11,D2,.TRUE.,IFAIL1)
            BB(I,J) = -D11
  120    CONTINUE
         IF (BBMAX.GT.D1*XMAX) D1 = BBMAX/XMAX
         IF (BBMAX.GT.2.0D0*EPS*XMAX) ID2 = 1
  140 CONTINUE
      IF ((D1.GT.0.5D0*D0) .AND. (L.NE.1)) GO TO 160
      D0 = D1
      IF (ID2.EQ.1) GO TO 60
      IFAIL = 0
      RETURN
  160 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F04AJF(N,IR,A,IA,P,B,IB)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. EXTENDED BLAS (JUNE 1986)
C
C     UNSYMSOL
C     SOLVES AX=B, WHERE A IS AN UNSYMMETRIC MATRIX AND B IS AN
C     N*IR
C     MATRIX OF IR RIGHT-HAND SIDES. THE SUBROUTINE F04AJF MUST BY
C     PRECEDED BY F03AFF IN WHICH L AND U ARE PRODUCED IN A(I,J),
C     FROM A, AND THE RECORD OF THE INTERCHANGES IS PRODUCED IN
C     P(I). AX=B IS SOLVED IN THREE STEPS, INTERCHANGE THE
C     ELEMENTS OF B, LY=B AND UX=Y. THE MATRICES Y AND THEN X ARE
C     OVERWRITTEN ON B.
C     1ST AUGUST 1971
C
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,IR), P(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  X
      INTEGER           I, I1, K
C     .. External Subroutines ..
      EXTERNAL          DTRSV
C     .. Executable Statements ..
C     INTERCHANGING ELEMENTS OF B
      DO 40 I = 1, N
         I1 = P(I) + 0.5D0
         IF (I1.EQ.I) GO TO 40
         DO 20 K = 1, IR
            X = B(I,K)
            B(I,K) = B(I1,K)
            B(I1,K) = X
   20    CONTINUE
   40 CONTINUE
      DO 60 K = 1, IR
C        SOLUTION OF LY= B
         CALL DTRSV('L','N','N',N,A,IA,B(1,K),1)
C        SOLUTION OF UX= Y
         CALL DTRSV('U','N','U',N,A,IA,B(1,K),1)
   60 CONTINUE
      RETURN
      END
      SUBROUTINE F04JAF(M,N,A,NRA,B,TOL,SIGMA,IRANK,WORK,LWORK,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (SVDLS1)
C
C     F04JAF RETURNS THE N ELEMENT VECTOR X, OF MINIMAL
C     LENGTH, THAT MINIMIZES THE EUCLIDEAN LENGTH OF THE M
C     ELEMENT VECTOR R GIVEN BY
C
C     R = B-A*X ,
C
C     WHERE A IS AN M*N (M.GE.N) MATRIX AND B IS AN M ELEMENT
C     VECTOR. X IS OVERWRITTEN ON B.
C
C     THE SOLUTION IS OBTAINED VIA A SINGULAR VALUE
C     DECOMPOSITION (SVD) OF THE MATRIX A GIVEN BY
C
C     A = Q*(D)*(P**T) ,
C           (0)
C
C     WHERE Q AND P ARE ORTHOGONAL AND D IS A DIAGONAL MATRIX WITH
C     NON-NEGATIVE DIAGONAL ELEMENTS, THESE BEING THE SINGULAR
C     VALUES OF A.
C
C     INPUT PARAMETERS.
C
C     M     - NUMBER OF ROWS OF A. M MUST BE AT LEAST N.
C
C     N     - NUMBER OF COLUMNS OF A. N MUST BE AT LEAST UNITY.
C
C     A     - AN M*N REAL MATRIX.
C
C     NRA   - ROW DIMENSION OF A AS DECLARED IN THE CALLING PROGRAM.
C             NRA MUST BE AT LEAST M.
C
C     B     - AN M ELEMENT REAL VECTOR.
C
C     TOL   - A RELATIVE TOLERANCE USED TO DETERMINE THE RANK OF A.
C             TOL SHOULD BE CHOSEN AS APPROXIMATELY THE
C             LARGEST RELATIVE ERROR IN THE ELEMENTS OF A.
C             FOR EXAMPLE IF THE ELEMENTS OF A ARE CORRECT
C             TO ABOUT 4 SIGNIFICANT FIGURES THEN TOL
C             SHOULD BE CHOSEN AS ABOUT 5.0*10.0**(-4).
C
C     IFAIL - THE USUAL FAILURE PARAMETER. IF IN DOUBT SET
C             IFAIL TO ZERO BEFORE CALLING THIS ROUTINE.
C
C     OUTPUT PARAMETERS.
C
C     A     - THE TOP N*N PART OF A WILL CONTAIN THE
C             ORTHOGONAL MATRIX P**T OF THE SVD.
C             THE REMAINDER OF A IS USED FOR INTERNAL WORKSPACE.
C
C     B     - THE FIRST N ELEMENTS OF B WILL CONTAIN THE
C             MINIMAL LEAST SQUARES SOLUTION VECTOR X.
C
C     SIGMA - IF M IS GREATER THAN IRANK THEN SIGMA WILL CONTAIN THE
C             STANDARD ERROR GIVEN BY
C             SIGMA=L(R)/SQRT(M-IRANK), WHERE L(R) DENOTES
C             THE EUCLIDEAN LENGTH OF THE RESIDUAL VECTOR
C             R. IF M=IRANK THEN SIGMA IS RETURNED AS ZERO.
C
C     IRANK - THE RANK OF THE MATRIX A.
C
C     IFAIL - ON NORMAL RETURN IFAIL WILL BE ZERO.
C             IN THE UNLIKELY EVENT THAT THE QR-ALGORITHM
C             FAILS TO FIND THE SINGULAR VALUES IN 50*N
C             ITERATIONS THEN IFAIL IS SET TO 2.
C             IF AN INPUT PARAMETER IS INCORRECTLY SUPPLIED
C             THEN IFAIL IS SET TO UNITY.
C
C     WORKSPACE PARAMETERS.
C
C     WORK  - A 4*N ELEMENT VECTOR.
C             ON RETURN THE FIRST N ELEMENTS OF WORK WILL
C             CONTAIN THE SINGULAR VALUES OF A ARRANGED IN
C             DESCENDING ORDER.
C             WORK(N+1) WILL CONTAIN THE TOTAL NUMBER OF ITERATIONS
C             TAKEN BY THE QR-ALGORITHM.
C
C     LWORK - THE LENGTH OF THE VECTOR WORK. LWORK MUST BE
C             AT LEAST 4*N.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04JAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SIGMA, TOL
      INTEGER           IFAIL, IRANK, LWORK, M, N, NRA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NRA,N), B(M), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           IERR, NNN, NP1, NP2
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           F02WDY, P01ABF
      EXTERNAL          F02WDY, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F02WAF, F04JAZ
C     .. Executable Statements ..
      IERR = IFAIL
      IF (IERR.EQ.0) IFAIL = 1
C
      IF (NRA.LT.M .OR. M.LT.N .OR. N.LT.1 .OR. LWORK.LT.4*N)
     *    GO TO 20
C
      NP1 = N + 1
      NP2 = NP1 + 1
      NNN = 3*N
C
      CALL F02WAF(M,N,A,NRA,.TRUE.,B,WORK,WORK(NP1),NNN,IFAIL)
C
      IF (IFAIL.NE.0) GO TO 20
C
      IRANK = F02WDY(N,WORK,TOL)
C
      CALL F04JAZ(M,N,IRANK,WORK,N,B,A,NRA,B,SIGMA,WORK(NP2))
C
      RETURN
C
   20 IFAIL = P01ABF(IERR,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F04JAY(N,IRANK,SV,LSV,B,PT,NRPT,X,WORK)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (SVDLSQ)
C
C     F04JAY RETURNS THE N ELEMENT VECTOR X GIVEN BY
C
C     X = P*(D**(-1))*B ,
C
C     WHERE D IS AN IRANK*IRANK NON-SINGULAR DIAGONAL MATRIX,
C     P CONTAINS THE FIRST IRANK COLUMNS OF AN N*N ORTHOGONAL
C     MATRIX AND B IS AN IRANK ELEMENT VECTOR.
C
C     THE ROUTINE MAY BE CALLED WITH IRANK=0 IN WHICH CASE X
C     IS RETURNED AS THE ZERO VECTOR.
C
C     INPUT PARAMETERS.
C
C     N     - NUMBER OF ROWS OF P. N MUST BE AT LEAST UNITY.
C
C     IRANK - ORDER OF THE MATRIX D.
C             IF IRANK=0 THEN SV, B, PT AND WORK ARE NOT REFERENCED.
C
C     SV    - AN IRANK ELEMENT VECTOR CONTAINING THE
C             DIAGONAL ELEMENTS OF D. SV MUST BE SUCH THAT
C             NO ELEMENT OF (D**(-1)*B WILL OVERFLOW.
C
C     LSV   - LSV MUST BE AT LEAST MAX(1,IRANK).
C
C     B     - AN IRANK ELEMENT VECTOR.
C
C     PT    - AN IRANK*N ELEMENT MATRIX CONTAINING THE MATRIX P**T.
C
C     NRPT  - ROW DIMENSION OF PT AS DECLARED IN THE
C             CALLING PROGRAM. NRPT MUST BE AT LEAST LSV.
C
C     OUTPUT PARAMETER.
C
C     X     - N ELEMENT VECTOR CONTAINING P*(D**(-1))*B.
C             IF IRANK=0 THEN X RETURNS THE ZERO VECTOR.
C             THE ROUTINE MAY BE CALLED WITH X=B OR WITH X=SV.
C
C     WORKSPACE PARAMETER.
C
C     WORK  - AN LSV ELEMENT VECTOR.
C             IF THE ROUTINE IS NOT CALLED WITH X=B THEN IT MAY BE
C             CALLED WITH WORK=B. SIMILARLY IF THE ROUTINE
C             IS NOT CALLED WITH X=SV THEN IT MAY BE CALLED
C             WITH WORK=SV.
C
C     Modified to call BLAS.
C     Jeremy Du Croz, NAG Central Office, October 1987.
C
C     .. Scalar Arguments ..
      INTEGER           IRANK, LSV, N, NRPT
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LSV), PT(NRPT,N), SV(LSV), WORK(LSV), X(N)
C     .. Local Scalars ..
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          DGEMV_
C     .. Executable Statements ..
      IF (IRANK.EQ.0) GO TO 40
C
      DO 20 I = 1, IRANK
         WORK(I) = B(I)/SV(I)
   20 CONTINUE
C
      CALL DGEMV_('Transpose',IRANK,N,1.0D0,PT,NRPT,WORK,1,0.0D0,X,1)
C
      RETURN
C
   40 DO 60 I = 1, N
         X(I) = 0.0D0
   60 CONTINUE
C
      RETURN
      END
      SUBROUTINE F04JAZ(M,N,IRANK,SV,LSV,B,PT,NRPT,X,SIGMA,WORK)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (SVDLS0)
C
C     F04JAZ RETURNS THE N ELEMENT VECTOR X, OF MINIMAL
C     LENGTH, THAT MINIMIZES THE EUCLIDEAN LENGTH OF THE M
C     ELEMENT VECTOR R GIVEN BY
C
C     R = B-A*X ,
C
C     WHERE B IS AN M ELEMENT VECTOR AND A IS AN M*N MATRIX,
C     FOLLOWING A SINGULAR VALUE DECOMPOSITION (SVD) OF A
C     GIVEN BY
C
C     A = Q*D*(P**T) ,
C
C     WHERE D IS A RECTANGULAR DIAGONAL MATRIX WHOSE DIAGONAL
C     ELEMENTS CONTAIN THE SINGULAR VALUES OF A IN DESCENDING
C     ORDER.
C
C     INPUT PARAMETERS.
C
C     M     - NUMBER OF ROWS OF A. M MUST BE AT LEAST UNITY.
C
C     N     - NUMBER OF COLUMNS OF A. N MUST BE AT LEAST UNITY.
C
C     IRANK - THE RANK OF THE MATRIX A. IRANK MUST BE SUCH THAT THE
C             ELEMENTS SV(I), I=1,2,...,IRANK ARE NON-NEGLIGIBLE.
C             IRANK MUST BE AT LEAST ZERO AND MUST NOT BE
C             LARGER THAN MIN(M,N).
C             ROUTINE F02WDY CAN BE USED TO DETERMINE RANK FOLLOWING
C             AN SVD.
C
C     SV    - AN LSV ELEMENT VECTOR CONTAINING THE POSITIVE
C             NON-NEGLIGIBLE SINGULAR VALUES OF A.
C
C     LSV   - LENGTH OF THE VECTOR SV.
C             LSV MUST BE AT LEAST MAX(1,IRANK).
C
C     B     - MUST CONTAIN THE M ELEMENT VECTOR (Q**T)*B, WHERE Q IS
C             THE LEFT-HAND ORTHOGONAL MATRIX OF THE SVD.
C
C     PT    - THE IRANK*N PART OF PT MUST CONTAIN THE FIRST
C             IRANK ROWS OF THE RIGHT-HAND ORTHOGONAL
C             MATRIX P**T OF THE SVD.
C
C     NRPT  - ROW DIMENSION OF PT AS DECLARED IN THE CALLING PROGRAM
C             NRPT MUST BE AT LEAST LSV.
C
C     OUTPUT PARAMETERS.
C
C     X     - THE N ELEMENT SOLUTION VECTOR.
C             THE ROUTINE MAY BE CALLED WITH X=B OR WITH X=SV.
C
C     SIGMA - IF M IS GREATER THAN IRANK THEN SIGMA WILL CONTAIN THE
C             STANDARD ERROR GIVEN BY
C             SIGMA=L(R)/SQRT(M-IRANK), WHERE L(R) DENOTES
C             THE EUCLIDEAN LENGTH OF THE RESIDUAL VECTOR
C             R. IF M=IRANK THEN SIGMA IS RETURNED AS ZERO.
C
C     WORKSPACE PARAMETER.
C
C     WORK  - AN LSV ELEMENT VECTOR.
C             IF THE ROUTINE IS NOT CALLED WITH X=B THEN IT MAY BE
C             CALLED WITH WORK=B. SIMILARLY IF THE ROUTINE
C             IS NOT CALLED WITH X=SV THEN IT MAY BE CALLED
C             WITH WORK=SV.
C
C     Modified to call BLAS.
C     Jeremy Du Croz, NAG Central Office, October 1987.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SIGMA
      INTEGER           IRANK, LSV, M, N, NRPT
C     .. Array Arguments ..
      DOUBLE PRECISION  B(M), PT(NRPT,N), SV(LSV), WORK(LSV), X(N)
C     .. Local Scalars ..
      INTEGER           IRP1, MMIR
C     .. External Functions ..
      DOUBLE PRECISION  F06EJF
      EXTERNAL          F06EJF
C     .. External Subroutines ..
      EXTERNAL          F04JAY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      SIGMA = 0.0D0
      IF (IRANK.EQ.M) GO TO 20
      IRP1 = IRANK + 1
      MMIR = M - IRANK
C
      SIGMA = F06EJF(MMIR,B(IRP1),1)/SQRT(DBLE(MMIR))
C
   20 CALL F04JAY(N,IRANK,SV,LSV,B,PT,NRPT,X,WORK)
C
      RETURN
      END
      SUBROUTINE F06AAZ ( SRNAME, INFO )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 15 REVISED. IER-915 (APR 1991).
C     .. Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*13       SRNAME
C     ..
C
C  Purpose
C  =======
C
C  F06AAZ  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*13.
C           On entry, SRNAME specifies the name of the routine which
C           called F06AAZ.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Local Scalars ..
      INTEGER            IERR, IFAIL
      CHARACTER*4        VARBNM
C     .. Local Arrays ..
      CHARACTER*80       REC (1)
C     .. External Functions ..
      INTEGER            P01ACF
      EXTERNAL           P01ACF
C     ..
C     .. Executable Statements ..
      WRITE (REC (1),99999) SRNAME, INFO
      IF (SRNAME(1:3).EQ.'F06') THEN
         IERR = -1
         VARBNM = '    '
      ELSE
         IERR = -INFO
         VARBNM = 'INFO'
      END IF
      IFAIL = 0
      IFAIL = P01ACF (IFAIL, IERR, SRNAME(1:6), VARBNM, 1, REC)
C
      RETURN
C
99999 FORMAT ( ' ** On entry to ', A13, ' parameter number ', I2,
     $         ' had an illegal value' )
C
C     End of F06AAZ.
C
      END
      DOUBLE PRECISION FUNCTION F06BMF( SCALE, SSQ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  SCALE, SSQ
C     ..
C
C  F06BMF returns the value norm given by
C
C     norm = ( scale*sqrt( ssq ), scale*sqrt( ssq ) .lt. flmax
C            (
C            ( flmax,             scale*sqrt( ssq ) .ge. flmax
C
C  via the function name.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION      FLMAX, FLMIN, NORM, SQT
      LOGICAL               FIRST
C     .. External Functions ..
      DOUBLE PRECISION      X02AMF
      EXTERNAL              X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC             SQRT
C     .. Save statement ..
      SAVE                  FIRST, FLMAX
C     .. Data statements ..
      DATA                  FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( FIRST )THEN
         FIRST = .FALSE.
         FLMIN =  X02AMF( )
         FLMAX =  1/FLMIN
      END IF
C
      SQT = SQRT( SSQ )
      IF( SCALE.LT.FLMAX/SQT )THEN
         NORM = SCALE*SQT
      ELSE
         NORM = FLMAX
      END IF
C
      F06BMF = NORM
      RETURN
C
C     End of F06BMF. ( SNORM )
C
      END
      SUBROUTINE F06EDF( N, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DSCAL_ ( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06EDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine DSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.EQ.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   10       CONTINUE
         ELSE IF( ALPHA.EQ.( -ONE ) )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = -X( IX )
   20       CONTINUE
         ELSE IF( ALPHA.NE.ONE )THEN
            DO 30, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   30       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06EDF. ( DSCAL )
C
      END
      SUBROUTINE F06EGF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DSWAP_ ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06EGF performs the operations
C
C     temp := x,   x := y,   y := temp.
C
C
C  Nag Fortran 77 version of the Blas routine DSWAP.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               TEMP    = X( IY )
               X( IY ) = Y( IY )
               Y( IY ) = TEMP
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06EGF. ( DSWAP )
C
      END
      DOUBLE PRECISION FUNCTION F06EJF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DNRM2_
      ENTRY                     DNRM2_ ( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                           INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  X( * )
C     ..
C
C  F06EJF returns the euclidean norm of a vector via the function
C  name, so that
C
C     F06EJF := sqrt( x'*x )
C
C
C  Nag Fortran 77 version of the Blas routine DNRM2.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 25-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      NORM, SCALE, SSQ
C     .. External Functions ..
      DOUBLE PRECISION      F06BMF
      EXTERNAL              F06BMF
C     .. External Subroutines ..
      EXTERNAL              F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC             ABS
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
         CALL F06FJF( N, X, INCX, SCALE, SSQ )
         NORM  = F06BMF( SCALE, SSQ )
      END IF
C
      F06EJF = NORM
      RETURN
C
C     End of F06EJF. ( DNRM2 )
C
      END
      SUBROUTINE F06FJF( N, X, INCX, SCALE, SUMSQ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   SCALE, SUMSQ
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FJF returns the values scl and smsq such that
C
C     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
C
C  where x( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is assumed
C  to be at least unity and the value of smsq will then satisfy
C
C     1.0 .le. smsq .le. ( sumsq + n ) .
C
C  scale is assumed to be non-negative and scl returns the value
C
C     scl = max( scale, abs( x( i ) ) ) .
C
C  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
C  scl and smsq are overwritten on SCALE and SUMSQ respectively.
C
C  The routine makes only one pass through the vector X.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   ABSXI
      INTEGER            IX
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SUMSQ = 1     + SUMSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ +       ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
C
C     End of F06FJF. ( SSSQ )
C
      END
      SUBROUTINE F06FRF( N, ALPHA, X, INCX, TOL, ZETA )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, TOL, ZETA
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FRF generates details of a generalized Householder reflection such
C  that
C
C     P*( alpha ) = ( beta ),   P'*P = I.
C       (   x   )   (   0  )
C
C  P is given in the form
C
C     P = I - ( zeta )*( zeta  z' ),
C             (   z  )
C
C  where z is an n element vector and zeta is a scalar that satisfies
C
C     1.0 .le. zeta .le. sqrt( 2.0 ).
C
C  zeta is returned in ZETA unless x is such that
C
C     max( abs( x( i ) ) ) .le. max( eps*abs( alpha ), tol )
C
C  where eps is the relative machine precision and tol is the user
C  supplied value TOL, in which case ZETA is returned as 0.0 and P can
C  be taken to be the unit matrix.
C
C  beta is overwritten on alpha and z is overwritten on x.
C  the routine may be called with  n = 0  and advantage is taken of the
C  case where  n = 1.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 30-August-1984.
C     Sven Hammarling, Nag Central Office.
C     This version dated 28-September-1984.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   BETA, EPS, SCALE, SSQ
      LOGICAL            FIRST
C     .. External Functions ..
      DOUBLE PRECISION   X02AJF
      EXTERNAL           X02AJF
C     .. External Subroutines ..
      EXTERNAL           F06FJF, DSCAL_
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
C     .. Save statement ..
      SAVE               EPS, FIRST
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         ZETA = ZERO
      ELSE IF( ( N.EQ.1 ).AND.( X( 1 ).EQ.ZERO ) )THEN
         ZETA = ZERO
      ELSE
C
         IF( FIRST )THEN
            FIRST = .FALSE.
            EPS   =  X02AJF( )
         END IF
C
C        Treat case where P is a 2 by 2 matrix specially.
C
         IF( N.EQ.1 )THEN
C
C           Deal with cases where  ALPHA = zero  and
C           abs( X( 1 ) ) .le. max( EPS*abs( ALPHA ), TOL )  first.
C
            IF( ALPHA.EQ.ZERO )THEN
               ZETA   =  ONE
               ALPHA  =  ABS ( X( 1 ) )
               X( 1 ) = -SIGN( ONE, X( 1 ) )
            ELSE IF( ABS( X( 1 ) ).LE.MAX( EPS*ABS( ALPHA ), TOL ) )THEN
               ZETA   =  ZERO
            ELSE
               IF( ABS( ALPHA ).GE.ABS( X( 1 ) ) )THEN
                  BETA = ABS( ALPHA ) *SQRT( 1 + ( X( 1 )/ALPHA )**2 )
               ELSE
                  BETA = ABS( X( 1 ) )*SQRT( 1 + ( ALPHA/X( 1 ) )**2 )
               END IF
               ZETA = SQRT( ( ABS( ALPHA ) + BETA )/BETA )
               IF( ALPHA.GE.ZERO )
     $            BETA = -BETA
               X( 1 ) = -X( 1 )/( ZETA*BETA )
               ALPHA  = BETA
            END IF
         ELSE
C
C           Now P is larger than 2 by 2.
C
            SSQ   = ONE
            SCALE = ZERO
            CALL F06FJF( N, X, INCX, SCALE, SSQ )
C
C           Treat cases where  SCALE = zero,
C           SCALE .le. max( EPS*abs( ALPHA ), TOL )  and
C           ALPHA = zero  specially.
C           Note that  SCALE = max( abs( X( i ) ) ).
C
            IF( ( SCALE.EQ.ZERO ).OR.
     $          ( SCALE.LE.MAX( EPS*ABS( ALPHA ), TOL ) ) )THEN
               ZETA  = ZERO
            ELSE IF( ALPHA.EQ.ZERO )THEN
               ZETA  = ONE
               ALPHA = SCALE*SQRT( SSQ )
               CALL DSCAL_( N, -1/ALPHA, X, INCX )
            ELSE
               IF( SCALE.LT.ABS( ALPHA ) )THEN
                  BETA = ABS( ALPHA )*SQRT( 1 + SSQ*( SCALE/ALPHA )**2 )
               ELSE
                  BETA = SCALE       *SQRT( SSQ +   ( ALPHA/SCALE )**2 )
               END IF
               ZETA = SQRT( ( BETA + ABS( ALPHA ) )/BETA )
               IF( ALPHA.GT.ZERO )
     $            BETA = -BETA
               CALL DSCAL_( N, -1/( ZETA*BETA ), X, INCX )
               ALPHA = BETA
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FRF. ( SGRFG )
C
      END
      SUBROUTINE F06PAF( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DGEMV_ ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGEMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry with BETA non-zero, the incremented array Y
C           must contain the vector y. On exit, Y is overwritten by the
C           updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PAF/DGEMV_ ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  y := alpha*A*x + y.
C
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y := alpha*A'*x + y.
C
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PAF (DGEMV ).
C
      END
      SUBROUTINE F06PJF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  DTRSV  solves one of the systems of equations
C
C     A*x = b,   or   A'*x = b,
C
C  where b and x are n element vectors and A is an n by n unit, or
C  non-unit, upper or lower triangular matrix.
C
C  No test for singularity or near-singularity is included in this
C  routine. Such tests must be performed before calling this routine.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the equations to be solved as
C           follows:
C
C              TRANS = 'N' or 'n'   A*x = b.
C
C              TRANS = 'T' or 't'   A'*x = b.
C
C              TRANS = 'C' or 'c'   A'*x = b.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element right-hand side vector b. On exit, X is overwritten
C           with the solution vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO .EQ.'U' .OR. UPLO .EQ.'u').AND.
     $         .NOT.(UPLO .EQ.'L' .OR. UPLO .EQ.'l')      )THEN
         INFO = 1
      ELSE IF( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 2
      ELSE IF( .NOT.(DIAG .EQ.'U' .OR. DIAG .EQ.'u').AND.
     $         .NOT.(DIAG .EQ.'N' .OR. DIAG .EQ.'n')      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PJF/DTRSV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  x := inv( A )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, 1, -1
                        X( I ) = X( I ) - TEMP*A( I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 40, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 30, I = J - 1, 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, N
                        X( I ) = X( I ) - TEMP*A( I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 70, I = J + 1, N
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := inv( A' )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  DO 90, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  DO 110, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  DO 130, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  DO 150, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06PJF (DTRSV ).
C
      END
      SUBROUTINE F06PMF( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DGER_  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGER   performs the rank 1 operation
C
C     A := alpha*x*y' + A,
C
C  where alpha is a scalar, x is an m element vector, y is an n element
C  vector and A is an m by n matrix.
C
C  Parameters
C  ==========
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( m - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the m
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients. On exit, A is
C           overwritten by the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PMF/DGER_  ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06PMF (DGER  ).
C
      END
      SUBROUTINE F06YAF(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,
     *                  LDC)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  Purpose
C  =======
C
C  DGEMM  performs one of the matrix-matrix operations
C
C     C := alpha*op( A )*op( B ) + beta*C,
C
C  where  op( X ) is one of
C
C     op( X ) = X   or   op( X ) = X',
C
C  alpha and beta are scalars, and A, B and C are matrices, with op( A )
C  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSA = 'N' or 'n',  op( A ) = A.
C
C              TRANSA = 'T' or 't',  op( A ) = A'.
C
C              TRANSA = 'C' or 'c',  op( A ) = A'.
C
C           Unchanged on exit.
C
C  TRANSB - CHARACTER*1.
C           On entry, TRANSB specifies the form of op( B ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSB = 'N' or 'n',  op( B ) = B.
C
C              TRANSB = 'T' or 't',  op( B ) = B'.
C
C              TRANSB = 'C' or 'c',  op( B ) = B'.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry,  M  specifies  the number  of rows  of the  matrix
C           op( A )  and of the  matrix  C.  M  must  be at least  zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry,  N  specifies the number  of columns of the matrix
C           op( B ) and the number of columns of the matrix C. N must be
C           at least zero.
C           Unchanged on exit.
C
C  K      - INTEGER.
C           On entry,  K  specifies  the number of columns of the matrix
C           op( A ) and the number of rows of the matrix op( B ). K must
C           be at least  zero.
C           Unchanged on exit.
C
C  ALPHA  - REAL            .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
C           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
C           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
C           part of the array  A  must contain the matrix  A,  otherwise
C           the leading  k by m  part of the array  A  must contain  the
C           matrix A.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
C           LDA must be at least  max( 1, m ), otherwise  LDA must be at
C           least  max( 1, k ).
C           Unchanged on exit.
C
C  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
C           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
C           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
C           part of the array  B  must contain the matrix  B,  otherwise
C           the leading  n by k  part of the array  B  must contain  the
C           matrix B.
C           Unchanged on exit.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
C           LDB must be at least  max( 1, k ), otherwise  LDB must be at
C           least  max( 1, n ).
C           Unchanged on exit.
C
C  BETA   - REAL            .
C           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C           supplied as zero then C need not be set on input.
C           Unchanged on exit.
C
C  C      - REAL             array of DIMENSION ( LDC, n ).
C           Before entry, the leading  m by n  part of the array  C must
C           contain the matrix  C,  except when  beta  is zero, in which
C           case C need not be set on entry.
C           On exit, the array  C  is overwritten by the  m by n  matrix
C           ( alpha*op( A )*op( B ) + beta*C ).
C
C  LDC    - INTEGER.
C           On entry, LDC specifies the first dimension of C as declared
C           in  the  calling  (sub)  program.   LDC  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 3 Blas routine.
C
C  -- Written on 8-February-1989.
C     Jack Dongarra, Argonne National Laboratory.
C     Iain Duff, AERE Harwell.
C     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C     Sven Hammarling, Numerical Algorithms Group Ltd.
C
C
C     .. Entry Points ..
      ENTRY             DGEMM_(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     *                  BETA,C,LDC)
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, BETA
      INTEGER           K, LDA, LDB, LDC, M, N
      CHARACTER*1       TRANSA, TRANSB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, INFO, J, L, NCOLA, NROWA, NROWB
      LOGICAL           NOTA, NOTB
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
C     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
C     and  columns of  A  and the  number of  rows  of  B  respectively.
C
      NOTA = (TRANSA.EQ.'N' .OR. TRANSA.EQ.'n')
      NOTB = (TRANSB.EQ.'N' .OR. TRANSB.EQ.'n')
      IF (NOTA) THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF (NOTB) THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
C
C     Test the input parameters.
C
      INFO = 0
      IF (( .NOT. NOTA) .AND. ( .NOT. (TRANSA.EQ.'C' .OR. TRANSA.EQ.'c')
     *    ) .AND. ( .NOT. (TRANSA.EQ.'T' .OR. TRANSA.EQ.'t'))) THEN
         INFO = 1
      ELSE IF (( .NOT. NOTB) .AND. ( .NOT. (TRANSB.EQ.'C' .OR.
     *         TRANSB.EQ.'c')) .AND. ( .NOT. (TRANSB.EQ.'T' .OR.
     *         TRANSB.EQ.'t'))) THEN
         INFO = 2
      ELSE IF (M.LT.0) THEN
         INFO = 3
      ELSE IF (N.LT.0) THEN
         INFO = 4
      ELSE IF (K.LT.0) THEN
         INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
         INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = 13
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F06YAF/DGEMM_ ',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (((ALPHA.EQ.ZERO) .OR. (K.EQ.0))
     *     .AND. (BETA.EQ.ONE))) RETURN
C
C     And if  alpha.eq.zero.
C
      IF (ALPHA.EQ.ZERO) THEN
         IF (BETA.EQ.ZERO) THEN
            DO 40 J = 1, N
               DO 20 I = 1, M
                  C(I,J) = ZERO
   20          CONTINUE
   40       CONTINUE
         ELSE
            DO 80 J = 1, N
               DO 60 I = 1, M
                  C(I,J) = BETA*C(I,J)
   60          CONTINUE
   80       CONTINUE
         END IF
         RETURN
      END IF
C
C     Start the operations.
C
      IF (NOTB) THEN
         IF (NOTA) THEN
C
C           Form  C := alpha*A*B + beta*C.
C
            DO 180 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 100 I = 1, M
                     C(I,J) = ZERO
  100             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 120 I = 1, M
                     C(I,J) = BETA*C(I,J)
  120             CONTINUE
               END IF
               DO 160 L = 1, K
                  IF (B(L,J).NE.ZERO) THEN
                     TEMP = ALPHA*B(L,J)
                     DO 140 I = 1, M
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  140                CONTINUE
                  END IF
  160          CONTINUE
  180       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B + beta*C
C
            DO 240 J = 1, N
               DO 220 I = 1, M
                  TEMP = ZERO
                  DO 200 L = 1, K
                     TEMP = TEMP + A(L,I)*B(L,J)
  200             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  220          CONTINUE
  240       CONTINUE
         END IF
      ELSE
         IF (NOTA) THEN
C
C           Form  C := alpha*A*B' + beta*C
C
            DO 340 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 260 I = 1, M
                     C(I,J) = ZERO
  260             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 280 I = 1, M
                     C(I,J) = BETA*C(I,J)
  280             CONTINUE
               END IF
               DO 320 L = 1, K
                  IF (B(J,L).NE.ZERO) THEN
                     TEMP = ALPHA*B(J,L)
                     DO 300 I = 1, M
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  300                CONTINUE
                  END IF
  320          CONTINUE
  340       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B' + beta*C
C
            DO 400 J = 1, N
               DO 380 I = 1, M
                  TEMP = ZERO
                  DO 360 L = 1, K
                     TEMP = TEMP + A(L,I)*B(J,L)
  360             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  380          CONTINUE
  400       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06YAF (DGEMM ).
C
      END
      SUBROUTINE F06YJF(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  Purpose
C  =======
C
C  DTRSM  solves one of the matrix equations
C
C     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
C
C  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
C  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
C
C     op( A ) = A   or   op( A ) = A'.
C
C  The matrix X is overwritten on B.
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C           On entry, SIDE specifies whether op( A ) appears on the left
C           or right of X as follows:
C
C              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
C
C              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
C
C           Unchanged on exit.
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix A is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSA = 'N' or 'n'   op( A ) = A.
C
C              TRANSA = 'T' or 't'   op( A ) = A'.
C
C              TRANSA = 'C' or 'c'   op( A ) = A'.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit triangular
C           as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of B. M must be at
C           least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of B.  N must be
C           at least zero.
C           Unchanged on exit.
C
C  ALPHA  - REAL            .
C           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
C           zero then  A is not referenced and  B need not be set before
C           entry.
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
C           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
C           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
C           upper triangular part of the array  A must contain the upper
C           triangular matrix  and the strictly lower triangular part of
C           A is not referenced.
C           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
C           lower triangular part of the array  A must contain the lower
C           triangular matrix  and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
C           A  are not referenced either,  but are assumed to be  unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
C           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
C           then LDA must be at least max( 1, n ).
C           Unchanged on exit.
C
C  B      - REAL             array of DIMENSION ( LDB, n ).
C           Before entry,  the leading  m by n part of the array  B must
C           contain  the  right-hand  side  matrix  B,  and  on exit  is
C           overwritten by the solution matrix  X.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in  the  calling  (sub)  program.   LDB  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 3 Blas routine.
C
C
C  -- Written on 8-February-1989.
C     Jack Dongarra, Argonne National Laboratory.
C     Iain Duff, AERE Harwell.
C     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C     Sven Hammarling, Numerical Algorithms Group Ltd.
C
C
C     .. Entry Points ..
      ENTRY             DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,
     *                  LDB)
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA
      INTEGER           LDA, LDB, M, N
      CHARACTER*1       DIAG, SIDE, TRANSA, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, INFO, J, K, NROWA
      LOGICAL           LSIDE, NOUNIT, UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      LSIDE = (SIDE.EQ.'L' .OR. SIDE.EQ.'l')
      IF (LSIDE) THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
C
      INFO = 0
      IF (( .NOT. LSIDE) .AND. ( .NOT. (SIDE.EQ.'R' .OR. SIDE.EQ.'r')))
     *    THEN
         INFO = 1
      ELSE IF (( .NOT. UPPER) .AND. ( .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.
     *         'l'))) THEN
         INFO = 2
      ELSE IF (( .NOT. (TRANSA.EQ.'N' .OR. TRANSA.EQ.'n'))
     *         .AND. ( .NOT. (TRANSA.EQ.'T' .OR. TRANSA.EQ.'t'))
     *         .AND. ( .NOT. (TRANSA.EQ.'C' .OR. TRANSA.EQ.'c'))) THEN
         INFO = 3
      ELSE IF (( .NOT. (DIAG.EQ.'U' .OR. DIAG.EQ.'u'))
     *         .AND. ( .NOT. (DIAG.EQ.'N' .OR. DIAG.EQ.'n'))) THEN
         INFO = 4
      ELSE IF (M.LT.0) THEN
         INFO = 5
      ELSE IF (N.LT.0) THEN
         INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
         INFO = 11
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F06YJF/DTRSM ',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF (N.EQ.0) RETURN
C
C     And when  alpha.eq.zero.
C
      IF (ALPHA.EQ.ZERO) THEN
         DO 40 J = 1, N
            DO 20 I = 1, M
               B(I,J) = ZERO
   20       CONTINUE
   40    CONTINUE
         RETURN
      END IF
C
C     Start the operations.
C
      IF (LSIDE) THEN
         IF ((TRANSA.EQ.'N' .OR. TRANSA.EQ.'n')) THEN
C
C           Form  B := alpha*inv( A )*B.
C
            IF (UPPER) THEN
               DO 120 J = 1, N
                  IF (ALPHA.NE.ONE) THEN
                     DO 60 I = 1, M
                        B(I,J) = ALPHA*B(I,J)
   60                CONTINUE
                  END IF
                  DO 100 K = M, 1, -1
                     IF (B(K,J).NE.ZERO) THEN
                        IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                        DO 80 I = 1, K - 1
                           B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                   CONTINUE
                     END IF
  100             CONTINUE
  120          CONTINUE
            ELSE
               DO 200 J = 1, N
                  IF (ALPHA.NE.ONE) THEN
                     DO 140 I = 1, M
                        B(I,J) = ALPHA*B(I,J)
  140                CONTINUE
                  END IF
                  DO 180 K = 1, M
                     IF (B(K,J).NE.ZERO) THEN
                        IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                        DO 160 I = K + 1, M
                           B(I,J) = B(I,J) - B(K,J)*A(I,K)
  160                   CONTINUE
                     END IF
  180             CONTINUE
  200          CONTINUE
            END IF
         ELSE
C
C           Form  B := alpha*inv( A' )*B.
C
            IF (UPPER) THEN
               DO 260 J = 1, N
                  DO 240 I = 1, M
                     TEMP = ALPHA*B(I,J)
                     DO 220 K = 1, I - 1
                        TEMP = TEMP - A(K,I)*B(K,J)
  220                CONTINUE
                     IF (NOUNIT) TEMP = TEMP/A(I,I)
                     B(I,J) = TEMP
  240             CONTINUE
  260          CONTINUE
            ELSE
               DO 320 J = 1, N
                  DO 300 I = M, 1, -1
                     TEMP = ALPHA*B(I,J)
                     DO 280 K = I + 1, M
                        TEMP = TEMP - A(K,I)*B(K,J)
  280                CONTINUE
                     IF (NOUNIT) TEMP = TEMP/A(I,I)
                     B(I,J) = TEMP
  300             CONTINUE
  320          CONTINUE
            END IF
         END IF
      ELSE
         IF ((TRANSA.EQ.'N' .OR. TRANSA.EQ.'n')) THEN
C
C           Form  B := alpha*B*inv( A ).
C
            IF (UPPER) THEN
               DO 420 J = 1, N
                  IF (ALPHA.NE.ONE) THEN
                     DO 340 I = 1, M
                        B(I,J) = ALPHA*B(I,J)
  340                CONTINUE
                  END IF
                  DO 380 K = 1, J - 1
                     IF (A(K,J).NE.ZERO) THEN
                        DO 360 I = 1, M
                           B(I,J) = B(I,J) - A(K,J)*B(I,K)
  360                   CONTINUE
                     END IF
  380             CONTINUE
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(J,J)
                     DO 400 I = 1, M
                        B(I,J) = TEMP*B(I,J)
  400                CONTINUE
                  END IF
  420          CONTINUE
            ELSE
               DO 520 J = N, 1, -1
                  IF (ALPHA.NE.ONE) THEN
                     DO 440 I = 1, M
                        B(I,J) = ALPHA*B(I,J)
  440                CONTINUE
                  END IF
                  DO 480 K = J + 1, N
                     IF (A(K,J).NE.ZERO) THEN
                        DO 460 I = 1, M
                           B(I,J) = B(I,J) - A(K,J)*B(I,K)
  460                   CONTINUE
                     END IF
  480             CONTINUE
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(J,J)
                     DO 500 I = 1, M
                        B(I,J) = TEMP*B(I,J)
  500                CONTINUE
                  END IF
  520          CONTINUE
            END IF
         ELSE
C
C           Form  B := alpha*B*inv( A' ).
C
            IF (UPPER) THEN
               DO 620 K = N, 1, -1
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(K,K)
                     DO 540 I = 1, M
                        B(I,K) = TEMP*B(I,K)
  540                CONTINUE
                  END IF
                  DO 580 J = 1, K - 1
                     IF (A(J,K).NE.ZERO) THEN
                        TEMP = A(J,K)
                        DO 560 I = 1, M
                           B(I,J) = B(I,J) - TEMP*B(I,K)
  560                   CONTINUE
                     END IF
  580             CONTINUE
                  IF (ALPHA.NE.ONE) THEN
                     DO 600 I = 1, M
                        B(I,K) = ALPHA*B(I,K)
  600                CONTINUE
                  END IF
  620          CONTINUE
            ELSE
               DO 720 K = 1, N
                  IF (NOUNIT) THEN
                     TEMP = ONE/A(K,K)
                     DO 640 I = 1, M
                        B(I,K) = TEMP*B(I,K)
  640                CONTINUE
                  END IF
                  DO 680 J = K + 1, N
                     IF (A(J,K).NE.ZERO) THEN
                        TEMP = A(J,K)
                        DO 660 I = 1, M
                           B(I,J) = B(I,J) - TEMP*B(I,K)
  660                   CONTINUE
                     END IF
  680             CONTINUE
                  IF (ALPHA.NE.ONE) THEN
                     DO 700 I = 1, M
                        B(I,K) = ALPHA*B(I,K)
  700                CONTINUE
                  END IF
  720          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06YJF (DTRSM ).
C
      END
      SUBROUTINE F07ADG(M,N,A,LDA,PIV,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C  Purpose
C  =======
C
C  F07ADG computes an LU factorization of a general m-by-n matrix A
C  using partial pivoting with row interchanges.
C
C  The factorization has the form
C     A = P * L * U
C  where P is a permutation matrix, L is lower triangular (lower
C  trapezoidal if m > n), and U is upper triangular with unit diagonal
C  elements (upper trapezoidal if m < n).
C
C  This is the Level 3 BLAS version of the Crout algorithm.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix A.  N >= 0.
C
C  A       (input/output) REAL array, dimension (LDA,N)
C          On entry, the m by n matrix to be factored.
C          On exit, the factors L and U from the factorization
C          A = P*L*U; the unit diagonal elements of U are not stored.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,M).
C
C  PIV     (output) REAL array, dimension (M)
C          The pivot indices; for 1 <= i <= min(M,N), row i of the
C          matrix was interchanged with row PIV(i). The rest of PIV is
C          used for workspace.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, approximate singularity has been detected
C               at the k-th stage; the factorization has not been
C               completed.
C
C  This is a modified version of the LAPACK routine F07ADF/DGETRF, in
C  which the INTEGER array IPIV has been replaced by a REAL array PIV,
C  row-equilibration is used in the choice of pivot, U has unit diagonal
C  elements, and the routine exits immediately if approximate
C  singularity is detected.
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), PIV(*)
C     .. Local Scalars ..
      INTEGER           I, IINFO, J, JB, NB
C     .. External Subroutines ..
      EXTERNAL          DGEMM_, DTRSM, F06AAZ, F07ADH, F07ADJ, F07ZAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, SQRT
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -4
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07ADG       ',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
C
C     Compute 2-norm of each row and store reciprocal in PIV.
C
      DO 20 I = 1, M
         PIV(I) = ZERO
   20 CONTINUE
      DO 60 J = 1, N
         DO 40 I = 1, M
            PIV(I) = PIV(I) + A(I,J)**2
   40    CONTINUE
   60 CONTINUE
      DO 80 I = 1, M
         IF (PIV(I).LE.ZERO) THEN
            INFO = I
            RETURN
         ELSE
            PIV(I) = ONE/SQRT(PIV(I))
         END IF
   80 CONTINUE
C
C     Determine the block size for this environment.
C
      CALL F07ZAZ(1,'F07ADG',NB,0)
      IF (NB.LE.1) THEN
C
C        Use unblocked code.
C
         CALL F07ADH(M,N,A,LDA,PIV,INFO)
      ELSE
C
C        Use blocked code.
C
         DO 120 J = 1, MIN(M,N), NB
            JB = MIN(MIN(M,N)-J+1,NB)
C
C           Update diagonal and subdiagonal blocks.
C
            CALL DGEMM_('No transpose','No transpose',M-J+1,JB,J-1,-ONE,
     *                 A(J,1),LDA,A(1,J),LDA,ONE,A(J,J),LDA)
C
C           Factorize diagonal and subdiagonal blocks and test for
C           approximate singularity.
C
            CALL F07ADH(M-J+1,JB,A(J,J),LDA,PIV(J),IINFO)
C
            IF (IINFO.GT.0) THEN
               INFO = IINFO + J - 1
               RETURN
            END IF
C
C           Update pivot indices and apply the interchanges to columns
C           1:J-1.
C
            DO 100 I = J, MIN(M,J+JB-1)
               PIV(I) = J - 1 + PIV(I)
  100       CONTINUE
            CALL F07ADJ(J-1,A,LDA,J,J+JB-1,PIV,1)
C
            IF (J+JB.LE.N) THEN
C
C              Apply the interchanges to columns J+JB:N.
C
               CALL F07ADJ(N-J-JB+1,A(1,J+JB),LDA,J,J+JB-1,PIV,1)
C
C              Compute block row of U.
C
            CALL DGEMM_('No transpose','No transpose',JB,N-J-JB+1,J-1,
     *                    -ONE,A(J,1),LDA,A(1,J+JB),LDA,ONE,A(J,J+JB),
     *                    LDA)
            CALL DTRSM('Left','Lower','No transpose','Non-unit',JB,
     *                    N-J-JB+1,ONE,A(J,J),LDA,A(J,J+JB),LDA)
            END IF
  120    CONTINUE
      END IF
      RETURN
C
C     End of F07ADG
C
      END
      SUBROUTINE F07ADH(M,N,A,LDA,PIV,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C  Purpose
C  =======
C
C  F07ADH computes an LU factorization of a general m-by-n matrix A
C  using partial pivoting with row interchanges.
C
C  The factorization has the form
C     A = P * L * U
C  where P is a permutation matrix, L is lower triangular (lower
C  trapezoidal if m > n), and U is upper triangular with unit diagonal
C  elements (upper trapezoidal if m < n).
C
C  This is the Level 2 BLAS version of the Crout algorithm.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix A.  N >= 0.
C
C  A       (input/output) REAL array, dimension (LDA,N)
C          On entry, the m by n matrix to be factored.
C          On exit, the factors L and U from the factorization
C          A = P*L*U; the unit diagonal elements of U are not stored.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,M).
C
C  PIV     (input/output) REAL array, dimension (M)
C          On entry, M scale factors for equilibrating the rows of A.
C          On exit, the pivot indices; for 1 <= i <= min(M,N), row i of
C          the matrix was interchanged with row PIV(i).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, approximate singularity has been detected a
C               the k-th stage; the factorization has not been
C               completed.
C
C  This is a modified version of the LAPACK routine F07ADZ/DGETF2, in
C  which the INTEGER array IPIV has been replaced by a REAL array PIV,
C  row-equilibration is used in the choice of pivot, U has unit diagonal
C  elements, and the routine exits immediately if approximate
C  singularity is detected.
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO, EIGHT
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0,EIGHT=8.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), PIV(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  THRESH, X, Y
      INTEGER           I, J, JP
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          DGEMV_, DSCAL_, DSWAP_, F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -4
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07ADH       ',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
C
C     Set threshold for test for singularity
C
      THRESH = EIGHT*X02AJF()
C
      DO 40 J = 1, MIN(M,N)
C
C        Update diagonal and subdiagonal elements in column J.
C
         CALL DGEMV_('No transpose',M-J+1,J-1,-ONE,A(J,1),LDA,A(1,J),1,
     *              ONE,A(J,J),1)
C
C        Find pivot and test for singularity.
C
         JP = J
         X = ZERO
         DO 20 I = J, M
            Y = ABS(A(I,J))*PIV(I)
            IF (Y.GT.X) THEN
               JP = I
               X = Y
            END IF
   20    CONTINUE
         PIV(JP) = PIV(J)
         PIV(J) = JP
         IF (X.GE.THRESH) THEN
C
C           Apply interchange to columns 1:N.
C
            IF (JP.NE.J) CALL DSWAP_(N,A(J,1),LDA,A(JP,1),LDA)
C
            IF (J+1.LE.N) THEN
C
C              Compute row of U.
C
               CALL DGEMV_('Transpose',J-1,N-J,-ONE,A(1,J+1),LDA,A(J,1),
     *                    LDA,ONE,A(J,J+1),LDA)
C
               CALL DSCAL_(N-J,ONE/A(J,J),A(J,J+1),LDA)
            END IF
C
         ELSE
C
C           If A( JP, J ) is small, set INFO to indicate that a small
C           pivot has been found.
C
            INFO = J
            RETURN
         END IF
   40 CONTINUE
      RETURN
C
C     End of F07ADH
C
      END
      SUBROUTINE F07ADJ(N,A,LDA,K1,K2,PIV,INCX)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C  Purpose
C  =======
C
C  F07ADJ performs a series of row interchanges on the matrix A.
C  One row interchange is initiated for each of rows K1 through K2 of A.
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The number of columns of the matrix A.
C
C  A       (input/output) REAL array, dimension (LDA,N)
C          On entry, the matrix of column dimension N to which the row
C          interchanges will be applied.
C          On exit, the permuted matrix.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.
C
C  K1      (input) INTEGER
C          The first element of PIV for which a row interchange will
C          be done.
C
C  K2      (input) INTEGER
C          The last element of PIV for which a row interchange will
C          be done.
C
C  PIV     (input) REAL array, dimension( M*abs(INCX) )
C          The vector of pivot indices.  Only the elements in positions
C          K1 through K2 of PIV are accessed.
C          PIV(K) = L implies rows K and L are to be interchanged.
C
C  INCX    (input) INTEGER
C          The increment between succesive values of PIV.  If PIV
C          is negative, the pivots are applied in reverse order.
C
C  This is a modified version of the LAPACK auxiliary routine
C  F07ADY/DLASWP, in which the INTEGER array IPIV has been replaced by
C  a REAL array PIV.
C
C     .. Scalar Arguments ..
      INTEGER           INCX, K1, K2, LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), PIV(*)
C     .. Local Scalars ..
      INTEGER           I, IP, IX
C     .. External Subroutines ..
      EXTERNAL          DSWAP_
C     .. Intrinsic Functions ..
      INTRINSIC         NINT
C     .. Executable Statements ..
C
C     Interchange row I with row PIV(I) for each of rows K1 through K2.
C
      IF (INCX.EQ.0) RETURN
      IF (INCX.GT.0) THEN
         IX = K1
      ELSE
         IX = 1 + (1-K2)*INCX
      END IF
      IF (INCX.EQ.1) THEN
         DO 20 I = K1, K2
            IP = NINT(PIV(I))
            IF (IP.NE.I) CALL DSWAP_(N,A(I,1),LDA,A(IP,1),LDA)
   20    CONTINUE
      ELSE IF (INCX.GT.1) THEN
         DO 40 I = K1, K2
            IP = NINT(PIV(IX))
            IF (IP.NE.I) CALL DSWAP_(N,A(I,1),LDA,A(IP,1),LDA)
            IX = IX + INCX
   40    CONTINUE
      ELSE IF (INCX.LT.0) THEN
         DO 60 I = K2, K1, -1
            IP = NINT(PIV(IX))
            IF (IP.NE.I) CALL DSWAP_(N,A(I,1),LDA,A(IP,1),LDA)
            IX = IX + INCX
   60    CONTINUE
      END IF
C
      RETURN
C
C     End of F07ADJ
C
      END
      INTEGER FUNCTION F07ZAY(NAME)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     F07ZAY returns a unique positive integer code
C     corresponding to a six-letter NAG routine name
C     given in NAME. If NAME is not recognised, 0
C     is returned.
C
C     .. Scalar Arguments ..
      CHARACTER*6             NAME
C     .. Local Scalars ..
      INTEGER                 J, K
      CHARACTER               NAME4, NAME5
C     .. Executable Statements ..
C
      IF (NAME(3:3).EQ.'7') THEN
         NAME4 = NAME(4:4)
         NAME5 = NAME(5:5)
      ELSE
         NAME4 = NAME(1:1)
         NAME5 = NAME(2:2)
      END IF
C
      IF (NAME4.EQ.'A') THEN
         J = 0
      ELSE IF (NAME4.EQ.'B') THEN
         J = 1
      ELSE IF (NAME4.EQ.'F') THEN
         J = 2
      ELSE IF (NAME4.EQ.'H') THEN
         J = 3
      ELSE IF (NAME4.EQ.'M') THEN
         J = 4
      ELSE IF (NAME4.EQ.'N') THEN
         J = 5
      ELSE IF (NAME4.EQ.'T') THEN
         J = 6
      ELSE
         J = -1
      END IF
C
      IF (NAME5.EQ.'D') THEN
         K = 0
      ELSE IF (NAME5.EQ.'J') THEN
         K = 1
      ELSE IF (NAME5.EQ.'R') THEN
         K = 2
      ELSE IF (NAME5.EQ.'W') THEN
         K = 3
      ELSE
         K = -1
      END IF
C
      IF (J.LT.0 .OR. K.LT.0) THEN
         F07ZAY = 0
      ELSE
C        F07ZAY is in the range 1-28.
         F07ZAY = 1 + 4*J + K
      END IF
C
      RETURN
C
      END
      SUBROUTINE F07ZAZ(ISPEC,NAME,IVAL,RWFLAG)
*
*     Mark 15 Release.  NAG Copyright 1991
*  -- NAG version of LAPACK auxiliary routine ILAENV
*     This version generated by program GENZAZ
*
*  Purpose
*  =======
*
*  F07ZAZ sets or returns problem-dependent
*  parameters for the local environment. See
*  ISPEC for a description of the parameters.
*
*  The problem-dependent parameters are contained
*  in the integer array IPARMS, and the value with
*  index ISPEC is set or copied to IVAL.
*
*  Arguments
*  =========
*
*  ISPEC (input) INTEGER
*     Specifies the parameter to be set or
*     returned by F07ZAZ.
*     = 1: the optimal blocksize; if this value
*          is 1, an unblocked algorithm will give
*          the best performance.
*     = 2: the minimum block size for which the
*          block routine should be used; if the
*          usable block size is less than this
*          value, an unblocked routine should be
*          used.
*     = 3: the crossover point (for N less than
*          this value, an unblocked routine should
*          be used)
*
*  NAME  (input) CHARACTER*(*)
*     The name of the calling subroutine.
*
*  IVAL  (input/output) INTEGER
*     the value of the parameter set or returned.
*
*  FLAG  (input) INTEGER
*     = 0: F07ZAZ returns in IVAL the value of
*          the parameter specified by ISPEC.
*     = 1: F07ZAZ sets the parameter specified
*          by ISPEC to the value in IVAL.
*
*  ==============================================
*
*     .. Parameters ..
      INTEGER           NSPECS, NCODES, MAXIC
      PARAMETER        (NSPECS=3,NCODES=  17,MAXIC=  28)
*     .. Scalar Arguments ..
      INTEGER           ISPEC, IVAL, RWFLAG
      CHARACTER*(*)     NAME
*     .. Local Scalars ..
      INTEGER           ICODE
*     .. Local Arrays ..
      INTEGER           IPARMS(NSPECS,NCODES),
     +                  POINT(MAXIC)
*     .. External Functions ..
      INTEGER           F07ZAY
      EXTERNAL          F07ZAY
*     .. Save statement ..
      SAVE              IPARMS, POINT
*     .. Data statements ..
      DATA              IPARMS(1,  1), IPARMS(2,  1), IPARMS(3,  1)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1,  2), IPARMS(2,  2), IPARMS(3,  2)
     +                  /   1,   1,   0 /
      DATA              IPARMS(1,  3), IPARMS(2,  3), IPARMS(3,  3)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1,  4), IPARMS(2,  4), IPARMS(3,  4)
     +                  /   1,   1,   0 /
      DATA              IPARMS(1,  5), IPARMS(2,  5), IPARMS(3,  5)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1,  6), IPARMS(2,  6), IPARMS(3,  6)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1,  7), IPARMS(2,  7), IPARMS(3,  7)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1,  8), IPARMS(2,  8), IPARMS(3,  8)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1,  9), IPARMS(2,  9), IPARMS(3,  9)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1, 10), IPARMS(2, 10), IPARMS(3, 10)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1, 11), IPARMS(2, 11), IPARMS(3, 11)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1, 12), IPARMS(2, 12), IPARMS(3, 12)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1, 13), IPARMS(2, 13), IPARMS(3, 13)
     +                  /  80,  25,   0 /
      DATA              IPARMS(1, 14), IPARMS(2, 14), IPARMS(3, 14)
     +                  /  48,  10,   0 /
      DATA              IPARMS(1, 15), IPARMS(2, 15), IPARMS(3, 15)
     +                  /   1,   1,   0 /
      DATA              IPARMS(1, 16), IPARMS(2, 16), IPARMS(3, 16)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1, 17), IPARMS(2, 17), IPARMS(3, 17)
     +                  /   1,   0,   0 /
      DATA              POINT /
     +                   1,   2,   3,   4,   5,   0,   6,   0,
     +                   7,   8,   9,  10,  11,   0,  12,   0,
     +                  13,   0,  14,   0,   0,   0,  15,   0,
     +                   0,  16,   0,  17
     +                  /
*     .. Executable Statements ..
*
*     Convert the NAG name to an integer code.
      ICODE = F07ZAY(NAME)
*
      IF (ISPEC.LT.1 .OR. ISPEC.GT.NSPECS) THEN
*        Invalid value for ISPEC
         IVAL = -1
      ELSE IF (ICODE.EQ.0) THEN
*        Invalid value for NAME
         IVAL = -2
      ELSE IF (POINT(ICODE).EQ.0) THEN
*        Invalid value for NAME
         IVAL = -2
      ELSE IF (RWFLAG.EQ.0) THEN
*        Read the value of a parameter
         IVAL = IPARMS(ISPEC,POINT(ICODE))
      ELSE
*        Set the value of a parameter
         IPARMS(ISPEC,POINT(ICODE)) = IVAL
      END IF
*
      RETURN
*
*     End of F07ZAZ
*
      END
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END
      SUBROUTINE P01ABW(N,NAME,INFORM,IERR,SRNAME)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     P01ABW increases the value of IERR by 1 and, if
C
C        ( mod( INFORM, 10 ).ne.1 ).or.( mod( INFORM/10, 10 ).ne.0 )
C
C     writes a message on the current error message channel giving the
C     value of N, a message to say that N is invalid and the strings
C     NAME and SRNAME.
C
C     NAME must be the name of the actual argument for N and SRNAME must
C     be the name of the calling routine.
C
C     This routine is intended for use when N is an invalid input
C     parameter to routine SRNAME. For example
C
C        IERR = 0
C        IF( N.NE.'Valid value' )
C     $     CALL P01ABW( N, 'N', IDIAG, IERR, SRNAME )
C
C  -- Written on 15-November-1984.
C     Sven Hammarling, Nag Central Office.
C
C     .. Scalar Arguments ..
      INTEGER           IERR, INFORM
      CHARACTER*(*)     N
      CHARACTER*(*)     NAME, SRNAME
C     .. Local Scalars ..
      INTEGER           NERR
C     .. Local Arrays ..
      CHARACTER*65      REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      IERR = IERR + 1
      IF ((MOD(INFORM,10).NE.1) .OR. (MOD(INFORM/10,10).NE.0)) THEN
         CALL X04AAF(0,NERR)
         WRITE (REC,FMT=99999) NAME, SRNAME, N
         CALL X04BAF(NERR,' ')
         CALL X04BAF(NERR,REC(1))
         CALL X04BAF(NERR,REC(2))
         CALL X04BAF(NERR,REC(3))
      END IF
      RETURN
C
C
C     End of P01ABW.
C
99999 FORMAT (' *****  Parameter  ',A,'  is invalid in routine  ',A,
     *  '  ***** ',/8X,'Value supplied is',/8X,A)
      END
      SUBROUTINE P01ABY(N,NAME,INFORM,IERR,SRNAME)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     P01ABY increases the value of IERR by 1 and, if
C
C        ( mod( INFORM, 10 ).ne.1 ).or.( mod( INFORM/10, 10 ).ne.0 )
C
C     writes a message on the current error message channel giving the
C     value of N, a message to say that N is invalid and the strings
C     NAME and SRNAME.
C
C     NAME must be the name of the actual argument for N and SRNAME must
C     be the name of the calling routine.
C
C     This routine is intended for use when N is an invalid input
C     parameter to routine SRNAME. For example
C
C        IERR = 0
C        IF( N.LT.1 )CALL P01ABY( N, 'N', IDIAG, IERR, SRNAME )
C
C  -- Written on 23-February-1984.  Sven.
C
C     .. Scalar Arguments ..
      INTEGER           IERR, INFORM, N
      CHARACTER*(*)     NAME, SRNAME
C     .. Local Scalars ..
      INTEGER           NERR
C     .. Local Arrays ..
      CHARACTER*65      REC(2)
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      IERR = IERR + 1
      IF ((MOD(INFORM,10).NE.1) .OR. (MOD(INFORM/10,10).NE.0)) THEN
         CALL X04AAF(0,NERR)
         WRITE (REC,FMT=99999) NAME, SRNAME, N
         CALL X04BAF(NERR,' ')
         CALL X04BAF(NERR,REC(1))
         CALL X04BAF(NERR,REC(2))
      END IF
      RETURN
C
C
C     End of P01ABY.
C
99999 FORMAT (' *****  Parameter  ',A,'  is invalid in routine  ',A,
     *  '  ***** ',/8X,'Value supplied is ',I6)
      END
      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
      STOP
      END
      INTEGER FUNCTION P01ACF(IFAIL,IERROR,SRNAME,VARBNM,NREC,REC)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     P01ACF is the error-handling routine for the F06 AND F07
C     Chapters of the NAG Fortran Library. It is a slightly modified
C     version of P01ABF.
C
C     P01ACF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ACF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ACF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME, VARBNM
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR, VARLEN
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, LEN, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
         VARLEN = 0
         DO 20 I = LEN(VARBNM), 1, -1
            IF (VARBNM(I:I).NE.' ') THEN
               VARLEN = I
               GO TO 40
            END IF
   20    CONTINUE
   40    CONTINUE
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 60 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   60       CONTINUE
            IF (IFAIL.NE.-13) THEN
               IF (VARLEN.NE.0) THEN
                  WRITE (MESS,FMT=99999) SRNAME, VARBNM(1:VARLEN),
     *              IERROR
               ELSE
                  WRITE (MESS,FMT=99998) SRNAME
               END IF
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ACF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': ',A,
     *       ' =',I6)
99998 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A)
      END
      DOUBLE PRECISION FUNCTION X02AJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
      DOUBLE PRECISION X02CON
      DATA X02CON /1.11022302462516D-16 /
C     .. Executable Statements ..
      X02AJF = X02CON
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AKF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  B**(EMIN-1)  (THE SMALLEST POSITIVE MODEL NUMBER)
C
      DOUBLE PRECISION X02CON
      DATA X02CON /2.22507385850721D-308 /
C     .. Executable Statements ..
      X02AKF = X02CON
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AMF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE 'SAFE RANGE' PARAMETER
C     I.E. THE SMALLEST POSITIVE MODEL NUMBER Z SUCH THAT
C     FOR ANY X WHICH SATISFIES X.GE.Z AND X.LE.1/Z
C     THE FOLLOWING CAN BE COMPUTED WITHOUT OVERFLOW, UNDERFLOW OR OTHER
C     ERROR
C
C        -X
C        1.0/X
C        SQRT(X)
C        LOG(X)
C        EXP(LOG(X))
C        Y**(LOG(X)/LOG(Y)) FOR ANY Y
C
      DOUBLE PRECISION X02CON
      DATA X02CON /2.22507385850739D-308 /
C     .. Executable Statements ..
      X02AMF = X02CON
      RETURN
      END
      INTEGER FUNCTION X02BHF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE MODEL PARAMETER, B.
C
C     .. Executable Statements ..
      X02BHF =     2
      RETURN
      END
      SUBROUTINE X03AAF(A,ISIZEA,B,ISIZEB,N,ISTEPA,ISTEPB,C1,C2,D1,D2,
     *                  SW,IFAIL)
C     MARK 10 RE-ISSUE. NAG COPYRIGHT 1982
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. IER-524 (AUG 1986).
C     DOUBLE PRECISION BASE VERSION
C
C     CALCULATES THE VALUE OF A SCALAR PRODUCT USING BASIC OR
C     ADDITIONAL PRECISION AND ADDS IT TO A BASIC OR ADDITIONAL
C     PRECISION INITIAL VALUE.
C
C     FOR THIS DOUBLE PRECISION VERSION, ALL ADDITIONAL (I.E.
C     QUADRUPLE) PRECISION COMPUTATION IS PERFORMED BY THE AUXILIARY
C     ROUTINE X03AAY.  SEE THE COMMENTS AT THE HEAD OF THAT ROUTINE
C     CONCERNING IMPLEMENTATION.
C
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='X03AAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C1, C2, D1, D2
      INTEGER           IFAIL, ISIZEA, ISIZEB, ISTEPA, ISTEPB, N
      LOGICAL           SW
C     .. Array Arguments ..
      DOUBLE PRECISION  A(ISIZEA), B(ISIZEB)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, IERR, IS, IT
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          X03AAY
C     .. Executable Statements ..
      IERR = 0
      IF (ISTEPA.LE.0 .OR. ISTEPB.LE.0) IERR = 1
      IF (ISIZEA.LT.(N-1)*ISTEPA+1 .OR. ISIZEB.LT.(N-1)*ISTEPB+1)
     *    IERR = 2
      IF (IERR.EQ.0) GO TO 20
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
   20 IFAIL = 0
      IF (SW) GO TO 80
C
C     BASIC PRECISION CALCULATION
C
      SUM = C1
      IF (N.LT.1) GO TO 60
      IS = 1
      IT = 1
      DO 40 I = 1, N
         SUM = SUM + A(IS)*B(IT)
         IS = IS + ISTEPA
         IT = IT + ISTEPB
   40 CONTINUE
   60 D1 = SUM
      D2 = 0.0D0
      RETURN
C
C     ADDITIONAL PRECISION COMPUTATION
C
   80 CALL X03AAY(A,ISIZEA,B,ISIZEB,N,ISTEPA,ISTEPB,C1,C2,D1,D2)
      RETURN
      END
      SUBROUTINE X03AAY(A,ISIZEA,B,ISIZEB,N,ISTEPA,ISTEPB,C1,C2,D1,D2)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DOUBLE PRECISION BASE VERSION
C
C     PERFORMS QUADRUPLE PRECISION COMPUTATION FOR X03AAF
C
C     **************************************************************
C     THIS FORTRAN CODE WILL NOT WORK ON ALL MACHINES. IT IS
C     PRESUMED TO WORK IF THE MACHINE SATISFIES ONE OF THE FOLLOWING
C     ASSUMPTIONS.
C
C     A.) THERE ARE AN EVEN NUMBER OF B-ARY DIGITS IN THE MANTISSA
C     -   OF A DOUBLE PRECISION NUMBER (WHERE B IS THE BASE FOR THE
C     -   REPRESENTATION OF FLOATING-POINT NUMBERS), AND THE
C     -   COMPUTED RESULT OF A DOUBLE PRECISION ADDITION,
C     -   SUBTRACTION OR MULTIPLICATION IS EITHER CORRECTLY ROUNDED
C     -   OR CORRECTLY CHOPPED.
C
C     B.) FLOATING-POINT NUMBERS ARE REPRESENTED TO THE BASE 2 (WITH
C     -   ANY NUMBER OF BITS IN THE MANTISSA OF A DOUBLE PRECISION
C     -   NUMBER), AND THE COMPUTED RESULT OF A DOUBLE PRECISION
C     -   ADDITION, SUBTRACTION OR MULTIPLICATION IS CORRECTLY
C     -   ROUNDED.
C
C     REFERENCES-
C
C     T.J. DEKKER  A FLOATING-POINT TECHNIQUE FOR EXTENDING THE
C     AVAILABLE PRECISION. NUMER. MATH. 18, 224-242 (1971)
C
C     S. LINNAINMAA  SOFTWARE FOR DOUBLED-PRECISION FLOATING-POINT
C     COMPUTATIONS. ACM TRANS. MATH. SOFTWARE 7, 272-283 (1981)
C
C     IF THE ABOVE ASSUMPTIONS ARE NOT SATISFIED, THIS ROUTINE MUST
C     BE IMPLEMENTED IN ASSEMBLY LANGUAGE. IN ANY CASE ASSEMBLY
C     LANGUAGE MAY BE PREFERABLE FOR GREATER EFFICIENCY.  CONSULT
C     NAG CENTRAL OFFICE.
C
C     THE ROUTINE MUST SIMULATE THE FOLLOWING QUADRUPLE PRECISION
C     CODING IN PSEUDO-FORTRAN, WHERE
C     - QEXTD CONVERTS FROM DOUBLE TO QUADRUPLE PRECISION
C     - DBLEQ CONVERTS FROM QUADRUPLE TO DOUBLE PRECISION.
C
C     QUADRUPLE PRECISION SUM
C     SUM = QEXTD(C1) + QEXTD(C2)
C     IF (N.LT.1) GO TO 80
C     IS = 1
C     IT = 1
C     DO 60 I = 1, N
C        SUM = SUM + QEXTD(A(IS))*QEXTD(B(IT))
C        IS = IS + ISTEPA
C        IT = IT + ISTEPB
C  60 CONTINUE
C  80 D1 = SUM -- ROUNDED TO DOUBLE PRECISION
C     D2 = DBLEQ(SUM-QEXTD(D1))
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C1, C2, D1, D2
      INTEGER           ISIZEA, ISIZEB, ISTEPA, ISTEPB, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(ISIZEA), B(ISIZEB)
C     .. Local Scalars ..
      DOUBLE PRECISION  AA, BB, CONS, HA, HB, R, S, SUM, SUMM, TA, TB,
     *                  Z, ZZ
      INTEGER           I, IS, IT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
C     ************* IMPLEMENTATION-DEPENDENT CONSTANT **************
C     CONS MUST BE SET TO  B**(T - INT(T/2)) + 1 , WHERE T IS THE
C     NUMBER OF B-ARY DIGITS IN THE MANTISSA OF A DOUBLE PRECISION
C     NUMBER.
C     FOR B = 16 AND T = 14 (E.G. IBM 370) OR
C     FOR B = 2 AND T = 56 (E.G. DEC VAX-11)
C      DATA CONS /268435457.0D0/
      DATA CONS /1.34217729000000D+8/ 
C     **************************************************************
C     .. Executable Statements ..
      SUM = C1 + C2
      SUMM = (C1-SUM) + C2
      IF (N.LT.1) GO TO 80
      IS = 1
      IT = 1
      DO 60 I = 1, N
         AA = A(IS)
         BB = B(IT)
         Z = AA*CONS
         HA = (AA-Z) + Z
         TA = AA - HA
         Z = BB*CONS
         HB = (BB-Z) + Z
         TB = BB - HB
         Z = AA*BB
         ZZ = (((HA*HB-Z)+HA*TB)+TA*HB) + TA*TB
         R = Z + SUM
         IF (ABS(Z).GT.ABS(SUM)) GO TO 20
         S = (((SUM-R)+Z)+ZZ) + SUMM
         GO TO 40
   20    S = (((Z-R)+SUM)+SUMM) + ZZ
   40    SUM = R + S
         SUMM = (R-SUM) + S
         IS = IS + ISTEPA
         IT = IT + ISTEPB
   60 CONTINUE
C  80 D1 = SUM + (SUMM+SUMM)
C     *************** IMPLEMENTATION DEPENDENT CODE ****************
C     THE PREVIOUS STATEMENT ASSUMES THAT THE RESULT OF A DOUBLE
C     PRECISION ADDITION IS TRUNCATED.  IF IT IS ROUNDED, THEN
C     THE STATEMENT MUST BE CHANGED TO
   80 D1 = SUM + SUMM
C     **************************************************************
      D2 = (SUM-D1) + SUMM
      RETURN
      END
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-829 (DEC 1989).
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/6/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END
