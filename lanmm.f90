!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   jrnmm : Jansen-Rit Neural Mass Model for a cortical neuron.
!       Here we analyze the homogeneous component of properly normalized
!       network, which boils down to study a self-coupled JR-NMM    
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION sigm(v,vr)
       DOUBLE PRECISION  e0,r
       DOUBLE PRECISION :: v
       DOUBLE PRECISION :: vr
       r=0.56
       e0 = 2.5
       sigm = 2.0 * e0 / (1.0 + EXP(r * (vr - v)))
      END FUNCTION sigm
      DOUBLE PRECISION FUNCTION dsigm(v,vr)
       DOUBLE PRECISION  e0,r
       DOUBLE PRECISION :: v
       DOUBLE PRECISION :: vr
       r=0.56
       e0 = 2.5
       dsigm =2*e0*r*exp(r*(vr-v))/((1+exp(r*(vr-v)))**2)
      END FUNCTION dsigm


      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)


      DOUBLE PRECISION :: C(13) = (/108.0, 33.7, 1.0, 135.0, 33.75, 70.0, 550.0, 1.0, 200.0, 100.0, 80.0, 200.0, 30.0/)
      DOUBLE PRECISION :: A(5) = (/3.25, 3.25, -22.0, 3.25, -30.0/)
      DOUBLE PRECISION :: b(5) = (/100.0, 100.0, 50.0, 100.0, 220.0/) * 1.0
      DOUBLE PRECISION :: incoming(5) = (/0.0, 0.0, 0.0, 0.0, 0.0/)
      DOUBLE PRECISION :: w(5) = (/6.0, 6.0, 6.0, 1.0, 6.0/)
      DOUBLE PRECISION y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, p1, p2
      DOUBLE PRECISION sigm,dsigm,s
      INTEGER i,j



       y0=U(1)
       y1=U(2)
       y2=U(3)
       y3=U(4)
       y4=U(5)
       y5=U(6)
       y6=U(7)
       y7=U(8)
       y8=U(9)
       y9=U(10)


       p1 =PAR(1) 
       p2 = PAR(2)

!       incoming(1) = C(1) * y1 + C(2) * y2 + C(11) * y3 + A(1) * p1 / b(1) 
!       incoming(2) = C(4) * y0
!       incoming(3) = C(5) * y0
!       incoming(4) = C(6) * y3 +  C(7) * y4 + C(12) * y0 + A(4) * p2 / b(4) 
!       incoming(5) = C(9) * y3 + C(10) * y4 + C(13) * y0

       F(1)= y5
       F(2)= y6
       F(3)= y7
       F(4)= y8
       F(5)= y9
       
       F(6)= A(1) * b(1) * sigm(C(1) * y1 + C(2) * y2 + C(11) * y3 + A(1) * p1 / b(1), w(1)) - 2 * b(1) * y5 - b(1)**2 * y0
       F(7)= A(2) * b(2) * sigm(C(4) * y0, w(2)) - 2 * b(2) * y6 - b(2)**2 * y1
       F(8)= A(3) * b(3) * sigm(C(5) * y0, w(3)) - 2 * b(3) * y7 - b(3)**2 * y2
       F(9)= A(4) * b(4) * sigm(C(6) * y3 +  C(7) * y4 + C(12) * y0 + A(4) * p2 / b(4) , w(4)) - 2 * b(4) * y8 - b(4)**2 * y3
       F(10)=A(5) * b(5) * sigm(C(9) * y3 + C(10) * y4 + C(13) * y0,w(5)) - 2 * b(5) * y9 - b(5)**2 * y4


        IF(IJAC.EQ.0)RETURN
        do i = 1, 10
                do j = 1, 10
                        DFDU(i,j)=0;
                end do
        end do


        DFDU(1,6)=1;
        DFDU(2,7)=1;
        DFDU(3,8)=1;
        DFDU(4,9)=1;
        DFDU(5,10)=1;

        DFDU(6,6)=-2*b(1);
        DFDU(7,7)=-2*b(2);
        DFDU(8,8)=-2*b(3);
        DFDU(9,9) = -2*b(4);
        DFDU(10,10) = -2*b(5);

        DFDU(6,1)=-b(1)**2;
        DFDU(6,2)=A(1)*b(1)*C(1)*dsigm(C(1) * y1 + C(2) * y2 + C(11) * y3 + A(1) * p1 / b(1),w(1));
        DFDU(6,3)=A(1)*b(1)*C(2)*dsigm(C(1) * y1 + C(2) * y2 + C(11) * y3 + A(1) * p1 / b(1),w(1));
        DFDU(6,4)=A(1)*b(1)*C(11)*dsigm(C(1) * y1 + C(2) * y2 + C(11) * y3 + A(1) * p1 / b(1),w(1));

        DFDU(7,1)=A(2)*b(2)*C(4)*dsigm(C(4) * y0,w(2));
        DFDU(7,2)=-b(2)**2;

        DFDU(8,1)=A(3)*b(3)*C(5)*dsigm(C(5) * y0,w(3));
        DFDU(8,3)=-b(3)**2;

        DFDU(9,1) = A(4)*b(4)*dsigm(C(6) * y3 +  C(7) * y4 + C(12) * y0 + A(4) * p2 / b(4) ,w(4))*C(12); 
        DFDU(9,4) = A(4)*b(4)*dsigm(C(6) * y3 +  C(7) * y4 + C(12) * y0 + A(4) * p2 / b(4) ,w(4))*C(6) - b(4)**2;
        DFDU(9,5) = A(4)*b(4)*dsigm(C(6) * y3 +  C(7) * y4 + C(12) * y0 + A(4) * p2 / b(4) ,w(4))*C(7); 

        DFDU(10,1) = A(5)*b(5)*dsigm(C(9) * y3 + C(10) * y4 + C(13) * y0,w(5))*C(13);
        DFDU(10,4) = A(5)*b(5)*dsigm(C(9) * y3 + C(10) * y4 + C(13) * y0,w(5))*C(9);
        DFDU(10,5) = A(5)*b(5)*dsigm(C(9) * y3 + C(10) * y4 + C(13) * y0,w(5))*C(10) - b(5)**2;


        IF(IJAC.EQ.1)RETURN

        do i = 1, 10
                do j = 1, 2
                        DFDP(i,j)=0;
                end do
        end do

        DFDP(6,1) = A(1)*dsigm(C(1) * y1 + C(2) * y2 + C(11) * y3 + A(1) * p1 / b(1),w(1))*A(1)
        DFDP(9,2) = A(4)*dsigm(C(6) * y3 +  C(7) * y4 + C(12) * y0 + A(4) * p2 / b(4) ,w(4))*A(4)

      END SUBROUTINE FUNC
!-----------------------------------------------------------
!-----------------------------------------------------------
      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=00.0 ! P
      PAR(2)=00.0 ! EPSILON

      U(1)=0.01
      U(2)=0.01
      U(3)=0.01
      U(4)=0.01
      U(5)=0.01
      U(6)=0.01
      U(7)=0.01
      U(8)=0.01
      U(9)=0.01
      U(10)=0.01
      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT


!      AUX_Y=C(1)*U(2,:) + C(2)*U(3,:) + C(11)*U(4,:)!!vP1
!      AUX_Y=C(6)*U(4,:) + C(7)*U(5,:) + C(12)*U(1,:)!!vP2


!     ------ --------- -------- ----- Max and min of P1
      DOUBLE PRECISION FUNCTION GETUMN1(U,NDX,NTST,NCOL)
!     ------ --------- -------- -----
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION :: C(13) = (/108.0, 33.75, 1.0, 135.0, 33.75, 70.0, 550.0, 1.0, 200.0, 100.0, 80.0, 200.0, 30.0/)
      GETUMN1 = MINVAL(C(1)*U(2,:) + C(2)*U(3,:) + C(11)*U(4,:))! + C(2)*U(3,:) + C(11)*U(4,:) + C(6)*U(11,:))

      END FUNCTION GETUMN1

      DOUBLE PRECISION FUNCTION GETUMX1(U,NDX,NTST,NCOL)
!     ------ --------- -------- -----
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION :: C(13) = (/108.0, 33.75, 1.0, 135.0, 33.75, 70.0, 550.0, 1.0, 200.0, 100.0, 80.0, 200.0, 30.0/)
      GETUMX1 = MAXVAL(C(1)*U(2,:) + C(2)*U(3,:) + C(11)*U(4,:))! + C(2)*U(3,:) + C(11)*U(4,:) + C(6)*U(11,:))

      END FUNCTION GETUMX1
!     ------ --------- -------- ----- Max and min of P2
      DOUBLE PRECISION FUNCTION GETUMN2(U,NDX,NTST,NCOL)
!     ------ --------- -------- -----
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION :: C(13) = (/108.0, 33.75, 1.0, 135.0, 33.75, 70.0, 550.0, 1.0, 200.0, 100.0, 80.0, 200.0, 30.0/)
      GETUMN2 = MINVAL(C(6)*U(4,:) + C(7)*U(5,:) + C(12)*U(1,:))! + C(2)*U(3,:) + C(11)*U(4,:) + C(6)*U(11,:))

      END FUNCTION GETUMN2

      DOUBLE PRECISION FUNCTION GETUMX2(U,NDX,NTST,NCOL)
!     ------ --------- -------- -----
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION :: C(13) = (/108.0, 33.75, 1.0, 135.0, 33.75, 70.0, 550.0, 1.0, 200.0, 100.0, 80.0, 200.0, 30.0/)
      GETUMX2 = MAXVAL(C(6)*U(4,:) + C(7)*U(5,:) + C(12)*U(1,:))! + C(2)*U(3,:) + C(11)*U(4,:) + C(6)*U(11,:))

      END FUNCTION GETUMX2
!     ------ --------- -------- -----
!     ------ --------- -------- -----
!     ------ --------- -------- -----

      SUBROUTINE PVLS(NDIM,U,PAR)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      DOUBLE PRECISION, EXTERNAL :: GETP,GETUMN1, GETUMX1, GETUMN2, GETUMX2
      DOUBLE PRECISION :: C(13) = (/108.0, 33.75, 1.0, 135.0, 33.75, 70.0, 550.0, 1.0, 200.0, 100.0, 80.0, 200.0, 30.0/)
      INTEGER NDX,NCOL,NTST
       NDX=NINT(GETP('NDX',0,U))
       NTST=NINT(GETP('NTST',0,U))
       NCOL=NINT(GETP('NCOL',0,U))
       PAR(3)=GETUMN1(U,NDX,NTST,NCOL) 
       PAR(4)=GETUMX1(U,NDX,NTST,NCOL)
       PAR(5)=GETUMN2(U,NDX,NTST,NCOL) 
       PAR(6)=GETUMX2(U,NDX,NTST,NCOL)
       PAR(7)= C(1)*GETP('MIN',2,U) + C(2)*GETP('MIN',3,U) + C(11)*GETP('MIN',4,U) 
       PAR(8)= C(6)*GETP('MIN',4,U) +  C(7)*GETP('MIN',5,U) + C(12)*GETP('MIN',1,U)

      END SUBROUTINE PVLS

