!===============================================================================
! Copyright 2005-2018 Intel Corporation All Rights Reserved.
!
! The source code,  information  and material  ("Material") contained  herein is
! owned by Intel Corporation or its  suppliers or licensors,  and  title to such
! Material remains with Intel  Corporation or its  suppliers or  licensors.  The
! Material  contains  proprietary  information  of  Intel or  its suppliers  and
! licensors.  The Material is protected by  worldwide copyright  laws and treaty
! provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
! modified, published,  uploaded, posted, transmitted,  distributed or disclosed
! in any way without Intel's prior express written permission.  No license under
! any patent,  copyright or other  intellectual property rights  in the Material
! is granted to  or  conferred  upon  you,  either   expressly,  by implication,
! inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
! property rights must be express and approved by Intel in writing.
!
! Unless otherwise agreed by Intel in writing,  you may not remove or alter this
! notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
! suppliers or licensors in any way.
!===============================================================================

!  Content:
!      F95 interface for BLAS routines
!*******************************************************************************
! This file was generated automatically!
!*******************************************************************************

MODULE F95_PRECISION
    INTEGER, PARAMETER :: SP = KIND(1.0E0)
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
END MODULE F95_PRECISION

MODULE BLAS95

INTERFACE ASUM
    PURE FUNCTION SASUM_F95(X)
        ! Fortran77 call:
        ! SASUM(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP) :: SASUM_F95
        REAL(WP), INTENT(IN) :: X(:)
    END FUNCTION SASUM_F95
    PURE FUNCTION SCASUM_F95(X)
        ! Fortran77 call:
        ! SCASUM(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP) :: SCASUM_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
    END FUNCTION SCASUM_F95
    PURE FUNCTION DASUM_F95(X)
        ! Fortran77 call:
        ! DASUM(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP) :: DASUM_F95
        REAL(WP), INTENT(IN) :: X(:)
    END FUNCTION DASUM_F95
    PURE FUNCTION DZASUM_F95(X)
        ! Fortran77 call:
        ! DZASUM(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP) :: DZASUM_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
    END FUNCTION DZASUM_F95
END INTERFACE ASUM

INTERFACE AXPY
        ! Default A=1
    PURE SUBROUTINE SAXPY_F95(X,Y,A)
        ! Fortran77 call:
        ! SAXPY(N,A,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN), OPTIONAL :: A
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SAXPY_F95
    PURE SUBROUTINE DAXPY_F95(X,Y,A)
        ! Fortran77 call:
        ! DAXPY(N,A,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN), OPTIONAL :: A
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DAXPY_F95
    PURE SUBROUTINE CAXPY_F95(X,Y,A)
        ! Fortran77 call:
        ! CAXPY(N,A,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: A
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE CAXPY_F95
    PURE SUBROUTINE ZAXPY_F95(X,Y,A)
        ! Fortran77 call:
        ! ZAXPY(N,A,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: A
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE ZAXPY_F95
END INTERFACE AXPY

INTERFACE COPY
    PURE SUBROUTINE SCOPY_F95(X,Y)
        ! Fortran77 call:
        ! SCOPY(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SCOPY_F95
    PURE SUBROUTINE DCOPY_F95(X,Y)
        ! Fortran77 call:
        ! DCOPY(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DCOPY_F95
    PURE SUBROUTINE CCOPY_F95(X,Y)
        ! Fortran77 call:
        ! CCOPY(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE CCOPY_F95
    PURE SUBROUTINE ZCOPY_F95(X,Y)
        ! Fortran77 call:
        ! ZCOPY(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE ZCOPY_F95
END INTERFACE COPY

INTERFACE DOT
    PURE FUNCTION SDOT_F95(X,Y)
        ! Fortran77 call:
        ! SDOT(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP) :: SDOT_F95
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END FUNCTION SDOT_F95
    PURE FUNCTION DDOT_F95(X,Y)
        ! Fortran77 call:
        ! DDOT(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP) :: DDOT_F95
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END FUNCTION DDOT_F95
END INTERFACE DOT

INTERFACE SDOT
    PURE FUNCTION SDSDOT_F95(SX,SY,SB)
        ! Fortran77 call:
        ! SDSDOT(N,SB,SX,INCX,SY,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP) :: SDSDOT_F95
        REAL(WP), INTENT(IN) :: SB
        REAL(WP), INTENT(IN) :: SX(:)
        REAL(WP), INTENT(IN) :: SY(:)
    END FUNCTION SDSDOT_F95
    PURE FUNCTION DSDOT_F95(SX,SY)
        ! Fortran77 call:
        ! DSDOT(N,SX,INCX,SY,INCY)
        USE F95_PRECISION, ONLY: WP => DP, SP
        REAL(WP) :: DSDOT_F95
        REAL(SP), INTENT(IN) :: SX(:)
        REAL(SP), INTENT(IN) :: SY(:)
    END FUNCTION DSDOT_F95
END INTERFACE SDOT

INTERFACE DOTC
    PURE FUNCTION CDOTC_F95(X,Y)
        ! Fortran77 call:
        ! CDOTC(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP) :: CDOTC_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END FUNCTION CDOTC_F95
    PURE FUNCTION ZDOTC_F95(X,Y)
        ! Fortran77 call:
        ! ZDOTC(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP) :: ZDOTC_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END FUNCTION ZDOTC_F95
END INTERFACE DOTC

INTERFACE DOTU
    PURE FUNCTION CDOTU_F95(X,Y)
        ! Fortran77 call:
        ! CDOTU(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP) :: CDOTU_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END FUNCTION CDOTU_F95
    PURE FUNCTION ZDOTU_F95(X,Y)
        ! Fortran77 call:
        ! ZDOTU(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP) :: ZDOTU_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END FUNCTION ZDOTU_F95
END INTERFACE DOTU

INTERFACE NRM2
    PURE FUNCTION SNRM2_F95(X)
        ! Fortran77 call:
        ! SNRM2(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP) :: SNRM2_F95
        REAL(WP), INTENT(IN) :: X(:)
    END FUNCTION SNRM2_F95
    PURE FUNCTION DNRM2_F95(X)
        ! Fortran77 call:
        ! DNRM2(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP) :: DNRM2_F95
        REAL(WP), INTENT(IN) :: X(:)
    END FUNCTION DNRM2_F95
    PURE FUNCTION SCNRM2_F95(X)
        ! Fortran77 call:
        ! SCNRM2(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP) :: SCNRM2_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
    END FUNCTION SCNRM2_F95
    PURE FUNCTION DZNRM2_F95(X)
        ! Fortran77 call:
        ! DZNRM2(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP) :: DZNRM2_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
    END FUNCTION DZNRM2_F95
END INTERFACE NRM2

INTERFACE ROT
    PURE SUBROUTINE SROT_F95(X,Y,C,S)
        ! Fortran77 call:
        ! SROT(N,X,INCX,Y,INCY,C,S)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN) :: C
        REAL(WP), INTENT(IN) :: S
        REAL(WP), INTENT(INOUT) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SROT_F95
    PURE SUBROUTINE DROT_F95(X,Y,C,S)
        ! Fortran77 call:
        ! DROT(N,X,INCX,Y,INCY,C,S)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN) :: C
        REAL(WP), INTENT(IN) :: S
        REAL(WP), INTENT(INOUT) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DROT_F95
    PURE SUBROUTINE CSROT_F95(X,Y,C,S)
        ! Fortran77 call:
        ! CSROT(N,X,INCX,Y,INCY,C,S)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN) :: C
        REAL(WP), INTENT(IN) :: S
        COMPLEX(WP), INTENT(INOUT) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE CSROT_F95
    PURE SUBROUTINE ZDROT_F95(X,Y,C,S)
        ! Fortran77 call:
        ! ZDROT(N,X,INCX,Y,INCY,C,S)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN) :: C
        REAL(WP), INTENT(IN) :: S
        COMPLEX(WP), INTENT(INOUT) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE ZDROT_F95
END INTERFACE ROT

INTERFACE ROTG
    PURE SUBROUTINE SROTG_F95(A,B,C,S)
        ! Fortran77 call:
        ! SROTG(A,B,C,S)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(INOUT) :: A
        REAL(WP), INTENT(INOUT) :: B
        REAL(WP), INTENT(OUT) :: C
        REAL(WP), INTENT(OUT) :: S
    END SUBROUTINE SROTG_F95
    PURE SUBROUTINE DROTG_F95(A,B,C,S)
        ! Fortran77 call:
        ! DROTG(A,B,C,S)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(INOUT) :: A
        REAL(WP), INTENT(INOUT) :: B
        REAL(WP), INTENT(OUT) :: C
        REAL(WP), INTENT(OUT) :: S
    END SUBROUTINE DROTG_F95
    PURE SUBROUTINE CROTG_F95(A,B,C,S)
        ! Fortran77 call:
        ! CROTG(A,B,C,S)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(INOUT) :: A
        COMPLEX(WP), INTENT(INOUT) :: B
        REAL(WP), INTENT(OUT) :: C
        COMPLEX(WP), INTENT(OUT) :: S
    END SUBROUTINE CROTG_F95
    PURE SUBROUTINE ZROTG_F95(A,B,C,S)
        ! Fortran77 call:
        ! ZROTG(A,B,C,S)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(INOUT) :: A
        COMPLEX(WP), INTENT(INOUT) :: B
        REAL(WP), INTENT(OUT) :: C
        COMPLEX(WP), INTENT(OUT) :: S
    END SUBROUTINE ZROTG_F95
END INTERFACE ROTG

INTERFACE ROTM
    PURE SUBROUTINE SROTM_F95(X,Y,PARAM)
        ! Fortran77 call:
        ! SROTM(N,X,INCX,Y,INCY,PARAM)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(INOUT) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
        REAL(WP), INTENT(IN) :: PARAM(5)
    END SUBROUTINE SROTM_F95
    PURE SUBROUTINE DROTM_F95(X,Y,PARAM)
        ! Fortran77 call:
        ! DROTM(N,X,INCX,Y,INCY,PARAM)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(INOUT) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
        REAL(WP), INTENT(IN) :: PARAM(5)
    END SUBROUTINE DROTM_F95
END INTERFACE ROTM

INTERFACE ROTMG
    PURE SUBROUTINE SROTMG_F95(D1,D2,X1,Y1,PARAM)
        ! Fortran77 call:
        ! SROTMG(D1,D2,X1,Y1,PARAM)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(INOUT) :: D1
        REAL(WP), INTENT(INOUT) :: D2
        REAL(WP), INTENT(INOUT) :: X1
        REAL(WP), INTENT(IN) :: Y1
        REAL(WP), INTENT(OUT) :: PARAM(5)
    END SUBROUTINE SROTMG_F95
    PURE SUBROUTINE DROTMG_F95(D1,D2,X1,Y1,PARAM)
        ! Fortran77 call:
        ! DROTMG(D1,D2,X1,Y1,PARAM)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(INOUT) :: D1
        REAL(WP), INTENT(INOUT) :: D2
        REAL(WP), INTENT(INOUT) :: X1
        REAL(WP), INTENT(IN) :: Y1
        REAL(WP), INTENT(OUT) :: PARAM(5)
    END SUBROUTINE DROTMG_F95
END INTERFACE ROTMG

INTERFACE SCAL
    PURE SUBROUTINE SSCAL_F95(X,A)
        ! Fortran77 call:
        ! SSCAL(N,A,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN) :: A
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE SSCAL_F95
    PURE SUBROUTINE DSCAL_F95(X,A)
        ! Fortran77 call:
        ! DSCAL(N,A,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN) :: A
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE DSCAL_F95
    PURE SUBROUTINE CSCAL_F95(X,A)
        ! Fortran77 call:
        ! CSCAL(N,A,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(IN) :: A
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE CSCAL_F95
    PURE SUBROUTINE ZSCAL_F95(X,A)
        ! Fortran77 call:
        ! ZSCAL(N,A,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(IN) :: A
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE ZSCAL_F95
    PURE SUBROUTINE CSSCAL_F95(X,A)
        ! Fortran77 call:
        ! CSSCAL(N,A,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN) :: A
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE CSSCAL_F95
    PURE SUBROUTINE ZDSCAL_F95(X,A)
        ! Fortran77 call:
        ! ZDSCAL(N,A,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN) :: A
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE ZDSCAL_F95
END INTERFACE SCAL

INTERFACE SWAP
    PURE SUBROUTINE SSWAP_F95(X,Y)
        ! Fortran77 call:
        ! SSWAP(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(INOUT) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SSWAP_F95
    PURE SUBROUTINE DSWAP_F95(X,Y)
        ! Fortran77 call:
        ! DSWAP(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(INOUT) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DSWAP_F95
    PURE SUBROUTINE CSWAP_F95(X,Y)
        ! Fortran77 call:
        ! CSWAP(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(INOUT) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE CSWAP_F95
    PURE SUBROUTINE ZSWAP_F95(X,Y)
        ! Fortran77 call:
        ! ZSWAP(N,X,INCX,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(INOUT) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE ZSWAP_F95
END INTERFACE SWAP

INTERFACE IAMAX
    PURE FUNCTION ISAMAX_F95(X)
        ! Fortran77 call:
        ! ISAMAX(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        INTEGER :: ISAMAX_F95
        REAL(WP), INTENT(IN) :: X(:)
    END FUNCTION ISAMAX_F95
    PURE FUNCTION IDAMAX_F95(X)
        ! Fortran77 call:
        ! IDAMAX(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        INTEGER :: IDAMAX_F95
        REAL(WP), INTENT(IN) :: X(:)
    END FUNCTION IDAMAX_F95
    PURE FUNCTION ICAMAX_F95(X)
        ! Fortran77 call:
        ! ICAMAX(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        INTEGER :: ICAMAX_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
    END FUNCTION ICAMAX_F95
    PURE FUNCTION IZAMAX_F95(X)
        ! Fortran77 call:
        ! IZAMAX(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        INTEGER :: IZAMAX_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
    END FUNCTION IZAMAX_F95
END INTERFACE IAMAX

INTERFACE IAMIN
    PURE FUNCTION ISAMIN_F95(X)
        ! Fortran77 call:
        ! ISAMIN(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        INTEGER :: ISAMIN_F95
        REAL(WP), INTENT(IN) :: X(:)
    END FUNCTION ISAMIN_F95
    PURE FUNCTION IDAMIN_F95(X)
        ! Fortran77 call:
        ! IDAMIN(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        INTEGER :: IDAMIN_F95
        REAL(WP), INTENT(IN) :: X(:)
    END FUNCTION IDAMIN_F95
    PURE FUNCTION ICAMIN_F95(X)
        ! Fortran77 call:
        ! ICAMIN(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        INTEGER :: ICAMIN_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
    END FUNCTION ICAMIN_F95
    PURE FUNCTION IZAMIN_F95(X)
        ! Fortran77 call:
        ! IZAMIN(N,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        INTEGER :: IZAMIN_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
    END FUNCTION IZAMIN_F95
END INTERFACE IAMIN

INTERFACE CABS1
    PURE FUNCTION SCABS1_F95(C)
        ! Fortran77 call:
        ! SCABS1(C)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP) :: SCABS1_F95
        COMPLEX(WP), INTENT(IN) :: C
    END FUNCTION SCABS1_F95
    PURE FUNCTION DCABS1_F95(Z)
        ! Fortran77 call:
        ! DCABS1(Z)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP) :: DCABS1_F95
        COMPLEX(WP), INTENT(IN) :: Z
    END FUNCTION DCABS1_F95
END INTERFACE CABS1

INTERFACE GBMV
        ! TRANS='N','C','T'; default: 'N'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE SGBMV_F95(A,X,Y,KL,M,ALPHA,BETA,TRANS)
        ! Fortran77 call:
        ! SGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        INTEGER, INTENT(IN), OPTIONAL :: KL
        INTEGER, INTENT(IN), OPTIONAL :: M
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SGBMV_F95
    PURE SUBROUTINE DGBMV_F95(A,X,Y,KL,M,ALPHA,BETA,TRANS)
        ! Fortran77 call:
        ! DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        INTEGER, INTENT(IN), OPTIONAL :: KL
        INTEGER, INTENT(IN), OPTIONAL :: M
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DGBMV_F95
    PURE SUBROUTINE CGBMV_F95(A,X,Y,KL,M,ALPHA,BETA,TRANS)
        ! Fortran77 call:
        ! CGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        INTEGER, INTENT(IN), OPTIONAL :: KL
        INTEGER, INTENT(IN), OPTIONAL :: M
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE CGBMV_F95
    PURE SUBROUTINE ZGBMV_F95(A,X,Y,KL,M,ALPHA,BETA,TRANS)
        ! Fortran77 call:
        ! ZGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        INTEGER, INTENT(IN), OPTIONAL :: KL
        INTEGER, INTENT(IN), OPTIONAL :: M
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE ZGBMV_F95
END INTERFACE GBMV

INTERFACE GEMV
        ! TRANS='N','C','T'; default: 'N'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE SGEMV_F95(A,X,Y,ALPHA,BETA,TRANS)
        ! Fortran77 call:
        ! SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SGEMV_F95
    PURE SUBROUTINE DGEMV_F95(A,X,Y,ALPHA,BETA,TRANS)
        ! Fortran77 call:
        ! DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DGEMV_F95
    PURE SUBROUTINE CGEMV_F95(A,X,Y,ALPHA,BETA,TRANS)
        ! Fortran77 call:
        ! CGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE CGEMV_F95
    PURE SUBROUTINE ZGEMV_F95(A,X,Y,ALPHA,BETA,TRANS)
        ! Fortran77 call:
        ! ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE ZGEMV_F95
    PURE SUBROUTINE SCGEMV_F95(A,X,Y,ALPHA,BETA,TRANS)
        ! Fortran77 call:
        ! SCGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        REAL(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SCGEMV_F95
    PURE SUBROUTINE DZGEMV_F95(A,X,Y,ALPHA,BETA,TRANS)
        ! Fortran77 call:
        ! DZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        REAL(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DZGEMV_F95
END INTERFACE GEMV

INTERFACE GER
        ! Default ALPHA=1
    PURE SUBROUTINE SGER_F95(A,X,Y,ALPHA)
        ! Fortran77 call:
        ! SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(INOUT) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE SGER_F95
    PURE SUBROUTINE DGER_F95(A,X,Y,ALPHA)
        ! Fortran77 call:
        ! DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(INOUT) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE DGER_F95
END INTERFACE GER

INTERFACE GERC
        ! Default ALPHA=1
    PURE SUBROUTINE CGERC_F95(A,X,Y,ALPHA)
        ! Fortran77 call:
        ! CGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(INOUT) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE CGERC_F95
    PURE SUBROUTINE ZGERC_F95(A,X,Y,ALPHA)
        ! Fortran77 call:
        ! ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(INOUT) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE ZGERC_F95
END INTERFACE GERC

INTERFACE GERU
        ! Default ALPHA=1
    PURE SUBROUTINE CGERU_F95(A,X,Y,ALPHA)
        ! Fortran77 call:
        ! CGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(INOUT) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE CGERU_F95
    PURE SUBROUTINE ZGERU_F95(A,X,Y,ALPHA)
        ! Fortran77 call:
        ! ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(INOUT) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE ZGERU_F95
END INTERFACE GERU

INTERFACE HBMV
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE CHBMV_F95(A,X,Y,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! CHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE CHBMV_F95
    PURE SUBROUTINE ZHBMV_F95(A,X,Y,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! ZHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE ZHBMV_F95
END INTERFACE HBMV

INTERFACE HEMV
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE CHEMV_F95(A,X,Y,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! CHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE CHEMV_F95
    PURE SUBROUTINE ZHEMV_F95(A,X,Y,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! ZHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE ZHEMV_F95
END INTERFACE HEMV

INTERFACE HER
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
    PURE SUBROUTINE CHER_F95(A,X,UPLO,ALPHA)
        ! Fortran77 call:
        ! CHER(UPLO,N,ALPHA,X,INCX,A,LDA)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(INOUT) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
    END SUBROUTINE CHER_F95
    PURE SUBROUTINE ZHER_F95(A,X,UPLO,ALPHA)
        ! Fortran77 call:
        ! ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(INOUT) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
    END SUBROUTINE ZHER_F95
END INTERFACE HER

INTERFACE HER2
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
    PURE SUBROUTINE CHER2_F95(A,X,Y,UPLO,ALPHA)
        ! Fortran77 call:
        ! CHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(INOUT) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE CHER2_F95
    PURE SUBROUTINE ZHER2_F95(A,X,Y,UPLO,ALPHA)
        ! Fortran77 call:
        ! ZHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(INOUT) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE ZHER2_F95
END INTERFACE HER2

INTERFACE HPMV
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE CHPMV_F95(AP,X,Y,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! CHPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: AP(:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE CHPMV_F95
    PURE SUBROUTINE ZHPMV_F95(AP,X,Y,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! ZHPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: AP(:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE ZHPMV_F95
END INTERFACE HPMV

INTERFACE HPR
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
    PURE SUBROUTINE CHPR_F95(AP,X,UPLO,ALPHA)
        ! Fortran77 call:
        ! CHPR(UPLO,N,ALPHA,X,INCX,AP)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(INOUT) :: AP(:)
        COMPLEX(WP), INTENT(IN) :: X(:)
    END SUBROUTINE CHPR_F95
    PURE SUBROUTINE ZHPR_F95(AP,X,UPLO,ALPHA)
        ! Fortran77 call:
        ! ZHPR(UPLO,N,ALPHA,X,INCX,AP)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(INOUT) :: AP(:)
        COMPLEX(WP), INTENT(IN) :: X(:)
    END SUBROUTINE ZHPR_F95
END INTERFACE HPR

INTERFACE HPR2
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
    PURE SUBROUTINE CHPR2_F95(AP,X,Y,UPLO,ALPHA)
        ! Fortran77 call:
        ! CHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(INOUT) :: AP(:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE CHPR2_F95
    PURE SUBROUTINE ZHPR2_F95(AP,X,Y,UPLO,ALPHA)
        ! Fortran77 call:
        ! ZHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(INOUT) :: AP(:)
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE ZHPR2_F95
END INTERFACE HPR2

INTERFACE SBMV
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE SSBMV_F95(A,X,Y,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! SSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SSBMV_F95
    PURE SUBROUTINE DSBMV_F95(A,X,Y,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! DSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DSBMV_F95
END INTERFACE SBMV

INTERFACE SPMV
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE SSPMV_F95(AP,X,Y,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! SSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: AP(:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SSPMV_F95
    PURE SUBROUTINE DSPMV_F95(AP,X,Y,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: AP(:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DSPMV_F95
END INTERFACE SPMV

INTERFACE SPR
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
    PURE SUBROUTINE SSPR_F95(AP,X,UPLO,ALPHA)
        ! Fortran77 call:
        ! SSPR(UPLO,N,ALPHA,X,INCX,AP)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(INOUT) :: AP(:)
        REAL(WP), INTENT(IN) :: X(:)
    END SUBROUTINE SSPR_F95
    PURE SUBROUTINE DSPR_F95(AP,X,UPLO,ALPHA)
        ! Fortran77 call:
        ! DSPR(UPLO,N,ALPHA,X,INCX,AP)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(INOUT) :: AP(:)
        REAL(WP), INTENT(IN) :: X(:)
    END SUBROUTINE DSPR_F95
END INTERFACE SPR

INTERFACE SPR2
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
    PURE SUBROUTINE SSPR2_F95(AP,X,Y,UPLO,ALPHA)
        ! Fortran77 call:
        ! SSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(INOUT) :: AP(:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE SSPR2_F95
    PURE SUBROUTINE DSPR2_F95(AP,X,Y,UPLO,ALPHA)
        ! Fortran77 call:
        ! DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(INOUT) :: AP(:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE DSPR2_F95
END INTERFACE SPR2

INTERFACE SYMV
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE SSYMV_F95(A,X,Y,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! SSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SSYMV_F95
    PURE SUBROUTINE DSYMV_F95(A,X,Y,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DSYMV_F95
END INTERFACE SYMV

INTERFACE SYR
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
    PURE SUBROUTINE SSYR_F95(A,X,UPLO,ALPHA)
        ! Fortran77 call:
        ! SSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(INOUT) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
    END SUBROUTINE SSYR_F95
    PURE SUBROUTINE DSYR_F95(A,X,UPLO,ALPHA)
        ! Fortran77 call:
        ! DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(INOUT) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
    END SUBROUTINE DSYR_F95
END INTERFACE SYR

INTERFACE SYR2
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
    PURE SUBROUTINE SSYR2_F95(A,X,Y,UPLO,ALPHA)
        ! Fortran77 call:
        ! SSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(INOUT) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE SSYR2_F95
    PURE SUBROUTINE DSYR2_F95(A,X,Y,UPLO,ALPHA)
        ! Fortran77 call:
        ! DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(INOUT) :: A(:,:)
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE DSYR2_F95
END INTERFACE SYR2

INTERFACE TBMV
        ! UPLO='U','L'; default: 'U'
        ! TRANS='N','C','T'; default: 'N'
        ! DIAG='N','U'; default: 'N'
    PURE SUBROUTINE STBMV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! STBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE STBMV_F95
    PURE SUBROUTINE DTBMV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! DTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE DTBMV_F95
    PURE SUBROUTINE CTBMV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! CTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE CTBMV_F95
    PURE SUBROUTINE ZTBMV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE ZTBMV_F95
END INTERFACE TBMV

INTERFACE TBSV
        ! UPLO='U','L'; default: 'U'
        ! TRANS='N','C','T'; default: 'N'
        ! DIAG='N','U'; default: 'N'
    PURE SUBROUTINE STBSV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! STBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE STBSV_F95
    PURE SUBROUTINE DTBSV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE DTBSV_F95
    PURE SUBROUTINE CTBSV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! CTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE CTBSV_F95
    PURE SUBROUTINE ZTBSV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! ZTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE ZTBSV_F95
END INTERFACE TBSV

INTERFACE TPMV
        ! UPLO='U','L'; default: 'U'
        ! TRANS='N','C','T'; default: 'N'
        ! DIAG='N','U'; default: 'N'
    PURE SUBROUTINE STPMV_F95(AP,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! STPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN) :: AP(:)
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE STPMV_F95
    PURE SUBROUTINE DTPMV_F95(AP,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! DTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN) :: AP(:)
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE DTPMV_F95
    PURE SUBROUTINE CTPMV_F95(AP,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! CTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN) :: AP(:)
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE CTPMV_F95
    PURE SUBROUTINE ZTPMV_F95(AP,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! ZTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN) :: AP(:)
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE ZTPMV_F95
END INTERFACE TPMV

INTERFACE TPSV
        ! UPLO='U','L'; default: 'U'
        ! TRANS='N','C','T'; default: 'N'
        ! DIAG='N','U'; default: 'N'
    PURE SUBROUTINE STPSV_F95(AP,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! STPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN) :: AP(:)
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE STPSV_F95
    PURE SUBROUTINE DTPSV_F95(AP,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN) :: AP(:)
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE DTPSV_F95
    PURE SUBROUTINE CTPSV_F95(AP,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! CTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN) :: AP(:)
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE CTPSV_F95
    PURE SUBROUTINE ZTPSV_F95(AP,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! ZTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN) :: AP(:)
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE ZTPSV_F95
END INTERFACE TPSV

INTERFACE TRMV
        ! UPLO='U','L'; default: 'U'
        ! TRANS='N','C','T'; default: 'N'
        ! DIAG='N','U'; default: 'N'
    PURE SUBROUTINE STRMV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE STRMV_F95
    PURE SUBROUTINE DTRMV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE DTRMV_F95
    PURE SUBROUTINE CTRMV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! CTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE CTRMV_F95
    PURE SUBROUTINE ZTRMV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE ZTRMV_F95
END INTERFACE TRMV

INTERFACE TRSV
        ! UPLO='U','L'; default: 'U'
        ! TRANS='N','C','T'; default: 'N'
        ! DIAG='N','U'; default: 'N'
    PURE SUBROUTINE STRSV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! STRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE STRSV_F95
    PURE SUBROUTINE DTRSV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE DTRSV_F95
    PURE SUBROUTINE CTRSV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! CTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE CTRSV_F95
    PURE SUBROUTINE ZTRSV_F95(A,X,UPLO,TRANS,DIAG)
        ! Fortran77 call:
        ! ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: X(:)
    END SUBROUTINE ZTRSV_F95
END INTERFACE TRSV

INTERFACE GEMM
        ! TRANSA='N','C','T'; default: 'N'
        ! TRANSB='N','C','T'; default: 'N'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE SGEMM_F95(A,B,C,TRANSA,TRANSB,ALPHA,BETA)
        ! Fortran77 call:
        ! SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSB
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: B(:,:)
        REAL(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE SGEMM_F95
    PURE SUBROUTINE DGEMM_F95(A,B,C,TRANSA,TRANSB,ALPHA,BETA)
        ! Fortran77 call:
        ! DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSB
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: B(:,:)
        REAL(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE DGEMM_F95
    PURE SUBROUTINE CGEMM_F95(A,B,C,TRANSA,TRANSB,ALPHA,BETA)
        ! Fortran77 call:
        ! CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSB
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE CGEMM_F95
    PURE SUBROUTINE ZGEMM_F95(A,B,C,TRANSA,TRANSB,ALPHA,BETA)
        ! Fortran77 call:
        ! ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSB
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE ZGEMM_F95
    PURE SUBROUTINE SCGEMM_F95(A,B,C,TRANSA,TRANSB,ALPHA,BETA)
        ! Fortran77 call:
        ! SCGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSB
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE SCGEMM_F95
    PURE SUBROUTINE DZGEMM_F95(A,B,C,TRANSA,TRANSB,ALPHA,BETA)
        ! Fortran77 call:
        ! DZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSB
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE DZGEMM_F95
END INTERFACE GEMM

INTERFACE HEMM
        ! SIDE='L','R'; default: 'L'
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE CHEMM_F95(A,B,C,SIDE,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE CHEMM_F95
    PURE SUBROUTINE ZHEMM_F95(A,B,C,SIDE,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE ZHEMM_F95
END INTERFACE HEMM

INTERFACE HERK
        ! UPLO='U','L'; default: 'U'
        ! TRANS='N','C'; default: 'N'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE CHERK_F95(A,C,UPLO,TRANS,ALPHA,BETA)
        ! Fortran77 call:
        ! CHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE CHERK_F95
    PURE SUBROUTINE ZHERK_F95(A,C,UPLO,TRANS,ALPHA,BETA)
        ! Fortran77 call:
        ! ZHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE ZHERK_F95
END INTERFACE HERK

INTERFACE HER2K
        ! UPLO='U','L'; default: 'U'
        ! TRANS='N','C'; default: 'N'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE CHER2K_F95(A,B,C,UPLO,TRANS,ALPHA,BETA)
        ! Fortran77 call:
        ! CHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE CHER2K_F95
    PURE SUBROUTINE ZHER2K_F95(A,B,C,UPLO,TRANS,ALPHA,BETA)
        ! Fortran77 call:
        ! ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE ZHER2K_F95
END INTERFACE HER2K

INTERFACE SYMM
        ! SIDE='L','R'; default: 'L'
        ! UPLO='U','L'; default: 'U'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE SSYMM_F95(A,B,C,SIDE,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! SSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: B(:,:)
        REAL(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE SSYMM_F95
    PURE SUBROUTINE DSYMM_F95(A,B,C,SIDE,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: B(:,:)
        REAL(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE DSYMM_F95
    PURE SUBROUTINE CSYMM_F95(A,B,C,SIDE,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! CSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE CSYMM_F95
    PURE SUBROUTINE ZSYMM_F95(A,B,C,SIDE,UPLO,ALPHA,BETA)
        ! Fortran77 call:
        ! ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE ZSYMM_F95
END INTERFACE SYMM

INTERFACE SYRK
        ! UPLO='U','L'; default: 'U'
        ! TRANS='N','C','T'; default: 'N'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE SSYRK_F95(A,C,UPLO,TRANS,ALPHA,BETA)
        ! Fortran77 call:
        ! SSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE SSYRK_F95
    PURE SUBROUTINE DSYRK_F95(A,C,UPLO,TRANS,ALPHA,BETA)
        ! Fortran77 call:
        ! DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE DSYRK_F95
    PURE SUBROUTINE CSYRK_F95(A,C,UPLO,TRANS,ALPHA,BETA)
        ! Fortran77 call:
        ! CSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE CSYRK_F95
    PURE SUBROUTINE ZSYRK_F95(A,C,UPLO,TRANS,ALPHA,BETA)
        ! Fortran77 call:
        ! ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE ZSYRK_F95
END INTERFACE SYRK

INTERFACE SYR2K
        ! UPLO='U','L'; default: 'U'
        ! TRANS='N','C','T'; default: 'N'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE SSYR2K_F95(A,B,C,UPLO,TRANS,ALPHA,BETA)
        ! Fortran77 call:
        ! SSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: B(:,:)
        REAL(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE SSYR2K_F95
    PURE SUBROUTINE DSYR2K_F95(A,B,C,UPLO,TRANS,ALPHA,BETA)
        ! Fortran77 call:
        ! DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: B(:,:)
        REAL(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE DSYR2K_F95
    PURE SUBROUTINE CSYR2K_F95(A,B,C,UPLO,TRANS,ALPHA,BETA)
        ! Fortran77 call:
        ! CSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE CSYR2K_F95
    PURE SUBROUTINE ZSYR2K_F95(A,B,C,UPLO,TRANS,ALPHA,BETA)
        ! Fortran77 call:
        ! ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE ZSYR2K_F95
END INTERFACE SYR2K

INTERFACE TRMM
        ! SIDE='L','R'; default: 'L'
        ! UPLO='U','L'; default: 'U'
        ! TRANSA='N','C','T'; default: 'N'
        ! DIAG='N','U'; default: 'N'
        ! Default ALPHA=1
    PURE SUBROUTINE STRMM_F95(A,B,SIDE,UPLO,TRANSA,DIAG,ALPHA)
        ! Fortran77 call:
        ! STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: B(:,:)
    END SUBROUTINE STRMM_F95
    PURE SUBROUTINE DTRMM_F95(A,B,SIDE,UPLO,TRANSA,DIAG,ALPHA)
        ! Fortran77 call:
        ! DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: B(:,:)
    END SUBROUTINE DTRMM_F95
    PURE SUBROUTINE CTRMM_F95(A,B,SIDE,UPLO,TRANSA,DIAG,ALPHA)
        ! Fortran77 call:
        ! CTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: B(:,:)
    END SUBROUTINE CTRMM_F95
    PURE SUBROUTINE ZTRMM_F95(A,B,SIDE,UPLO,TRANSA,DIAG,ALPHA)
        ! Fortran77 call:
        ! ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: B(:,:)
    END SUBROUTINE ZTRMM_F95
END INTERFACE TRMM

INTERFACE TRSM
        ! SIDE='L','R'; default: 'L'
        ! UPLO='U','L'; default: 'U'
        ! TRANSA='N','C','T'; default: 'N'
        ! DIAG='N','U'; default: 'N'
        ! Default ALPHA=1
    PURE SUBROUTINE STRSM_F95(A,B,SIDE,UPLO,TRANSA,DIAG,ALPHA)
        ! Fortran77 call:
        ! STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: B(:,:)
    END SUBROUTINE STRSM_F95
    PURE SUBROUTINE DTRSM_F95(A,B,SIDE,UPLO,TRANSA,DIAG,ALPHA)
        ! Fortran77 call:
        ! DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(INOUT) :: B(:,:)
    END SUBROUTINE DTRSM_F95
    PURE SUBROUTINE CTRSM_F95(A,B,SIDE,UPLO,TRANSA,DIAG,ALPHA)
        ! Fortran77 call:
        ! CTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: B(:,:)
    END SUBROUTINE CTRSM_F95
    PURE SUBROUTINE ZTRSM_F95(A,B,SIDE,UPLO,TRANSA,DIAG,ALPHA)
        ! Fortran77 call:
        ! ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: SIDE
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: DIAG
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(INOUT) :: B(:,:)
    END SUBROUTINE ZTRSM_F95
END INTERFACE TRSM

INTERFACE GEMMT
        ! UPLO='U','L'; default: 'U'
        ! TRANSA='N','C','T'; default: 'N'
        ! TRANSB='N','C','T'; default: 'N'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE SGEMMT_F95(A,B,C,UPLO,TRANSA,TRANSB,ALPHA,BETA)
        ! Fortran77 call:
        ! SGEMMT(UPLO,TRANSA,TRANSB,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSB
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: B(:,:)
        REAL(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE SGEMMT_F95
    PURE SUBROUTINE DGEMMT_F95(A,B,C,UPLO,TRANSA,TRANSB,ALPHA,BETA)
        ! Fortran77 call:
        ! DGEMMT(UPLO,TRANSA,TRANSB,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSB
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: B(:,:)
        REAL(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE DGEMMT_F95
    PURE SUBROUTINE CGEMMT_F95(A,B,C,UPLO,TRANSA,TRANSB,ALPHA,BETA)
        ! Fortran77 call:
        ! CGEMMT(UPLO,TRANSA,TRANSB,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSB
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE CGEMMT_F95
    PURE SUBROUTINE ZGEMMT_F95(A,B,C,UPLO,TRANSA,TRANSB,ALPHA,BETA)
        ! Fortran77 call:
        ! ZGEMMT(UPLO,TRANSA,TRANSB,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSB
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE ZGEMMT_F95
END INTERFACE GEMMT

INTERFACE AXPYI
        ! Default A=1
    PURE SUBROUTINE SAXPYI_F95(X,INDX,Y,A)
        ! Fortran77 call:
        ! SAXPYI(NZ,A,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN), OPTIONAL :: A
        REAL(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SAXPYI_F95
    PURE SUBROUTINE DAXPYI_F95(X,INDX,Y,A)
        ! Fortran77 call:
        ! DAXPYI(NZ,A,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN), OPTIONAL :: A
        REAL(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DAXPYI_F95
    PURE SUBROUTINE CAXPYI_F95(X,INDX,Y,A)
        ! Fortran77 call:
        ! CAXPYI(NZ,A,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: A
        COMPLEX(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE CAXPYI_F95
    PURE SUBROUTINE ZAXPYI_F95(X,INDX,Y,A)
        ! Fortran77 call:
        ! ZAXPYI(NZ,A,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: A
        COMPLEX(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE ZAXPYI_F95
END INTERFACE AXPYI

INTERFACE DOTI
    PURE FUNCTION SDOTI_F95(X,INDX,Y)
        ! Fortran77 call:
        ! SDOTI(NZ,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP) :: SDOTI_F95
        REAL(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END FUNCTION SDOTI_F95
    PURE FUNCTION DDOTI_F95(X,INDX,Y)
        ! Fortran77 call:
        ! DDOTI(NZ,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP) :: DDOTI_F95
        REAL(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END FUNCTION DDOTI_F95
END INTERFACE DOTI

INTERFACE DOTCI
    PURE FUNCTION CDOTCI_F95(X,INDX,Y)
        ! Fortran77 call:
        ! CDOTCI(NZ,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP) :: CDOTCI_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END FUNCTION CDOTCI_F95
    PURE FUNCTION ZDOTCI_F95(X,INDX,Y)
        ! Fortran77 call:
        ! ZDOTCI(NZ,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP) :: ZDOTCI_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END FUNCTION ZDOTCI_F95
END INTERFACE DOTCI

INTERFACE DOTUI
    PURE FUNCTION CDOTUI_F95(X,INDX,Y)
        ! Fortran77 call:
        ! CDOTUI(NZ,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP) :: CDOTUI_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END FUNCTION CDOTUI_F95
    PURE FUNCTION ZDOTUI_F95(X,INDX,Y)
        ! Fortran77 call:
        ! ZDOTUI(NZ,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP) :: ZDOTUI_F95
        COMPLEX(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END FUNCTION ZDOTUI_F95
END INTERFACE DOTUI

INTERFACE GTHR
    PURE SUBROUTINE SGTHR_F95(X,INDX,Y)
        ! Fortran77 call:
        ! SGTHR(NZ,Y,X,INDX)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(OUT) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE SGTHR_F95
    PURE SUBROUTINE DGTHR_F95(X,INDX,Y)
        ! Fortran77 call:
        ! DGTHR(NZ,Y,X,INDX)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(OUT) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE DGTHR_F95
    PURE SUBROUTINE CGTHR_F95(X,INDX,Y)
        ! Fortran77 call:
        ! CGTHR(NZ,Y,X,INDX)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(OUT) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE CGTHR_F95
    PURE SUBROUTINE ZGTHR_F95(X,INDX,Y)
        ! Fortran77 call:
        ! ZGTHR(NZ,Y,X,INDX)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(OUT) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        COMPLEX(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE ZGTHR_F95
END INTERFACE GTHR

INTERFACE GTHRZ
    PURE SUBROUTINE SGTHRZ_F95(X,INDX,Y)
        ! Fortran77 call:
        ! SGTHRZ(NZ,Y,X,INDX)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(OUT) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SGTHRZ_F95
    PURE SUBROUTINE DGTHRZ_F95(X,INDX,Y)
        ! Fortran77 call:
        ! DGTHRZ(NZ,Y,X,INDX)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(OUT) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DGTHRZ_F95
    PURE SUBROUTINE CGTHRZ_F95(X,INDX,Y)
        ! Fortran77 call:
        ! CGTHRZ(NZ,Y,X,INDX)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(OUT) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE CGTHRZ_F95
    PURE SUBROUTINE ZGTHRZ_F95(X,INDX,Y)
        ! Fortran77 call:
        ! ZGTHRZ(NZ,Y,X,INDX)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(OUT) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE ZGTHRZ_F95
END INTERFACE GTHRZ

INTERFACE ROTI
        ! Default C=1
        ! Default S=1
    PURE SUBROUTINE SROTI_F95(X,INDX,Y,C,S)
        ! Fortran77 call:
        ! SROTI(NZ,X,INDX,Y,C,S)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN) :: C
        REAL(WP), INTENT(IN) :: S
        REAL(WP), INTENT(INOUT) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE SROTI_F95
    PURE SUBROUTINE DROTI_F95(X,INDX,Y,C,S)
        ! Fortran77 call:
        ! DROTI(NZ,X,INDX,Y,C,S)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN) :: C
        REAL(WP), INTENT(IN) :: S
        REAL(WP), INTENT(INOUT) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        REAL(WP), INTENT(IN) :: Y(:)
    END SUBROUTINE DROTI_F95
END INTERFACE ROTI

INTERFACE SCTR
    PURE SUBROUTINE SSCTR_F95(X,INDX,Y)
        ! Fortran77 call:
        ! SSCTR(NZ,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        REAL(WP), INTENT(OUT) :: Y(:)
    END SUBROUTINE SSCTR_F95
    PURE SUBROUTINE DSCTR_F95(X,INDX,Y)
        ! Fortran77 call:
        ! DSCTR(NZ,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        REAL(WP), INTENT(OUT) :: Y(:)
    END SUBROUTINE DSCTR_F95
    PURE SUBROUTINE CSCTR_F95(X,INDX,Y)
        ! Fortran77 call:
        ! CSCTR(NZ,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        COMPLEX(WP), INTENT(OUT) :: Y(:)
    END SUBROUTINE CSCTR_F95
    PURE SUBROUTINE ZSCTR_F95(X,INDX,Y)
        ! Fortran77 call:
        ! ZSCTR(NZ,X,INDX,Y)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(IN) :: X(:)
        INTEGER, INTENT(IN) :: INDX(:)
        COMPLEX(WP), INTENT(OUT) :: Y(:)
    END SUBROUTINE ZSCTR_F95
END INTERFACE SCTR

INTERFACE GEMM3M
        ! TRANSA='N','C','T'; default: 'N'
        ! TRANSB='N','C','T'; default: 'N'
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE CGEMM3M_F95(A,B,C,TRANSA,TRANSB,ALPHA,BETA)
        ! Fortran77 call:
        ! CGEMM3M(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => SP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSB
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE CGEMM3M_F95
    PURE SUBROUTINE ZGEMM3M_F95(A,B,C,TRANSA,TRANSB,ALPHA,BETA)
        ! Fortran77 call:
        ! ZGEMM3M(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        USE F95_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSA
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANSB
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: B(:,:)
        COMPLEX(WP), INTENT(INOUT) :: C(:,:)
    END SUBROUTINE ZGEMM3M_F95
END INTERFACE GEMM3M

INTERFACE AXPBY
        ! Default ALPHA=1
        ! Default BETA=1
    PURE SUBROUTINE SAXPBY_F95(X,Y,ALPHA,BETA)
        ! Fortran77 call:
        ! SAXPBY(N,ALPHA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE SAXPBY_F95
    PURE SUBROUTINE DAXPBY_F95(X,Y,ALPHA,BETA)
        ! Fortran77 call:
        ! DAXPBY(N,ALPHA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: X(:)
        REAL(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE DAXPBY_F95
    PURE SUBROUTINE CAXPBY_F95(X,Y,ALPHA,BETA)
        ! Fortran77 call:
        ! CAXPBY(N,ALPHA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE CAXPBY_F95
    PURE SUBROUTINE ZAXPBY_F95(X,Y,ALPHA,BETA)
        ! Fortran77 call:
        ! ZAXPBY(N,ALPHA,X,INCX,BETA,Y,INCY)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: X(:)
        COMPLEX(WP), INTENT(INOUT) :: Y(:)
    END SUBROUTINE ZAXPBY_F95
END INTERFACE AXPBY

INTERFACE GEM2V
        ! Default ALPHA=1
        ! Default BETA=0
    PURE SUBROUTINE SGEM2VU_F95(A,X1,X2,Y1,Y2,ALPHA,BETA)
        ! Fortran77 call:
        ! SGEM2VU(M,N,ALPHA,A,LDA,X1,INCX1,X2,INCX2,BETA,Y1,INCY1,Y2,
        !   INCY2)
        USE F95_PRECISION, ONLY: WP => SP
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: X1(:)
        REAL(WP), INTENT(IN) :: X2(:)
        REAL(WP), INTENT(INOUT) :: Y1(:)
        REAL(WP), INTENT(INOUT) :: Y2(:)
    END SUBROUTINE SGEM2VU_F95
    PURE SUBROUTINE DGEM2VU_F95(A,X1,X2,Y1,Y2,ALPHA,BETA)
        ! Fortran77 call:
        ! DGEM2VU(M,N,ALPHA,A,LDA,X1,INCX1,X2,INCX2,BETA,Y1,INCY1,Y2,
        !   INCY2)
        USE F95_PRECISION, ONLY: WP => DP
        REAL(WP), INTENT(IN), OPTIONAL :: ALPHA
        REAL(WP), INTENT(IN), OPTIONAL :: BETA
        REAL(WP), INTENT(IN) :: A(:,:)
        REAL(WP), INTENT(IN) :: X1(:)
        REAL(WP), INTENT(IN) :: X2(:)
        REAL(WP), INTENT(INOUT) :: Y1(:)
        REAL(WP), INTENT(INOUT) :: Y2(:)
    END SUBROUTINE DGEM2VU_F95
    PURE SUBROUTINE CGEM2VC_F95(A,X1,X2,Y1,Y2,ALPHA,BETA)
        ! Fortran77 call:
        ! CGEM2VC(M,N,ALPHA,A,LDA,X1,INCX1,X2,INCX2,BETA,Y1,INCY1,Y2,
        !   INCY2)
        USE F95_PRECISION, ONLY: WP => SP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X1(:)
        COMPLEX(WP), INTENT(IN) :: X2(:)
        COMPLEX(WP), INTENT(INOUT) :: Y1(:)
        COMPLEX(WP), INTENT(INOUT) :: Y2(:)
    END SUBROUTINE CGEM2VC_F95
    PURE SUBROUTINE ZGEM2VC_F95(A,X1,X2,Y1,Y2,ALPHA,BETA)
        ! Fortran77 call:
        ! ZGEM2VC(M,N,ALPHA,A,LDA,X1,INCX1,X2,INCX2,BETA,Y1,INCY1,Y2,
        !   INCY2)
        USE F95_PRECISION, ONLY: WP => DP
        COMPLEX(WP), INTENT(IN), OPTIONAL :: ALPHA
        COMPLEX(WP), INTENT(IN), OPTIONAL :: BETA
        COMPLEX(WP), INTENT(IN) :: A(:,:)
        COMPLEX(WP), INTENT(IN) :: X1(:)
        COMPLEX(WP), INTENT(IN) :: X2(:)
        COMPLEX(WP), INTENT(INOUT) :: Y1(:)
        COMPLEX(WP), INTENT(INOUT) :: Y2(:)
    END SUBROUTINE ZGEM2VC_F95
END INTERFACE GEM2V

INTERFACE SGEMM_BATCH
        ! TRANSA_ARRAY=Array where each element is one of 'N','C','T'; default: 'N'
        ! TRANSB_ARRAY=Array where each element is one of 'N','C','T'; default: 'N'
        ! ALPHA_ARRAY=Array of alpha values; default: 1
        ! BETA_ARRAY=Array of beta values; default: 0
    PURE SUBROUTINE SGEMM_BATCH_F95(A_ARRAY,B_ARRAY,C_ARRAY,M_ARRAY,N_ARRAY,K_ARRAY, &
         GROUP_SIZE,TRANSA_ARRAY,TRANSB_ARRAY,ALPHA_ARRAY,BETA_ARRAY)
        ! Fortran77 call:
        ! SGEMM_BATCH(TRANSA_ARRAY,TRANSB_ARRAY,M_ARRAY,N_ARRAY,K_ARRAY,ALPHA_ARRAY,A_ARRAY,LDA_ARRAY,B_ARRAY,LDB_ARRAY,BETA_ARRAY,C_ARRAY,LDC_ARRAY,GROUP_COUNT,GROUP_SIZE)
        USE F95_PRECISION, ONLY: WP => SP
        USE, INTRINSIC :: ISO_C_BINDING
        ! TRANSA_ARRAY: INOUT intent instead of IN because PURE.
        CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: TRANSA_ARRAY(:)
        ! TRANSB_ARRAY: INOUT intent instead of IN because PURE.
        CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: TRANSB_ARRAY(:)
        ! ALPHA_ARRAY: INOUT intent instead of IN because PURE.
        REAL(WP), INTENT(INOUT), OPTIONAL :: ALPHA_ARRAY(:)
        ! BETA_ARRAY: INOUT intent instead of IN because PURE.
        REAL(WP), INTENT(INOUT), OPTIONAL :: BETA_ARRAY(:)
        INTEGER, INTENT(IN) :: M_ARRAY(:)
        INTEGER, INTENT(IN) :: N_ARRAY(:)
        INTEGER, INTENT(IN) :: K_ARRAY(:)
        INTEGER, INTENT(IN) :: GROUP_SIZE(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(IN) :: A_ARRAY(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(IN) :: B_ARRAY(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(INOUT) :: C_ARRAY(:)
    END SUBROUTINE SGEMM_BATCH_F95
END INTERFACE SGEMM_BATCH

INTERFACE DGEMM_BATCH
        ! TRANSA_ARRAY=Array where each element is one of 'N','C','T'; default: 'N'
        ! TRANSB_ARRAY=Array where each element is one of 'N','C','T'; default: 'N'
        ! ALPHA_ARRAY=Array of alpha values; default: 1
        ! BETA_ARRAY=Array of beta values; default: 0
    PURE SUBROUTINE DGEMM_BATCH_F95(A_ARRAY,B_ARRAY,C_ARRAY,M_ARRAY,N_ARRAY,K_ARRAY, &
         GROUP_SIZE,TRANSA_ARRAY,TRANSB_ARRAY,ALPHA_ARRAY,BETA_ARRAY)
        ! Fortran77 call:
        ! DGEMM_BATCH(TRANSA_ARRAY,TRANSB_ARRAY,M_ARRAY,N_ARRAY,K_ARRAY,ALPHA_ARRAY,A_ARRAY,LDA_ARRAY,B_ARRAY,LDB_ARRAY,BETA_ARRAY,C_ARRAY,LDC_ARRAY,GROUP_COUNT,GROUP_SIZE)
        USE F95_PRECISION, ONLY: WP => DP
        USE, INTRINSIC :: ISO_C_BINDING
        ! TRANSA_ARRAY: INOUT intent instead of IN because PURE.
        CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: TRANSA_ARRAY(:)
        ! TRANSB_ARRAY: INOUT intent instead of IN because PURE.
        CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: TRANSB_ARRAY(:)
        ! ALPHA_ARRAY: INOUT intent instead of IN because PURE.
        REAL(WP), INTENT(INOUT), OPTIONAL :: ALPHA_ARRAY(:)
        ! BETA_ARRAY: INOUT intent instead of IN because PURE.
        REAL(WP), INTENT(INOUT), OPTIONAL :: BETA_ARRAY(:)
        INTEGER, INTENT(IN) :: M_ARRAY(:)
        INTEGER, INTENT(IN) :: N_ARRAY(:)
        INTEGER, INTENT(IN) :: K_ARRAY(:)
        INTEGER, INTENT(IN) :: GROUP_SIZE(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(IN) :: A_ARRAY(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(IN) :: B_ARRAY(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(INOUT) :: C_ARRAY(:)
    END SUBROUTINE DGEMM_BATCH_F95
END INTERFACE DGEMM_BATCH

INTERFACE CGEMM_BATCH
        ! TRANSA_ARRAY=Array where each element is one of 'N','C','T'; default: 'N'
        ! TRANSB_ARRAY=Array where each element is one of 'N','C','T'; default: 'N'
        ! ALPHA_ARRAY=Array of alpha values; default: 1
        ! BETA_ARRAY=Array of beta values; default: 0
    PURE SUBROUTINE CGEMM_BATCH_F95(A_ARRAY,B_ARRAY,C_ARRAY,M_ARRAY,N_ARRAY,K_ARRAY, &
         GROUP_SIZE,TRANSA_ARRAY,TRANSB_ARRAY,ALPHA_ARRAY,BETA_ARRAY)
        ! Fortran77 call:
        ! CGEMM_BATCH(TRANSA_ARRAY,TRANSB_ARRAY,M_ARRAY,N_ARRAY,K_ARRAY,ALPHA_ARRAY,A_ARRAY,LDA_ARRAY,B_ARRAY,LDB_ARRAY,BETA_ARRAY,C_ARRAY,LDC_ARRAY,GROUP_COUNT,GROUP_SIZE)
        USE F95_PRECISION, ONLY: WP => SP
        USE, INTRINSIC :: ISO_C_BINDING
        ! TRANSA_ARRAY: INOUT intent instead of IN because PURE.
        CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: TRANSA_ARRAY(:)
        ! TRANSB_ARRAY: INOUT intent instead of IN because PURE.
        CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: TRANSB_ARRAY(:)
        ! ALPHA_ARRAY: INOUT intent instead of IN because PURE.
        COMPLEX(WP), INTENT(INOUT), OPTIONAL :: ALPHA_ARRAY(:)
        ! BETA_ARRAY: INOUT intent instead of IN because PURE.
        COMPLEX(WP), INTENT(INOUT), OPTIONAL :: BETA_ARRAY(:)
        INTEGER, INTENT(IN) :: M_ARRAY(:)
        INTEGER, INTENT(IN) :: N_ARRAY(:)
        INTEGER, INTENT(IN) :: K_ARRAY(:)
        INTEGER, INTENT(IN) :: GROUP_SIZE(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(IN) :: A_ARRAY(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(IN) :: B_ARRAY(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(INOUT) :: C_ARRAY(:)
    END SUBROUTINE CGEMM_BATCH_F95
END INTERFACE CGEMM_BATCH

INTERFACE ZGEMM_BATCH
        ! TRANSA_ARRAY=Array where each element is one of 'N','C','T'; default: 'N'
        ! TRANSB_ARRAY=Array where each element is one of 'N','C','T'; default: 'N'
        ! ALPHA_ARRAY=Array of alpha values; default: 1
        ! BETA_ARRAY=Array of beta values; default: 0
    PURE SUBROUTINE ZGEMM_BATCH_F95(A_ARRAY,B_ARRAY,C_ARRAY,M_ARRAY,N_ARRAY,K_ARRAY, &
         GROUP_SIZE,TRANSA_ARRAY,TRANSB_ARRAY,ALPHA_ARRAY,BETA_ARRAY)
        ! Fortran77 call:
        ! ZGEMM_BATCH(TRANSA_ARRAY,TRANSB_ARRAY,M_ARRAY,N_ARRAY,K_ARRAY,ALPHA_ARRAY,A_ARRAY,LDA_ARRAY,B_ARRAY,LDB_ARRAY,BETA_ARRAY,C_ARRAY,LDC_ARRAY,GROUP_COUNT,GROUP_SIZE)
        USE F95_PRECISION, ONLY: WP => DP
        USE, INTRINSIC :: ISO_C_BINDING
        ! TRANSA_ARRAY: INOUT intent instead of IN because PURE.
        CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: TRANSA_ARRAY(:)
        ! TRANSB_ARRAY: INOUT intent instead of IN because PURE.
        CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: TRANSB_ARRAY(:)
        ! ALPHA_ARRAY: INOUT intent instead of IN because PURE.
        COMPLEX(WP), INTENT(INOUT), OPTIONAL :: ALPHA_ARRAY(:)
        ! BETA_ARRAY: INOUT intent instead of IN because PURE.
        COMPLEX(WP), INTENT(INOUT), OPTIONAL :: BETA_ARRAY(:)
        INTEGER, INTENT(IN) :: M_ARRAY(:)
        INTEGER, INTENT(IN) :: N_ARRAY(:)
        INTEGER, INTENT(IN) :: K_ARRAY(:)
        INTEGER, INTENT(IN) :: GROUP_SIZE(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(IN) :: A_ARRAY(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(IN) :: B_ARRAY(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(INOUT) :: C_ARRAY(:)
    END SUBROUTINE ZGEMM_BATCH_F95
END INTERFACE ZGEMM_BATCH

INTERFACE CGEMM3M_BATCH
        ! TRANSA_ARRAY=Array where each element is one of 'N','C','T'; default: 'N'
        ! TRANSB_ARRAY=Array where each element is one of 'N','C','T'; default: 'N'
        ! ALPHA_ARRAY=Array of alpha values; default: 1
        ! BETA_ARRAY=Array of beta values; default: 0
    PURE SUBROUTINE CGEMM3M_BATCH_F95(A_ARRAY,B_ARRAY,C_ARRAY,M_ARRAY,N_ARRAY,K_ARRAY, &
         GROUP_SIZE,TRANSA_ARRAY,TRANSB_ARRAY,ALPHA_ARRAY,BETA_ARRAY)
        ! Fortran77 call:
        ! CGEMM3M_BATCH(TRANSA_ARRAY,TRANSB_ARRAY,M_ARRAY,N_ARRAY,K_ARRAY,ALPHA_ARRAY,A_ARRAY,LDA_ARRAY,B_ARRAY,LDB_ARRAY,BETA_ARRAY,C_ARRAY,LDC_ARRAY,GROUP_COUNT,GROUP_SIZE)
        USE F95_PRECISION, ONLY: WP => SP
        USE, INTRINSIC :: ISO_C_BINDING
        ! TRANSA_ARRAY: INOUT intent instead of IN because PURE.
        CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: TRANSA_ARRAY(:)
        ! TRANSB_ARRAY: INOUT intent instead of IN because PURE.
        CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: TRANSB_ARRAY(:)
        ! ALPHA_ARRAY: INOUT intent instead of IN because PURE.
        COMPLEX(WP), INTENT(INOUT), OPTIONAL :: ALPHA_ARRAY(:)
        ! BETA_ARRAY: INOUT intent instead of IN because PURE.
        COMPLEX(WP), INTENT(INOUT), OPTIONAL :: BETA_ARRAY(:)
        INTEGER, INTENT(IN) :: M_ARRAY(:)
        INTEGER, INTENT(IN) :: N_ARRAY(:)
        INTEGER, INTENT(IN) :: K_ARRAY(:)
        INTEGER, INTENT(IN) :: GROUP_SIZE(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(IN) :: A_ARRAY(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(IN) :: B_ARRAY(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(INOUT) :: C_ARRAY(:)
    END SUBROUTINE CGEMM3M_BATCH_F95
END INTERFACE CGEMM3M_BATCH

INTERFACE ZGEMM3M_BATCH
        ! TRANSA_ARRAY=Array where each element is one of 'N','C','T'; default: 'N'
        ! TRANSB_ARRAY=Array where each element is one of 'N','C','T'; default: 'N'
        ! ALPHA_ARRAY=Array of alpha values; default: 1
        ! BETA_ARRAY=Array of beta values; default: 0
    PURE SUBROUTINE ZGEMM3M_BATCH_F95(A_ARRAY,B_ARRAY,C_ARRAY,M_ARRAY,N_ARRAY,K_ARRAY, &
         GROUP_SIZE,TRANSA_ARRAY,TRANSB_ARRAY,ALPHA_ARRAY,BETA_ARRAY)
        ! Fortran77 call:
        ! ZGEMM3M_BATCH(TRANSA_ARRAY,TRANSB_ARRAY,M_ARRAY,N_ARRAY,K_ARRAY,ALPHA_ARRAY,A_ARRAY,LDA_ARRAY,B_ARRAY,LDB_ARRAY,BETA_ARRAY,C_ARRAY,LDC_ARRAY,GROUP_COUNT,GROUP_SIZE)
        USE F95_PRECISION, ONLY: WP => DP
        USE, INTRINSIC :: ISO_C_BINDING
        ! TRANSA_ARRAY: INOUT intent instead of IN because PURE.
        CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: TRANSA_ARRAY(:)
        ! TRANSB_ARRAY: INOUT intent instead of IN because PURE.
        CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: TRANSB_ARRAY(:)
        ! ALPHA_ARRAY: INOUT intent instead of IN because PURE.
        COMPLEX(WP), INTENT(INOUT), OPTIONAL :: ALPHA_ARRAY(:)
        ! BETA_ARRAY: INOUT intent instead of IN because PURE.
        COMPLEX(WP), INTENT(INOUT), OPTIONAL :: BETA_ARRAY(:)
        INTEGER, INTENT(IN) :: M_ARRAY(:)
        INTEGER, INTENT(IN) :: N_ARRAY(:)
        INTEGER, INTENT(IN) :: K_ARRAY(:)
        INTEGER, INTENT(IN) :: GROUP_SIZE(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(IN) :: A_ARRAY(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(IN) :: B_ARRAY(:)
        INTEGER(KIND=C_INTPTR_T), INTENT(INOUT) :: C_ARRAY(:)
    END SUBROUTINE ZGEMM3M_BATCH_F95
END INTERFACE ZGEMM3M_BATCH

END MODULE BLAS95
