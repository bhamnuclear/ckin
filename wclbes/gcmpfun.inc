*
* $Id: gcmpfun.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: gcmpfun.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
*
* gcmpfun.inc
*
#if (defined(CERNLIB_DOUBLE))&&(!defined(CERNLIB_F90))&&(defined(CERNLIB_LINUX)||!defined(CERNLIB_QF2C))
#include "def64.inc"
     +      GREAL,GIMAG,XARG,YARG
#include "defc64.inc"
     +      ZARG,GCONJG,GCMPLX
      GREAL( ZARG)=DREAL( ZARG)
      GIMAG( ZARG)=DIMAG( ZARG)
      GCONJG(ZARG)=DCONJG(ZARG)
      GCMPLX(XARG,YARG)=DCMPLX(XARG,YARG)

#elif (defined(CERNLIB_DOUBLE))&&(defined(CERNLIB_QF2C))&&(!defined(CERNLIB_LINUX))
#include "def64.inc"
     +      GREAL,GIMAG,XARG,YARG
#include "gen/defc64.inc"
     +      ZARG,GCONJG,GCMPLX

#elif (defined(CERNLIB_DOUBLE))&&(defined(CERNLIB_F90))
#include "def64.inc"
     +      GREAL,GIMAG,XARG,YARG
#include "gen/defc64.inc"
     +      ZARG,GCONJG,GCMPLX
      GREAL( ZARG)= REAL( ZARG,2)
      GIMAG( ZARG)=AIMAG( ZARG)
      GCONJG(ZARG)= CONJG(ZARG)
      GCMPLX(XARG,YARG)= CMPLX(XARG,YARG,2)

#elif 1
#include "def64.inc"
     +      GREAL,GIMAG,XARG,YARG
#include "defc64.inc"
     +      ZARG,GCONJG,GCMPLX
      GREAL( ZARG)=REAL( ZARG)
      GIMAG( ZARG)=AIMAG(ZARG)
      GCONJG(ZARG)=CONJG(ZARG)
      GCMPLX(XARG,YARG)= CMPLX(XARG,YARG)
#endif
