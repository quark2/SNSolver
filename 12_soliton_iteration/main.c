#include <stdio.h>
#include <complex.h>
#include <math.h>


#define MY_PI   3.14159265359


typedef double complex scalar;
//typedef double scalar;


double g_dLambda = 0.5;


///////////////////////////////////////////////////////////////////////
// Beginning of Thomas tridiagonal algorithm
///////////////////////////////////////////////////////////////////////


// Note : psA[ 0 ] should be occupied and the actual value of psA be from psA[ 1 ]
int myCalcMtxTridiagonal(scalar *psA, scalar *psB, scalar *psC, int nN, 
                         scalar *psP, scalar *psQ, scalar *psR)
{
  int i;
  
  scalar sFrac;
  
  psP[ 1 ] = psC[ 1 ] / psB[ 1 ];
  psQ[ 1 ] = 1.0 / psB[ 1 ];
  psR[ 1 ] = 0.0;
  
  for ( i = 2 ; i <= nN - 1 ; i++ ) {
    sFrac = 1.0 / ( psB[ i ] - psA[ i ] * psP[ i - 1 ] );
    psP[ i ] = psC[ i ] * sFrac;
    psQ[ i ] = sFrac;
    psR[ i ] = -psA[ i ] * sFrac;
  }
  
  /*i = nN - 1; // although it already has this value...
  sFrac = 1.0 / ( psB[ i ] - psA[ i ] * psP[ i - 1 ] );
  psQ[ i ] = sFrac;
  psR[ i ] = -psA[ i ] * sFrac;*/
  
  return 0;
}


int myCalcSolTridiagonal(scalar *psD, scalar *psP, scalar *psQ, scalar *psR, int nN, scalar *psX) {
  int i;
  
  scalar *psDP;
  
  // Calculating the array of d'
  
  psDP = psX;
  psDP[ 1 ] = psD[ 1 ] * psQ[ 1 ];
  
  for ( i = 2 ; i <= nN - 1 ; i++ ) {
    psDP[ i ] = psD[ i ] * psQ[ i ] + psDP[ i - 1 ] * psR[ i ];
  }
  
  // Calculating the solution
  
  psX[ nN - 1 ] = psDP[ nN - 1 ];
  
  for ( i = nN - 2 ; i >= 1 ; i-- ) {
    psX[ i ] = psDP[ i ] - psP[ i ] * psX[ i + 1 ];
  }
  
  return 0;
}


///////////////////////////////////////////////////////////////////////
// End of Thomas tridiagonal algorithm
///////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////
// Beginning of Crank-Nicolson method
///////////////////////////////////////////////////////////////////////


int myCreateMtxTridiagonal(scalar sDX, scalar sDT, scalar sAlpha, int nN, 
                           scalar *psV, scalar *psA, scalar *psB, scalar *psC)
{
  int i;
  
  scalar sQ, sP, sR;
  scalar sCoeffPart1, sCoeffPart2;
  
  for ( i = 1 ; i <= nN - 1 ; i++ ) {
    sQ = 1.0 * i / nN * 2.0 - 1.0;
    sP = 1 - sQ * sQ;
    sR = 0.50 * sDT;
    
    sCoeffPart1 = sAlpha * I * sR * sP * sP * sP / sDX / sDX;
    sCoeffPart2 = 1.5 * sAlpha * I * sR * sQ * sP * sP / sDX;
    
    psA[ i ] = -sCoeffPart1 - sCoeffPart2;
    psB[ i ] = 1 + 2.0 * sCoeffPart1 + I * sR * psV[ i ];
    psC[ i ] = -sCoeffPart1 + sCoeffPart2;
  }
  
  psA[ 0 ] = 0.0;
  
  return 0;
}


int myCreateVecInput(scalar sDX, scalar sDT, scalar sAlpha, 
                     scalar *psV, scalar *psU, int nN, scalar *psD)
{
  int i;
  
  //scalar sDiag, sOffDiag;
  scalar sDiag, sUpperDiag, sLowerDiag;
  
  scalar sQ, sP, sR;
  scalar sCoeffPart1, sCoeffPart2;
  
  // The boundary condition; at infinity point
  psU[ 0 ] = psU[ nN ] = 0.0;
  
  //psD[ 0 ] = sDiag * psU[ 0 ] + sOffDiag * psU[ 1 ];
  
  for ( i = 1 ; i <= nN - 1 ; i++ ) {
    sQ = 1.0 * i / nN * 2.0 - 1.0;
    sP = 1 - sQ * sQ;
    sR = 0.50 * sDT;
    
    sCoeffPart1 = sR * sAlpha * I * sP * sP * sP / sDX / sDX;
    sCoeffPart2 = 1.5 * sAlpha * I * sR * sQ * sP * sP / sDX;
    
    sUpperDiag = sCoeffPart1 + sCoeffPart2;
    sDiag = 1 - 2.0 * sCoeffPart1 - I * sR * psV[ i ];
    sLowerDiag = sCoeffPart1 - sCoeffPart2;
    
    psD[ i ] = sUpperDiag * psU[ i - 1 ] + sDiag * psU[ i ] + sLowerDiag * psU[ i + 1 ];
  }
  
  //psD[ nN - 1 ] = sOffDiag * psU[ 1 ] + sDiag * psU[ nN - 1 ];
  
  return 0;
}


///////////////////////////////////////////////////////////////////////
// End of Crank-Nicolson method
///////////////////////////////////////////////////////////////////////


// An exact solution of Schroedinger equation : (1 / sqrt(1 + 4at)) exp{-x^2/(1 + 4at)}, 
// where a = ihbar / 2m


double complex myExactSolution(double complex dD, double dX, double dT) {
  double dEta = 1.0;
  double dV = 0.5;
  //return 1.0 / csqrt(1 + 4 * dAlpha * I * dT) * cexp(-dX * dX / ( 1 + 4 * dAlpha * I * dT ));
  return ( sqrt(2 * dD) / g_dLambda ) * dEta / cosh(dEta * ( dX - dV * dT )) * 
    cexp(I / ( 2 * dD ) * dV * ( dX - 0.5 * dV * dT ) + I * dEta * dEta * dD * dT );
}


int myCalcPotential(int nNX, scalar *psU, scalar *psV) {
  int i;
  
  for ( i = 0 ; i <= nNX ; i++ ) {
    psV[ i ] = -g_dLambda * g_dLambda * cabs(psU[ i ]) * cabs(psU[ i ]);
  }
  
  return 0;
}


double complex myFuncInit(double complex dAlpha, double dQ) {
  return myExactSolution(dAlpha, dQ / sqrt(1.0 - dQ * dQ), 0.0);
}


int main() {
  int i, j;
  
  double complex arrdA[ 4096 ], arrdB[ 4096 ], arrdC[ 4096 ], arrdD[ 4096 ];
  double complex arrdP[ 4096 ], arrdQ[ 4096 ], arrdR[ 4096 ];
  double complex arrdU[ 4096 ], arrdV[ 4096 ];
  double complex arrdUN1[ 4096 ], arrdUN2[ 4096 ], arrdVN[ 4096 ];
  
  double complex *pdUN1, *pdUN2;
  double complex *pdDummySwap;
  
  double dXMax, dTMax;
  double dIntX, dIntT;
  int nNX, nNT;
  
  double complex dAlpha;
  
  int nLoop;
  
  // Initialization: setting some variables
  
  dXMax = 2.0;
  dTMax = 3.00;
  nNX = 200;
  nNT = 1000;
  dIntX = dXMax / nNX;
  dIntT = dTMax / nNT;
  
  dAlpha = 0.5 / 4.0 / 1.0 / 1.0;
  
  // Initialization: making the initial state
  
  arrdU[ 0   ] = 0.0;
  arrdU[ nNX ] = 0.0;
  
  for ( i = 1 ; i <= nNX - 1 ; i++ ) {
    arrdU[ i ] = myFuncInit(dAlpha, 1.0 * i / nNX * 2.0 - 1.0);
    arrdV[ i ] = -g_dLambda * g_dLambda * cabs(arrdU[ i ]) * cabs(arrdU[ i ]);
  }
  
  pdUN1 = arrdUN1;
  pdUN2 = arrdUN2;
  
  // Calculation of the numerical method
  
  for ( i = 0 ; i < nNT ; i++ ) {
    for ( j = 0 ; j <= nNX ; j++ ) {
      arrdVN[ j ] = arrdV[ j ];
    }
    
    myCreateMtxTridiagonal(dIntX, dIntT, dAlpha, nNX, arrdVN, arrdA, arrdB, arrdC);
    myCalcMtxTridiagonal(arrdA, arrdB, arrdC, nNX, arrdP, arrdQ, arrdR);
    myCreateVecInput(dIntX, dIntT, dAlpha, arrdV, arrdU, nNX, arrdD);
    myCalcSolTridiagonal(arrdD, arrdP, arrdQ, arrdR, nNX, pdUN1);
    pdDummySwap = pdUN1; pdUN1 = pdUN2; pdUN2 = pdDummySwap;
    
    nLoop = 0;
    while ( 1 ) {
      myCalcPotential(nNX, pdUN2, arrdVN);
      
      myCreateMtxTridiagonal(dIntX, dIntT, dAlpha, nNX, arrdVN, arrdA, arrdB, arrdC);
      myCalcMtxTridiagonal(arrdA, arrdB, arrdC, nNX, arrdP, arrdQ, arrdR);
      myCreateVecInput(dIntX, dIntT, dAlpha, arrdV, arrdU, nNX, arrdD);
      myCalcSolTridiagonal(arrdD, arrdP, arrdQ, arrdR, nNX, pdUN1);
      pdDummySwap = pdUN1; pdUN1 = pdUN2; pdUN2 = pdDummySwap;
      
      //for ( j = 0 ; j <= nNX ; j++ ) if ( cabs(pdUN2[ j ] - pdUN1[ j ]) > 0.000001 ) break;
      for ( j = 0 ; j <= nNX ; j++ ) {
        if ( creal(pdUN1[ j ]) != 0 ) {
          if ( fabs(( creal(pdUN2[ j ]) - creal(pdUN1[ j ]) ) / creal(pdUN1[ j ])) > 0.000001 ) break;
        } else {
          if ( fabs(creal(pdUN2[ j ]) - creal(pdUN1[ j ])) > 0.000001 ) break;
        }
        
        if ( cimag(pdUN1[ j ]) != 0 ) {
          if ( fabs(( cimag(pdUN2[ j ]) - cimag(pdUN1[ j ]) ) / cimag(pdUN1[ j ])) > 0.000001 ) break;
        } else {
          if ( fabs(cimag(pdUN2[ j ]) - cimag(pdUN1[ j ])) > 0.000001 ) break;
        }
      }
      
      if ( j > nNX ) break;
      
      if ( nLoop++ >= 100 ) break;
    }
    //printf("%i - %i\n", i, nLoop);
    
    for ( j = 0 ; j <= nNX ; j++ ) {
      arrdU[ j ] = pdUN2[ j ];
      arrdV[ j ] = arrdVN[ j ];
    }
  }
  
  // Output
  
  for ( i = 1 ; i <= nNX - 1 ; i++ ) {
    double dQ, dX;
    double complex cYCal, cYExt;
    
    dQ = 1.0 * i / nNX * 2.0 - 1.0;
    dX = dQ / sqrt(1.0 - dQ * dQ);
    
    cYCal = arrdU[ i ];
    cYExt = myExactSolution(dAlpha, dX, 0.0 + dTMax);
    
    printf("%3i %13.6f     %11.6f %11.6f   %11.6f     %11.6f %11.6f   %11.6f\n", 
            i, dX, 
            //creal(cYCal), creal(cYExt), ( creal(cYCal) - creal(cYExt) ) / creal(cYExt), 
            //cimag(cYCal), cimag(cYExt), ( cimag(cYCal) - cimag(cYExt) ) / cimag(cYExt));
            cabs(cYCal), cabs(cYExt), ( cabs(cYCal) - cabs(cYExt) ) / cabs(cYExt), 
            carg(cYCal), carg(cYExt), ( carg(cYCal) - carg(cYExt) ) / carg(cYExt));
  }
  
  return 0;
}


