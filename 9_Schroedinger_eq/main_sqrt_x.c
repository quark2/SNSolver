#include <stdio.h>
#include <complex.h>
#include <math.h>


#define MY_PI   3.14159265359


//typedef double complex scalar;
typedef double scalar;


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


int myCalcSolTridiagonal(scalar *psD, scalar *psP, scalar *psQ, scalar *psR, int nN, scalar *psDP, scalar *psX) {
  int i;
  
  // Calculating the array of d'
  
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
                           scalar *psA, scalar *psB, scalar *psC)
{
  int i;
  
  scalar sQ, sP;
  scalar sCoeffPart1, sCoeffPart2;
  
  for ( i = 1 ; i <= nN - 1 ; i++ ) {
    sQ = 1.0 * i / nN * 2.0 - 1.0;
    sP = 1 - sQ * sQ;
    
    sCoeffPart1 = 0.50 * sAlpha * sDT * sP * sP * sP / sDX / sDX;
    sCoeffPart2 = 0.75 * sAlpha * sDT * sQ * sP * sP / sDX;
    
    psA[ i ] = -sCoeffPart1 - sCoeffPart2;
    psB[ i ] = 1 + 2.0 * sCoeffPart1;
    psC[ i ] = -sCoeffPart1 + sCoeffPart2;
  }
  
  psA[ 0 ] = 0.0;
  
  return 0;
}


int myCreateVecInput(scalar sDX, scalar sDT, scalar sAlpha, 
                     scalar *psU, int nN, scalar *psD)
{
  int i;
  
  //scalar sDiag, sOffDiag;
  scalar sDiag, sUpperDiag, sLowerDiag;
  
  scalar sQ, sP;
  scalar sCoeffPart1, sCoeffPart2;
  
  // The boundary condition; at infinity point
  psU[ 0 ] = psU[ nN ] = 0.0;
  
  //psD[ 0 ] = sDiag * psU[ 0 ] + sOffDiag * psU[ 1 ];
  
  for ( i = 1 ; i <= nN - 1 ; i++ ) {
    sQ = 1.0 * i / nN * 2.0 - 1.0;
    sP = 1 - sQ * sQ;
    
    sCoeffPart1 = 0.50 * sAlpha * sDT * sP * sP * sP / sDX / sDX;
    sCoeffPart2 = 0.75 * sAlpha * sDT * sQ * sP * sP / sDX;
    
    sUpperDiag = sCoeffPart1 + sCoeffPart2;
    sDiag = 1 - 2.0 * sCoeffPart1;
    sLowerDiag = sCoeffPart1 - sCoeffPart2;
    
    psD[ i ] = sUpperDiag * psU[ i - 1 ] + sDiag * psU[ i ] + sLowerDiag * psU[ i + 1 ];
  }
  
  //psD[ nN - 1 ] = sOffDiag * psU[ 1 ] + sDiag * psU[ nN - 1 ];
  
  return 0;
}


///////////////////////////////////////////////////////////////////////
// End of Crank-Nicolson method
///////////////////////////////////////////////////////////////////////


// An exact solution of diffusion equation : t^{-1/2} exp{-x^2/4at}


double myExactSolution(double dAlpha, double dX, double dT) {
  return 1.0 / sqrt(dT) * exp(-dX * dX / ( 4 * dAlpha * dT ));
}


double myFuncInit(double dAlpha, double dQ) {
  return myExactSolution(dAlpha, dQ / sqrt(1.0 - dQ * dQ), 1.0);
}


int main() {
  int i;
  
  double pdA[ 4096 ], pdB[ 4096 ], pdC[ 4096 ];
  double pdP[ 4096 ], pdQ[ 4096 ], pdR[ 4096 ];
  double pdU[ 4096 ], pdD[ 4096 ], pdDP[ 4096 ];
  
  double dXMax, dTMax;
  double dIntX, dIntT;
  int nNX,nNT;
  
  double dAlpha;
  
  // Initialization: setting some variables
  
  dXMax = 2.0;
  dTMax = 50.00;
  nNX = 2000;
  nNT = 25000;
  dIntX = dXMax / nNX;
  dIntT = dTMax / nNT;
  
  dAlpha = 1.0;
  
  // Initialization: making the initial state
  
  pdU[ 0   ] = 0.0;
  pdU[ nNX ] = 0.0;
  
  for ( i = 1 ; i <= nNX - 1 ; i++ ) {
    pdU[ i ] = myFuncInit(dAlpha, 1.0 * i / nNX * 2.0 - 1.0);
  }
  
  // Calculation of the numerical method
  
  myCreateMtxTridiagonal(dIntX, dIntT, dAlpha, nNX, pdA, pdB, pdC);
  myCalcMtxTridiagonal(pdA, pdB, pdC, nNX, pdP, pdQ, pdR);
  
  for ( i = 0 ; i < nNT ; i++ ) {
    myCreateVecInput(dIntX, dIntT, dAlpha, pdU, nNX, pdD);
    myCalcSolTridiagonal(pdD, pdP, pdQ, pdR, nNX, pdDP, pdU);
  }
  
  // Output
  
  for ( i = 1 ; i <= nNX - 1 ; i++ ) {
    double dQ, dX;
    double dYCal, dYExt;
    
    dQ = 1.0 * i / nNX * 2.0 - 1.0;
    dX = dQ / sqrt(1.0 - dQ * dQ);
    
    dYCal = pdU[ i ];
    dYExt = myExactSolution(dAlpha, dX, 1.0 + dTMax);
    
    printf("%3i %15.8f %13.8f %13.8f %lf\n", i, dX, dYCal, dYExt, ( dYCal - dYExt ) / dYExt);
  }
  
  return 0;
}


