#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>


#define MY_PI   3.14159265359


#define MIDX(N, I, J) ( ( N ) * ( I ) + ( J ) )


typedef double complex scalar;
//typedef double scalar;


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
// Beginning of Newton-Raphson method
///////////////////////////////////////////////////////////////////////


// Swapping I-th row and J-th row
int mySwapRow(int nN, scalar *psA, int nI, int nJ) {
  int i;
  
  scalar sDummy;
  int nIdxI, nIdxJ;
  
  for ( i = 0 ; i < nN ; i++ ) {
    nIdxI = MIDX(nN, nI, i);
    nIdxJ = MIDX(nN, nJ, i);
    
    sDummy = psA[ nIdxI ];
    psA[ nIdxI ] = psA[ nIdxJ ];
    psA[ nIdxJ ] = sDummy;
  }
  
  return 0;
}


// Multiplying a scalar on a row
int myMulScalarRow(int nN, scalar *psA, scalar sC, int nI) {
  int i;
  
  for ( i = 0 ; i < nN ; i++ ) {
    psA[ MIDX(nN, nI, i) ] *= sC;
  }
  
  return 0;
}


// Adding C * J-th row into the I-th row
int myAddRowMul(int nN, scalar *psA, int nI, scalar sC, int nJ) {
  int i;
  
  for ( i = 0 ; i < nN ; i++ ) {
    psA[ MIDX(nN, nI, i) ] += sC * psA[ MIDX(nN, nJ, i) ];
  }
  
  return 0;
}


// Multiplying N x N matrix A^{-1} to X on the left by Gaussian elimination
int myCalcMulInvMtx(int nN, scalar *psA, scalar *psX) {
  int i, j, k;
  
  scalar sC;
  
  for ( i = 0 ; i < nN ; i++ ) {
    // Step1 : If the pivot is 0, then make it non-zero!
    if ( psA[ MIDX(nN, i, i) ] == 0 ) {
      for ( j = i + 1 ; j < nN ; j++ ) {
        if ( psA[ MIDX(nN, j, j) ] != 0 ) break;
      }
      
      // In this case we found that A is non-invertible...
      if ( j >= nN ) return -1;
      
      mySwapRow(nN, psA, i, j);
      
      //mySwapRow(nN, psX, i, j);
      sC = psX[ i ];
      psX[ i ] = psX[ j ];
      psX[ j ] = sC;
    }
    
    // Step2 : Making the pivot 1
    sC = 1.0 / psA[ MIDX(nN, i, i) ];
    
    myMulScalarRow(nN, psA, sC, i);
    //myMulScalarRow(nN, psX, sC, i);
    psX[ i ] *= sC;
    
    // Step3 : Eliminating all off-diagonal components on the current column
    for ( j = 0 ; j < nN ; j++ ) {
      if ( j == i ) continue;
      
      sC = -psA[ MIDX(nN, j, i) ];
      if ( sC == 0 ) continue;
      
      myAddRowMul(nN, psA, j, sC, i);
      //myAddRowMul(nN, psX, j, sC, i);
      psX[ j ] += sC * psX[ i ];
    }
  }
  
  return 0;
}


// One step for Newton-Raphson method
int myCalcOneStepNewtonRaphson(int nN, scalar *psU, scalar *psUN, 
                               scalar sR, scalar sLambda, scalar sDT, 
                               scalar *psGrad, scalar *psOut)
{
  int i, j;
  
  int nIdx;
  
  scalar sSumU, sSumUN;
  scalar sAbsU, sAbsUN;
  
  int nRes;
  
  // Calculating the gradient
  for ( i = 0 ; i < nN ; i++ ) {
    for ( j = 0 ; j < nN ; j++ ) {
      nIdx = MIDX(nN, i, j);
      
      if ( i == j ) {
        sAbsUN = cabs(psUN[ i ]);
        
        psGrad[ nIdx ] = 2.0 * I - 2.0 * sR;
        psGrad[ nIdx ] -= sDT * sLambda * sAbsUN * ( 2 * psUN[ i ] + sAbsUN );
      } else if ( i - j == 1 || i - j == -1 ) {
        psGrad[ nIdx ] = sR;
      } else {
        psGrad[ nIdx ] = 0;
      }
    }
  }
  
  //sSumU = 0;
  //for ( i = 0 ; i < nN ; i++ ) sSumU += psU[ i ];
  
  //sSumUN = 0;
  //for ( i = 0 ; i < nN ; i++ ) sSumUN += psUN[ i ];
  
  // Calculating -F(U^{n+1})
  for ( i = 0 ; i < nN ; i++ ) {
    psOut[ i ]  = ( 2.0 * I - 2.0 * sR ) * psUN[ i ];
    psOut[ i ] -= ( 2.0 * I + 2.0 * sR ) *  psU[ i ] - sR * sSumU;
    
    if ( i > 0 ) {
      psOut[ i ] +=  sR * psUN[ i - 1 ];
      psOut[ i ] -= -sR * psU[ i - 1 ];
    }
    
    if ( i < nN - 1 ) {
      psOut[ i ] +=  sR * psUN[ i + 1 ];
      psOut[ i ] -= -sR * psU[ i + 1 ];
    }
    
    sAbsU  = cabs(psU[ i ]);
    sAbsUN = cabs(psUN[ i ]);
    psOut[ i ] -= sDT * sLambda * ( sAbsUN * sAbsUN * psUN[ i ] + sAbsU * sAbsU * psU[ i ] );
  }
  
  // Multiplying the inverse of the gradient
  nRes = myCalcMulInvMtx(nN, psGrad, psOut);
  if ( nRes != 0 ) return nRes; // Do not forget to check the insanity
  
  // Finisher: Add the above 'difference' into the previous result
  for ( i = 0 ; i < nN ; i++ ) {
    psOut[ i ] = psUN[ i ] - psOut[ i ];
  }
  
  return 0;
}


int myCalcNewtonRaphsonMethod(int nN, scalar *psU, 
                              scalar sR, scalar sLambda, scalar sDT, 
                              scalar *psUN)
{
  int i;
  
  int nIter;
  
  scalar *psGrad;
  scalar *psUDummy1, *psUDummy2;
  
  scalar *psDummy;
  
  int nRes;
  
  psGrad = (scalar *)malloc(sizeof(scalar) * nN * nN);
  psUDummy1 = (scalar *)malloc(sizeof(scalar) * nN);
  psUDummy2 = (scalar *)malloc(sizeof(scalar) * nN);
  
  for ( i = 0 ; i < nN ; i++ ) {
    psUDummy1[ i ] = psU[ i ];
  }
  
  for ( nIter = 0 ; nIter < 20 ; nIter++ ) {
    nRes = myCalcOneStepNewtonRaphson(nN, psU, psUDummy1, sR, sLambda, sDT, psGrad, psUDummy2);
    if ( nRes != 0 ) return -1;
    
    printf("  Iter %i\n", nIter);
    for ( i = 0 ; i < nN ; i++ ) {
      printf("%13.6f\t%13.6f\n", cabs(psUDummy1[ i ]), cabs(psUDummy2[ i ]));
    }
    
    for ( i = 0 ; i < nN ; i++ ) { 
      double sDiff = cabs(psUDummy1[ i ] - psUDummy2[ i ]);
      if ( sDiff > 0.00001 * cabs(psUDummy1[ i ]) && sDiff > 0.00001 * cabs(psUDummy2[ i ]) ) break;
    }
    
    if ( i >= nN ) break;
    
    psDummy = psUDummy1;
    psUDummy1 = psUDummy2;
    psUDummy2 = psDummy;
  }
  
  //if ( nIter >= 1000 ) return -2;
  
  for ( i = 0 ; i < nN ; i++ ) {
    psUN[ i ] = psUDummy1[ i ];
  }
  
  free(psGrad);
  free(psUDummy1);
  free(psUDummy2);
  
  return 0;
}


///////////////////////////////////////////////////////////////////////
// End of Newton-Raphson method
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


// An exact solution of Schroedinger equation : (1 / sqrt(1 + 4at)) exp{-x^2/(1 + 4at)}, 
// where a = ihbar / 2m


double complex myExactSolution(double complex dLambda, double dX, double dT) {
  //return 1.0 / csqrt(1 + 4 * dAlpha * dT) * cexp(-dX * dX / ( 1 + 4 * dAlpha * dT ));
  double dV = 0.001;
  double dEta = 1.0;
  return sqrt(2 / -dLambda) * dEta / cosh(dX - dV * dT) * cexp(I * ( 0.5 * dV * (dX - 0.5 * dV * dT) + dEta * dEta * dT ));
}


double complex myFuncInit(double complex dLambda, double dQ) {
  //return myExactSolution(dAlpha, dQ / sqrt(1.0 - dQ * dQ), 0.0);
  return myExactSolution(dLambda, dQ, 0.0);
}


int main() {
  int i;
  
  double complex pdA[ 4096 ], pdB[ 4096 ], pdC[ 4096 ];
  double complex pdP[ 4096 ], pdQ[ 4096 ], pdR[ 4096 ];
  double complex pdU[ 4096 ], pdD[ 4096 ], pdDP[ 4096 ];
  
  double dXMax, dTMax;
  double dIntX, dIntT;
  int nNX, nNT;
  
  double complex dAlpha;
  double complex dLambda;
  
  // Initialization: setting some variables
  
  dXMax = 10.0;
  dTMax = 10.00;
  nNX = 200;
  nNT = 100;
  dIntX = dXMax / nNX;
  dIntT = dTMax / nNT;
  
  dAlpha = 0.01 * I;
  dLambda = -0.01;
  
  // Initialization: making the initial state
  
  pdU[ 0   ] = 0.0;
  pdU[ nNX ] = 0.0;
  
  //for ( i = 1 ; i <= nNX - 1 ; i++ ) {
  for ( i = 0 ; i <= nNX ; i++ ) {
    pdU[ i ] = myFuncInit(dLambda, 1.0 * i / nNX * 2.0 - 1.0);
  }
  
  // Calculation of the numerical method
  
  //myCreateMtxTridiagonal(dIntX, dIntT, dAlpha, nNX, pdA, pdB, pdC);
  //myCalcMtxTridiagonal(pdA, pdB, pdC, nNX, pdP, pdQ, pdR);
  
  for ( i = 0 ; i < nNT ; i++ ) {
    printf("Step %i...\n", i);
    if ( myCalcNewtonRaphsonMethod(nNX + 1, pdU, dIntT / ( dIntX * dIntX ), dLambda, dIntT, pdU) != 0 )
      printf("Error!\n");
    //myCreateVecInput(dIntX, dIntT, dAlpha, pdU, nNX, pdD);
    //myCalcSolTridiagonal(pdD, pdP, pdQ, pdR, nNX, pdDP, pdU);
    if ( i == 0 ) break;
  }
  
  // Output
  
  for ( i = 1 ; i <= nNX - 1 ; i++ ) {
    double dQ, dX;
    double complex cYCal, cYExt;
    
    dQ = 1.0 * i / nNX * 2.0 - 1.0;
    dX = dQ / sqrt(1.0 - dQ * dQ);
    
    cYCal = pdU[ i ];
    cYExt = myExactSolution(dLambda, dX, 0.0 + dTMax);
    
    printf("%3i %13.6f     %11.6f %11.6f (r: %11.6f)     %11.6f %11.6f (r: %11.6f)\n", 
            i, dX, 
            creal(cYCal), creal(cYExt), ( creal(cYCal) - creal(cYExt) ) / creal(cYExt), 
            cimag(cYCal), cimag(cYExt), ( cimag(cYCal) - cimag(cYExt) ) / cimag(cYExt));
  }
  
  return 0;
}


