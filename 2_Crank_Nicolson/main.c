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
    
    psP[ 0 ] = psC[ 0 ] / psB[ 0 ];
    psQ[ 0 ] = 1.0 / psB[ 0 ];
    psR[ 0 ] = 0.0;
    
    for ( i = 1 ; i < nN - 1 ; i++ ) {
        sFrac = 1.0 / ( psB[ i ] - psA[ i ] * psP[ i - 1 ] );
        psP[ i ] = psC[ i ] * sFrac;
        psQ[ i ] = sFrac;
        psR[ i ] = -psA[ i ] * sFrac;
    }
    
    i = nN - 1; // although it already has this value...
    sFrac = 1.0 / ( psB[ i ] - psA[ i ] * psP[ i - 1 ] );
    psQ[ i ] = sFrac;
    psR[ i ] = -psA[ i ] * sFrac;
    
    return 0;
}


int myCalcSolTridiagonal(scalar *psD, scalar *psP, scalar *psQ, scalar *psR, int nN, scalar *psDP, scalar *psX) {
    int i;
    
    // Calculating the array of d'
    
    psDP[ 0 ] = psD[ 0 ] * psQ[ 0 ];
    
    for ( i = 1 ; i < nN ; i++ ) {
        psDP[ i ] = psD[ i ] * psQ[ i ] + psDP[ i - 1 ] * psR[ i ];
    }
    
    // Calculating the solution
    
    psX[ nN - 1 ] = psDP[ nN - 1 ];
    
    for ( i = nN - 2 ; i >= 0 ; i-- ) {
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
    
    scalar sR;
    scalar sDiag, sOffDiag;
    
    sR = 0.5 * sAlpha * sDT / sDX / sDX;
    sDiag = 1 + 2.0 * sR;
    sOffDiag = -sR;
    
    for ( i = 0 ; i < nN ; i++ ) {
        psA[ i ] = sOffDiag;
        psB[ i ] = sDiag;
        psC[ i ] = sOffDiag;
    }
    
    psA[ 0 ] = 0.0;
    
    return 0;
}


int myCreateVecInput(scalar sDX, scalar sDT, scalar sAlpha, 
                     scalar *psU, int nN, scalar *psD)
{
    int i;
    
    scalar sR;
    scalar sDiag, sOffDiag;
    
    sR = 0.5 * sAlpha * sDT / sDX / sDX;
    sDiag = 1 - 2.0 * sR;
    sOffDiag = sR;
    
    psD[ 0 ] = sDiag * psU[ 0 ] + sOffDiag * psU[ 1 ];
    
    for ( i = 1 ; i < nN - 1 ; i++ ) {
        psD[ i ] = sOffDiag * psU[ i - 1 ] + sDiag * psU[ i ] + sOffDiag * psU[ i + 1 ];
    }
    
    psD[ nN - 1 ] = sOffDiag * psU[ 1 ] + sDiag * psU[ nN - 1 ];
    
    return 0;
}


///////////////////////////////////////////////////////////////////////
// End of Crank-Nicolson method
///////////////////////////////////////////////////////////////////////


double myFuncInit(double dX) {
    return sin(MY_PI * dX);
}


int main() {
    int i;
    
    double pdA[ 512 ], pdB[ 512 ], pdC[ 512 ];
    double pdP[ 512 ], pdQ[ 512 ], pdR[ 512 ];
    double pdU[ 512 ], pdD[ 512 ], pdDP[ 512 ];
    
    double dXMax, dTMax;
    double dIntX, dIntT;
    int nNX,nNT;
    
    double dAlpha;
    double dC;
    
    // Initialization: setting some variables
    
    dXMax = 1.0;
    dTMax = 0.5;
    nNX = 20;
    nNT = 3000;
    dIntX = dXMax / nNX;
    dIntT = dTMax / nNT;
    
    dAlpha = 1.0;
    
    dC = dAlpha * dAlpha * dIntT / ( dIntX * dIntX );
    
    // Initialization: making the initial state
    
    pdU[ nNX - 1 ] = 0.0;
    
    for ( i = 0 ; i < nNX - 1 ; i++ ) {
        pdU[ i ] = myFuncInit( ( i + 1 ) * dIntX );
    }
    
    // Calculation of the numerical method
    
    myCreateMtxTridiagonal(dIntX, dIntT, dAlpha, nNX, pdA, pdB, pdC);
    myCalcMtxTridiagonal(pdA, pdB, pdC, nNX - 1, pdP, pdQ, pdR);
    
    for ( i = 0 ; i < nNT ; i++ ) {
        myCreateVecInput(dIntX, dIntT, dAlpha, pdU, nNX - 1, pdD);
        myCalcSolTridiagonal(pdD, pdP, pdQ, pdR, nNX - 1, pdDP, pdU);
    }
    
    // Output
    
    for ( i = 0 ; i < nNX - 1 ; i++ ) {
        printf("%3i %11.8f %13.8f\n", i + 1, dIntX * ( i + 1 ), pdU[ i ]);
    }
    
    return 0;
}


