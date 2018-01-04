#include <stdio.h>
#include <math.h>


#define MY_PI   3.14159265359


int myInitTridiagonal(double *pdA, double *pdB, double dCRe, double dCIm, int nNX) {
    int i;
    
    pdA[ 0 ] = 1.0 + dCRe;
    pdB[ 0 ] = -dCRe / ( 2.0 * pdA[ 0 ] );
    
    for ( i = 1 ; i < nNX - 2 ; i++) {
        pdA[ i ] = 1.0 + dCRe + dCRe * pdB[ i - 1 ] / 2.0;
        pdB[ i ] = -dCRe / ( 2.0 * pdA[ i ] );
    }  
    
    pdA[ nNX - 2 ] = 1.0 + dCRe + 0.5 * dCRe * pdB[ nNX - 3 ];
    
    return 0;
}


int myCalcNextTimeStep(double *pdU, double *pdBuf, double *pdA, double *pdB, 
                       double dCRe, double dCIm, int nNX)
{
    int i;
    
    pdBuf[ 0 ] = ( ( 1.0 - dCRe ) * pdU[ 0 ] + dCRe * pdU[ 1 ] / 2.0 ) / pdA[ 0 ];
    
    for ( i = 1 ; i < nNX - 1 ; i++ ) {
        pdBuf[ i ] = ( ( 1.0 - dCRe ) * pdU[ i ] + 0.5 * dCRe * 
            ( pdU[ i + 1 ] + pdU[ i - 1 ] + pdBuf[ i - 1 ] ) ) / pdA[ i ];
    }
    
    pdU[ nNX - 2 ] = pdBuf[ nNX - 2 ];
    
    for ( i = ( nNX - 2 ) - 1 ; i >= 0 ; i-- ) {
        pdU[ i ] = pdBuf[ i ] - pdB[ i ] * pdU[ i + 1 ];
    }
    
    return 0;
}


double myFuncInit(double dX) {
    return sin(MY_PI * dX);
}


int main() {
    int i;
    
    double pdU[ 512 ], pdA[ 512 ], pdB[ 512 ], pdBuf[ 512 ];
    
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
    
    myInitTridiagonal(pdA, pdB, dC, 0, nNX);
    
    for ( i = 0 ; i < nNT ; i++ ) {
        myCalcNextTimeStep(pdU, pdBuf, pdA, pdB, dC, 0, nNX);
    }
    
    // Output
    
    for ( i = 0 ; i < nNX - 1 ; i++ ) {
        printf("%3i %11.8f %13.8f\n", i + 1, dIntX * ( i + 1 ), pdU[ i ]);
    }
    
    return 0;
}


