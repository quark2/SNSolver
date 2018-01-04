#include <stdio.h>
#include <complex.h>
#include <math.h>


#define MY_PI   3.14159265359

#define MY_MAXB 4


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
                     scalar *psU, scalar sBCAdvL, scalar sBCAdvR, int nN, scalar *psD)
{
    int i;
    
    scalar sR;
    scalar sDiag, sOffDiag;
    
    sR = 0.5 * sAlpha * sDT / sDX / sDX;
    sDiag = 1 - 2.0 * sR;
    sOffDiag = sR;
    
    for ( i = 1 ; i <= nN - 1 ; i++ ) {
        psD[ i - 1 ] = sOffDiag * psU[ i - 1 ] + sDiag * psU[ i ] + sOffDiag * psU[ i + 1 ];
    }
    
    psD[ 0      ] += sOffDiag * sBCAdvL;
    psD[ nN - 2 ] += sOffDiag * sBCAdvR;
    
    return 0;
}


double g_dTTest = 0.0;


double myFuncExactSol(double dX, double dT);


int myCalcEffFreeBC(scalar sDX, scalar sDT, scalar sAlpha, scalar *psP, scalar *psQ, scalar *psR, 
                    scalar *psU, scalar *psD, scalar *psTPrev, scalar *psTRes, int nNX)
{
    int i;
    
    scalar arrsCoeffX[] = {6.0, -15.0, 20.0, -15.0, 6.0, -1.0};
    scalar arrsCoeffT[] = {3.0, -3.0, 1.0};
    //scalar arrsCoeff[] = {3.0, -3.0, 1.0};
    //scalar arrsCoeff[] = {5.0, -10.0, 10.0, -5.0, 1.0};
    //scalar arrsCoeff[] = {10.0, -45.0, 120.0, -210.0, 252.0, -210.0, 120.0, -45.0, 10.0, -1.0};
    
    int nNumLvX = (int)( sizeof(arrsCoeffX) / sizeof(arrsCoeffX[ 0 ]) );
    int nNumLvT = (int)( sizeof(arrsCoeffT) / sizeof(arrsCoeffT[ 0 ]) );
    
    scalar sDiv, sNum;
    scalar sValBound;
    scalar sSumViaX;
    
    for ( i = 0 ; i < nNX - 1 ; i++ ) {
        psD[ i ] = 0.0;
    }
    
    psD[ nNX - 2 ] = -0.5 * sAlpha * sDT / sDX / sDX;
    
    myCalcSolTridiagonal(psD, psP, psQ, psR, nNX - 1, psD, psD);
    //printf("<< %13.8f, %13.8f\n", psD[ nNX - 2 ], psD[ nNX - 3 ]);
    
    sDiv = 1.0;
    sNum = 0.0;
    
    for ( i = 0 ; i < nNumLvX ; i++ ) {
        sDiv += arrsCoeffX[ i ] * psD[ nNX - 2 - i ];
        sNum += arrsCoeffX[ i ] * psU[ nNX - 1 - i ];
        //printf("%13.8f, %13.8f, %13.8f\n", arrsCoeffX[ i ], psD[ nNX - 2 - i ], psU[ nNX - 1 - i ]);
    }
    
    sValBound = sNum / sDiv; // This is the u^{n+1}_{nX+1}, one of which we seek
    //printf("--- %13.8f, %13.8f, %13.8f, %13.8f\n", sNum, sDiv, sValBound, myFuncExactSol(5.0, g_dTTest));
    
    for ( i = 0 ; i < nNX - 1 ; i++ ) {
        psU[ i + 1 ] -= psD[ i ] * sValBound;
    }
    
    //psU[ nNX ] = sValBound;
    
    if ( psTPrev != NULL ) {
        sSumViaX = 0.0;
        
        for ( i = 0 ; i < nNumLvT ; i++ ) {
            sSumViaX += arrsCoeffT[ i ] * psU[ nNX - 1 - i ];
        }
        
        
    }
    
    return 0;
}


///////////////////////////////////////////////////////////////////////
// End of Crank-Nicolson method
///////////////////////////////////////////////////////////////////////


double g_dTInit = 1.5;
double g_dXMin = -5.0, g_dXMax = 5.0;
double g_dDiffusionConst = 1.0;


double myFuncExactSol(double dX, double dT) {
    return sqrt(g_dTInit / dT) * exp(-( dX * dX ) / ( 4.0 * g_dDiffusionConst * dT ));
}


double myFuncInit(double dX) {
    return myFuncExactSol(dX, g_dTInit);
}


double myBCLeft(double dT) {
    return myFuncExactSol(g_dXMin, dT);
}


double myBCRight(double dT) {
    return myFuncExactSol(g_dXMax, dT);
}


#define RP_MAXN 4096


int main() {
    int i;
    
    double pdA[ RP_MAXN ], pdB[ RP_MAXN ], pdC[ RP_MAXN ];
    double pdP[ RP_MAXN ], pdQ[ RP_MAXN ], pdR[ RP_MAXN ];
    double pdU[ RP_MAXN ], pdD[ RP_MAXN ], pdDP[ RP_MAXN ];
    
    double dXMin, dXMax, dXCurr;
    double dIntX;
    int nNT;
    
    double dTInit, dTFin, dTCurr;
    double dIntT;
    int nNX;
    
    double dAlpha;
    
    // Initialization: setting some variables for PDE itself
    
    dXMin = g_dXMin;
    dXMax = g_dXMax;
    nNX = 512;
    dIntX = ( dXMax - dXMin ) / nNX;
    
    dTInit = g_dTInit;
    dTFin  = 3.0;
    nNT = 32768;
    dIntT = ( dTFin - dTInit ) / nNT;
    
    dAlpha = g_dDiffusionConst;
    
    // Initialization: making the initial state
    
    for ( i = 0 ; i <= nNX ; i++ ) {
        dXCurr = dXMin + i * dIntX;
        pdU[ i ] = myFuncInit(dXCurr);
    }
    
    // Calculation of the numerical method
    
    myCreateMtxTridiagonal(dIntX, dIntT, dAlpha, nNX, pdA, pdB, pdC);
    myCalcMtxTridiagonal(pdA, pdB, pdC, nNX - 1, pdP, pdQ, pdR);
    
    for ( i = 0 ; i < nNT ; i++ ) {
        double dBCAdvL, dBCAdvR;
        
        // Setting BC
        dTCurr = dTInit + i * dIntT; // That's why I starts at 1, not 0 (and ends at nNT, not nNT - 1)
        g_dTTest = dTCurr + dIntT;
        
        dBCAdvL = myBCLeft(dTCurr + dIntT);
        //dBCAdvR = myBCRight(dTCurr + dIntT);
        dBCAdvR = 0; // There is no Adv BC in this time; it will be handled in CalcEffFreeBC()
        
        // Calculating the function by tridiagonal matrix algorithm
        myCreateVecInput(dIntX, dIntT, dAlpha, pdU, dBCAdvL, dBCAdvR, nNX, pdD);
        myCalcSolTridiagonal(pdD, pdP, pdQ, pdR, nNX - 1, pdDP, pdU + 1);
        //for ( int j = 0 ; i == 0 && j <= nNX ; j++ ) printf("%13.8f\n", pdU[ j ]);
        myCalcEffFreeBC(dIntX, dIntT, dAlpha, pdP, pdQ, pdR, pdU, pdD, nNX);
        
        pdU[ 0 ] = dBCAdvL;
        //pdU[ nNX ] = dBCAdvR; // It is already calculated
    }
    
    // Output
    dTFin = g_dTInit + i * dIntT;
    
    for ( i = 0 ; i <= nNX ; i++ ) {
        double dResNum;
        double dResExa;
        
        dXCurr = dXMin + i * dIntX;
        
        dResNum = pdU[ i ];
        dResExa = myFuncExactSol(dXCurr, dTFin);
        printf("%4i %11.8f %13.8f %13.8f %13.8f\n", i, dXCurr, dResNum, dResExa, ( dResExa - dResNum ) / dResExa);
    }
    
    return 0;
}


