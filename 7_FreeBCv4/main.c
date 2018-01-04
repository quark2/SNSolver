#include <stdio.h>
#include <complex.h>
#include <math.h>


#define MY_PI   3.14159265359


//typedef double complex scalar;
typedef double scalar;


///////////////////////////////////////////////////////////////////////
// Beginning of corrected free boundary condition
///////////////////////////////////////////////////////////////////////


int myApplyFreeBCToMtx(int nN, scalar *psA, scalar *psB, scalar *psC, 
                       scalar *psRatioL, scalar *psRatioR)
{
    int nIdxL, nIdxR;
    
    nIdxL = 1;
    nIdxR = nN - 1;
    
    *psRatioL = psA[ nIdxL ] / psC[ nIdxL + 1 ];
    *psRatioR = psC[ nIdxR ] / psA[ nIdxR - 1 ];
    
    psB[ nIdxL ] +=  3 * psA[ nIdxL ] - ( *psRatioL ) * psA[ nIdxL + 1 ];
    psC[ nIdxL ] += -3 * psA[ nIdxL ] - ( *psRatioL ) * psB[ nIdxL + 1 ];
    
    psB[ nIdxR ] +=  3 * psC[ nIdxR ] - ( *psRatioR ) * psC[ nIdxR - 1 ];
    psA[ nIdxR ] += -3 * psC[ nIdxR ] - ( *psRatioR ) * psB[ nIdxR - 1 ];
    
    //printf("%lf, %lf, %lf, %lf, %lf, %lf\n", *psRatioL, *psRatioR, psB[ nIdxL ], psC[ nIdxL ], psB[ nIdxR ], psA[ nIdxR ]);
    
    return 0;
}


int myApplyFreeBCToVecInput(int nN, scalar sRatioL, scalar sRatioR, scalar *psY) {
    psY[ 1      ] -= sRatioL * psY[ 2 ];
    psY[ nN - 1 ] -= sRatioR * psY[ nN - 2 ];
    
    return 0;
}


int myGetBoundaryByFreeBC(int nN, scalar *psU) {
    int nIdxL, nIdxR;
    
    nIdxL = 0;
    nIdxR = nN;
    
    psU[ nIdxL ] = 3 * psU[ nIdxL + 1 ] - 3 * psU[ nIdxL + 2 ] * psU[ nIdxL + 3 ];
    psU[ nIdxR ] = 3 * psU[ nIdxR - 1 ] - 3 * psU[ nIdxR - 2 ] * psU[ nIdxR - 3 ];
    
    return 0;
}


///////////////////////////////////////////////////////////////////////
// End of corrected free boundary condition
///////////////////////////////////////////////////////////////////////


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
    
    psDP[ 1 ] = psD[ 1 ] * psQ[ 1 ];
    
    for ( i = 2 ; i < nN ; i++ ) {
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
    
    scalar sR;
    scalar sDiag, sOffDiag;
    
    sR = 0.5 * sAlpha * sDT / sDX / sDX;
    sDiag = 1 + 2.0 * sR;
    sOffDiag = -sR;
    
    for ( i = 0 ; i <= nN ; i++ ) {
        psA[ i ] = sOffDiag;
        psB[ i ] = sDiag;
        psC[ i ] = sOffDiag;
    }
    
    //psA[ 0 ] = 0.0;
    
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
        psD[ i ] = sOffDiag * psU[ i - 1 ] + sDiag * psU[ i ] + sOffDiag * psU[ i + 1 ];
    }
    
    //psD[ 1      ] += sOffDiag * sBCAdvL;
    //psD[ nN - 1 ] += sOffDiag * sBCAdvR;
    
    return 0;
}


///////////////////////////////////////////////////////////////////////
// End of Crank-Nicolson method
///////////////////////////////////////////////////////////////////////


double g_dTInit = 0.1;
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
    
    double dBCRatioL, dBCRatioR;
    double dBCAdvL, dBCAdvR;
    
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
    nNX = 1024;
    dIntX = ( dXMax - dXMin ) / nNX;
    
    nNT =  1;
    dTInit = g_dTInit;
    dTFin  = 0.1 * 0.1 * nNT;
    dIntT = ( dTFin - dTInit ) / nNT;
    
    dAlpha = g_dDiffusionConst;
    
    // Initialization: making the initial state
    
    for ( i = 0 ; i <= nNX ; i++ ) {
        dXCurr = dXMin + i * dIntX;
        pdU[ i ] = myFuncInit(dXCurr);
    }
    
    // Calculation of the numerical method
    
    myCreateMtxTridiagonal(dIntX, dIntT, dAlpha, nNX, pdA, pdB, pdC);
    myApplyFreeBCToMtx(nNX, pdA, pdB, pdC, &dBCRatioL, &dBCRatioR);
    myCalcMtxTridiagonal(pdA, pdB, pdC, nNX, pdP, pdQ, pdR);
    
    for ( i = 0 ; i < nNT ; i++ ) {
        // Setting BC
        dTCurr = dTInit + i * dIntT; // That's why I starts at 1, not 0 (and ends at nNT, not nNT - 1)
        
        //dBCAdvL = myBCLeft(dTCurr + dIntT);
        //dBCAdvR = myBCRight(dTCurr + dIntT);
        
        // Calculating the function by tridiagonal matrix algorithm
        myCreateVecInput(dIntX, dIntT, dAlpha, pdU, dBCAdvL, dBCAdvR, nNX, pdD);
        myApplyFreeBCToVecInput(nNX, dBCRatioL, dBCRatioR, pdD);
        myCalcSolTridiagonal(pdD, pdP, pdQ, pdR, nNX, pdDP, pdU);
        
        //pdU[ 0   ] = dBCAdvL;
        //pdU[ nNX ] = dBCAdvR;
        myGetBoundaryByFreeBC(nNX, pdU);
    }
    
    // Output
    
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


