#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NGRID 1000

// ICs
const double x0 = 0.5;
const double pLeft = 1.0;
const double pRight = 0.1;
const double rhoLeft = 1.0;
const double rhoRight = 0.125;
const double uLeft = 0.0;
const double uRight = 0.0;
const double Gamma = 1.4;
const double a1 = 1.0;

// helper vars and speed of sound
double G1, G2, G3, G4, G5, G6, AR, BR, cLeft, cRight;

// star region variables
double pStar, uStar, uShock, cStarLeft, rhoStarLeft, rhoStarRight;

void get_star_region();

void main() {
    double xi;
    G1 = (Gamma - 1.0)/(2.0 * Gamma);
    G2 = (Gamma + 1.0)/(2.0 * Gamma);
    G3 = 2.0 * Gamma / (Gamma - 1.0);
    G4 = 2.0 / (Gamma - 1.0);
    G5 = 2.0 / (Gamma + 1.0);
    G6 = (Gamma - 1.0)/(Gamma + 1.0);
    AR = G5/rhoRight;
    BR = G6*pRight;
    cLeft = sqrt(Gamma * pLeft / rhoLeft);
    cRight = sqrt(Gamma * pRight / rhoRight);

    double x[NGRID];
    double rho[NGRID];
    for(int i = 0; i < NGRID; i++) {
        x[i] = i * 1.0/(NGRID - 1.0);
    }
    printf("Initial x = %f\n Final x = %f\n", x[0], x[NGRID-1]);

    get_star_region();
    int region;
    for(int i = 0; i < NGRID; i++) {
        xi = (x[i] - x0) / 0.2;
        if(xi <= uLeft - cLeft) rho[i] = rhoLeft;
        else if(xi <= uStar - cStarLeft) {
            double u = G6 * uLeft + G5 * (cLeft + xi);
            rho[i] = rhoLeft * pow((u-xi)/cLeft, G4);
        }
        else if(xi <= uStar) rho[i] = rhoStarLeft;
        else if(xi <= uShock) rho[i] = rhoStarRight;
        else rho[i] = rhoRight;
    }

    FILE* file_ptr;
    file_ptr = fopen("sod.dat", "w");
    for(int i = 0; i < NGRID; i++) fprintf(file_ptr, "%f %f\n", x[i],rho[i]);

    fclose(file_ptr);
}



void get_star_region() {
    double pPrev, f, dfdp;
    double TOL = 1e-10;
    int NRUN = 10e3;

    pStar = 0.5*(pLeft + pRight);
    pPrev = 0.0;
    while((fabs(pStar - pPrev)/(0.5 * (pStar + pPrev)) >= TOL) && (NRUN >= 0)) {
        pPrev = pStar;
        f = G4 * cLeft * (pow(pStar/pLeft, G1) - 1.0)                                  // fL
            + (pStar - pRight) * sqrt(AR/(pStar + BR));                                // fR
        dfdp = sqrt(AR/(pStar + BR)) * (1.0 - 0.5 * (pStar - pRight)/(BR + pStar))     // fL'
            + pow(pStar/pLeft, -G2) / (rhoLeft * cLeft);                               // fR'
        pStar = pStar - f/dfdp;
        NRUN--;
    }
    if(NRUN<0) printf("Error: number runs exceeded max NRUN=%d\n",NRUN);
\

    /* --- Compute necessary star-region values --- */
    uStar = 0.5 * (uLeft + uRight) + 0.5 * ((pStar - pRight) * sqrt(AR/(pStar + BR)) - G4 * cLeft * (pow(pStar/pLeft, G1) - 1.0));
    cStarLeft = cLeft * pow(pStar/pLeft, 1.0/G3);
    uShock = uRight + cRight * sqrt(G2 * (pStar/pRight - 1.0) + 1.0);

    rhoStarLeft = rhoLeft * pow(pStar/pLeft, 1.0/Gamma);
    rhoStarRight = rhoRight * (pStar/pRight + G6) / (G6 * pStar/pRight + 1.0); // Rankine-Hugoniot Shock Condition
}