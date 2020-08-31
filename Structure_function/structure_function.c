//
//  structure_function.c
//  Astrophysical_Turbulence
//
//  Created by Mrinal Jetti on 29/08/20.
//  Copyright © 2020 Mrinal Jetti. All rights reserved.
//

#include "structure_function.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>




int main(int argc, char *argv[])
{
    //Variable declaration
    int i,j;
    double x[GRID_SIZE],y[GRID_SIZE];
    double f[GRID_SIZE][GRID_SIZE];
    double **strucfunc;
    double *strucfunc_rbin;

    //Making the Meshed Grid
    for(i=0;i<GRID_SIZE;i++)
    {
        x[i] = i*(UNIT_LENGTH);
    }

    for(j=0;j<GRID_SIZE;j++)
    {
        y[j] = j*(UNIT_LENGTH);
    }

    //Declaration of Variation in the box
    for(i=0;i<GRID_SIZE;i++)
    {
        for(j=0;j<GRID_SIZE;j++)
        {
            f[i][j] = sin(x[i] + y[j]);
        }
    }

    strucfunc = structure_function(f);

    strucfunc_rbin = structurefunction_rbin(strucfunc);
    
    
    for(j=0;j<GRID_SIZE;j++)
    {
        printf("%f\n",strucfunc_rbin[j]);
    }

    return 0;

}













