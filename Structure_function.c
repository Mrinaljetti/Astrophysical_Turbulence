//
//  Structure_function.c
//  Multiphase_Turbulence_C
//
//  Created by Mrinal Jetti on 26/08/20.
//  Copyright © 2020 Mrinal Jetti. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>





int main()
{
    const int G = 32;   //Grid size
    const double PI = 3.141592653589793238;
    const int BOX_LENTGH = 2*PI;
    double x[G],y[G],function[G][G];
    int i,j,mesh_width;
    double strucfunc[2*G][2*G];
    
    // MAKING THE MESHED GRID
    for(i=0;i<G;i++)
    {
        mesh_width = BOX_LENTGH/(G);
        x[i] = i*mesh_width;
    }
    
    for(j=0;j<G;j++)
    {
        mesh_width = BOX_LENTGH/(G);
        y[j] = j*mesh_width;
    }
    
    //FUNCTION DECLERATION
    for(i=0;i<G;i++)
    {
        for(j=0;j<G;j++)
        {
            function[i][j] = sin(x[i] + y[j] );
        }
    }
    
    strucfunc = structure_function(function);
    
    
}


