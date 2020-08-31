//
//  structure_function.h
//  Astrophysical_Turbulence
//
//  Created by Mrinal Jetti on 29/08/20.
//  Copyright © 2020 Mrinal Jetti. All rights reserved.
//

#ifndef structure_function_h
#define structure_function_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//constant delaration
#define GRID_SIZE 32
#define PI  3.141592653589793238
#define BOX_LENGTH  2*PI
#define UNIT_LENGTH  BOX_LENGTH/GRID_SIZE


double ** structure_function(double f[GRID_SIZE][GRID_SIZE])
{
    //strucfunc is the grid which stores structure function values
    //with each point corresponding to rx,ry,rz, the x,y,z components of r vector
    //the centre of the grid [GRID_SIZE,GRID_SIZE] corresponds to rx=ry=rz=0
    
    //variable declaration
    int i,j,rx,ry;
    double q1_sum,q2_sum,q3_sum,q4_sum;
    int rx1,rx2,rx3,rx4,ry1,ry2,ry3,ry4;
    
    
    //Allocating 2D/3D array memory to strucfunc grid using malloc
    double ** strucfunc = (double **)malloc(sizeof(double *)*(2*GRID_SIZE)) ;
    for(i=0;i<2*GRID_SIZE;i++)
    {
        strucfunc[i] = (double *)malloc(sizeof(double)*(2*GRID_SIZE));
    }

    //Calculating values of structure functions
    for(rx=0;rx<GRID_SIZE;rx++)
    {
        for(ry=0;ry<GRID_SIZE;ry++)
        {
            q1_sum =0;q2_sum=0;q3_sum=0;q4_sum=0;
            for(i=0;i<GRID_SIZE;i++)
            {
                for(j=0;j<GRID_SIZE;j++)
                {
                    rx1 = rx; rx2 = rx; rx3 = rx; rx4 = rx; ry1 = ry; ry2 = ry; ry3 = ry; ry4 = ry;
                    
                    //When any index goes out of bound using periodic boundary
                    //conditions to appropriately set the index value
                    if(rx1+i>=GRID_SIZE){ rx1 = rx1 - GRID_SIZE; }
                    
                    if(ry1+j>=GRID_SIZE){ ry1 = ry1 - GRID_SIZE; }
                        
                    if(rx2>i)           { rx2 = rx2 - GRID_SIZE; }
                    
                    if(ry2+j>=GRID_SIZE){ ry2 = ry2 - GRID_SIZE; }
                    
                    if(rx3>i)           { rx3 = rx3 - GRID_SIZE; }
                    
                    if(ry3>j)           { ry3 = ry3 - GRID_SIZE; }
                    
                    if(rx4+i>=GRID_SIZE){ rx4 = rx4 - GRID_SIZE; }
                    
                    if(ry4>j)           { ry4 = ry4 - GRID_SIZE; }
                        
                    //Calculation of the structure function
                    q1_sum = q1_sum + pow((f[i+rx1][j+ry1]-f[i][j]),2);
                    q2_sum = q2_sum + pow((f[i-rx2][j+ry2]-f[i][j]),2);
                    q3_sum = q3_sum + pow((f[i-rx3][j-ry3]-f[i][j]),2);
                    q4_sum = q4_sum + pow((f[i+rx4][j-ry4]-f[i][j]),2);
                }
            }
            
            //Allocating the Structure function values appropriately
            //to the meshed grid of rx ,ry explained before
            strucfunc[GRID_SIZE+rx][GRID_SIZE+ry] = q1_sum/(pow(GRID_SIZE,2));
            strucfunc[GRID_SIZE-rx][GRID_SIZE+ry] = q2_sum/(pow(GRID_SIZE,2));
            strucfunc[GRID_SIZE-rx][GRID_SIZE-ry] = q3_sum/(pow(GRID_SIZE,2));
            strucfunc[GRID_SIZE+rx][GRID_SIZE-ry] = q4_sum/(pow(GRID_SIZE,2));
        }
    }
    return strucfunc;
}


double * structurefunction_rbin(double **strucfunc)
{
    //Allocating 1Darray memomry to strucfunc_rbin array
    //strucfunc_rbin stores the values of the structure function wrt to mod(r vector) bin
    double *strucfunc_rbin = (double *)malloc(sizeof(double)*GRID_SIZE);
    
    //variable declaration
    int i,j,rx,ry;
    
    for(i=0;i<GRID_SIZE;i++)
    {
        int indices[GRID_SIZE*GRID_SIZE][2]={ };
        int count = 0;
        
        for(rx=0;rx<2*GRID_SIZE;rx++)
        {
            for(ry=0;ry<2*GRID_SIZE;ry++)
            {
                //Creating shells(of unit depth) of increasing radius
                //and noting down the (rx,ry) indices that fall in the shell
                //indices array stores these indices such that each row of indices
                //has (rx,ry) which fall into the shell
                //count notes down the number of points that fall into the shell
                if(((i*i) <= ( (rx-GRID_SIZE)*(rx-GRID_SIZE) + (ry-GRID_SIZE)*(ry-GRID_SIZE) )) &&  (((rx-GRID_SIZE)*(rx-GRID_SIZE) + (ry-GRID_SIZE)*(ry-GRID_SIZE) )< (i+1)*(i+1)))
                {
                    indices[count][0]=rx;
                    indices[count][1]=ry;
                    count = count+1;
                    
                }
            }
        }
        //Normalizing the structure function wrt to rbin
        double sum = 0;
        if(count==0)
        {
            strucfunc_rbin[i] = 0;
        }
        else
        {
            for(j=0;j<count;j++)
            {
                sum = sum + strucfunc[indices[j][0]][indices[j][1]];
            }
            strucfunc_rbin[i] = sum/count;
        }
    }
    
    return strucfunc_rbin;
}



#endif /* structure_function_h */
