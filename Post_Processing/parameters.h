//
//  parameters.c
//  target_0
//
//  Created by Mrinal Jetti on 22/09/20.
//  Copyright © 2020 Mrinal Jetti. All rights reserved.
//

#include <stdio.h>
#include <math.h>

#define nx 64
#define ny 64
#define nz 64
#define nv 5
#define f1 5
#define f2 25
#define fstep 1
#define prefix "data."
#define suffix ".dbl"
#define datdir "/Users/mj/Documents/Pluto/PLUTO/my_problems/turbulence/"
#define data_prefix "dataarray_step."
#define data_suffix ".txt"
#define sb_data_dir "/Users/mj/Documents/Pluto/PLUTO/my_problems/turbulence/sb_data.txt"
#define temp_data_dir "/Users/mj/Documents/Pluto/PLUTO/my_problems/turbulence/temp_data.txt"
#define mrms_data_dir "/Users/mj/Documents/Pluto/PLUTO/my_problems/turbulence/mrms_data.txt"



#define CONST_PI      3.14159265358979
#define CONST_mp      1.67262171e-24     //Proton mass
#define CONST_kB      1.3806505e-16      // Boltzmann constant
#define CONST_mu      0.5                //average mass in CONST_mp units
#define CONST_pc      3.0856775807e18    //Parsec

#define UNIT_VELOCITY (1.e8)
#define UNIT_DENSITY  (CONST_mp*.1)
#define UNIT_LENGTH   (CONST_pc*40.e3)
#define gamma         5./3.                //Value of gamma

#define lambda(T)     2e-27*sqrt(T)
