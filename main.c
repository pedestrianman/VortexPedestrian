//===================================================
// Author: Iñaki Echeverria Huarte
// Project: COVID19 SIMULATIONS
//===================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>  

//#define POSITIONS // Store Positions
//#define ANGULARMOM_TIMESERIES // Store Angular Momentum Time Series
//#define MEAN_ANGULARMOM // Store Mean Angular Momentum Simulation

#define pi 3.141592653589
#define sod sizeof(double)

// Main Routines
void simulate(int rep,int iters, int N, int N_CCW, int N_CW,int TP, double Ped_Exponents[3][2], 
            double Wall_Exponents[3][3], double drift_coeff, double vel_happy, double gamma_wall, double* Ang_Mom_NormPos);
void   ran_seed(long j);
double ran_ran2();
double sign(double x);
inline double rad2deg(double rad);
inline double deg2rad(double deg);
inline double dot_product2D(double *v, double *u);
inline double cross_product2D(double *v, double *u);
inline double mod_2D(double *v);
inline double angular_difference(double *vel_part, double *n);
void correct_posWall(double *x, double *Pos_wall, double *n_wall,double *NewPos_wall);
void calc_force_wall(int rep, int TP, double *pos_part,double *vel_part,int charge_part,
                double* NewPos_wall,double gamma_wall,double Wall_Exponents[3][3],
                double *fnwall_x,double *fnwall_y,double *ftwall_x, double*ftwall_y);
void angularmom(double *x, double *v, int N, double *SS, double *L, double *LnormPos, double *LnormTot);
unsigned long long int vseed;
unsigned long long int vran;

// Characteristic Lengths
double SSX = (11.37 / 2.0); // System Size X
double SSY = (6.7 / 2.0);   // System Size Y
double r = 5.;        //Neighbour radius detections
double r_wall = 3.;   // Wall radius detection
double r_ped = 0.25; // Particle radius

// Simulation Time
double t = 0.0; 
double time_end = 400; // Simulation Time in Seconds
double dt = 0.1e-1; // Time step
double wfPos_s; // Writing frequency positions in secs
double wfL_s; // Calculation frequency of L in secs
int wfPos_it; // Writing frequency positions in iterations
int wfL_it; // Calculation frequency of L in iterations
int reps = 1; // Number of repetitions to be performed


//===================================================
// the main function
//===================================================
int main(int argc, char **argv){

    time_t seed;
    ran_seed(time(&seed)); // Initial random seed with process time
    srand((unsigned) time(&seed+1));
    
    // Time simulation
    wfPos_s = 0.5;
    wfPos_it = round(wfPos_s/dt);

    wfL_s = 0.5;
    wfL_it = round(wfL_s/dt);
    int iters = round(time_end/wfL_s);
    int iters_rm = round(200./wfL_s);
    
    // Force Weights
    double drift_coeff = 4.0; // Weight parameter autopropulsion force
    double vel_happy = 1.5; // Desired Velocity
    double gamma_wall; // Damping Coefficient
    double Ped_Exponents[3][2] = { // Repulsion Force Ped-Ped. First col: Short Range/ Second col: Long Range
        {200.,13.0},
        {2*r_ped, 0.85},
        {13.0, 2*r_ped}
    };
    double Wall_Exponents[3][3] = { // Repulsion Force Wall-Ped. First col: Short Range/ Second col: Long Range/ Third Col: Turning Force
        {200, 15.0, 9},
        {r_ped, 0.4, 0.4},
        {15., r_ped, r_ped}
    };

    double MeanL = 0;

    //int people[12] = {8,10,12,15,18,20,22,25,28,30,32,35}; //For Turning Preference with 60% CW and 40% CCW
    int people[14] = {8,10,12,14,16,18,20,22,24,26,28,30,32,34}; //For normal simulation
    int films = sizeof(people) / sizeof(people[0]);

    for (gamma_wall = 1.5; gamma_wall <= 1.5; gamma_wall += 0.05){
        for(int film = 0; film < films; film++){
            
            int N = people[film];
            double CCW_P = 0.6;        // CCW Percentage
            double CW_P = 1.0 - CCW_P; // CW Percentage
            int N_CCW = round(N * CCW_P);
            int N_CW = round(N * CW_P);
            int TP = 0; // 0 for no turning preference, 1 otherwise

            printf("Simulating %d reps with N = %d  - Gamma: %.3f ¿TP? = %d\n",reps,N,gamma_wall, TP);

            double *Ang_Mom_NormPos = (double*)malloc(sod*iters*reps);
            for (int i=0; i<iters*reps; i++){Ang_Mom_NormPos[i]  = 0.0;} // All zeros
            
            // SIMULATION CALL //
            for(int rep=0; rep < reps; rep++){    
                simulate(rep,iters,N,N_CCW,N_CW,TP,Ped_Exponents, Wall_Exponents, drift_coeff, vel_happy, gamma_wall, Ang_Mom_NormPos);
            }

            #ifdef ANGULARMOM_TIMESERIES
            char filename1[100];
            sprintf(filename1,"YourFileName",N,gamma_wall);
            FILE *file1 = fopen(filename1, "w");
        
            for (int i=0; i<iters; i++){
                for (int rep=0; rep < reps; rep++){
                    fprintf(file1,"%lf\t", Ang_Mom_NormPos[rep*iters + i]);
                }
                fprintf(file1,"\n");
            }
            fclose(file1);
            #endif

            #ifdef MEAN_ANGULARMOM
            char filename2[100];
            sprintf(filename2,"YourFileName",N,gamma_wall);
            FILE *file2 = fopen(filename2, "w");

            for (int rep=0; rep < reps; rep++){
                MeanL = 0.0;
                for (int i=iters_rm; i<iters; i++){
                    MeanL += Ang_Mom_NormPos[rep*iters + i];
                }
                MeanL /= (iters-iters_rm);
                fprintf(file2,"%f\n", MeanL);
            }
            fclose(file2);
            #endif

            free(Ang_Mom_NormPos);

        }// loop over people
    } // loop over gamma values
    
    return 0;
}


//==================================================
// simulation
//==================================================
void simulate(int rep,int iters, int N, int N_CCW, int N_CW,int TP, double Ped_Exponents[3][2], double Wall_Exponents[3][3], double drift_coeff, double vel_happy, double gamma_wall, double* Ang_Mom_NormPos){
    
    int i, j, k; // Index Variables
    double SS[2] = {2. * SSX, 2. * SSY};
    double *x = (double*)malloc(sod*2*N); // Positions
    double *v = (double*)malloc(sod*2*N); // Velocities
    double *f = (double*)malloc(sod*2*N); // Forces
    int *type   = (int*)malloc(sizeof(int)*N); // CCW or CW turn preference
    
    // Initialization
    for (i=0; i<N; i++){ type[i] = 0;} // All zeros
    for (i=0; i<2*N; i++){x[i] = v[i] = f[i]  = 0.0;} // All zeros
    
    // Reading HexGrid Positions
    FILE *file3 = fopen("HexGrid_InitialPos.txt", "rt");
    double *hexgrid = (double*)malloc(sod*2); // Hexagonal Grid Positions Vector 
    double xgrid,ygrid;
    int cont = 0;
    while(fscanf(file3,"%lf %lf",&xgrid,&ygrid)!= EOF){
        hexgrid[2*cont] = xgrid; hexgrid[2*cont + 1] = ygrid;
        hexgrid = realloc(hexgrid,2*sod*(cont + 2));
        cont ++;
    }
    fclose(file3);

    // Random Permutation
    int perm[cont];
    for (int i = 0; i < cont; i++){
        perm[i] = i;
    }

    // Random permutation the order
    for (int i = 0; i < cont; i++) {
        int j, t;
        j = rand() % (cont-i) + i;
        t = perm[j]; perm[j] = perm[i]; perm[i] = t; // Swap i and j
    }

    //-------------------------------------------------
    // initialize
    
    double vcmx = 0.0;
    double vcmy = 0.0;

    for (i=0; i<N; i++){

        int idx = perm[i];
        double t = 2*pi*ran_ran2();

        x[2*i+0] = hexgrid[2*idx];    
        x[2*i+1] = hexgrid[2*idx + 1];  
        v[2*i+0] = vel_happy * cos(t);
        v[2*i+1] = vel_happy * sin(t);

        vcmx+= v[2*i+0];
        vcmy+= v[2*i+1];

        if(i < N_CCW){
            type[i] = 1;
        }else{
            type[i] = -1;
        }
    }
    free(hexgrid);

    // The center of mass has zero velocity initially
    vcmx /= (double)N; vcmy /= (double)N;
    
    for(i=0; i<N; i++){
        v[2*i+0] -=  vcmx;
        v[2*i+1] -=  vcmy;
    }

    #ifdef POSITIONS
    FILE *file4;
    if(rep%1 == 0){
        char filename[100];
        sprintf(filename,"YourFileName",N,rep);
        file4 = fopen(filename, "w");   
    }
    #endif

    int iter = 0;

    for (t=0.0; t<time_end; t+=dt){                 /////////////// TEMPORAL LOOP ///////////////
        
        #ifdef POSITIONS
        if(rep%1 == 0 && iter%wfPos_it == 0){ // Print positions every wf iters
            for (i=0; i<N; i++){
                fprintf(file4,"%lf\t%lf\t%lf\t%lf\t%d\n", x[2*i],x[2*i+1], v[2*i],v[2*i+1],type[i]);
            }   
        }
        #endif
        
        // Neighbours detection
        double *Pdist = (double*)malloc(sod*N*N); //Distances
        double delta,dist,l;
        for (i=0; i<N*N; i++){Pdist[i] = -1.0;} 

        for(i=0; i<N; i++){
            for(j=i+1; j<N; j++){
                dist = 0;

                for(k = 0; k< 2;k++){
                    delta = x[2*j+k] - x[2*i+k];
                    dist += delta*delta;
                }
                l = sqrt(dist);
                
                if(l < r) {
                    Pdist[i*N + j] = l;
                    Pdist[j*N + i] = l; 
                }
            }
        }

        for (i=0; i<N; i++){ // LOOP OVER EACH PARTICLE
            
            int ind;
            double vlen;
            double co, co1, norm;
            double pos_part[2];
            double vel_part[2];
            int charge_part;    

            f[2*i+0] = 0.0;
            f[2*i+1] = 0.0;
            pos_part[0] = x[2*i],
            pos_part[1] = x[2*i+1];
            vel_part[0] = v[2*i],
            vel_part[1] = v[2*i+1];
            charge_part = type[i];

            //===============================================
            // Repulsion Force Against Pedestrians

            for(j=0; j<N ; j++){

                ind = i*N + j;

                if(i != j && Pdist[ind] > 1e-6){
            
                    l = Pdist[ind];
                    co = 0.0;
                   
                    if (l > Ped_Exponents[1][0]){   
                        co = Ped_Exponents[0][1] * exp(-(l - Ped_Exponents[2][1]) / Ped_Exponents[1][1]);
                    }
                    else{   
                        co1 = (1.0 - l/ Ped_Exponents[1][0]);
                        co = Ped_Exponents[0][0] * co1 * sqrt(co1) + Ped_Exponents[2][0];
                    }
                    for (k = 0; k < 2; k++){
                        norm = x[2*j+k]-x[2*i+k];
                        f[2*i+k] += -norm/l*co;
                    }
                }
            }

            //===============================================
            // Repulsion Force Against Walls

            double fnwall_x;
            double fnwall_y;
            double ftwall_x;
            double ftwall_y;
            double pos_wall[2];
            double n_wall[2];
            double *NewPos_wall = malloc(sod*2);

            fnwall_x = 0.0;
            fnwall_y = 0.0;
            ftwall_x = 0.0;
            ftwall_y = 0.0;


            if(x[2*i] > -r_wall + SSX){ // Right
                pos_wall[0] = SSX; pos_wall[1] = 0.0; 
                n_wall[0] = -1. ; n_wall[1] = 0.;
                correct_posWall(pos_part,pos_wall,n_wall,NewPos_wall);
                calc_force_wall(rep,TP,pos_part,vel_part,charge_part,NewPos_wall,gamma_wall,Wall_Exponents,&fnwall_x,&fnwall_y,&ftwall_x,&ftwall_y);
            }
            if(x[2*i] < r_wall - SSX){ // Left
                pos_wall[0] = -SSX; pos_wall[1] = 0.0;
                n_wall[0] = 1. ; n_wall[1] = 0.;
                correct_posWall(pos_part,pos_wall,n_wall,NewPos_wall);
                calc_force_wall(rep,TP,pos_part,vel_part,charge_part,NewPos_wall,gamma_wall,Wall_Exponents,&fnwall_x,&fnwall_y,&ftwall_x,&ftwall_y);
            }            
            if(x[2*i+1] < r_wall - SSY){ // Bottom
                pos_wall[0] = 0.0; pos_wall[1] = -SSY; 
                n_wall[0] = 0. ; n_wall[1] = 1.;
                correct_posWall(pos_part,pos_wall,n_wall,NewPos_wall);
                calc_force_wall(rep,TP,pos_part,vel_part,charge_part,NewPos_wall,gamma_wall,Wall_Exponents,&fnwall_x,&fnwall_y,&ftwall_x,&ftwall_y);
            }

            if(x[2*i+1] > -r_wall + SSY){ // Top
                pos_wall[0] = 0.0; pos_wall[1] = SSY; 
                n_wall[0] = 0. ; n_wall[1] = -1.;
                correct_posWall(pos_part,pos_wall,n_wall,NewPos_wall);
                calc_force_wall(rep,TP,pos_part,vel_part,charge_part,NewPos_wall,gamma_wall,Wall_Exponents,&fnwall_x,&fnwall_y,&ftwall_x,&ftwall_y);
            }

            free(NewPos_wall);
            f[2*i+0] += fnwall_x + ftwall_x;
            f[2*i+1] += fnwall_y + ftwall_y;
            

            //====================================
            // self-propulsion
            vlen = sqrt(v[2*i+0]*v[2*i+0] + v[2*i+1]*v[2*i+1]);
            
            if (vlen > 1e-6){
                f[2*i+0] += drift_coeff*(vel_happy - vlen)*v[2*i+0]/vlen;
                f[2*i+1] += drift_coeff*(vel_happy - vlen)*v[2*i+1]/vlen;
            }
            
        }
        // Now integrate the forces since we have found them

        for (i=0;i<N;i++){

            // Euler Method 
            v[2*i+0] += f[2*i+0] * dt;
            v[2*i+1] += f[2*i+1] * dt;

            x[2*i+0] += v[2*i+0] * dt;
            x[2*i+1] += v[2*i+1] * dt;

            // just check for errors
            if (x[2*i+0] > SSX || x[2*i+0] < -SSX ||
                x[2*i+1] > SSY || x[2*i+1] < -SSY)
            printf("out of bounds\n");
        }

        // Angular Momentum measure

        if(iter%wfL_it == 0){
            int indx = round(iter/wfL_it);
            double L = 0.0; 
            double LnormPos = 0.0;
            double LnormTot = 0.0;
            angularmom(x, v, N, SS,&L,&LnormPos,&LnormTot);
            Ang_Mom_NormPos[rep*iters + indx] = LnormPos;
        }

        iter ++;

        free(Pdist);
    }

    // end of the program, cleanup
    //----------------------------------------------

    #ifdef POSITIONS
    if(rep%1 == 0){
        fclose(file4);
    }
    #endif

    free(x);
    free(v);
    free(f);
    free(type);
    
}

//=======================================
// walls interaction functions
//=======================================
void correct_posWall(double *x, double *Pos_wall, double *n_wall,double *NewPos_wall) {

    int k;
    double A,B,D,distance;
    
    double wall2part[2] = {x[0]-Pos_wall[0],x[1]-Pos_wall[1]};

    double where = dot_product2D(wall2part,n_wall);

    if(where > 0){ n_wall[0] *= -1.; n_wall[1] *= -1.;}
    D = -1.*dot_product2D(n_wall,Pos_wall);
    A = n_wall[0];
    B = n_wall[1];
    distance = A*x[0] + B*x[1] + D;
    distance = fabs(distance);

    for(k = 0; k<2 ; k++) NewPos_wall[k] = x[k] + distance*n_wall[k];
    
}

void calc_force_wall(int rep, int TP, double *pos_part,double *vel_part,int charge_part,
                double* NewPos_wall,double gamma_wall,double Wall_Exponents[3][3],
                double *fnwall_x,double *fnwall_y,double *ftwall_x, double*ftwall_y){

    double diff[2] = {NewPos_wall[0]-pos_part[0],NewPos_wall[1]-pos_part[1]};
    double l = sqrt(diff[0]*diff[0] + diff[1]*diff[1]);
    double n[2] = {diff[0]/l,diff[1]/l};
    double Ang_n = rad2deg(atan2(n[1],n[0]));
    double Ang_t = Ang_n + 90.*(double)charge_part;
    double t[2] = {cos(deg2rad(Ang_t)),sin(deg2rad(Ang_t))};
    double ang_diff = angular_difference(vel_part, n);
    ang_diff = fabs(ang_diff);

    double fact = 1.0;
    double damp = gamma_wall;
    
    // Normal repulsive force
    double co, co1;
    if(l > Wall_Exponents[1][0]){
        co = Wall_Exponents[0][1] * exp(-(l - Wall_Exponents[2][1]) / Wall_Exponents[1][1]);
        if(co < 0.01*Wall_Exponents[0][1]) damp = 0.0; // Damping is active only when the exponential part starts   
    }else{
        co1 = (1.0 - l/Wall_Exponents[1][0]);
        co = Wall_Exponents[0][0] * co1 * sqrt(co1) + Wall_Exponents[2][0];
    }
    
    *fnwall_x += (-co - fact*damp*dot_product2D(vel_part,n))*n[0];
    *fnwall_y += (-co - fact*damp*dot_product2D(vel_part,n))*n[1]; 

    // Tangential repulsive force
    if(TP){
        if(ang_diff > 90.) fact = 0.0; // pedestrian moves away from the wall
        double c1 =  Wall_Exponents[0][2] * exp(-(l - Wall_Exponents[2][2]) / Wall_Exponents[1][2]);
        *ftwall_x += fact * c1 * t[0] * cos(deg2rad(ang_diff));
        *ftwall_y += fact * c1 * t[1] * cos(deg2rad(ang_diff));
    }
}

//==========================================
// measurement functions
//=========================================

void angularmom(double *x, double *v, int N, double *SS, double *L, double *LnormPos, double *LnormTot){
    
    int i=0;
    double CP;

    for (i=0; i<N; i++){
        double pos_part[2] = {x[2*i],x[2*i+1]};
        double vel_part[2] = {v[2*i],v[2*i+1]};
        CP = cross_product2D(pos_part, vel_part);
        *L += CP;
        *LnormPos += CP/mod_2D(pos_part);
        *LnormTot += CP/(mod_2D(pos_part)*mod_2D(vel_part));
    }
    *L /= (double)N;
    *LnormPos /= (double)N;
    *LnormTot /= (double)N;
}

//=================================================
// extra stuff
//=================================================

void ran_seed(long j){
  vseed = j;  vran = 4101842887655102017LL;
  vran ^= vseed; 
  vran ^= vran >> 21; vran ^= vran << 35; vran ^= vran >> 4;
  vran = vran * 2685821657736338717LL;
}

double ran_ran2(){
    vran ^= vran >> 21; vran ^= vran << 35; vran ^= vran >> 4;
    unsigned long long int t = vran * 2685821657736338717LL;
    return 5.42101086242752217e-20*t;
}

double sign(double x){
    
    if (x > 0) return 1.0;
    if (x < 0) return -1.0;
    return 0.0;
}

inline double dot_product2D(double *v, double *u)
{
    double result = 0.0;
    for (int i = 0; i < 2; i++)
        result += v[i]*u[i];
    return result;
}

inline double cross_product2D(double *v, double *u)
{   
    double result = 0.0;
    result = v[0]*u[1]-v[1]*u[0]; 
    return result;
}

inline double angular_difference(double *vel_part, double *n)
{   

    double c = cross_product2D(vel_part,n);
    double d = dot_product2D(vel_part,n); 
    double thetaInRad = atan2(c,d);
    double thetaInDeg = rad2deg(thetaInRad);
    return thetaInDeg;

}

inline double mod_2D(double *v)
{
    return sqrt(v[0]*v[0] + v[1]*v[1]);
}

inline double rad2deg(double rad)
{
    return rad*180.0/pi;
}

inline double deg2rad(double deg)
{
    return deg*pi/180.0;
}

