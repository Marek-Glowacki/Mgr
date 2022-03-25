  
#include <math.h>
#include <stdio.h>
#include <windows.h>
#include <stdlib.h>

#define samples 32000
#define TrainRange 500
#define M 20

void DataGenerator(double *x, double *d){
    double H[M] = {0.0625, 0.125, 0.25, 0.5, 1}; //Primary path

    for(int i = 0; i<samples;i++){ //random data
    x[i] = rand()/(double)RAND_MAX;
    }
    
    for(int i = 0;i<samples;i++){
        for(int j =0;j< M; j++){
            if((i - j) >= 0){
                d[i] += x[i-j] * H[j];
            }
        }
    }
}


double Aweighting(int sampleN, double *x){
    #define FL 7 // Filter lenght
    double y = 0; // output
    double coef[FL] = {0.630620946823873,-1.261241893647745,-0.630620946823875,2.522483787295494,-0.630620946823872,-1.261241893647748,0.630620946823874};
    double buff[FL] = {0.0};

    // Past FL samples
    for(int i =0;i<FL;i++){
        if((sampleN - i) >= 0){
        buff[FL - i] = x[sampleN - i];
        }
    }


    //Convolve
     for(int j =0;j< FL; j++){
            
            y += 0.2*( buff[FL-j] * coef[j]);

            
        }   

    return y;
}



void PANC(double *d, double *x, double *error, double *Yout){
    double X[M]= {0.0};
    double Xh[M]= {0.0};
    double W[M]= {0.0};
    //double H[M] = {0.0125, 0.025, 0.05, 0.1, 0.2}; //Secondary path
    double H[M] = {0.5, 0.12, 1.0, 1.8, 1.69}; 
    double Sh[M] = {0.0};
    double out[samples]= {0.0};
    double xs[samples]= {0.0};
    double xh[samples]= {0.0};
    double eh[samples]= {0.0};

    /* Debugging variables */
    double y = 0.0;
    double e = 0.0;
    double t = 0.0;

FILE *F;
F = fopen("data.csv", "w+");


for(int i =0;i<samples;i++){
    //debug
    t = x[i];

    //Estimate(Sh, H);
    /* Sh filtration */
    for(int j =0;j< M; j++){
            if((i - j) >= 0){
                /*Propagation - secondary path model*/
               // xs[i] += x[i-j] * Sh[j]; //with Estimate();
                xs[i] += x[i-j] * H[j]; //without Estimate(); 
            }
        }
    

    /* Weighting xs*/
    xh[i] = Aweighting(i,xs);

    /* Preparing adaptation imput*/
    for(int j = i;j>(i - M);j--){
        if(j >= 0){
            X[M + (j - i) -1 ] = x[j];
            Xh[M + (j - i) -1 ] = xh[j];
            }else{ break;}
    }

    /* Calculating output */
    for(int k = 0;k < M; k++){
        Yout[i] +=(W[k] * X[k]); 
    }

    /* Output propagation - secondary path */
    for(int j =0;j< M; j++){
            if((i - j) >= 0){
                y += Yout[i-j] * H[j];
                /* Simple output limits  for debugging */
                /*
                if(y > 5){
                    out[i] = 5;
                }else if(y < - 5){
                    out[i] = -5;
                }else{
                    out[i] = y;
                }
                */
                out[i] = y;
            }
    }

    /* Calculating error */
    e = d[i] - out[i];
    error[i] = e;
    
    /* Weighting error */
    eh[i] = Aweighting(i,error);

    /* Filter coefficients update */
    for(int k = 0;k<M;k++){
        W[k] = W[k] + (2*0.2* eh[i]) * Xh[k];
        //fprintf(F, "%d,%f,%f\n", k, eh[i], Xh[k]);
    } 

    /* Printing variables */
    fprintf(F, "%d,%f,%f,%f,%f,%f,%f\n",i,xh[i], x[i], xs[i], error[i], eh[i], out[i]);
    //fprintf(F, "%d,%f,%f,%f\n",i,xh[i], x[i], xs[i]);
}

fclose(F);

}



void main(){
    double noise[samples] = {0.0};
    double error[samples]= {0.0};
    double desired[samples]= {0.0};
    double Anoise[samples]= {0.0};
   
    /* Funcian calling */
    DataGenerator(noise, desired);
    PANC(desired,noise, error, Anoise);

}