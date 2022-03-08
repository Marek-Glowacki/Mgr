#include <math.h>
#include <stdio.h>
#include <windows.h>
#include <stdlib.h>

#define samples 32000
#define TrainRange 500
#define M 5

void DataGenerator(double *x){
    for(int i = 0; i<samples;i++){ //random data
    x[i] = rand()/(double)RAND_MAX;
    //x[i] = 1;

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
            
            y += buff[FL-j] * coef[j];

            
        }   

    return y;
}
    



void main(){
    double noise[samples] = {0.0};
    double output[samples] = {0.0};

    
    FILE *D;
    D = fopen("Aweight.csv", "w+") ;
   
    DataGenerator(noise);
     
    for(int i=0;i<samples;i++){

    output[i] = Aweighting(i,noise);
    
    }


    //Kontrolnie do pliku
    for (int i = 0; i< samples; i++){
        fprintf(D, "%d,%f,%f\n",i,noise[i],output[i]);
        
    }

    fclose(D);
    

    //return 0;
}