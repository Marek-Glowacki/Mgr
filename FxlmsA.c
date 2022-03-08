  
#include <math.h>
#include <stdio.h>
#include <windows.h>
#include <stdlib.h>

#define samples 32000
#define TrainRange 500
#define M 5

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

void Estimate(double *Sh, double *H){
double Input[TrainRange] = {0.0} ;
double Desired[TrainRange] = {0.0};
double EstimateData[TrainRange] = {0.0};
double NewInput[M] = {0.0};
double Temp[M] = {0.0};
double EstimateError = 0;
double Y = 0;

for(int i = 0;i<TrainRange;i++){
    //Sinus
    Input[i] = rand()/(double)RAND_MAX; 
    
 //Propagacja
    for(int j =0;j< M; j++){
            if((i - j) >= 0){
                Desired[i] += Input[i-j] * H[j];
            }
        }
}
for(int i = 0;i<TrainRange;i++){

    for(int j = i;j>(i - M);j--){
        if(j >= 0){
            NewInput[M + (j - i) - 1] = Input[j];
        }
        else{ break;}
    }
   
    Y = 0;
    for(int k = 0;k < M; k++){
        Y +=(Temp[k] * NewInput[k]); 
    }
    
    EstimateError = Desired[i] - Y;

    for(int k = 0;k<M;k++){
        Temp[k] = Temp[k] + (0.2 * EstimateError * NewInput[k]);
      
    }
}
for(int i=0;i<M;i++){
        Sh[i] = Temp[M - i -1];
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









void FxLMS(double *d, double *x, double *error, double *Yout){
    double X[M]= {0.0};
    double X2[M]= {0.0};
    double Xh[M]= {0.0};
    double W[M]= {0.0};
    double H[M] = {0.0125, 0.025, 0.05, 0.1, 0.2}; //Secondary path
    double Sh[M] = {0.0};
    double out[samples]= {0.0};
    double xs[samples]= {0.0};
    Estimate(Sh, H);
    for(int i = 0;i<samples;i++){

    //Sh(z) / filtracja
    for(int j =0;j< M; j++){
            if((i - j) >= 0){
                xs[i] += x[i-j] * Sh[j];

            }
        }    

    //Wagowanie xs
    xs[i] = Aweighting(i,xs);

    //Przygotowanie wejscia 
    for(int j = i;j>(i - M);j--){
        if(j >= 0){
            X[M + (j - i) - 1] = x[j];
            X2[M + (j - i) - 1] = xs[j]; 
            Xh[M+ (j - i) -1] = xs[j]/(xs[j] * xs[j]);   
        }
        else{ break;}
    }



    //Obliczenie wartosci wyjscia
    for(int k = 0;k < M; k++){
        Yout[i] +=(W[k] * X[k]); 
    }

    //Propagacja fali - wyjscie filtra
    for(int j =0;j< M; j++){
            if((i - j) >= 0){
                out[i] += Yout[i-j] * H[j];
            }
        }


    //error
    error[i] = d[i] - out[i];

    //Wagowanie bledu
    error[i] = Aweighting(i,error);

    //aktualizacja wagi
    for(int k = 0;k<M;k++){
        W[k] = W[k] + 2*(0.2 * error[i] * Xh[k]);
    } 
    
    }

}

void main(){
    double noise[samples] = {0.0};
    double error[samples]= {0.0};
    double desired[samples]= {0.0};
    double Anoise[samples]= {0.0};

    
    FILE *D;
    D = fopen("Data.csv", "w+") ;
   
    
    DataGenerator(noise, desired);
    FxLMS(desired,noise, error, Anoise);

    //Kontrolnie do pliku
    for (int i = 0; i< samples; i++){
        fprintf(D, "%d,%f,%f,%f,%f\n",i,noise[i], desired[i], Anoise[i], error[i]);
        
    }

    fclose(D);
    

    //return 0;
}