//
//  ViewController.m
//  numericalMethods
//
//  Created by Ilya on 02.10.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import "ViewController.h"

@interface ViewController ()

@end

@implementation ViewController

//tridiagonal algo usage
//    double x[5] = {-14, -55, 49, 86,8};
//    double a[5] = {0, 7, -4, 7, 4};
//    double b[5] = {8, -19, 21, -23, -7};
//    double c[5] = {-2, 9, -8, 9, 0};
//    size_t N = 5;
//    solve_tridiagonal_in_place_destructive(x, N, a, b, c);
//    //should be -1.0000	3.0000	1.0000	-5.0000	-4.0000
//    NSLog(@"answer: %f %f %f %f %f", x[0], x[1], x[2], x[3], x[4]);


void solve_tridiagonal_in_place_destructive(double x[], const size_t N, const double a[], const double b[], double c[]) {
    /* unsigned integer of same size as pointer */
    size_t i;
    
    /*
     solves Ax = v where A is a tridiagonal matrix consisting of vectors a, b, c
     note that contents of input vector c will be modified, making this a one-time-use function
     x[] - initially contains the input vector v, and returns the solution x. indexed from [0, ..., N - 1]
     N — number of equations
     a[] - subdiagonal (means it is the diagonal below the main diagonal) -- indexed from [1, ..., N - 1]
     b[] - the main diagonal, indexed from [0, ..., N - 1]
     c[] - superdiagonal (means it is the diagonal above the main diagonal) -- indexed from [0, ..., N - 2]
     */
    
    c[0] = c[0] / b[0];
    x[0] = x[0] / b[0];
    
    /* loop from 1 to N - 1 inclusive */
    for (i = 1; i < N; i++) {
        double m = 1.0 / (b[i] - a[i] * c[i - 1]);
        c[i] = c[i] * m;
        x[i] = (x[i] - a[i] * x[i - 1]) * m;
    }
    
    /* loop from N - 2 to 0 inclusive, safely testing loop end condition */
    for (i = N - 1; i-- > 0; )
        x[i] = x[i] - c[i] * x[i + 1];
}

double phi(double x){
    return sin(x);//x+sin(M_PI*x);
}

double phi_0(double t){
    return exp(-t);//0;
}
double phi_l(double t){
    return -exp(-t);//1;
}

double f(double x, double t){
    return 0;
}

double **explicitScheme(int K, int N, double a, double b, double c, double tau, double h, double alpha, double betta, double gamma, double delta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    for (int i = 0; i <= N; ++i)
        U[0][i] = phi(i*h);
    
    for (int k = 0; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            U[k+1][i] = U[k][i] + a*tau/h/h*(U[k][i+1] - 2*U[k][i] + U[k][i-1]) + b*tau/2/h*(U[k][i+1] - U[k][i-1]) + c*tau*U[k][i] + tau*f(i * h, k * tau);
        }
        
      //  двухточечная первого
      //  U[k+1][0] = -alpha/h/(betta - alpha/h)*U[k+1][1]  + phi_0((k + 1) * tau)/(betta - alpha/h);
      //  U[k+1][N] = gamma/h/(delta + gamma/h)*U[k+1][N-1] + phi_l((k + 1) * tau)/(delta + gamma/h);

      //  трехточечная второго
      //  U[k+1][0] = alpha/2/h/(betta - 3*alpha/2/h)*(U[k+1][2] - 4*U[k+1][1]) + phi_0((k+1) * tau)/(betta - 3*alpha/2/h);
      //  U[k+1][N] = gamma/2/h/(delta - 3*gamma/2/h)*(4*U[k+1][N-1] - U[k+1][N-2]) + phi_l((k+1) * tau)/(delta - 3*gamma/2/h);

      //  Двухточечная второго
          U[k+1][0] = 1/(2*a*alpha/h + h*alpha/tau - c*h*alpha - betta*(2*a - b*h))*(U[k+1][1]*2*a*alpha/h + U[k][0]*h*alpha/tau - phi_0((k+1) * tau)*(2*a-b*h)+f(0, (k+1)*tau)*h*alpha);
          U[k+1][N] = 1/(2*a*gamma/h + h*gamma/tau - c*h*gamma + delta*(2*a + b*h))*(U[k+1][N-1]*2*a*gamma/h + U[k][N]*h*gamma/tau + phi_l((k+1) * tau)*(2*a+b*h)+f(N, (k+1)*tau)*h*gamma);
    }
    return U;
}


- (void)viewDidLoad
{
    [super viewDidLoad];
    
    int K = 1000;
    int N = 10;
    
    double l = 3.14; //1
    double T = 10;
    
    double tau = T / K, h = l / N;
    
    double a = 1, b = 0, c = 0;
    double alpha = 1.0, betta = 0.0, gamma = 1.0, delta = 0.0;//double alpha = 0.0, betta = 1.0, gamma = 0.0, delta = 1.0;
    
    if(a*tau/h/h > 0.5){
        printf("bad sigma");
        exit(0);
    }
    
    double **U = explicitScheme(K, N, a, b, c, tau, h, alpha, betta, gamma, delta);
    
    for (int i = 0; i < K; ++i) {
        printf("K = %f\n", i*tau);
        for (int j = 0; j <= N; ++j) {
            printf("%f	%f\n", j*h, U[i][j]);
        }
        printf("\n\n\n");
    }
}

- (void)didReceiveMemoryWarning
{
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

@end
