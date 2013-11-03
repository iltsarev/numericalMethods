//
//  hyperbolicViewController.m
//  numericalMethods
//
//  Created by Ilya on 03.11.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import "hyperbolicViewController.h"
#import "CorePlotViewController.h"

#define NUMBERS_ONLY @"1234567890.-"
#define CHARACTER_LIMIT 5

double a,b,c,e,alpha,betta,gama,delta,l;
long startApproximation;

@interface hyperbolicViewController ()

@end

@implementation hyperbolicViewController
@synthesize HUD;
@synthesize system;

double * processTridiagonalMatrixH(double *x, const size_t N, const double *a, const double *b, double *c) {
    
    size_t i;
    c[0] = c[0] / b[0];
    x[0] = x[0] / b[0];
    
    for (i = 1; i < N; i++) {
        double m = 1.0 / (b[i] - a[i] * c[i - 1]);
        c[i] = c[i] * m;
        x[i] = (x[i] - a[i] * x[i - 1]) * m;
    }
    
    
    for (i = N - 1; i-- > 0; )
        x[i] = x[i] - c[i] * x[i + 1];
    
    double *answer = (double *)malloc((N)*sizeof(double));
    for (int j = 0; j < N; ++j)
        answer[j] = x[j];
    return answer;
}

#pragma mark - init funcs
- (double)functionH: (double)x withTime:(double)t  {
    return exp(-t)*cos(x);
    //return exp(-a*t)*sin(x);
}

double phiH(double x){
    return sin(x);//x+sin(M_PI*x);
}

double phiH_0(double t){
    return exp(-t);
    //return exp(-a*t);//0;
}
double phiH_l(double t){
    return -exp(-t);
    //    return -exp(-a*t);//1;
}

double psiH_1(double x){
    return cos(x);
    //return 0;
}

double psiH_1d(double x){
    return -sin(x);
}

double psiH_1dd(double x){
    return -cos(x);
}

double psiH_2(double x){
    return -cos(x);
    //return 0;
}


double fH(double x, double t){
    return sin(x)*exp(-t);
    //   return 0;
}

#pragma mark - hyperbolic_explicitScheme

double **hyperbolic_explicitScheme_TwoPoint_FirstOrder(int K, int N, double a, double b, double c, double e, double tau, double h, double alpha, double betta, double gama, double delta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    
    switch (startApproximation) {
        case 0:
        {
            //Первый порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + tau*psiH_2(x_i);
            }
            break;
        }
        case 1:
        {
            //Второй порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + (tau - e*tau*tau/2)*psiH_2(x_i) + a*tau*tau/2*psiH_1dd(x_i) + b*tau*tau/2*psiH_1d(x_i) + c*tau*tau/2*psiH_1(x_i) + tau*tau/2*fH(x_i, 0);
            }
            break;
        }
        default:
            break;
    }
    
    for (int k = 1; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            U[k+1][i] = 1/(e*tau+2)*(2*a*tau*tau/h/h*(U[k][i+1] - 2*U[k][i] + U[k][i-1]) + b*tau*tau/h*(U[k][i+1] - U[k][i-1]) + U[k][i]*(2*c*tau*tau + 4) + U[k-1][i]*(e*tau -2) + 2*tau*tau*fH(i * h, k * tau));
        }
        //  двухточечная первого
        U[k+1][0] = -alpha/h/(betta - alpha/h)*U[k+1][1]  + phiH_0((k + 1) * tau)/(betta - alpha/h);
        U[k+1][N] = gama/h/(delta + gama/h)*U[k+1][N-1] + phiH_l((k + 1) * tau)/(delta + gama/h);
    }
    
    return U;
}

double **hyperbolic_explicitScheme_TwoPoint_SecondOrder(int K, int N, double a, double b, double c, double e, double tau, double h, double alpha, double betta, double gama, double delta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    
    switch (startApproximation) {
        case 0:
        {
            //Первый порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + tau*psiH_2(x_i);
            }
            break;
        }
        case 1:
        {
            //Второй порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + (tau - e*tau*tau/2)*psiH_2(x_i) + a*tau*tau/2*psiH_1dd(x_i) + b*tau*tau/2*psiH_1d(x_i) + c*tau*tau/2*psiH_1(x_i) + tau*tau/2*fH(x_i, 0);
            }
            break;
        }
        default:
            break;
    }
    
    
    for (int k = 1; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            U[k+1][i] = 1/(e*tau+2)*(2*a*tau*tau/h/h*(U[k][i+1] - 2*U[k][i] + U[k][i-1]) + b*tau*tau/h*(U[k][i+1] - U[k][i-1]) + U[k][i]*(2*c*tau*tau + 4) + U[k-1][i]*(e*tau -2) + 2*tau*tau*fH(i * h, k * tau));
        }
        //  Двухточечная второго
        U[k+1][0] = (2*h*alpha/tau/tau * U[k][0] - h*alpha/tau/tau*U[k-1][0] + h*e*alpha/2/tau*U[k-1][0] + h*alpha*fH(0, (k + 1)*tau) - phiH_0((k+1)*tau)*(2*a-b*h) + 2*a*alpha/h*U[k+1][1]) / (2*a*alpha/h + h*alpha/tau/tau + h*e*alpha/2/tau - c*h*alpha - betta*(2*a-b*h));
        U[k+1][N] = (2*h*gama/tau/tau*U[k][N] - h*gama/tau/tau*U[k-1][N] + h*e*gama/2/tau*U[k-1][N] + h*gama*fH(N, (k + 1)*tau) + phiH_l((k+1)*tau)*(2*a + b*h)  + 2*a*gama/h*U[k+1][N-1]) / (2*a*gama/h + h*gama/tau/tau + h*e*gama/2/tau - c*h*gama + delta*(2*a+b*h));
    }
    
    return U;
}


double **hyperbolic_explicitScheme_ThreePoint_SecondOrder(int K, int N, double a, double b, double c, double e, double tau, double h, double alpha, double betta, double gama, double delta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    switch (startApproximation) {
        case 0:
        {
            //Первый порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + tau*psiH_2(x_i);
            }
            break;
        }
        case 1:
        {
            //Второй порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + (tau - e*tau*tau/2)*psiH_2(x_i) + a*tau*tau/2*psiH_1dd(x_i) + b*tau*tau/2*psiH_1d(x_i) + c*tau*tau/2*psiH_1(x_i) + tau*tau/2*fH(x_i, 0);
            }
            break;
        }
        default:
            break;
    }
    
    for (int k = 1; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            U[k+1][i] = 1/(e*tau+2)*(2*a*tau*tau/h/h*(U[k][i+1] - 2*U[k][i] + U[k][i-1]) + b*tau*tau/h*(U[k][i+1] - U[k][i-1]) + U[k][i]*(2*c*tau*tau + 4) + U[k-1][i]*(e*tau -2) + 2*tau*tau*fH(i * h, k * tau));
        }
        //  трехточечная второго
                U[k+1][0] = alpha/2/h/(betta - 3*alpha/2/h)*(U[k+1][2] - 4*U[k+1][1]) + phiH_0((k+1) * tau)/(betta - 3*alpha/2/h);
                U[k+1][N] = gama/2/h/(delta - 3*gama/2/h)*(4*U[k+1][N-1] - U[k+1][N-2]) + phiH_l((k+1) * tau)/(delta - 3*gama/2/h);
    }
    
    return U;
}


#pragma mark - hyperbolic_unstableScheme

double **hyperbolic_unstableScheme_TwoPoint_FirstOrder(int K, int N, double a, double b, double c, double e, double tau, double h, double alpha, double betta, double gama, double delta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    switch (startApproximation) {
        case 0:
        {
            //Первый порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + tau*psiH_2(x_i);
            }
            break;
        }
        case 1:
        {
            //Второй порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + (tau - e*tau*tau/2)*psiH_2(x_i) + a*tau*tau/2*psiH_1dd(x_i) + b*tau*tau/2*psiH_1d(x_i) + c*tau*tau/2*psiH_1(x_i) + tau*tau/2*fH(x_i, 0);
            }
            break;
        }
        default:
            break;
    }
    
    for (int k = 1; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            U[k+1][i] = 1/(e*tau+2)*(2*a*tau*tau/h/h*(U[k-1][i+1]- 2*U[k-1][i] + U[k-1][i-1]) + b*tau*tau/h*(U[k-1][i+1] - U[k-1][i-1]) + U[k-1][i]*(2*c*tau*tau + e*tau - 2) + 4*U[k][i] + 2*tau*tau*fH(i * h, (k - 1) * tau));
        }
        //  двухточечная первого
        U[k+1][0] = -alpha/h/(betta - alpha/h)*U[k+1][1]  + phiH_0((k + 1) * tau)/(betta - alpha/h);
        U[k+1][N] = gama/h/(delta + gama/h)*U[k+1][N-1] + phiH_l((k + 1) * tau)/(delta + gama/h);
    }
    
    return U;
}

double **hyperbolic_unstableScheme_TwoPoint_SecondOrder(int K, int N, double a, double b, double c, double e, double tau, double h, double alpha, double betta, double gama, double delta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    switch (startApproximation) {
        case 0:
        {
            //Первый порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + tau*psiH_2(x_i);
            }
            break;
        }
        case 1:
        {
            //Второй порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + (tau - e*tau*tau/2)*psiH_2(x_i) + a*tau*tau/2*psiH_1dd(x_i) + b*tau*tau/2*psiH_1d(x_i) + c*tau*tau/2*psiH_1(x_i) + tau*tau/2*fH(x_i, 0);
            }
            break;
        }
        default:
            break;
    }
    
    for (int k = 1; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            U[k+1][i] = 1/(e*tau+2)*(2*a*tau*tau/h/h*(U[k-1][i+1]- 2*U[k-1][i] + U[k-1][i-1]) + b*tau*tau/h*(U[k-1][i+1] - U[k-1][i-1]) + U[k-1][i]*(2*c*tau*tau + e*tau - 2) + 4*U[k][i] + 2*tau*tau*fH(i * h, (k - 1) * tau));
        }
        //  Двухточечная второго
        U[k+1][0] = (2*h*alpha/tau/tau * U[k][0] - h*alpha/tau/tau*U[k-1][0] + h*e*alpha/2/tau*U[k-1][0] + h*alpha*fH(0, (k + 1)*tau) - phiH_0((k+1)*tau)*(2*a-b*h) + 2*a*alpha/h*U[k+1][1]) / (2*a*alpha/h + h*alpha/tau/tau + h*e*alpha/2/tau - c*h*alpha - betta*(2*a-b*h));
        U[k+1][N] = (2*h*gama/tau/tau*U[k][N] - h*gama/tau/tau*U[k-1][N] + h*e*gama/2/tau*U[k-1][N] + h*gama*fH(N, (k + 1)*tau) + phiH_l((k+1)*tau)*(2*a + b*h)  + 2*a*gama/h*U[k+1][N-1]) / (2*a*gama/h + h*gama/tau/tau + h*e*gama/2/tau - c*h*gama + delta*(2*a+b*h));
    }
    
    return U;
}


double **hyperbolic_unstableScheme_ThreePoint_SecondOrder(int K, int N, double a, double b, double c, double e, double tau, double h, double alpha, double betta, double gama, double delta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    switch (startApproximation) {
        case 0:
        {
            //Первый порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + tau*psiH_2(x_i);
            }
            break;
        }
        case 1:
        {
            //Второй порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + (tau - e*tau*tau/2)*psiH_2(x_i) + a*tau*tau/2*psiH_1dd(x_i) + b*tau*tau/2*psiH_1d(x_i) + c*tau*tau/2*psiH_1(x_i) + tau*tau/2*fH(x_i, 0);
            }
            break;
        }
        default:
            break;
    }
    
    for (int k = 1; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            U[k+1][i] = 1/(e*tau+2)*(2*a*tau*tau/h/h*(U[k-1][i+1]- 2*U[k-1][i] + U[k-1][i-1]) + b*tau*tau/h*(U[k-1][i+1] - U[k-1][i-1]) + U[k-1][i]*(2*c*tau*tau + e*tau - 2) + 4*U[k][i] + 2*tau*tau*fH(i * h, (k - 1) * tau));
        }
        //  трехточечная второго
        U[k+1][0] = alpha/2/h/(betta - 3*alpha/2/h)*(U[k+1][2] - 4*U[k+1][1]) + phiH_0((k+1) * tau)/(betta - 3*alpha/2/h);
        U[k+1][N] = gama/2/h/(delta - 3*gama/2/h)*(4*U[k+1][N-1] - U[k+1][N-2]) + phiH_l((k+1) * tau)/(delta - 3*gama/2/h);
    }
    
    return U;
}


#pragma mark - hyperbolic_implicitScheme

double **hyperbolic_implicitScheme_TwoPoint_FirstOrder(int K, int N, double a, double b, double c, double e, double tau, double h, double alpha, double betta, double gama, double delta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    double *lower = (double *)malloc((N+1) * sizeof(double));
    double *mid = (double *)malloc((N+1) * sizeof(double));
    double *upper = (double *)malloc((N+1) * sizeof(double));
    double *answer = (double *)malloc((N+1) * sizeof(double));
    
    memset(lower, 0, (N+1) * sizeof(double));
    memset(upper, 0, (N+1) * sizeof(double));
    memset(mid, 0, (N+1) * sizeof(double));
    memset(answer, 0, (N+1) * sizeof(double));
    
    
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    switch (startApproximation) {
        case 0:
        {
            //Первый порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + tau*psiH_2(x_i);
            }
            break;
        }
        case 1:
        {
            //Второй порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + (tau - e*tau*tau/2)*psiH_2(x_i) + a*tau*tau/2*psiH_1dd(x_i) + b*tau*tau/2*psiH_1d(x_i) + c*tau*tau/2*psiH_1(x_i) + tau*tau/2*fH(x_i, 0);
            }
            break;
        }
        default:
            break;
    }
    for (int k = 1; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            lower[i] = 2*a*tau*tau/h/h - b*tau*tau/h;
            mid[i] = -2 - e*tau - 4*a*tau*tau/h/h + 2*c*tau*tau;
            upper[i] = 2*a*tau*tau/h/h + b*tau*tau/h;
            answer[i] = -4*U[k][i] + U[k-1][i]*(2-e*tau) - 2*tau*tau*fH(i * h, (k+1) * tau);
        }
        //Двухточечная первого порядка
        lower[0] = 0.0;
        mid[0] = (betta - alpha/h);
        upper[0] = alpha/h;
        answer[0] = phiH_0((k+1)*tau);

        lower[N] = -gama/h;
        mid[N] = delta + gama/h;
        upper[N] = 0.0;
        answer[N] = phiH_l((k+1)*tau);
        
        U[k+1] = processTridiagonalMatrixH(answer, N+1, lower, mid, upper);
    }
    
    return U;
}

double **hyperbolic_implicitScheme_TwoPoint_SecondOrder(int K, int N, double a, double b, double c, double e, double tau, double h, double alpha, double betta, double gama, double delta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    double *lower = (double *)malloc((N+1) * sizeof(double));
    double *mid = (double *)malloc((N+1) * sizeof(double));
    double *upper = (double *)malloc((N+1) * sizeof(double));
    double *answer = (double *)malloc((N+1) * sizeof(double));
    
    memset(lower, 0, (N+1) * sizeof(double));
    memset(upper, 0, (N+1) * sizeof(double));
    memset(mid, 0, (N+1) * sizeof(double));
    memset(answer, 0, (N+1) * sizeof(double));
    
    
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    switch (startApproximation) {
        case 0:
        {
            //Первый порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + tau*psiH_2(x_i);
            }
            break;
        }
        case 1:
        {
            //Второй порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + (tau - e*tau*tau/2)*psiH_2(x_i) + a*tau*tau/2*psiH_1dd(x_i) + b*tau*tau/2*psiH_1d(x_i) + c*tau*tau/2*psiH_1(x_i) + tau*tau/2*fH(x_i, 0);
            }
            break;
        }
        default:
            break;
    }
    for (int k = 1; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            lower[i] = 2*a*tau*tau/h/h - b*tau*tau/h;
            mid[i] = -2 - e*tau - 4*a*tau*tau/h/h + 2*c*tau*tau;
            upper[i] = 2*a*tau*tau/h/h + b*tau*tau/h;
            answer[i] = -4*U[k][i] + U[k-1][i]*(2-e*tau) - 2*tau*tau*fH(i * h, (k+1) * tau);
        }
        //Двухточечная второго порядка
        lower[0] = 0.0;
        mid[0] = 2*a*alpha/h + h*alpha/tau/tau + h*e*alpha/2/tau - c*h*alpha - betta*(2*a-b*h);
        upper[0] = -2*a*alpha/h;
        answer[0] = 2*h*alpha/tau/tau * U[k][0] - h*alpha/tau/tau*U[k-1][0] + h*e*alpha/2/tau*U[k-1][0] + h*alpha*fH(0, (k + 1)*tau) - phiH_0((k+1)*tau)*(2*a-b*h);
        
        lower[N] = -2*a*gama/h;
        mid[N] = 2*a*gama/h + h*gama/tau/tau + h*e*gama/2/tau - c*h*gama + delta*(2*a+b*h);
        upper[N] = 0.0;
        answer[N] = 2*h*gama/tau/tau*U[k][N] - h*gama/tau/tau*U[k-1][N] + h*e*gama/2/tau*U[k-1][N] + h*gama*fH(N, (k + 1)*tau) + phiH_l((k+1)*tau)*(2*a + b*h);
        
        
        U[k+1] = processTridiagonalMatrixH(answer, N+1, lower, mid, upper);
    }
    
    return U;
}


double **hyperbolic_implicitScheme_ThreePoint_SecondOrder(int K, int N, double a, double b, double c, double e, double tau, double h, double alpha, double betta, double gama, double delta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    double *lower = (double *)malloc((N+1) * sizeof(double));
    double *mid = (double *)malloc((N+1) * sizeof(double));
    double *upper = (double *)malloc((N+1) * sizeof(double));
    double *answer = (double *)malloc((N+1) * sizeof(double));
    
    memset(lower, 0, (N+1) * sizeof(double));
    memset(upper, 0, (N+1) * sizeof(double));
    memset(mid, 0, (N+1) * sizeof(double));
    memset(answer, 0, (N+1) * sizeof(double));
    
    
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    switch (startApproximation) {
        case 0:
        {
            //Первый порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + tau*psiH_2(x_i);
            }
            break;
        }
        case 1:
        {
            //Второй порядок
            for (int i = 0; i <= N; ++i){
                double x_i = i * h;
                U[0][i] = psiH_1(x_i);
                U[1][i] = psiH_1(x_i) + (tau - e*tau*tau/2)*psiH_2(x_i) + a*tau*tau/2*psiH_1dd(x_i) + b*tau*tau/2*psiH_1d(x_i) + c*tau*tau/2*psiH_1(x_i) + tau*tau/2*fH(x_i, 0);
            }
            break;
        }
        default:
            break;
    }
    for (int k = 1; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            lower[i] = 2*a*tau*tau/h/h - b*tau*tau/h;
            mid[i] = -2 - e*tau - 4*a*tau*tau/h/h + 2*c*tau*tau;
            upper[i] = 2*a*tau*tau/h/h + b*tau*tau/h;
            answer[i] = -4*U[k][i] + U[k-1][i]*(2-e*tau) - 2*tau*tau*fH(i * h, (k+1) * tau);
        }
        //Трехточечная второго порядка
        lower[0] = 0.0;
        mid[0] = (betta - 3*alpha/2/h +lower[1]/upper[1]*alpha/2/h);
        upper[0] = 2*alpha/h + mid[1]/upper[1]*alpha/2/h;
        answer[0] = phiH_0((k+1)*tau) + answer[1]/upper[1]*alpha/2/h;

        lower[N] = -2*gama/h - mid[N-1]/lower[N-1]*gama/2/h;
        mid[N] = delta + 3*gama/2/h - upper[N-1]/lower[N-1]*gama/2/h;
        upper[N] = 0.0;
        answer[N] = phiH_l((k+1)*tau) - answer[N-1]/lower[N-1]*gama/2/h;
        
        U[k+1] = processTridiagonalMatrixH(answer, N+1, lower, mid, upper);
    }
    
    return U;
}


#pragma mark - navigation funcs

-(void)mainMenu:(id)sender{
    [self dismissViewControllerAnimated:YES completion:^{}];
}

-(void)nextScreen:(id)sender{
    [UIView animateWithDuration:0.3 animations:^{
        CGRect frame = bg2.frame;
        frame.origin.x = 0;
        bg2.frame = frame;
    }];
}

-(void)backScreen:(id)sender{
    [UIView animateWithDuration:0.3 animations:^{
        CGRect frame = bg2.frame;
        frame.origin.x = 330;
        bg2.frame = frame;
    }];
}

-(void)proceedScheme:(id)sender{
    int K = [KField.text intValue];
    int N = [NField.text intValue];
    
    l = [lField.text doubleValue];
    double T = [TField.text doubleValue];
    
    double tau = T / K, h = l / N;
    
    a = [aField.text doubleValue], b = [bField.text doubleValue], c = [cField.text doubleValue], e = [eField.text doubleValue];
    alpha = [alphaField.text doubleValue], betta = [bettaField.text doubleValue], gama = [gammaField.text doubleValue], delta = [deltaField.text doubleValue];
    long scheme = [schemePicker selectedRowInComponent:0]; // 0 -- явная, 1 -- неявная, 2 -- нестабильная
    long order = [orderPicker selectedRowInComponent:0];   // 0 -- 2-1, 1 -- 3-2, 2 -- 2-2
    startApproximation = [approxPicker selectedRowInComponent:0]; // 0 -- первого порядка, 2 -- второго порядка
    
    HUD = [[MBProgressHUD alloc] initWithView:bg2];
    [bg2 addSubview:HUD];
    HUD.removeFromSuperViewOnHide = YES;
    
    //Для гиперболических
    float sigma = a*tau/h;
    if(sigma > 1.0 && scheme == 0){
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = [NSString stringWithFormat:@"Не устойчива! σ = %.1f (σ > 1.0)", sigma];
        HUD.detailsLabelText =  @"Измените параметры сетки или схему";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    
    
    HUD.mode = MBProgressHUDAnimationFade;
    HUD.labelText = @"Идет расчет";
    //HUD.detailsLabelText =  @"Измените параметры сетки или схему";
    [HUD show:YES];
    dispatch_async(dispatch_get_global_queue( DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^(void){
        double **U;
       
        switch (scheme) {
            case 0:{
                switch (order) {
                    case 0:{
                        U = hyperbolic_explicitScheme_TwoPoint_FirstOrder(K, N, a, b, c, e, tau, h, alpha, betta, gama, delta);
                        break;
                    }
                    case 1:{
                        U = hyperbolic_explicitScheme_ThreePoint_SecondOrder(K, N, a, b, c, e, tau, h, alpha, betta, gama, delta);
                        break;
                    }
                    case 2:{
                        U = hyperbolic_explicitScheme_TwoPoint_SecondOrder(K, N, a, b, c, e, tau, h, alpha, betta, gama, delta);
                        break;
                    }
                    default:
                        break;
                }
                break;
            }
            case 1:{
                switch (order) {
                    case 0:{
                        U = hyperbolic_implicitScheme_TwoPoint_FirstOrder(K, N, a, b, c, e, tau, h, alpha, betta, gama, delta);
                        break;
                    }
                    case 1:{
                        U = hyperbolic_implicitScheme_ThreePoint_SecondOrder(K, N, a, b, c, e, tau, h, alpha, betta, gama, delta);
                        break;
                    }
                    case 2:{
                        U = hyperbolic_implicitScheme_TwoPoint_SecondOrder(K, N, a, b, c, e, tau, h, alpha, betta, gama, delta);
                        break;
                    }
                    default:
                        break;
                }
                break;
            }
            case 2:{
                switch (order) {
                    case 0:{
                        U = hyperbolic_unstableScheme_TwoPoint_FirstOrder(K, N, a, b, c, e, tau, h, alpha, betta, gama, delta);
                        break;
                    }
                    case 1:{
                        U = hyperbolic_unstableScheme_ThreePoint_SecondOrder(K, N, a, b, c, e, tau, h, alpha, betta, gama, delta);
                        break;
                    }
                    case 2:{
                        U = hyperbolic_unstableScheme_TwoPoint_SecondOrder(K, N, a, b, c, e, tau, h, alpha, betta, gama, delta);
                        break;
                    }
                    default:
                        break;
                }
                break;
            }
            default:
                break;
        }
        
        NSMutableDictionary *dataDict = [[NSMutableDictionary alloc] initWithCapacity:K];
        NSMutableArray *contentArrayErr = [[NSMutableArray alloc] init];
        NSMutableDictionary *dataDictAnalytic = [[NSMutableDictionary alloc] init];
        
        for (int i = 0; i < K; ++i) {
            NSMutableArray *contentArray = [[NSMutableArray alloc] init];
            NSMutableArray *contentArrayAnalytic = [[NSMutableArray alloc] init];
            
            for (int j = 0; j <= N; ++j) {
                id xAnalytic = [NSNumber numberWithDouble:(j)*h];
                id yAnalytic = [NSNumber numberWithDouble:[self functionH:((j)*h) withTime:i*tau]];
                [contentArrayAnalytic addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:xAnalytic, @"x", yAnalytic, @"y", nil]];
            }
            //epsiHlon[k][n] = max_i(u[k][i] - reshenie[k][i])
            id err = [NSNumber numberWithDouble:0.0];
            for (int j = 0; j <= N; ++j) {
                id x = [NSNumber numberWithDouble:j*h];
                id y = [NSNumber numberWithDouble:U[i][j]];
                
                err = [NSNumber numberWithDouble:([err doubleValue] > fabs([self functionH:((j)*h) withTime:i*tau] - U[i][j])) ? [err doubleValue] : fabs([self functionH:((j)*h) withTime:i*tau] - U[i][j])];
                printf("%f	%f\n", [x doubleValue], [y doubleValue]);
                [contentArray addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:x, @"x", y, @"y", nil]];
            }
            id time = [NSNumber numberWithDouble:i*tau];
            printf("K = %f\n", [time doubleValue]);
            [dataDict setObject:contentArray forKey:time];
            [dataDictAnalytic setObject:contentArrayAnalytic forKey:time];
            [contentArrayErr addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:time, @"x", err, @"y", nil]];
        }
        dispatch_async(dispatch_get_main_queue(), ^(void){
            [HUD hide:YES];
            [UIView animateWithDuration:0.3 animations:^{
                CorePlotViewController *viewControllerToPresent = [[CorePlotViewController alloc] initWithNibName:@"CorePlotViewController" bundle:nil];
                viewControllerToPresent.dictForPlot = dataDict;
                viewControllerToPresent.dictForPlotAnalytic = dataDictAnalytic;
                viewControllerToPresent.dictForPlotErr = contentArrayErr;
                [self presentViewController:viewControllerToPresent animated:YES completion:^{}];
                viewControllerToPresent.view.backgroundColor = [UIColor grayColor];
            }];
        });
    });
}

-(void)keyboardShown:(NSNotification *)note{
    [UIView animateWithDuration:0.3 animations:^{
        //CGRect frame = bg.frame;
        //frame.origin.y = 30;
        //bg.frame = frame;
        scrollBg.contentOffset = CGPointMake(0, 100);
    }];
}

-(void)keyboardHidden:(NSNotification *)note{
    [UIView animateWithDuration:0.3 animations:^{
        scrollBg.contentOffset = CGPointMake(0, 0);
        
    }];
}

#pragma mark - viewController delegate
- (void)viewDidLoad
{
    [super viewDidLoad];
    
    system = [[UIImageView alloc] initWithFrame:CGRectMake(32, 22, 256, 134)];
    system.image = [UIImage imageNamed:@"NM_sys2.png"];
    [self.view addSubview:system];
    
    [self.view setBackgroundColor:[UIColor whiteColor]];
    
    bg = [[UIView alloc] initWithFrame:CGRectMake(10, 160, 300, self.view.frame.size.height - 180)];
    bg.backgroundColor = [UIColor colorWithRed:141./255. green:158./255. blue:143./255. alpha:1.0f];
    
    CALayer * imgLayer1 = bg.layer;
    [imgLayer1 setBorderColor: [[UIColor blackColor] CGColor]];
    [imgLayer1 setBorderWidth:0.5f];
    [imgLayer1 setShadowColor: [[UIColor blackColor] CGColor]];
    [imgLayer1 setShadowOpacity:0.9f];
    [imgLayer1 setShadowOffset: CGSizeMake(0, 1)];
    [imgLayer1 setShadowRadius:3.0];
    imgLayer1.shouldRasterize = NO;
    
    scrollBg = [[UIScrollView alloc] initWithFrame:CGRectMake(0, 0, self.view.frame.size.width, self.view.frame.size.height)];
    
    if(self.view.frame.size.height < 500){
        [[NSNotificationCenter defaultCenter] addObserver:self selector:@selector(keyboardShown:) name:UIKeyboardDidShowNotification object:nil];
        [[NSNotificationCenter defaultCenter] addObserver:self selector:@selector(keyboardHidden:) name:UIKeyboardWillHideNotification object:nil];
        
        scrollBg.contentSize = CGSizeMake(320, 600);
    }
    
    [self.view addSubview:scrollBg];
    
    [scrollBg addSubview:bg];
    
    UILabel *aLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 10, 30, 20)];
    aLabel.text = @"a =";
    [bg addSubview:aLabel];
    aField = [[UITextField alloc] initWithFrame:CGRectMake(50, 12, 50, 20)];
    aField.delegate = self;
    aField.text = @"1";
    [bg addSubview:aField];
    
    UILabel *bLabel = [[UILabel alloc] initWithFrame:CGRectMake(130, 10, 30, 20)];
    bLabel.text = @"b =";
    [bg addSubview:bLabel];
    bField = [[UITextField alloc] initWithFrame:CGRectMake(160, 12, 50, 20)];
    bField.delegate = self;
    bField.text = @"1";
    [bg addSubview:bField];
    
    UILabel *cLabel = [[UILabel alloc] initWithFrame:CGRectMake(220, 10, 30, 20)];
    cLabel.text = @"c =";
    [bg addSubview:cLabel];
    cField = [[UITextField alloc] initWithFrame:CGRectMake(250, 12, 50, 20)];
    cField.delegate = self;
    cField.text = @"-1";
    [bg addSubview:cField];
    
    UILabel *eLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 30, 30, 20)];
    eLabel.text = @"e =";
    [bg addSubview:eLabel];
    eField = [[UITextField alloc] initWithFrame:CGRectMake(50, 32, 50, 20)];
    eField.delegate = self;
    eField.text = @"3";
    [bg addSubview:eField];
    
    UILabel *lLabel = [[UILabel alloc] initWithFrame:CGRectMake(130, 30, 30, 20)];
    lLabel.text = @"l =";
    [bg addSubview:lLabel];
    lField = [[UITextField alloc] initWithFrame:CGRectMake(160, 32, 50, 20)];
    lField.delegate = self;
    lField.text = @"3.14";
    [bg addSubview:lField];
    
    UILabel *alphaLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 60, 60, 20)];
    alphaLabel.text = @"α =";
    [bg addSubview:alphaLabel];
    alphaField = [[UITextField alloc] initWithFrame:CGRectMake(80, 62, 50, 20)];
    alphaField.delegate = self;
    alphaField.text = @"0";
    [bg addSubview:alphaField];
    
    UILabel *bettaLabel = [[UILabel alloc] initWithFrame:CGRectMake(150, 60, 60, 20)];
    bettaLabel.text = @"β =";
    [bg addSubview:bettaLabel];
    bettaField = [[UITextField alloc] initWithFrame:CGRectMake(210, 62, 50, 20)];
    bettaField.delegate = self;
    bettaField.text = @"1";
    [bg addSubview:bettaField];
    
    UILabel *gamaLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 80, 60, 20)];
    gamaLabel.text = @"γ =";
    [bg addSubview:gamaLabel];
    gammaField = [[UITextField alloc] initWithFrame:CGRectMake(80, 82, 50, 20)];
    gammaField.delegate = self;
    gammaField.text = @"0";
    [bg addSubview:gammaField];
    
    UILabel *deltaLabel = [[UILabel alloc] initWithFrame:CGRectMake(150, 80, 60, 20)];
    deltaLabel.text = @"δ =";
    [bg addSubview:deltaLabel];
    deltaField = [[UITextField alloc] initWithFrame:CGRectMake(210, 82, 50, 20)];
    deltaField.delegate = self;
    deltaField.text = @"1";
    [bg addSubview:deltaField];
    
    UILabel *TLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 130, 200, 20)];
    TLabel.text = @"Финальное время (T):";
    [bg addSubview:TLabel];
    TField = [[UITextField alloc] initWithFrame:CGRectMake(250, 132, 50, 20)];
    TField.delegate = self;
    TField.text = @"10";
    [bg addSubview:TField];
    
    UILabel *KLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 150, 230, 20)];
    KLabel.text = @"Разбиений по времени (K):";
    [bg addSubview:KLabel];
    KField = [[UITextField alloc] initWithFrame:CGRectMake(250, 150, 50, 20)];
    KField.delegate = self;
    KField.text = @"2000";
    [bg addSubview:KField];
    
    UILabel *NLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 170, 230, 20)];
    NLabel.text = @"Разбиений по OX (N):";
    [bg addSubview:NLabel];
    NField = [[UITextField alloc] initWithFrame:CGRectMake(250, 170, 50, 20)];
    NField.delegate = self;
    NField.text = @"30";
    [bg addSubview:NField];
    
    UIButton *next;
    
    if (self.view.bounds.size.height < 500) {
        next = [[UIButton alloc] initWithFrame:CGRectMake(230, 260, 60, 40)];
    }
    else {
        next = [[UIButton alloc] initWithFrame:CGRectMake(230, 320, 60, 40)];
    }
    
    UILabel *blab = [[UILabel alloc] initWithFrame:CGRectMake(0, 0, 60, 20)];
    blab.text = @"Далее";
    [next addSubview:blab];
    [next addTarget:self action:@selector(nextScreen:) forControlEvents:UIControlEventTouchUpInside];
    
    [bg addSubview:next];
    
    UIButton *mainMenu;
    if ([[UIScreen mainScreen] bounds].size.height < 500) {
        mainMenu = [[UIButton alloc] initWithFrame:CGRectMake(40, 260, 60, 40)];
    }
    else {
        mainMenu = [[UIButton alloc] initWithFrame:CGRectMake(40, 320, 60, 40)];
    }
    
    UILabel *blab4 = [[UILabel alloc] initWithFrame:CGRectMake(0, 0, 60, 20)];
    blab4.text = @"Назад";
    [mainMenu addSubview:blab4];
    [mainMenu addTarget:self action:@selector(mainMenu:) forControlEvents:UIControlEventTouchUpInside];
    [bg addSubview:mainMenu];
    
    //bg 2
    bg2 = [[UIView alloc] initWithFrame:CGRectMake(330, 22, 320, self.view.frame.size.height - 42)];
    bg2.backgroundColor = [UIColor colorWithRed:136./255. green:135./255. blue:148./255. alpha:1.0];
    
    CALayer * imgLayer2 = bg2.layer;
    [imgLayer2 setBorderColor: [[UIColor blackColor] CGColor]];
    [imgLayer2 setBorderWidth:0.5f];
    [imgLayer2 setShadowColor: [[UIColor blackColor] CGColor]];
    [imgLayer2 setShadowOpacity:0.9f];
    [imgLayer2 setShadowOffset: CGSizeMake(0, 1)];
    [imgLayer2 setShadowRadius:3.0];
    imgLayer2.shouldRasterize = NO;
    
    [self.view addSubview:bg2];
    
    UILabel *order = [[UILabel alloc] initWithFrame:CGRectMake(0, 5, 320, 15)];
    order.text = @"Аппроксимация:";
    order.textAlignment = NSTextAlignmentCenter;
    
    UILabel *scheme = [[UILabel alloc] initWithFrame:CGRectMake(0, 160, 320, 15)];
    scheme.text = @"Схема:";
    scheme.textAlignment = NSTextAlignmentCenter;
    
    UILabel *approximation = [[UILabel alloc] initWithFrame:CGRectMake(0, 310, 320, 15)];
    approximation.text = @"Аппрокс. начальных значений:";
    approximation.textAlignment = NSTextAlignmentCenter;
    
    orderPicker = [[UIPIckerViewOrder alloc] initWithFrame:CGRectMake(0, 10, 320, 60)];
    schemePicker = [[UIPickerViewSchemeHyperbolic alloc] initWithFrame:CGRectMake(0, 160, 320, 60)];
    approxPicker = [[UIPickerViewApproximation alloc] initWithFrame:CGRectMake(0, 285, 320, 40)];
    
    [bg2 addSubview:order];
    [bg2 addSubview:scheme];
    [bg2 addSubview:approximation];
    [bg2 addSubview:orderPicker];
    [bg2 addSubview:schemePicker];
    [bg2 addSubview:approxPicker];
    
    UIButton *nextButton;
    
    if ([[UIScreen mainScreen] bounds].size.height < 500) {
        nextButton = [[UIButton alloc] initWithFrame:CGRectMake(230, 390, 60, 40)];
    }
    else {
        nextButton = [[UIButton alloc] initWithFrame:CGRectMake(230, 450, 60, 40)];
    }
    
    
    UILabel *blab2 = [[UILabel alloc] initWithFrame:CGRectMake(0, 0, 60, 20)];
    blab2.text = @"Далее";
    [nextButton addSubview:blab2];
    [nextButton addTarget:self action:@selector(proceedScheme:) forControlEvents:UIControlEventTouchUpInside];
    
    [bg2 addSubview:nextButton];
    
    
    UIButton *backButton;
    if ([[UIScreen mainScreen] bounds].size.height < 500) {
        backButton = [[UIButton alloc] initWithFrame:CGRectMake(40, 390, 60, 40)];
    }
    else {
        backButton = [[UIButton alloc] initWithFrame:CGRectMake(40, 450, 60, 40)];
    }
    
    UILabel *blab3 = [[UILabel alloc] initWithFrame:CGRectMake(0, 0, 60, 20)];
    blab3.text = @"Назад";
    [backButton addSubview:blab3];
    [backButton addTarget:self action:@selector(backScreen:) forControlEvents:UIControlEventTouchUpInside];
    [bg2 addSubview:backButton];
}

- (void)didReceiveMemoryWarning
{
    [super didReceiveMemoryWarning];
}

#pragma mark - TextField delegate

-(BOOL) textFieldShouldReturn:(UITextField *)textField{
    [textField resignFirstResponder];
    return YES;
}

- (BOOL)textField:(UITextField *)textField shouldChangeCharactersInRange:(NSRange)range replacementString:(NSString *)string  {
    NSUInteger newLength = [textField.text length] + [string length] - range.length;
    NSCharacterSet *cs = [[NSCharacterSet characterSetWithCharactersInString:NUMBERS_ONLY] invertedSet];
    NSString *filtered = [[string componentsSeparatedByCharactersInSet:cs] componentsJoinedByString:@""];
    return (([string isEqualToString:filtered])&&(newLength <= CHARACTER_LIMIT));
}


@end
