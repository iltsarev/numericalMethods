//
//  ellipticViewController.m
//  numericalMethods
//
//  Created by Ilya on 23.11.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import "ellipticViewController.h"
#import "CorePlotViewController.h"

#define ALPHABET2 @"qwertyuiopasdfghjklzxcvbnm1234567890-+*/()$"


double a,b,bx,by,c,e,alpha,betta,gama,delta,l,lx,ly;
long startApproximation;

@interface ellipticViewController ()

@end

@implementation ellipticViewController
@synthesize HUD;
@synthesize system;

double * processTridiagonalMatrixE(double *x, const size_t N, const double *a, const double *b, double *c) {
    
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
    //return exp(-t)*cos(x);
    NSArray *values = @[[NSNumber numberWithDouble:x], [NSNumber numberWithDouble:t], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"t", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionRealFunc evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) phiH_0:(double)t{
    //    return exp(-t);
    NSArray *values = @[ [NSNumber numberWithDouble:t], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"t", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPhi_0 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}
-(double) phiH_l:(double)t{
    //return -exp(-t);
    NSArray *values = @[ [NSNumber numberWithDouble:t], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"t", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPhi_l evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
    
}

-(double) psiH_1:(double)x{
    // return cos(x);
    NSArray *values = @[ [NSNumber numberWithDouble:x], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPsi_1 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
    
}

-(double) psiH_1d:(double) x{
    //return -sin(x);
    NSArray *values = @[ [NSNumber numberWithDouble:x], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPsi_1d evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) psiH_1dd:(double) x{
    //return -cos(x);
    NSArray *values = @[ [NSNumber numberWithDouble:x], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPsi_1dd evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) psiH_2:(double) x{
    //return -cos(x);
    NSArray *values = @[ [NSNumber numberWithDouble:x], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPsi_2 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}


-(double) fH:(double) x :(double) t{
    //return sin(x)*exp(-t);
    NSArray *values = @[[NSNumber numberWithDouble:x], [NSNumber numberWithDouble:t], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"t", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionF evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

#pragma mark - elliptic Liebmann
-(double **)elliptic_liebmann:(int) K :(int) Nx :(int) Ny :(double) a :(double) b :(double) c :(double) e :(double) tau :(double) hx :(double) hy :(double) alpha1 :(double) betta1 :(double) alpha2 :(double) betta2 :(double) alpha3 :(double) betta3 :(double) alpha4 :(double) betta4{
    
    double alphax, alphay, bettax, bettay, gammax,gammay, deltax, deltay;
    
    double ***U = (double ***)malloc(K * sizeof(double **));
    for (int i = 0; i < K; i++){
        U[i] = (double **)malloc((Nx+1) * sizeof(double *));
        for (int j = 0; j < Nx+1; j++)
            U[i][j] = (double *)malloc((Ny+1) * sizeof(double));
    }
    
    for (int i = 0; i < Nx+1; i++){
        for (int j = 0; j < Ny+1; j++){
            U[0][i][j] = 0;
        }
    }
    
    for (int k = 1; k < K - 1; ++k) {
        for (int i = 1; i < Nx; ++i) {
            for (int j = 1; j < Ny; ++j) {
                U[k][i][j] = (1/(2*hx*hx - c*hx*hx*hy*hy +2*hy*hy)) * (U[k-1][i+1][j]*hy*hy*(1 + bx*hx/2) + U[k][i-1][j]*hy*hy*(1-bx*hx/2) + U[k-1][i][j+1]*hx*hx*(1+by*hy/2) + U[k][i][j-1]*hx*hx*(1 - by*hy/2) + hx*hx*hy*hy*[self fH:i*hx :j*hy]);
            }
        }
    }
    
//    U[0][j] = (alphax/2/hx)/((bettax-3*alphax)/2/hx) * (U[2][j] = 4*U[1][j]) + phix0(y_j)/(bettax - 3*alphax/2/hx);
//    U[Nx][j] = (gammax/2/hx)/(deltax + 3*gammax/2/hx) * (4*U[nx-1][j] - U[nx-2][j]) + phixl(y_j)/(deltax + 3*gammax/2/hx);
//    U[i][0] = (alphay/2/hy)/(bettay - 3*alphay/2/hy) * (U[i][2] - 4*U[i][1]) + phiy0(x_i)/(bettay - 3*alphay/2/hy)
//    U[i][ny] = (gammay/2/hy)/(deltay + 3*gammay/2/hy) * (4*U[i][ny-1] - U[i][ny-2]) + phiyl(x_i)/(deltay + 3*gammay/2/hy);
    
    //relaxation
//    U[k][i][j] = w*U[p][i][j] + (1-w)*U[k-1][i][j];
    
#warning U[last]
    return U[0];
}

#pragma mark - elliptic Seidel
-(double ***)elliptic_seidel:(int) K :(int) Nx :(int) Ny :(double) a :(double) b :(double) c :(double) e :(double) tau :(double) hx :(double) hy :(double) alpha :(double) betta :(double) gama :(double) delta{
    
    double ***U = (double ***)malloc(K * sizeof(double **));
    for (int i = 0; i < K; i++){
        U[i] = (double **)malloc((Nx+1) * sizeof(double *));
        for (int j = 0; j < Nx+1; j++){
            U[i][j] = (double *)malloc((Ny+1) * sizeof(double));
        }
    }
    
    for (int i = 0; i < Nx+1; i++){
        for (int j = 0; j < Ny+1; j++){
            U[0][i][j] = 0;
        }
    }
    
    for (int k = 1; k < K - 1; ++k) {
        for (int i = 1; i < Nx; ++i) {
            for (int j = 1; j < Ny; ++j) {
                U[k][i][j] = (1/(2*hx*hx - c*hx*hx*hy*hy + 2*hy*hy)) * (U[k-1][i+1][j]*hy*hy*(1+bx*hx/2) + U[k][i-1][j]*hy*hy*(1-bx*hx/2) + U[k-1][i][j+1]*hx*hx*(1+by*hy/2) + U[k][i][j-1]*(1-by*hy/2) + hx*hx*hy*hy*[self fH:i*hx :j*hy]);
            }
        }
    }
    
    return U;
}

//#pragma mark - hyperbolic_explicitScheme
//
//-(double **)hyperbolic_explicitScheme_TwoPoint_FirstOrder:(int) K :(int) N :(double) a :(double) b :(double) c :(double) e :(double) tau :(double) h :(double) alpha :(double) betta :(double) gama :(double) delta{
//    double **U = (double **)malloc(K * sizeof(double *));
//    
//    for (int i = 0; i < K; i++)
//        U[i] = (double *)malloc((N+1) * sizeof(double));
//    
//    
//    switch (startApproximation) {
//        case 0:
//        {
//            //Первый порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + tau*[self psiH_2:x_i];
//            }
//            break;
//        }
//        case 1:
//        {
//            //Второй порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + (tau - e*tau*tau/2)*[self psiH_2:x_i] + a*tau*tau/2*[self psiH_1dd:x_i] + b*tau*tau/2*[self psiH_1d:x_i] + c*tau*tau/2*[self psiH_1:x_i] + tau*tau/2*[self fH:x_i :0];
//            }
//            break;
//        }
//        default:
//            break;
//    }
//    
//    for (int k = 1; k < K - 1; ++k) {
//        for (int i = 1; i < N; ++i) {
//            U[k+1][i] = 1/(e*tau+2)*(2*a*tau*tau/h/h*(U[k][i+1] - 2*U[k][i] + U[k][i-1]) + b*tau*tau/h*(U[k][i+1] - U[k][i-1]) + U[k][i]*(2*c*tau*tau + 4) + U[k-1][i]*(e*tau -2) + 2*tau*tau*[self fH:i*h :k*tau]);
//        }
//        //  двухточечная первого
//        U[k+1][0] = -alpha/h/(betta - alpha/h)*U[k+1][1]  + [self phiH_0:(k+1)*tau]/(betta - alpha/h);
//        U[k+1][N] = gama/h/(delta + gama/h)*U[k+1][N-1] + [self phiH_l:(k+1)*tau]/(delta + gama/h);
//    }
//    
//    return U;
//}
//
//-(double **)hyperbolic_explicitScheme_TwoPoint_SecondOrder:(int) K :(int) N :(double) a :(double) b :(double) c :(double) e :(double) tau :(double) h :(double) alpha :(double) betta :(double) gama :(double) delta{
//    double **U = (double **)malloc(K * sizeof(double *));
//    
//    for (int i = 0; i < K; i++)
//        U[i] = (double *)malloc((N+1) * sizeof(double));
//    
//    
//    switch (startApproximation) {
//        case 0:
//        {
//            //Первый порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + tau*[self psiH_2:x_i];
//            }
//            break;
//        }
//        case 1:
//        {
//            //Второй порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + (tau - e*tau*tau/2)*[self psiH_2:x_i] + a*tau*tau/2*[self psiH_1dd:x_i] + b*tau*tau/2*[self psiH_1d:x_i] + c*tau*tau/2*[self psiH_1:x_i] + tau*tau/2*[self fH:x_i :0];
//            }
//            break;
//        }
//        default:
//            break;
//    }
//    
//    for (int k = 1; k < K - 1; ++k) {
//        for (int i = 1; i < N; ++i) {
//            U[k+1][i] = 1/(e*tau+2)*(2*a*tau*tau/h/h*(U[k][i+1] - 2*U[k][i] + U[k][i-1]) + b*tau*tau/h*(U[k][i+1] - U[k][i-1]) + U[k][i]*(2*c*tau*tau + 4) + U[k-1][i]*(e*tau -2) + 2*tau*tau*[self fH:i*h :k*tau]);
//        }
//        //  Двухточечная второго
//        U[k+1][0] = (2*h*alpha/tau/tau * U[k][0] - h*alpha/tau/tau*U[k-1][0] + h*e*alpha/2/tau*U[k-1][0] + h*alpha*[self fH:0 :(k+1)*tau] - [self phiH_0:(k+1)*tau]*(2*a-b*h) + 2*a*alpha/h*U[k+1][1]) / (2*a*alpha/h + h*alpha/tau/tau + h*e*alpha/2/tau - c*h*alpha - betta*(2*a-b*h));
//        U[k+1][N] = (2*h*gama/tau/tau*U[k][N] - h*gama/tau/tau*U[k-1][N] + h*e*gama/2/tau*U[k-1][N] + h*gama*[self fH:N :(k+1)*tau] + [self phiH_l:(k+1)*tau]*(2*a + b*h)  + 2*a*gama/h*U[k+1][N-1]) / (2*a*gama/h + h*gama/tau/tau + h*e*gama/2/tau - c*h*gama + delta*(2*a+b*h));
//    }
//    
//    return U;
//}
//
//
//-(double **)hyperbolic_explicitScheme_ThreePoint_SecondOrder:(int) K :(int) N :(double) a :(double) b :(double) c :(double) e :(double) tau :(double) h :(double) alpha :(double) betta :(double) gama :(double) delta{
//    double **U = (double **)malloc(K * sizeof(double *));
//    
//    for (int i = 0; i < K; i++)
//        U[i] = (double *)malloc((N+1) * sizeof(double));
//    
//    switch (startApproximation) {
//        case 0:
//        {
//            //Первый порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + tau*[self psiH_2:x_i];
//            }
//            break;
//        }
//        case 1:
//        {
//            //Второй порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + (tau - e*tau*tau/2)*[self psiH_2:x_i] + a*tau*tau/2*[self psiH_1dd:x_i] + b*tau*tau/2*[self psiH_1d:x_i] + c*tau*tau/2*[self psiH_1:x_i] + tau*tau/2*[self fH:x_i :0];
//            }
//            break;
//        }
//        default:
//            break;
//    }
//    
//    for (int k = 1; k < K - 1; ++k) {
//        for (int i = 1; i < N; ++i) {
//            U[k+1][i] = 1/(e*tau+2)*(2*a*tau*tau/h/h*(U[k][i+1] - 2*U[k][i] + U[k][i-1]) + b*tau*tau/h*(U[k][i+1] - U[k][i-1]) + U[k][i]*(2*c*tau*tau + 4) + U[k-1][i]*(e*tau -2) + 2*tau*tau*[self fH:i*h :k*tau]);
//        }
//        //  трехточечная второго
//        U[k+1][0] = alpha/2/h/(betta - 3*alpha/2/h)*(U[k+1][2] - 4*U[k+1][1]) + [self phiH_0:(k+1)*tau]/(betta - 3*alpha/2/h);
//        U[k+1][N] = gama/2/h/(delta - 3*gama/2/h)*(4*U[k+1][N-1] - U[k+1][N-2]) + [self phiH_l:(k+1)*tau]/(delta - 3*gama/2/h);
//    }
//    
//    return U;
//}
//
//
//#pragma mark - hyperbolic_unstableScheme
//
//-(double **)hyperbolic_unstableScheme_TwoPoint_FirstOrder:(int) K :(int) N :(double) a :(double) b :(double) c :(double) e :(double) tau :(double) h :(double) alpha :(double) betta :(double) gama :(double) delta{
//    double **U = (double **)malloc(K * sizeof(double *));
//    
//    for (int i = 0; i < K; i++)
//        U[i] = (double *)malloc((N+1) * sizeof(double));
//    
//    switch (startApproximation) {
//        case 0:
//        {
//            //Первый порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + tau*[self psiH_2:x_i];
//            }
//            break;
//        }
//        case 1:
//        {
//            //Второй порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + (tau - e*tau*tau/2)*[self psiH_2:x_i] + a*tau*tau/2*[self psiH_1dd:x_i] + b*tau*tau/2*[self psiH_1d:x_i] + c*tau*tau/2*[self psiH_1:x_i] + tau*tau/2*[self fH:x_i :0];
//            }
//            break;
//        }
//        default:
//            break;
//    }
//    
//    for (int k = 1; k < K - 1; ++k) {
//        for (int i = 1; i < N; ++i) {
//            U[k+1][i] = 1/(e*tau+2)*(2*a*tau*tau/h/h*(U[k-1][i+1]- 2*U[k-1][i] + U[k-1][i-1]) + b*tau*tau/h*(U[k-1][i+1] - U[k-1][i-1]) + U[k-1][i]*(2*c*tau*tau + e*tau - 2) + 4*U[k][i] + 2*tau*tau*[self fH:i*h :(k-1)*tau]);
//        }
//        //  двухточечная первого
//        U[k+1][0] = -alpha/h/(betta - alpha/h)*U[k+1][1]  + [self phiH_0:(k+1)*tau]/(betta - alpha/h);
//        U[k+1][N] = gama/h/(delta + gama/h)*U[k+1][N-1] + [self phiH_l:(k+1)*tau]/(delta + gama/h);
//    }
//    
//    return U;
//}
//
//-(double **)hyperbolic_unstableScheme_TwoPoint_SecondOrder:(int) K :(int) N :(double) a :(double) b :(double) c :(double) e :(double) tau :(double) h :(double) alpha :(double) betta :(double) gama :(double) delta{
//    double **U = (double **)malloc(K * sizeof(double *));
//    
//    for (int i = 0; i < K; i++)
//        U[i] = (double *)malloc((N+1) * sizeof(double));
//    
//    switch (startApproximation) {
//        case 0:
//        {
//            //Первый порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + tau*[self psiH_2:x_i];
//            }
//            break;
//        }
//        case 1:
//        {
//            //Второй порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + (tau - e*tau*tau/2)*[self psiH_2:x_i] + a*tau*tau/2*[self psiH_1dd:x_i] + b*tau*tau/2*[self psiH_1d:x_i] + c*tau*tau/2*[self psiH_1:x_i] + tau*tau/2*[self fH:x_i :0];
//            }
//            break;
//        }
//        default:
//            break;
//    }
//    
//    for (int k = 1; k < K - 1; ++k) {
//        for (int i = 1; i < N; ++i) {
//            U[k+1][i] = 1/(e*tau+2)*(2*a*tau*tau/h/h*(U[k-1][i+1]- 2*U[k-1][i] + U[k-1][i-1]) + b*tau*tau/h*(U[k-1][i+1] - U[k-1][i-1]) + U[k-1][i]*(2*c*tau*tau + e*tau - 2) + 4*U[k][i] + 2*tau*tau*[self fH:i*h :(k-1)*tau]);
//        }
//        //  Двухточечная второго
//        U[k+1][0] = (2*h*alpha/tau/tau * U[k][0] - h*alpha/tau/tau*U[k-1][0] + h*e*alpha/2/tau*U[k-1][0] + h*alpha*[self fH:0 :(k+1)*tau] - [self phiH_0:(k+1)*tau]*(2*a-b*h) + 2*a*alpha/h*U[k+1][1]) / (2*a*alpha/h + h*alpha/tau/tau + h*e*alpha/2/tau - c*h*alpha - betta*(2*a-b*h));
//        U[k+1][N] = (2*h*gama/tau/tau*U[k][N] - h*gama/tau/tau*U[k-1][N] + h*e*gama/2/tau*U[k-1][N] + h*gama*[self fH:N :(k+1)*tau] + [self phiH_l:(k+1)*tau]*(2*a + b*h)  + 2*a*gama/h*U[k+1][N-1]) / (2*a*gama/h + h*gama/tau/tau + h*e*gama/2/tau - c*h*gama + delta*(2*a+b*h));
//    }
//    
//    return U;
//}
//
//
//-(double **)hyperbolic_unstableScheme_ThreePoint_SecondOrder:(int) K :(int) N :(double) a :(double) b :(double) c :(double) e :(double) tau :(double) h :(double) alpha :(double) betta :(double) gama :(double) delta{
//    
//    double **U = (double **)malloc(K * sizeof(double *));
//    
//    for (int i = 0; i < K; i++)
//        U[i] = (double *)malloc((N+1) * sizeof(double));
//    
//    switch (startApproximation) {
//        case 0:
//        {
//            //Первый порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + tau*[self psiH_2:x_i];
//            }
//            break;
//        }
//        case 1:
//        {
//            //Второй порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + (tau - e*tau*tau/2)*[self psiH_2:x_i] + a*tau*tau/2*[self psiH_1dd:x_i] + b*tau*tau/2*[self psiH_1d:x_i] + c*tau*tau/2*[self psiH_1:x_i] + tau*tau/2*[self fH:x_i :0];
//            }
//            break;
//        }
//        default:
//            break;
//    }
//    
//    
//    for (int k = 1; k < K - 1; ++k) {
//        for (int i = 1; i < N; ++i) {
//            U[k+1][i] = 1/(e*tau+2)*(2*a*tau*tau/h/h*(U[k-1][i+1]- 2*U[k-1][i] + U[k-1][i-1]) + b*tau*tau/h*(U[k-1][i+1] - U[k-1][i-1]) + U[k-1][i]*(2*c*tau*tau + e*tau - 2) + 4*U[k][i] + 2*tau*tau*[self fH:i*h :(k-1)*tau]);
//        }
//        //  трехточечная второго
//        U[k+1][0] = alpha/2/h/(betta - 3*alpha/2/h)*(U[k+1][2] - 4*U[k+1][1]) + [self phiH_0:(k+1)*tau]/(betta - 3*alpha/2/h);
//        U[k+1][N] = gama/2/h/(delta - 3*gama/2/h)*(4*U[k+1][N-1] - U[k+1][N-2]) + [self phiH_l:(k+1)*tau]/(delta - 3*gama/2/h);
//    }
//    
//    return U;
//}
//
//
//#pragma mark - hyperbolic_implicitScheme
//
//-(double **)hyperbolic_implicitScheme_TwoPoint_FirstOrder:(int) K :(int) N :(double) a :(double) b :(double) c :(double) e :(double) tau :(double) h :(double) alpha :(double) betta :(double) gama :(double) delta{
//    
//    double **U = (double **)malloc(K * sizeof(double *));
//    
//    double *lower = (double *)malloc((N+1) * sizeof(double));
//    double *mid = (double *)malloc((N+1) * sizeof(double));
//    double *upper = (double *)malloc((N+1) * sizeof(double));
//    double *answer = (double *)malloc((N+1) * sizeof(double));
//    
//    memset(lower, 0, (N+1) * sizeof(double));
//    memset(upper, 0, (N+1) * sizeof(double));
//    memset(mid, 0, (N+1) * sizeof(double));
//    memset(answer, 0, (N+1) * sizeof(double));
//    
//    
//    for (int i = 0; i < K; i++)
//        U[i] = (double *)malloc((N+1) * sizeof(double));
//    
//    switch (startApproximation) {
//        case 0:
//        {
//            //Первый порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + tau*[self psiH_2:x_i];
//            }
//            break;
//        }
//        case 1:
//        {
//            //Второй порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + (tau - e*tau*tau/2)*[self psiH_2:x_i] + a*tau*tau/2*[self psiH_1dd:x_i] + b*tau*tau/2*[self psiH_1d:x_i] + c*tau*tau/2*[self psiH_1:x_i] + tau*tau/2*[self fH:x_i :0];
//            }
//            break;
//        }
//        default:
//            break;
//    }
//    for (int k = 1; k < K - 1; ++k) {
//        for (int i = 1; i < N; ++i) {
//            lower[i] = 2*a*tau*tau/h/h - b*tau*tau/h;
//            mid[i] = -2 - e*tau - 4*a*tau*tau/h/h + 2*c*tau*tau;
//            upper[i] = 2*a*tau*tau/h/h + b*tau*tau/h;
//            answer[i] = -4*U[k][i] + U[k-1][i]*(2-e*tau) - 2*tau*tau*[self fH:i*h :(k+1)*tau];
//        }
//        //Двухточечная первого порядка
//        lower[0] = 0.0;
//        mid[0] = (betta - alpha/h);
//        upper[0] = alpha/h;
//        answer[0] = [self phiH_0:(k+1)*tau];
//        
//        lower[N] = -gama/h;
//        mid[N] = delta + gama/h;
//        upper[N] = 0.0;
//        answer[N] = [self phiH_l:(k+1)*tau];
//        
//        U[k+1] = processTridiagonalMatrixE(answer, N+1, lower, mid, upper);
//    }
//    free(lower);
//    free(mid);
//    free(upper);
//    return U;
//}
//
//-(double **)hyperbolic_implicitScheme_TwoPoint_SecondOrder:(int) K :(int) N :(double) a :(double) b :(double) c :(double) e :(double) tau :(double) h :(double) alpha :(double) betta :(double) gama :(double) delta{
//    
//    double **U = (double **)malloc(K * sizeof(double *));
//    
//    double *lower = (double *)malloc((N+1) * sizeof(double));
//    double *mid = (double *)malloc((N+1) * sizeof(double));
//    double *upper = (double *)malloc((N+1) * sizeof(double));
//    double *answer = (double *)malloc((N+1) * sizeof(double));
//    
//    memset(lower, 0, (N+1) * sizeof(double));
//    memset(upper, 0, (N+1) * sizeof(double));
//    memset(mid, 0, (N+1) * sizeof(double));
//    memset(answer, 0, (N+1) * sizeof(double));
//    
//    
//    for (int i = 0; i < K; i++)
//        U[i] = (double *)malloc((N+1) * sizeof(double));
//    
//    switch (startApproximation) {
//        case 0:
//        {
//            //Первый порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + tau*[self psiH_2:x_i];
//            }
//            break;
//        }
//        case 1:
//        {
//            //Второй порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + (tau - e*tau*tau/2)*[self psiH_2:x_i] + a*tau*tau/2*[self psiH_1dd:x_i] + b*tau*tau/2*[self psiH_1d:x_i] + c*tau*tau/2*[self psiH_1:x_i] + tau*tau/2*[self fH:x_i :0];
//            }
//            break;
//        }
//        default:
//            break;
//    }
//    for (int k = 1; k < K - 1; ++k) {
//        for (int i = 1; i < N; ++i) {
//            lower[i] = 2*a*tau*tau/h/h - b*tau*tau/h;
//            mid[i] = -2 - e*tau - 4*a*tau*tau/h/h + 2*c*tau*tau;
//            upper[i] = 2*a*tau*tau/h/h + b*tau*tau/h;
//            answer[i] = -4*U[k][i] + U[k-1][i]*(2-e*tau) - 2*tau*tau*[self fH:i*h :(k+1)*tau];
//        }
//        //Двухточечная второго порядка
//        lower[0] = 0.0;
//        mid[0] = 2*a*alpha/h + h*alpha/tau/tau + h*e*alpha/2/tau - c*h*alpha - betta*(2*a-b*h);
//        upper[0] = -2*a*alpha/h;
//        answer[0] = 2*h*alpha/tau/tau * U[k][0] - h*alpha/tau/tau*U[k-1][0] + h*e*alpha/2/tau*U[k-1][0] + h*alpha*[self fH:0 :(k+1)*tau] - [self phiH_0:(k+1)*tau]*(2*a-b*h);
//        
//        lower[N] = -2*a*gama/h;
//        mid[N] = 2*a*gama/h + h*gama/tau/tau + h*e*gama/2/tau - c*h*gama + delta*(2*a+b*h);
//        upper[N] = 0.0;
//        answer[N] = 2*h*gama/tau/tau*U[k][N] - h*gama/tau/tau*U[k-1][N] + h*e*gama/2/tau*U[k-1][N] + h*gama*[self fH:N :(k+1)*tau] + [self phiH_l:(k+1)*tau]*(2*a + b*h);
//        
//        
//        U[k+1] = processTridiagonalMatrixE(answer, N+1, lower, mid, upper);
//    }
//    free(lower);
//    free(mid);
//    free(upper);
//    return U;
//}
//
//
//-(double **)hyperbolic_implicitScheme_ThreePoint_SecondOrder:(int) K :(int) N :(double) a :(double) b :(double) c :(double) e :(double) tau :(double) h :(double) alpha :(double) betta :(double) gama :(double) delta{
//    
//    double **U = (double **)malloc(K * sizeof(double *));
//    
//    double *lower = (double *)malloc((N+1) * sizeof(double));
//    double *mid = (double *)malloc((N+1) * sizeof(double));
//    double *upper = (double *)malloc((N+1) * sizeof(double));
//    double *answer = (double *)malloc((N+1) * sizeof(double));
//    
//    memset(lower, 0, (N+1) * sizeof(double));
//    memset(upper, 0, (N+1) * sizeof(double));
//    memset(mid, 0, (N+1) * sizeof(double));
//    memset(answer, 0, (N+1) * sizeof(double));
//    
//    
//    for (int i = 0; i < K; i++)
//        U[i] = (double *)malloc((N+1) * sizeof(double));
//    
//    switch (startApproximation) {
//        case 0:
//        {
//            //Первый порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + tau*[self psiH_2:x_i];
//            }
//            break;
//        }
//        case 1:
//        {
//            //Второй порядок
//            for (int i = 0; i <= N; ++i){
//                double x_i = i * h;
//                U[0][i] = [self psiH_1:x_i];
//                U[1][i] = [self psiH_1:x_i] + (tau - e*tau*tau/2)*[self psiH_2:x_i] + a*tau*tau/2*[self psiH_1dd:x_i] + b*tau*tau/2*[self psiH_1d:x_i] + c*tau*tau/2*[self psiH_1:x_i] + tau*tau/2*[self fH:x_i :0];
//            }
//            break;
//        }
//        default:
//            break;
//    }
//    
//    for (int k = 1; k < K - 1; ++k) {
//        for (int i = 1; i < N; ++i) {
//            lower[i] = 2*a*tau*tau/h/h - b*tau*tau/h;
//            mid[i] = -2 - e*tau - 4*a*tau*tau/h/h + 2*c*tau*tau;
//            upper[i] = 2*a*tau*tau/h/h + b*tau*tau/h;
//            answer[i] = -4*U[k][i] + U[k-1][i]*(2-e*tau) - 2*tau*tau*[self fH:i*h :(k+1)*tau];
//        }
//        //Трехточечная второго порядка
//        lower[0] = 0.0;
//        mid[0] = (betta - 3*alpha/2/h +lower[1]/upper[1]*alpha/2/h);
//        upper[0] = 2*alpha/h + mid[1]/upper[1]*alpha/2/h;
//        answer[0] = [self phiH_0:(k+1)*tau] + answer[1]/upper[1]*alpha/2/h;
//        
//        lower[N] = -2*gama/h - mid[N-1]/lower[N-1]*gama/2/h;
//        mid[N] = delta + 3*gama/2/h - upper[N-1]/lower[N-1]*gama/2/h;
//        upper[N] = 0.0;
//        answer[N] = [self phiH_l:(k+1)*tau] - answer[N-1]/lower[N-1]*gama/2/h;
//        
//        U[k+1] = processTridiagonalMatrixE(answer, N+1, lower, mid, upper);
//    }
//    free(lower);
//    free(mid);
//    free(upper);
//    return U;
//}


#pragma mark - navigation funcs

-(void)mainMenu:(id)sender{
    [self dismissViewControllerAnimated:YES completion:^{}];
}

-(void)nextScreenFuncs:(id)sender{
    [UIView animateWithDuration:0.3 animations:^{
        CGRect frame = bg1.frame;
        frame.origin.x = 0;
        bg1.frame = frame;
    }];
}

-(void)backScreenFuncs:(id)sender{
    [UIView animateWithDuration:0.3 animations:^{
        CGRect frame = bg1.frame;
        frame.origin.x = 330;
        bg1.frame = frame;
    }];
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

-(NSString *)proceedString:(NSString*)string {
    const char *cstring = [string UTF8String];
    char *newcstr = (char*)malloc(2*sizeof(char)*(string.length));
    
    int j = 0;
    for (int i = 0; i < string.length; ++i) {
        if (islower(cstring[i]) && !islower(cstring[i+1]) && (cstring[i+1] != '(')){
            newcstr[j++] = '$';
        }
        newcstr[j++] = cstring[i];
    }
    newcstr[j++] = '\0';
    
    NSString * res = [NSString stringWithUTF8String:newcstr];
    free(newcstr);
    return  res;
}


-(void)proceedScheme:(id)sender{
    
    HUD = [[MBProgressHUD alloc] initWithView:bg2];
    [bg2 addSubview:HUD];
    HUD.removeFromSuperViewOnHide = YES;
    
    [DDParser setDefaultPowerAssociativity:DDOperatorAssociativityRight];
    evaluator = [[DDMathEvaluator alloc] init];
    NSError *error = nil;
    
    NSString* realFunc = [self proceedString:UField.text];
    NSString* phi0 = [self proceedString:Phi0Field.text];
    NSString* phil = [self proceedString:PhilField.text];
    NSString* psi1 = [self proceedString:Psi1Field.text];
    NSString* psi2 = [self proceedString:Psi2Field.text];
    NSString* psi1d = [self proceedString:Psi1dField.text];
    NSString* psi1dd = [self proceedString:Psi1ddField.text];
    NSString* f = [self proceedString:FField.text];
    
    expressionRealFunc = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:realFunc error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPhi_0 = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:phi0 error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPhi_l = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:phil error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPsi_1 = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:psi1 error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPsi_2 = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:psi2 error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPsi_1d = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:psi1d error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPsi_1dd = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:psi1dd error:&error] error:&error] parsedExpressionWithError:&error];
    expressionF = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:f error:&error] error:&error] parsedExpressionWithError:&error];
    
    NSArray *values = @[[NSNumber numberWithDouble:1], [NSNumber numberWithDouble:1], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"t", @"a", @"b", @"c"];
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    
    if ([expressionF evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция f(x,t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPhi_0 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция φ0(t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPhi_l evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция φl(t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPsi_1 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция ψ1(t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPsi_2 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция ψ2(t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPsi_1d evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция ψ1'(t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPsi_1dd evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция ψ1''(t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionRealFunc evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция U(x,t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    
    int K = [KField.text intValue];
    int Nx = [NField.text intValue];
#warning Ny!
    int Ny = [NField.text intValue];
    
    lx = [[lField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    double T = [[TField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    
    double tau = T / K;
    double hx = lx / Nx;
    double hy = ly / Ny;
    
    a = [[aField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    b = [[bField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    c = [[cField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    e = [[eField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    alpha = [[alphaField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    betta = [[bettaField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    gama = [[gammaField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    delta = [[deltaField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    long scheme = [schemePicker selectedRowInComponent:0]; // 0 -- явная, 1 -- неявная, 2 -- нестабильная
    long order = [orderPicker selectedRowInComponent:0];   // 0 -- 2-1, 1 -- 3-2, 2 -- 2-2
    startApproximation = [approxPicker selectedRowInComponent:0]; // 0 -- первого порядка, 2 -- второго порядка
    
//    //Для гиперболических
//    float sigma = a*tau/hx;
//    if(sigma > 1.0 && scheme == 0){
//        HUD.mode = MBProgressHUDModeText;
//        HUD.labelText = [NSString stringWithFormat:@"Не устойчива! σ = %.1f (σ > 1.0)", sigma];
//        HUD.detailsLabelText =  @"Измените параметры сетки или схему";
//        [HUD show:YES];
//        [HUD hide:YES afterDelay:3];
//        return;
//    }
    
    
    HUD.mode = MBProgressHUDAnimationFade;
    HUD.labelText = @"Идет расчет";
    [HUD show:YES];
    dispatch_async(dispatch_get_global_queue( DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^(void){
        double ***U;
        
        switch (scheme) {
            case 0:{
                switch (order) {
                    case 0:{
                     //   U = [self hyperbolic_explicitScheme_TwoPoint_FirstOrder:K :N :a :b :c :e :tau :h :alpha :betta :gama :delta];
                        break;
                    }
                    case 1:{
                    //    U = [ self hyperbolic_explicitScheme_ThreePoint_SecondOrder:K :N :a :b :c :e :tau :h :alpha :betta :gama :delta];
                        break;
                    }
                    case 2:{
                   //     U = [self hyperbolic_explicitScheme_TwoPoint_SecondOrder:K :N :a :b :c :e :tau :h :alpha :betta :gama :delta];
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
                  //      U = [self hyperbolic_implicitScheme_TwoPoint_FirstOrder:K :N :a :b :c :e :tau :h :alpha :betta :gama :delta];
                        break;
                    }
                    case 1:{
                   //     U = [self hyperbolic_implicitScheme_ThreePoint_SecondOrder:K :N :a :b :c :e :tau :h :alpha :betta :gama :delta];
                        break;
                    }
                    case 2:{
                 //       U = [self hyperbolic_implicitScheme_TwoPoint_SecondOrder:K :N :a :b :c :e :tau :h :alpha :betta :gama :delta];
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
                 //       U = [self hyperbolic_unstableScheme_TwoPoint_FirstOrder:K :N :a :b :c :e :tau :h :alpha :betta :gama :delta];
                        break;
                    }
                    case 1:{
                 //       U = [self hyperbolic_unstableScheme_ThreePoint_SecondOrder:K :N :a :b :c :e :tau :h :alpha :betta :gama :delta];
                        break;
                    }
                    case 2:{
                  //      U = [self hyperbolic_unstableScheme_TwoPoint_SecondOrder:K :N :a :b :c :e :tau :h :alpha :betta :gama :delta];
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
        
        
#warning сейчас строит только для x ( в момент y = 0)
        for (int i = 0; i < K; ++i) {
            NSMutableArray *contentArray = [[NSMutableArray alloc] init];
            NSMutableArray *contentArrayAnalytic = [[NSMutableArray alloc] init];
            
            for (int j = 0; j <= Nx; ++j) {
                id xAnalytic = [NSNumber numberWithDouble:(j)*hx];
                id yAnalytic = [NSNumber numberWithDouble:[self functionH:((j)*hx) withTime:i*tau]];
                [contentArrayAnalytic addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:xAnalytic, @"x", yAnalytic, @"y", nil]];
            }
            //epsiHlon[k][n] = max_i(u[k][i] - reshenie[k][i])
            id err = [NSNumber numberWithDouble:0.0];
            for (int j = 0; j <= Nx; ++j) {
                id x = [NSNumber numberWithDouble:j*hx];
                id y = [NSNumber numberWithDouble:U[i][j][0]];
                
                err = [NSNumber numberWithDouble:([err doubleValue] > fabs([self functionH:((j)*hx) withTime:i*tau] - U[i][j][0])) ? [err doubleValue] : fabs([self functionH:((j)*hx) withTime:i*tau] - U[i][j][0])];
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
    if (bg1.frame.origin.x > 0){
        [UIView animateWithDuration:0.3 animations:^{
            if(self.view.frame.size.height < 500)
                scrollBg.contentOffset = CGPointMake(0, 130);
            else
                scrollBg.contentOffset = CGPointMake(0, 30);
        }];
    }
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
    
    CALayer * imgLayer = bg.layer;
    [imgLayer setBorderColor: [[UIColor blackColor] CGColor]];
    [imgLayer setBorderWidth:0.5f];
    [imgLayer setShadowColor: [[UIColor blackColor] CGColor]];
    [imgLayer setShadowOpacity:0.9f];
    [imgLayer setShadowOffset: CGSizeMake(0, 1)];
    [imgLayer setShadowRadius:3.0];
    imgLayer.shouldRasterize = NO;
    
    scrollBg = [[UIScrollView alloc] initWithFrame:CGRectMake(0, 0, self.view.frame.size.width, self.view.frame.size.height)];
    
    [[NSNotificationCenter defaultCenter] addObserver:self selector:@selector(keyboardShown:) name:UIKeyboardDidShowNotification object:nil];
    [[NSNotificationCenter defaultCenter] addObserver:self selector:@selector(keyboardHidden:) name:UIKeyboardWillHideNotification object:nil];
    
    if(self.view.frame.size.height < 500)
        scrollBg.contentSize = CGSizeMake(320, 600);
    else
        scrollBg.contentSize = CGSizeMake(320, 500);
    
    [self.view addSubview:scrollBg];
    
    [scrollBg addSubview:bg];
    
    UILabel *aLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 10, 30, 20)];
    aLabel.text = @"a =";
    [bg addSubview:aLabel];
    aField = [[NMFTextField alloc] initWithFrame:CGRectMake(50, 12, 50, 20)];
    aField.text = @"1";
    [bg addSubview:aField];
    
    UILabel *bLabel = [[UILabel alloc] initWithFrame:CGRectMake(130, 10, 30, 20)];
    bLabel.text = @"b =";
    [bg addSubview:bLabel];
    bField = [[NMFTextField alloc] initWithFrame:CGRectMake(160, 12, 50, 20)];
    bField.text = @"1";
    [bg addSubview:bField];
    
    UILabel *cLabel = [[UILabel alloc] initWithFrame:CGRectMake(220, 10, 30, 20)];
    cLabel.text = @"c =";
    [bg addSubview:cLabel];
    cField = [[NMFTextField alloc] initWithFrame:CGRectMake(250, 12, 50, 20)];
    cField.text = @"-1";
    [bg addSubview:cField];
    
    UILabel *eLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 30, 30, 20)];
    eLabel.text = @"e =";
    [bg addSubview:eLabel];
    eField = [[NMFTextField alloc] initWithFrame:CGRectMake(50, 32, 50, 20)];
    eField.text = @"3";
    [bg addSubview:eField];
    
    UILabel *lLabel = [[UILabel alloc] initWithFrame:CGRectMake(130, 30, 30, 20)];
    lLabel.text = @"l =";
    [bg addSubview:lLabel];
    lField = [[NMFTextField alloc] initWithFrame:CGRectMake(160, 32, 50, 20)];
    lField.text = @"3.14";
    [bg addSubview:lField];
    
    UILabel *alphaLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 60, 60, 20)];
    alphaLabel.text = @"α =";
    [bg addSubview:alphaLabel];
    alphaField = [[NMFTextField alloc] initWithFrame:CGRectMake(80, 62, 50, 20)];
    alphaField.text = @"0";
    [bg addSubview:alphaField];
    
    UILabel *bettaLabel = [[UILabel alloc] initWithFrame:CGRectMake(150, 60, 60, 20)];
    bettaLabel.text = @"β =";
    [bg addSubview:bettaLabel];
    bettaField = [[NMFTextField alloc] initWithFrame:CGRectMake(210, 62, 50, 20)];
    bettaField.text = @"1";
    [bg addSubview:bettaField];
    
    UILabel *gamaLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 80, 60, 20)];
    gamaLabel.text = @"γ =";
    [bg addSubview:gamaLabel];
    gammaField = [[NMFTextField alloc] initWithFrame:CGRectMake(80, 82, 50, 20)];
    gammaField.text = @"0";
    [bg addSubview:gammaField];
    
    UILabel *deltaLabel = [[UILabel alloc] initWithFrame:CGRectMake(150, 80, 60, 20)];
    deltaLabel.text = @"δ =";
    [bg addSubview:deltaLabel];
    deltaField = [[NMFTextField alloc] initWithFrame:CGRectMake(210, 82, 50, 20)];
    deltaField.text = @"1";
    [bg addSubview:deltaField];
    
    UILabel *TLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 130, 200, 20)];
    TLabel.text = @"Финальное время (T):";
    [bg addSubview:TLabel];
    TField = [[NMFTextField alloc] initWithFrame:CGRectMake(250, 132, 50, 20)];
    TField.text = @"5";
    [bg addSubview:TField];
    
    UILabel *KLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 150, 230, 20)];
    KLabel.text = @"Разбиений по времени (K):";
    [bg addSubview:KLabel];
    KField = [[NMFTextField alloc] initWithFrame:CGRectMake(250, 150, 50, 20)];
    KField.text = @"400";
    [bg addSubview:KField];
    
    UILabel *NLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 170, 230, 20)];
    NLabel.text = @"Разбиений по OX (N):";
    [bg addSubview:NLabel];
    NField = [[NMFTextField alloc] initWithFrame:CGRectMake(250, 170, 50, 20)];
    NField.text = @"20";
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
    [next addTarget:self action:@selector(nextScreenFuncs:) forControlEvents:UIControlEventTouchUpInside];
    
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
    //bg1
    bg1 = [[UIView alloc] initWithFrame:CGRectMake(330, 160, 320, self.view.frame.size.height - 180)];
    bg1.backgroundColor = [UIColor colorWithRed:136./255. green:135./255. blue:148./255. alpha:1.0];
    CALayer * imgLayer1 = bg1.layer;
    [imgLayer1 setBorderColor: [[UIColor blackColor] CGColor]];
    [imgLayer1 setBorderWidth:0.5f];
    [imgLayer1 setShadowColor: [[UIColor blackColor] CGColor]];
    [imgLayer1 setShadowOpacity:0.9f];
    [imgLayer1 setShadowOffset: CGSizeMake(0, 1)];
    [imgLayer1 setShadowRadius:3.0];
    imgLayer1.shouldRasterize = NO;
    [self.view addSubview:bg1];
    
    UIButton *nextFuncs;
    if (self.view.bounds.size.height < 500)
        nextFuncs = [[UIButton alloc] initWithFrame:CGRectMake(230, 260, 60, 40)];
    else
        nextFuncs = [[UIButton alloc] initWithFrame:CGRectMake(230, 320, 60, 40)];
    UILabel *blabFuncs = [[UILabel alloc] initWithFrame:CGRectMake(0, 0, 60, 20)];
    blabFuncs.text = @"Далее";
    [nextFuncs addSubview:blabFuncs];
    [nextFuncs addTarget:self action:@selector(nextScreen:) forControlEvents:UIControlEventTouchUpInside];
    [bg1 addSubview:nextFuncs];
    
    UIButton *backFuncs;
    if ([[UIScreen mainScreen] bounds].size.height < 500)
        backFuncs = [[UIButton alloc] initWithFrame:CGRectMake(40, 260, 60, 40)];
    else
        backFuncs = [[UIButton alloc] initWithFrame:CGRectMake(40, 320, 60, 40)];
    
    UILabel *blabFuncsBack = [[UILabel alloc] initWithFrame:CGRectMake(0, 0, 60, 20)];
    blabFuncsBack.text = @"Назад";
    [backFuncs addSubview:blabFuncsBack];
    [backFuncs addTarget:self action:@selector(backScreenFuncs:) forControlEvents:UIControlEventTouchUpInside];
    [bg1 addSubview:backFuncs];
    
    UILabel *FLabel = [[UILabel alloc] initWithFrame:CGRectMake(10, 10, 80, 20)];
    FLabel.text = @"f(x,t) =";
    [bg1 addSubview:FLabel];
    FField = [[UITextField alloc] initWithFrame:CGRectMake(84, 12, 230, 20)];
    FField.text = @"sin(x)exp(-t)";
    FField.delegate = self;
    [bg1 addSubview:FField];
    
    UILabel *Phi0Label = [[UILabel alloc] initWithFrame:CGRectMake(10, 30, 80, 20)];
    Phi0Label.text = @"φ0(t) =";
    [bg1 addSubview:Phi0Label];
    Phi0Field = [[UITextField alloc] initWithFrame:CGRectMake(84, 32, 230, 20)];
    Phi0Field.text = @"exp(-t)";
    Phi0Field.delegate = self;
    [bg1 addSubview:Phi0Field];
    
    UILabel *PhilLabel = [[UILabel alloc] initWithFrame:CGRectMake(10, 50, 80, 20)];
    PhilLabel.text = @"φl(t) =";
    [bg1 addSubview:PhilLabel];
    PhilField = [[UITextField alloc] initWithFrame:CGRectMake(84, 52, 230, 20)];
    PhilField.text = @"-exp(-t)";
    PhilField.delegate = self;
    [bg1 addSubview:PhilField];
    
    UILabel *Psi1Label = [[UILabel alloc] initWithFrame:CGRectMake(10, 70, 80, 20)];
    Psi1Label.text = @"ψ1(x) =";
    [bg1 addSubview:Psi1Label];
    Psi1Field = [[UITextField alloc] initWithFrame:CGRectMake(84, 72, 230, 20)];
    Psi1Field.text = @"cos(x)";
    Psi1Field.delegate = self;
    [bg1 addSubview:Psi1Field];
    
    UILabel *Psi2Label = [[UILabel alloc] initWithFrame:CGRectMake(10, 90, 80, 20)];
    Psi2Label.text = @"ψ2(x) =";
    [bg1 addSubview:Psi2Label];
    Psi2Field = [[UITextField alloc] initWithFrame:CGRectMake(84, 92, 230, 20)];
    Psi2Field.text = @"-cos(x)";
    Psi2Field.delegate = self;
    [bg1 addSubview:Psi2Field];
    
    UILabel *Psi1dLabel = [[UILabel alloc] initWithFrame:CGRectMake(10, 120, 80, 20)];
    Psi1dLabel.text = @"ψ1'(x) =";
    [bg1 addSubview:Psi1dLabel];
    Psi1dField = [[UITextField alloc] initWithFrame:CGRectMake(84, 122, 230, 20)];
    Psi1dField.text = @"-sin(x)";
    Psi1dField.delegate = self;
    [bg1 addSubview:Psi1dField];
    
    UILabel *Psi1ddLabel = [[UILabel alloc] initWithFrame:CGRectMake(10, 140, 80, 20)];
    Psi1ddLabel.text = @"ψ1''(x) =";
    [bg1 addSubview:Psi1ddLabel];
    Psi1ddField = [[UITextField alloc] initWithFrame:CGRectMake(84, 142, 230, 20)];
    Psi1ddField.text = @"-cos(x)";
    Psi1ddField.delegate = self;
    [bg1 addSubview:Psi1ddField];
    
    UILabel *ULabel = [[UILabel alloc] initWithFrame:CGRectMake(10, 160, 80, 20)];
    ULabel.text = @"U(x,t) =";
    [bg1 addSubview:ULabel];
    UField = [[UITextField alloc] initWithFrame:CGRectMake(84, 162, 230, 20)];
    UField.text = @"exp(-t)cos(x)";
    UField.delegate = self;
    [bg1 addSubview:UField];
    
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
    NSCharacterSet *cs = [[NSCharacterSet characterSetWithCharactersInString:ALPHABET2] invertedSet];
    NSString *filtered = [[string componentsSeparatedByCharactersInSet:cs] componentsJoinedByString:@""];
    return (([string isEqualToString:filtered]));
}


@end
