//
//  ViewController.m
//  numericalMethods
//
//  Created by Ilya on 02.10.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import "ViewController.h"


#define NUMBERS_ONLY @"1234567890."
#define CHARACTER_LIMIT 5

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


double * solve_tridiagonal_in_place_destructive(double *x, const size_t N, const double *a, const double *b, double *c) {

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
    
    double *answer = (double *)malloc((N)*sizeof(double));
    for (int j = 0; j < N; ++j)
        answer[j] = x[j];
    return answer;
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

double **implicitScheme_TwoPoint_FirstOrder(int K, int N, double a, double b, double c, double tau, double h, double alpha, double betta, double gamma, double delta, double tetta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    double *lower = (double *)malloc((N+1) * sizeof(double));
    double *mid = (double *)malloc((N+1) * sizeof(double));
    double *upper = (double *)malloc((N+1) * sizeof(double));
    double *answer = (double *)malloc((N+1) * sizeof(double));
    
    memset(lower, 0, (N+1) * sizeof(double));
    memset(upper, 0, (N+1) * sizeof(double));
    memset(mid, 0, (N+1) * sizeof(double));
    memset(answer, 0, (N+1) * sizeof(double));
    
    //double tetta = 1.0;
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    for (int i = 0; i <= N; ++i)
        U[0][i] = phi(i*h);
    
    for (int k = 0; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            lower[i] = tetta*tau*(a/h/h - b/2/h);
            mid[i] = -1 + tetta*tau*(-2*a/h/h + c);
            upper[i] = tetta*tau*(a/h/h + b/2/h);
            answer[i] = -U[k][i] + (tetta - 1)*tau*(a*(U[k][i+1] - 2*U[k][i] + U[k][i-1])/h/h + b*(U[k][i+1] - U[k][i-1])/2/h + c*U[k][i]) - tau*f(i * h, k * tau);
        }
        
        //Двухточечная первого порядка
        lower[0] = 0.0;
        mid[0] = (betta - alpha/h);
        upper[0] = alpha/h;
        answer[0] = phi_0((k+1)*tau);
        
        lower[N] = -gamma/h;
        mid[N] = delta + gamma/h;
        upper[N] = 0.0;
        answer[N] = phi_l((k+1)*tau);
        
        U[k+1] = solve_tridiagonal_in_place_destructive(answer, N+1, lower, mid, upper);
    }
    return U;
}

double **implicitScheme_ThreePoint_SecondOrder(int K, int N, double a, double b, double c, double tau, double h, double alpha, double betta, double gamma, double delta, double tetta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    double *lower = (double *)malloc((N+1) * sizeof(double));
    double *mid = (double *)malloc((N+1) * sizeof(double));
    double *upper = (double *)malloc((N+1) * sizeof(double));
    double *answer = (double *)malloc((N+1) * sizeof(double));
    
    memset(lower, 0, (N+1) * sizeof(double));
    memset(upper, 0, (N+1) * sizeof(double));
    memset(mid, 0, (N+1) * sizeof(double));
    memset(answer, 0, (N+1) * sizeof(double));
    
    //double tetta = 1.0;
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    for (int i = 0; i <= N; ++i)
        U[0][i] = phi(i*h);
    
    for (int k = 0; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            lower[i] = tetta*tau*(a/h/h - b/2/h);
            mid[i] = -1 + tetta*tau*(-2*a/h/h + c);
            upper[i] = tetta*tau*(a/h/h + b/2/h);
            answer[i] = -U[k][i] + (tetta - 1)*tau*(a*(U[k][i+1] - 2*U[k][i] + U[k][i-1])/h/h + b*(U[k][i+1] - U[k][i-1])/2/h + c*U[k][i]) - tau*f(i * h, k * tau);
        }
        
        
        //Трехточечная второго порядка
        lower[0] = 0.0;
        mid[0] = (betta - 3*alpha/2/h +lower[1]/upper[1]*alpha/2/h);
        upper[0] = 2*alpha/h + mid[1]/upper[1]*alpha/2/h;
        answer[0] = phi_0((k+1)*tau) + answer[1]/upper[1]*alpha/2/h;
        
        lower[N] = -2*gamma/h - mid[N-1]/lower[N-1]*gamma/2/h;
        mid[N] = delta + 3*gamma/2/h - upper[N-1]/lower[N-1]*gamma/2/h;
        upper[N] = 0.0;
        answer[N] = phi_l((k+1)*tau) - answer[N-1]/lower[N-1]*gamma/2/h;
        
        U[k+1] = solve_tridiagonal_in_place_destructive(answer, N+1, lower, mid, upper);
    }
    return U;
}

double **implicitScheme_TwoPoint_SecondOrder(int K, int N, double a, double b, double c, double tau, double h, double alpha, double betta, double gamma, double delta, double tetta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    double *lower = (double *)malloc((N+1) * sizeof(double));
    double *mid = (double *)malloc((N+1) * sizeof(double));
    double *upper = (double *)malloc((N+1) * sizeof(double));
    double *answer = (double *)malloc((N+1) * sizeof(double));
    
    memset(lower, 0, (N+1) * sizeof(double));
    memset(upper, 0, (N+1) * sizeof(double));
    memset(mid, 0, (N+1) * sizeof(double));
    memset(answer, 0, (N+1) * sizeof(double));
    
    //double tetta = 1.0;
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    for (int i = 0; i <= N; ++i)
        U[0][i] = phi(i*h);
    
    for (int k = 0; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            lower[i] = tetta*tau*(a/h/h - b/2/h);
            mid[i] = -1 + tetta*tau*(-2*a/h/h + c);
            upper[i] = tetta*tau*(a/h/h + b/2/h);
            answer[i] = -U[k][i] + (tetta - 1)*tau*(a*(U[k][i+1] - 2*U[k][i] + U[k][i-1])/h/h + b*(U[k][i+1] - U[k][i-1])/2/h + c*U[k][i]) - tau*f(i * h, k * tau);
        }
        
        //Двухточечная второго порядка
        lower[0] = 0.0;
        mid[0] = 2*a*alpha/h + h*alpha/tau - c*h*alpha - betta*(2*a - b*h);
        upper[0] = -2*a*alpha/h;
        answer[0] = U[k][0]*h*alpha/tau - phi_0((k+1) * tau)*(2*a -b*h) + f(0, (k+1)*tau)*h*alpha;
        
        lower[N] = -2*a*gamma/h;
        mid[N] = 2*a*gamma/h + h*gamma/tau - c*h*gamma + delta*(2*a + b*h);
        upper[N] = 0.0;
        answer[N] = U[k][N]*h*gamma/tau + phi_l((k + 1) * tau) * (2*a + b*h) + f(N, (k+1)*tau)*h*alpha;
        
        U[k+1] = solve_tridiagonal_in_place_destructive(answer, N+1, lower, mid, upper);
    }
    return U;
}
double **explicitScheme_TwoPoint_FirstOrder(int K, int N, double a, double b, double c, double tau, double h, double alpha, double betta, double gamma, double delta){
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
        U[k+1][0] = -alpha/h/(betta - alpha/h)*U[k+1][1]  + phi_0((k + 1) * tau)/(betta - alpha/h);
        U[k+1][N] = gamma/h/(delta + gamma/h)*U[k+1][N-1] + phi_l((k + 1) * tau)/(delta + gamma/h);
    }
    return U;
}

double **explicitScheme_ThreePoint_SecondOrder(int K, int N, double a, double b, double c, double tau, double h, double alpha, double betta, double gamma, double delta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    for (int i = 0; i <= N; ++i)
        U[0][i] = phi(i*h);
    
    for (int k = 0; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            U[k+1][i] = U[k][i] + a*tau/h/h*(U[k][i+1] - 2*U[k][i] + U[k][i-1]) + b*tau/2/h*(U[k][i+1] - U[k][i-1]) + c*tau*U[k][i] + tau*f(i * h, k * tau);
        }
        
        //  трехточечная второго
          U[k+1][0] = alpha/2/h/(betta - 3*alpha/2/h)*(U[k+1][2] - 4*U[k+1][1]) + phi_0((k+1) * tau)/(betta - 3*alpha/2/h);
          U[k+1][N] = gamma/2/h/(delta - 3*gamma/2/h)*(4*U[k+1][N-1] - U[k+1][N-2]) + phi_l((k+1) * tau)/(delta - 3*gamma/2/h);
    }
    return U;
}

double **explicitScheme_TwoPoint_SecondOrder(int K, int N, double a, double b, double c, double tau, double h, double alpha, double betta, double gamma, double delta){
    double **U = (double **)malloc(K * sizeof(double *));
    
    for (int i = 0; i < K; i++)
        U[i] = (double *)malloc((N+1) * sizeof(double));
    
    for (int i = 0; i <= N; ++i)
        U[0][i] = phi(i*h);
    
    for (int k = 0; k < K - 1; ++k) {
        for (int i = 1; i < N; ++i) {
            U[k+1][i] = U[k][i] + a*tau/h/h*(U[k][i+1] - 2*U[k][i] + U[k][i-1]) + b*tau/2/h*(U[k][i+1] - U[k][i-1]) + c*tau*U[k][i] + tau*f(i * h, k * tau);
        }
        //  Двухточечная второго
            U[k+1][0] = 1/(2*a*alpha/h + h*alpha/tau - c*h*alpha - betta*(2*a - b*h))*(U[k+1][1]*2*a*alpha/h + U[k][0]*h*alpha/tau - phi_0((k+1) * tau)*(2*a-b*h)+f(0, (k+1)*tau)*h*alpha);
           U[k+1][N] = 1/(2*a*gamma/h + h*gamma/tau - c*h*gamma + delta*(2*a + b*h))*(U[k+1][N-1]*2*a*gamma/h + U[k][N]*h*gamma/tau + phi_l((k+1) * tau)*(2*a+b*h)+f(N, (k+1)*tau)*h*gamma);
    }
    return U;
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
        frame.origin.x = 320;
        bg2.frame = frame;
    }];
}

-(void)backFromGraph:(id)sender{
    [UIView animateWithDuration:0.3 animations:^{
        CGRect frame = bg3.frame;
        frame.origin.y = -1000;
        bg3.frame = frame;
    }];
}

-(void)proceedScheme:(id)sender{
    int K = [KField.text intValue];//1000;
    int N = [NField.text intValue];//10;
    
    l = [lField.text doubleValue];//3.14; //1
    double T = [TField.text doubleValue];//10;
    
    double tau = T / K, h = l / N;
    
    a = [aField.text doubleValue], b = [bField.text doubleValue], c = [cField.text doubleValue];
    alpha = [alphaField.text doubleValue], betta = [bettaField.text doubleValue], gamma = [gammaField.text doubleValue], delta = [deltaField.text doubleValue];//double alpha = 0.0, betta = 1.0, gamma = 0.0, delta = 1.0;
    
    int scheme = [schemePicker selectedRowInComponent:0]; // 0 -- явная, 1 -- неявная, 2 -- Кранка
    int order = [orderPicker selectedRowInComponent:0];   // 0 -- 2-1, 1 -- 3-2, 2 -- 2-2
    double tetta = 0;
    
    if(a*tau/h/h > 0.5 && scheme == 0){
        printf("bad sigma");
        exit(0);
    }
    double **U;
    switch (scheme) {
        case 0:{
            switch (order) {
                case 0:{
                    U = explicitScheme_TwoPoint_FirstOrder(K, N, a, b, c, tau, h, alpha, betta, gamma, delta);
                    break;
                }
                case 1:{
                    U = explicitScheme_ThreePoint_SecondOrder(K, N, a, b, c, tau, h, alpha, betta, gamma, delta);
                    break;
                }
                case 2:{
                    U = explicitScheme_TwoPoint_SecondOrder(K, N, a, b, c, tau, h, alpha, betta, gamma, delta);
                    break;
                }
                    
                default:
                    break;
            }
            break;
        }
        case 1:{
            tetta = 1.0;
            switch (order) {
                case 0:{
                    U = implicitScheme_TwoPoint_FirstOrder(K, N, a, b, c, tau, h, alpha, betta, gamma, delta, tetta);
                    break;
                }
                case 1:{
                    U = implicitScheme_ThreePoint_SecondOrder(K, N, a, b, c, tau, h, alpha, betta, gamma, delta, tetta);
                    break;
                }
                case 2:{
                    U = implicitScheme_TwoPoint_SecondOrder(K, N, a, b, c, tau, h, alpha, betta, gamma, delta, tetta);
                    break;
                }
                    
                default:
                    break;
            }
            break;
        }
        case 2:{
            tetta = 0.5;
            switch (order) {
                case 0:{
                    U = implicitScheme_TwoPoint_FirstOrder(K, N, a, b, c, tau, h, alpha, betta, gamma, delta, tetta);
                    break;
                }
                case 1:{
                    U = implicitScheme_ThreePoint_SecondOrder(K, N, a, b, c, tau, h, alpha, betta, gamma, delta, tetta);
                    break;
                }
                case 2:{
                    U = implicitScheme_TwoPoint_SecondOrder(K, N, a, b, c, tau, h, alpha, betta, gamma, delta, tetta);
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
    
    [UIView animateWithDuration:0.3 animations:^{
        CGRect frame = bg3.frame;
        frame.origin.y = 0;
        bg3.frame = frame;
    }];
    //printf for grapher
    for (int i = 0; i < K; ++i) {
        printf("K = %f\n", i*tau);
        for (int j = 0; j <= N; ++j) {
            printf("%f	%f\n", j*h, U[i][j]);
        }
        printf("\n\n\n");
    }
    
}

- (void)viewDidLoad
{
    [super viewDidLoad];
    [self.view setBackgroundColor:[UIColor whiteColor]];
    
    bg = [[UIView alloc] initWithFrame:CGRectMake(0, 160, 320, self.view.frame.size.height - 160)];
    bg.backgroundColor = [UIColor orangeColor];
    [self.view addSubview:bg];
    
    UILabel *aLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 10, 30, 20)];
    aLabel.text = @"a =";
    [bg addSubview:aLabel];
    aField = [[UITextField alloc] initWithFrame:CGRectMake(50, 12, 50, 20)];
    aField.delegate = self;
    aField.text = @"1";
    [bg addSubview:aField];
    
    UILabel *bLabel = [[UILabel alloc] initWithFrame:CGRectMake(100, 10, 30, 20)];
    bLabel.text = @"b =";
    [bg addSubview:bLabel];
    bField = [[UITextField alloc] initWithFrame:CGRectMake(130, 12, 50, 20)];
    bField.delegate = self;
    bField.text = @"0";
    [bg addSubview:bField];
    
    UILabel *cLabel = [[UILabel alloc] initWithFrame:CGRectMake(180, 10, 30, 20)];
    cLabel.text = @"c =";
    [bg addSubview:cLabel];
    cField = [[UITextField alloc] initWithFrame:CGRectMake(210, 12, 50, 20)];
    cField.delegate = self;
    cField.text = @"0";
    [bg addSubview:cField];
    
    UILabel *alphaLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 40, 60, 20)];
    alphaLabel.text = @"alpha =";
    [bg addSubview:alphaLabel];
    alphaField = [[UITextField alloc] initWithFrame:CGRectMake(80, 42, 50, 20)];
    alphaField.delegate = self;
    alphaField.text = @"1";
    [bg addSubview:alphaField];
    
    UILabel *bettaLabel = [[UILabel alloc] initWithFrame:CGRectMake(150, 40, 60, 20)];
    bettaLabel.text = @"betta =";
    [bg addSubview:bettaLabel];
    bettaField = [[UITextField alloc] initWithFrame:CGRectMake(210, 42, 50, 20)];
    bettaField.delegate = self;
    bettaField.text = @"0";
    [bg addSubview:bettaField];
    
    UILabel *gammaLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 80, 60, 20)];
    gammaLabel.text = @"gamma =";
    [bg addSubview:gammaLabel];
    gammaField = [[UITextField alloc] initWithFrame:CGRectMake(80, 82, 50, 20)];
    gammaField.delegate = self;
    gammaField.text = @"1";
    [bg addSubview:gammaField];
    
    UILabel *deltaLabel = [[UILabel alloc] initWithFrame:CGRectMake(150, 80, 60, 20)];
    deltaLabel.text = @"delta =";
    [bg addSubview:deltaLabel];
    deltaField = [[UITextField alloc] initWithFrame:CGRectMake(210, 82, 50, 20)];
    deltaField.delegate = self;
    deltaField.text = @"0";
    [bg addSubview:deltaField];
    
    UILabel *lLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 100, 30, 20)];
    lLabel.text = @"l =";
    [bg addSubview:lLabel];
    lField = [[UITextField alloc] initWithFrame:CGRectMake(50, 102, 50, 20)];
    lField.delegate = self;
    lField.text = @"3.14";
    [bg addSubview:lField];
    
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
    KField.text = @"1000";
    [bg addSubview:KField];
    
    UILabel *NLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 170, 230, 20)];
    NLabel.text = @"Разбиенией по OX (N):";
    [bg addSubview:NLabel];
    NField = [[UITextField alloc] initWithFrame:CGRectMake(250, 170, 50, 20)];
    NField.delegate = self;
    NField.text = @"10";
    [bg addSubview:NField];

    UIButton *next = [[UIButton alloc] initWithFrame:CGRectMake(230, 320, 60, 40)];
    UILabel *blab = [[UILabel alloc] initWithFrame:CGRectMake(0, 0, 60, 20)];
    blab.text = @"Далее";
    [next addSubview:blab];
    [next addTarget:self action:@selector(nextScreen:) forControlEvents:UIControlEventTouchUpInside];
    [bg addSubview:next];
//bg 2
    bg2 = [[UIView alloc] initWithFrame:CGRectMake(320, 160, 320, self.view.frame.size.height - 160)];
    bg2.backgroundColor = [UIColor lightGrayColor];
    [self.view addSubview:bg2];
    
    orderPicker = [[UIPIckerViewOrder alloc] initWithFrame:CGRectMake(0, 0, 320, 100)];
    schemePicker = [[UIPickerViewScheme alloc] initWithFrame:CGRectMake(0, 140, 320, 100)];
   
    [bg2 addSubview:orderPicker];
    [bg2 addSubview:schemePicker];

    UIButton *nextButton = [[UIButton alloc] initWithFrame:CGRectMake(230, 320, 60, 40)];
    UILabel *blab2 = [[UILabel alloc] initWithFrame:CGRectMake(0, 0, 60, 20)];
    blab2.text = @"Далее";
    [nextButton addSubview:blab2];
    [nextButton addTarget:self action:@selector(proceedScheme:) forControlEvents:UIControlEventTouchUpInside];
    [bg2 addSubview:nextButton];
    
    UIButton *backButton = [[UIButton alloc] initWithFrame:CGRectMake(40, 320, 60, 40)];
    UILabel *blab3 = [[UILabel alloc] initWithFrame:CGRectMake(0, 0, 60, 20)];
    blab3.text = @"Назад";
    [backButton addSubview:blab3];
    [backButton addTarget:self action:@selector(backScreen:) forControlEvents:UIControlEventTouchUpInside];
    [bg2 addSubview:backButton];
    
//bg3
    bg3 = [[UIView alloc] initWithFrame:CGRectMake(0, -1000, 320, self.view.frame.size.height)];
    bg3.backgroundColor = [UIColor grayColor];
    [self.view addSubview:bg3];
    
    UIButton *bb = [[UIButton alloc] initWithFrame:CGRectMake(40, 520, 60, 40)];
    UILabel *blab4 = [[UILabel alloc] initWithFrame:CGRectMake(0, 0, 60, 20)];
    blab4.text = @"Назад";
    [bb addSubview:blab4];
    [bb addTarget:self action:@selector(backFromGraph:) forControlEvents:UIControlEventTouchUpInside];
    [bg3 addSubview:bb];

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

#pragma mark - PickerView delegate



- (void)didReceiveMemoryWarning
{
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

@end
