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

- (void)viewDidLoad
{
    [super viewDidLoad];
    double x[5] = {-14, -55, 49, 86,8};
    double a[5] = {0, 7, -4, 7, 4};
    double b[5] = {8, -19, 21, -23, -7};
    double c[5] = {-2, 9, -8, 9, 0};
    size_t N = 5;
    solve_tridiagonal_in_place_destructive(x, N, a, b, c);
    //should be -1.0000	3.0000	1.0000	-5.0000	-4.0000
    NSLog(@"answer: %f %f %f %f %f", x[0], x[1], x[2], x[3], x[4]);

}

- (void)didReceiveMemoryWarning
{
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

@end
