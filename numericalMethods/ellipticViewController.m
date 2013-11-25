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


double bx,by,c,alpha1,alpha2,alpha3,alpha4,betta1,betta2,betta3,betta4,l,lx,ly, epsilon;
long startApproximation;

@interface ellipticViewController ()

@end

@implementation ellipticViewController
@synthesize HUD;
@synthesize system;

#pragma mark - init funcs
- (double)analyticFunction:(double)x y:(double)y{
    NSArray *values = @[[NSNumber numberWithDouble:x], [NSNumber numberWithDouble:y], [NSNumber numberWithDouble:bx], [NSNumber numberWithDouble:by], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"y", @"bx", @"by", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionRealFunc evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) f:(double) x y:(double) y{
    NSArray *values = @[[NSNumber numberWithDouble:x], [NSNumber numberWithDouble:y], [NSNumber numberWithDouble:bx], [NSNumber numberWithDouble:by], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"y", @"bx", @"by", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionF evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) phi1:(double) y{
    NSArray *values = @[[NSNumber numberWithDouble:y], [NSNumber numberWithDouble:bx], [NSNumber numberWithDouble:by], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"y", @"bx", @"by", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPhi1 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) phi2:(double) y{
    NSArray *values = @[[NSNumber numberWithDouble:y], [NSNumber numberWithDouble:bx], [NSNumber numberWithDouble:by], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"y", @"bx", @"by", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPhi2 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) phi3:(double) x{
    NSArray *values = @[[NSNumber numberWithDouble:x], [NSNumber numberWithDouble:bx], [NSNumber numberWithDouble:by], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"bx", @"by", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPhi3 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) phi4:(double) x{
    NSArray *values = @[[NSNumber numberWithDouble:x], [NSNumber numberWithDouble:bx], [NSNumber numberWithDouble:by], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"bx", @"by", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPhi4 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

#pragma mark - elliptic Liebmann
-(double **)elliptic_liebmann:(int) Nx Ny:(int) Ny bx:(double) bx by:(double) by c:(double) c hx:(double) hx hy:(double) hy alpha1:(double) alpha1 betta1:(double) betta1 alpha2:(double) alpha2 betta2:(double) betta2 alpha3:(double) alpha3 betta3:(double) betta3 alpha4:(double) alpha4 betta4:(double) betta4 epsilon:(double) epsilon{
    
    int iteration = 0;
    int MAX_ITERATION = 1000;
    
    double **U = (double **)malloc((Nx +1) * sizeof(double *));
    for (int i = 0; i < (Nx +1); i++){
        U[i] = (double *)malloc((Ny+1) * sizeof(double));
    }
    
    double **Uprev = (double **)malloc((Nx +1) * sizeof(double *));
    for (int i = 0; i < (Nx +1); i++){
        Uprev[i] = (double *)malloc((Ny+1) * sizeof(double));
    }
    
    for (int i = 0; i < Nx+1; i++){
        for (int j = 0; j < Ny+1; j++){
            Uprev[i][j] = 0;
        }
    }
    
    while (true){
        for (int j = 0; j <= Ny; ++j)
            U[0][j] = (alpha1/2/hx)/((betta1-3*alpha1)/2/hx) * (Uprev[2][j] - 4*Uprev[1][j]) + [self phi1:(hy * j)]/(betta1 - 3*alpha1/2/hx);
        for (int i = 0; i <= Nx; ++i)
            U[i][0] = (alpha3/2/hy)/(betta3 - 3*alpha3/2/hy) * (Uprev[i][2] - 4*Uprev[i][1]) + [self phi3:(hx * i)]/(betta3 - 3*alpha3/2/hy);
            
        for (int i = 1; i < Nx; ++i) {
            for (int j = 1; j < Ny; ++j) {
                U[i][j] = (1/(2*hx*hx - c*hx*hx*hy*hy +2*hy*hy)) * (Uprev[i+1][j]*hy*hy*(1 + bx*hx/2) + Uprev[i-1][j]*hy*hy*(1-bx*hx/2) + Uprev[i][j+1]*hx*hx*(1+by*hy/2) + Uprev[i][j-1]*hx*hx*(1 - by*hy/2) + hx*hx*hy*hy*[self f:i*hx y:j*hy]);
                
                double w = 1.0; // 0 < w <2
                U[i][j] = w*U[i][j] + (1-w)*Uprev[i][j];
            }
        }
        for (int j = 0; j <= Ny; ++j)
            U[Nx][j] = (alpha2/2/hx)/(betta2 + 3*alpha2/2/hx) * (4*Uprev[Nx-1][j] - Uprev[Nx-2][j]) + [self phi2:(hy * j)]/(betta2 + 3*alpha2/2/hx);
        for (int i = 0; i <= Nx; ++i)
            U[i][Ny] = (alpha4/2/hy)/(betta4 + 3*alpha4/2/hy) * (4*Uprev[i][Ny-1] - Uprev[i][Ny-2]) + [self phi4:(hx * i)]/(betta4 + 3*alpha4/2/hy);
        
        double cur_eps = 0.0;
        for (int i = 0; i <= Nx; ++i) {
            for (int j = 0; j <= Ny; ++j) {
                cur_eps = MAX(cur_eps, fabs(U[i][j] - Uprev[i][j]));
                Uprev[i][j] = U[i][j];
            }
        }
        
        if(cur_eps < epsilon || iteration >= MAX_ITERATION)
            break;
        iteration++;
    }
    
    return U;
}

#pragma mark - elliptic Seidel
-(double **)elliptic_seidel:(int) Nx Ny:(int) Ny bx:(double) bx by:(double) by c:(double) c hx:(double) hx hy:(double) hy alpha1:(double) alpha1 betta1:(double) betta1 alpha2:(double) alpha2 betta2:(double) betta2 alpha3:(double) alpha3 betta3:(double) betta3 alpha4:(double) alpha4 betta4:(double) betta4 epsilon:(double) epsilon{
    
    int iteration = 0;
    int MAX_ITERATION = 1000;
    
    double **U = (double **)malloc((Nx +1) * sizeof(double *));
    for (int i = 0; i < (Nx +1); i++){
        U[i] = (double *)malloc((Ny+1) * sizeof(double));
    }
    
    double **Uprev = (double **)malloc((Nx +1) * sizeof(double *));
    for (int i = 0; i < (Nx +1); i++){
        Uprev[i] = (double *)malloc((Ny+1) * sizeof(double));
    }
    
    for (int i = 0; i < Nx+1; i++){
        for (int j = 0; j < Ny+1; j++){
            Uprev[i][j] = 0;
        }
    }
    
    while (true){
        for (int j = 0; j <= Ny; ++j)
            U[0][j] = (alpha1/2/hx)/((betta1-3*alpha1)/2/hx) * (Uprev[2][j] - 4*Uprev[1][j]) + [self phi1:(hy * j)]/(betta1 - 3*alpha1/2/hx);
        for (int i = 0; i <= Nx; ++i)
            U[i][0] = (alpha3/2/hy)/(betta3 - 3*alpha3/2/hy) * (Uprev[i][2] - 4*Uprev[i][1]) + [self phi3:(hx * i)]/(betta3 - 3*alpha3/2/hy);
        
        for (int i = 1; i < Nx; ++i) {
            for (int j = 1; j < Ny; ++j) {
                U[i][j] = (1/(2*hx*hx - c*hx*hx*hy*hy +2*hy*hy)) * (Uprev[i+1][j]*hy*hy*(1 + bx*hx/2) + U[i-1][j]*hy*hy*(1-bx*hx/2) + Uprev[i][j+1]*hx*hx*(1+by*hy/2) + U[i][j-1]*hx*hx*(1 - by*hy/2) + hx*hx*hy*hy*[self f:i*hx y:j*hy]);
                
                double w = 1.0;
                U[i][j] = w*U[i][j] + (1-w)*Uprev[i][j];
            }
        }
        for (int j = 0; j <= Ny; ++j)
            U[Nx][j] = (alpha2/2/hx)/(betta2 + 3*alpha2/2/hx) * (4*U[Nx-1][j] - U[Nx-2][j]) + [self phi2:(hy * j)]/(betta2 + 3*alpha2/2/hx);
        for (int i = 0; i <= Nx; ++i)
            U[i][Ny] = (alpha4/2/hy)/(betta4 + 3*alpha4/2/hy) * (4*U[i][Ny-1] - U[i][Ny-2]) + [self phi4:(hx * i)]/(betta4 + 3*alpha4/2/hy);
        
        double cur_eps = 0.0;
        for (int i = 0; i <= Nx; ++i) {
            for (int j = 0; j <= Ny; ++j) {
                cur_eps = MAX(cur_eps, fabs(U[i][j] - Uprev[i][j]));
                Uprev[i][j] = U[i][j];
            }
        }
        if(cur_eps < epsilon || iteration >= MAX_ITERATION)
            break;
        iteration++;
    }
    return U;
}

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
    NSString* phi1 = [self proceedString:Phi1Field.text];
    NSString* phi2 = [self proceedString:Phi2Field.text];
    NSString* phi3 = [self proceedString:Phi3Field.text];
    NSString* phi4 = [self proceedString:Phi4Field.text];
    NSString* f = [self proceedString:FField.text];
    
    expressionRealFunc = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:realFunc error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPhi1 = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:phi1 error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPhi2 = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:phi2 error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPhi3 = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:phi3 error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPhi4 = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:phi4 error:&error] error:&error] parsedExpressionWithError:&error];
    expressionF = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:f error:&error] error:&error] parsedExpressionWithError:&error];
    
    NSArray *values = @[[NSNumber numberWithDouble:1], [NSNumber numberWithDouble:1], [NSNumber numberWithDouble:bx], [NSNumber numberWithDouble:by], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"y", @"bx", @"by", @"c"];
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    
    if ([expressionF evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция f(x,y).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPhi1 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция φ1(y).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPhi2 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция φ2(y).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPhi3 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция φ3(x).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }

    if ([expressionPhi4 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция φ4(x).";
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
    
    int Nx = [NxField.text intValue];
    int Ny = [NyField.text intValue];
    
    lx = [[lxField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    ly = [[lxField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];

    double hx = lx / Nx;
    double hy = ly / Ny;
    
    bx = [[bxField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    by = [[byField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    c = [[cField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    alpha1 = [[alpha1Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    betta1 = [[betta1Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    alpha2 = [[alpha2Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    betta2 = [[betta2Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    alpha3 = [[alpha3Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    betta3 = [[betta3Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    alpha4 = [[alpha4Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    betta4 = [[betta4Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    long scheme = [schemePicker selectedRowInComponent:0]; // 0 -- либмана, 1 -- зейделя
    long variable = [variablePicker selectedRowInComponent:0];   // 0 -- X, 1 -- Y
    
    HUD.mode = MBProgressHUDAnimationFade;
    HUD.labelText = @"Идет расчет";
    [HUD show:YES];
    dispatch_async(dispatch_get_global_queue( DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^(void){
        double **U;
//        :(int) Nx :(int) Ny :(double) bx :(double) by :(double) c :(double) hx :(double) hy :(double) alpha1 :(double) betta1 :(double) alpha2 :(double) betta2 :(double) alpha3 :(double) betta3 :(double) alpha4 :(double) betta4{
        U = [self elliptic_liebmann:Nx Ny:Ny bx:bx by:by c:c hx:hx hy:hy alpha1:alpha1 betta1:betta1 alpha2:alpha2 betta2:betta2 alpha3:alpha3 betta3:betta3 alpha4:alpha4 betta4:betta4 epsilon:epsilon];
//        U = [self elliptic_liebmann:Nx :Ny :2 :2 :4 :hx :hy :0 :1 :0 :1 :0 :1 :0 :1];
//        U = [self elliptic_seidel:Nx :Ny :2 :2 :4 :hx :hy :0 :1 :0 :1 :0 :1 :0 :1];
        
        switch (scheme) {
            case 0:{
                U = [self elliptic_liebmann:Nx Ny:Ny bx:bx by:by c:c hx:hx hy:hy alpha1:alpha1 betta1:betta1 alpha2:alpha2 betta2:betta2 alpha3:alpha3 betta3:betta3 alpha4:alpha4 betta4:betta4 epsilon:epsilon];
                break;
            }
            case 1:{
                U = [self elliptic_seidel:Nx Ny:Ny bx:bx by:by c:c hx:hx hy:hy alpha1:alpha1 betta1:betta1 alpha2:alpha2 betta2:betta2 alpha3:alpha3 betta3:betta3 alpha4:alpha4 betta4:betta4 epsilon:epsilon];
                break;
            }
            default:
                break;
        }
        
        NSMutableDictionary *dataDict;
        NSMutableArray *contentArrayErr;
        NSMutableDictionary *dataDictAnalytic;
        switch (variable) {
            case 0:{
                dataDict = [[NSMutableDictionary alloc] initWithCapacity:Nx];
                contentArrayErr = [[NSMutableArray alloc] init];
                dataDictAnalytic = [[NSMutableDictionary alloc] init];
                
                for (int i = 0; i <= Nx; ++i) {
                    NSMutableArray *contentArray = [[NSMutableArray alloc] init];
                    NSMutableArray *contentArrayAnalytic = [[NSMutableArray alloc] init];
                    
                    for (int j = 0; j <= Ny; ++j) {
                        id xAnalytic = [NSNumber numberWithDouble:(j)*hy];
                        id yAnalytic = [NSNumber numberWithDouble:[self analyticFunction:((i)*hx) y:j*hy]];
                        [contentArrayAnalytic addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:xAnalytic, @"x", yAnalytic, @"y", nil]];
                    }
                    //epsiHlon[k][n] = max_i(u[k][i] - reshenie[k][i])
                    id err = [NSNumber numberWithDouble:0.0];
                    for (int j = 0; j <= Ny; ++j) {
                        id x = [NSNumber numberWithDouble:j*hy];
                        id y = [NSNumber numberWithDouble:U[i][j]];
                        
                        err = [NSNumber numberWithDouble:([err doubleValue] > fabs([self analyticFunction:((i)*hx) y:j*hy] - U[i][j])) ? [err doubleValue] : fabs([self analyticFunction:((i)*hx) y:j*hy] - U[i][j])];
                        printf("%f	%f\n", [x doubleValue], [y doubleValue]);
                        [contentArray addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:x, @"x", y, @"y", nil]];
                    }
                    id time = [NSNumber numberWithDouble:i*hx];
                    printf("K = %f\n", [time doubleValue]);
                    [dataDict setObject:contentArray forKey:time];
                    [dataDictAnalytic setObject:contentArrayAnalytic forKey:time];
                    [contentArrayErr addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:time, @"x", err, @"y", nil]];
                }
                break;
            }
            case 1:{
                dataDict = [[NSMutableDictionary alloc] initWithCapacity:Ny];
                contentArrayErr = [[NSMutableArray alloc] init];
                dataDictAnalytic = [[NSMutableDictionary alloc] init];
                for (int i = 0; i <= Ny; ++i) {
                    NSMutableArray *contentArray = [[NSMutableArray alloc] init];
                    NSMutableArray *contentArrayAnalytic = [[NSMutableArray alloc] init];
                    
                    for (int j = 0; j <= Nx; ++j) {
                        id xAnalytic = [NSNumber numberWithDouble:j*hx];
                        id yAnalytic = [NSNumber numberWithDouble:[self analyticFunction:(j*hx) y:i*hy]];
                        [contentArrayAnalytic addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:xAnalytic, @"x", yAnalytic, @"y", nil]];
                    }
                    //epsiHlon[k][n] = max_i(u[k][i] - reshenie[k][i])
                    id err = [NSNumber numberWithDouble:0.0];
                    for (int j = 0; j <= Nx; ++j) {
                        id x = [NSNumber numberWithDouble:j*hx];
                        id y = [NSNumber numberWithDouble:U[j][i]];
                        
                        err = [NSNumber numberWithDouble:([err doubleValue] > fabs([self analyticFunction:(j*hx) y:i*hy] - U[j][i])) ? [err doubleValue] : fabs([self analyticFunction:(j*hx) y:i*hy] - U[j][i])];
                        printf("%f	%f\n", [x doubleValue], [y doubleValue]);
                        [contentArray addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:x, @"x", y, @"y", nil]];
                    }
                    id time = [NSNumber numberWithDouble:i*hy];
                    printf("K = %f\n", [time doubleValue]);
                    [dataDict setObject:contentArray forKey:time];
                    [dataDictAnalytic setObject:contentArrayAnalytic forKey:time];
                    [contentArrayErr addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:time, @"x", err, @"y", nil]];
                }

                break;
            }
            default:
                break;
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
    
    system = [[UIImageView alloc] initWithFrame:CGRectMake(36, 22, 248, 134)];
    system.image = [UIImage imageNamed:@"NM_sys3.png"];
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
    
    UILabel *bxLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 10, 35, 20)];
    bxLabel.text = @"bx =";
    [bg addSubview:bxLabel];
    bxField = [[NMFTextField alloc] initWithFrame:CGRectMake(55, 12, 50, 20)];
    bxField.text = @"2";
    [bg addSubview:bxField];
    
    UILabel *byLabel = [[UILabel alloc] initWithFrame:CGRectMake(130, 10, 35, 20)];
    byLabel.text = @"by =";
    [bg addSubview:byLabel];
    byField = [[NMFTextField alloc] initWithFrame:CGRectMake(165, 12, 50, 20)];
    byField.text = @"2";
    [bg addSubview:byField];
    
    UILabel *cLabel = [[UILabel alloc] initWithFrame:CGRectMake(220, 10, 30, 20)];
    cLabel.text = @"c =";
    [bg addSubview:cLabel];
    cField = [[NMFTextField alloc] initWithFrame:CGRectMake(250, 12, 50, 20)];
    cField.text = @"4";
    [bg addSubview:cField];
    
    UILabel *eLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 30, 43, 20)];
    eLabel.text = @"eps =";
    [bg addSubview:eLabel];
    epsilonField = [[NMFTextField alloc] initWithFrame:CGRectMake(63, 32, 50, 20)];
    epsilonField.text = @"0.0001";
    [bg addSubview:epsilonField];
    
    UILabel *lxLabel = [[UILabel alloc] initWithFrame:CGRectMake(130, 30, 30, 20)];
    lxLabel.text = @"lx =";
    [bg addSubview:lxLabel];
    lxField = [[NMFTextField alloc] initWithFrame:CGRectMake(160, 32, 50, 20)];
    lxField.text = @"1.57";
    [bg addSubview:lxField];
    
    UILabel *lyLabel = [[UILabel alloc] initWithFrame:CGRectMake(220, 30, 30, 20)];
    lyLabel.text = @"ly =";
    [bg addSubview:lyLabel];
    lyField = [[NMFTextField alloc] initWithFrame:CGRectMake(250, 32, 50, 20)];
    lyField.text = @"1.57";
    [bg addSubview:lyField];
    
    UILabel *alpha1Label = [[UILabel alloc] initWithFrame:CGRectMake(20, 60, 60, 20)];
    alpha1Label.text = @"α1 =";
    [bg addSubview:alpha1Label];
    alpha1Field = [[NMFTextField alloc] initWithFrame:CGRectMake(80, 62, 50, 20)];
    alpha1Field.text = @"0";
    [bg addSubview:alpha1Field];
    
    UILabel *betta1Label = [[UILabel alloc] initWithFrame:CGRectMake(150, 60, 60, 20)];
    betta1Label.text = @"β1 =";
    [bg addSubview:betta1Label];
    betta1Field = [[NMFTextField alloc] initWithFrame:CGRectMake(210, 62, 50, 20)];
    betta1Field.text = @"1";
    [bg addSubview:betta1Field];
    
    UILabel *alpha2Label = [[UILabel alloc] initWithFrame:CGRectMake(20, 80, 60, 20)];
    alpha2Label.text = @"α2 =";
    [bg addSubview:alpha2Label];
    alpha2Field = [[NMFTextField alloc] initWithFrame:CGRectMake(80, 82, 50, 20)];
    alpha2Field.text = @"0";
    [bg addSubview:alpha2Field];
    
    UILabel *betta2Label = [[UILabel alloc] initWithFrame:CGRectMake(150, 80, 60, 20)];
    betta2Label.text = @"β2 =";
    [bg addSubview:betta2Label];
    betta2Field = [[NMFTextField alloc] initWithFrame:CGRectMake(210, 82, 50, 20)];
    betta2Field.text = @"1";
    [bg addSubview:betta2Field];
    
    UILabel *alpha3Label = [[UILabel alloc] initWithFrame:CGRectMake(20, 100, 60, 20)];
    alpha3Label.text = @"α3 =";
    [bg addSubview:alpha3Label];
    alpha3Field = [[NMFTextField alloc] initWithFrame:CGRectMake(80, 102, 50, 20)];
    alpha3Field.text = @"0";
    [bg addSubview:alpha3Field];
    
    UILabel *betta3Label = [[UILabel alloc] initWithFrame:CGRectMake(150, 100, 60, 20)];
    betta3Label.text = @"β3 =";
    [bg addSubview:betta3Label];
    betta3Field = [[NMFTextField alloc] initWithFrame:CGRectMake(210, 102, 50, 20)];
    betta3Field.text = @"1";
    [bg addSubview:betta3Field];
   
    UILabel *alpha4Label = [[UILabel alloc] initWithFrame:CGRectMake(20, 120, 60, 20)];
    alpha4Label.text = @"α4 =";
    [bg addSubview:alpha4Label];
    alpha4Field = [[NMFTextField alloc] initWithFrame:CGRectMake(80, 122, 50, 20)];
    alpha4Field.text = @"0";
    [bg addSubview:alpha4Field];
    
    UILabel *betta4Label = [[UILabel alloc] initWithFrame:CGRectMake(150, 120, 60, 20)];
    betta4Label.text = @"β4 =";
    [bg addSubview:betta4Label];
    betta4Field = [[NMFTextField alloc] initWithFrame:CGRectMake(210, 122, 50, 20)];
    betta4Field.text = @"1";
    [bg addSubview:betta4Field];
    
    UILabel *NXLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 150, 230, 20)];
    NXLabel.text = @"Разбиений по OX:";
    [bg addSubview:NXLabel];
    NxField = [[NMFTextField alloc] initWithFrame:CGRectMake(250, 150, 50, 20)];
    NxField.text = @"10";
    [bg addSubview:NxField];
    
    UILabel *NYLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 170, 230, 20)];
    NYLabel.text = @"Разбиений по OY:";
    [bg addSubview:NYLabel];
    NyField = [[NMFTextField alloc] initWithFrame:CGRectMake(250, 170, 50, 20)];
    NyField.text = @"10";
    [bg addSubview:NyField];
    
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
    FLabel.text = @"f(x,y) =";
    [bg1 addSubview:FLabel];
    FField = [[UITextField alloc] initWithFrame:CGRectMake(84, 12, 230, 20)];
    FField.text = @"0";
    FField.delegate = self;
    [bg1 addSubview:FField];
    
    UILabel *Phi1Label = [[UILabel alloc] initWithFrame:CGRectMake(10, 30, 80, 20)];
    Phi1Label.text = @"φ1(y) =";
    [bg1 addSubview:Phi1Label];
    Phi1Field = [[UITextField alloc] initWithFrame:CGRectMake(84, 32, 230, 20)];
    Phi1Field.text = @"exp(-y)*cos(y)";
    Phi1Field.delegate = self;
    [bg1 addSubview:Phi1Field];
    
    UILabel *Phi2Label = [[UILabel alloc] initWithFrame:CGRectMake(10, 50, 80, 20)];
    Phi2Label.text = @"φ2(y) =";
    [bg1 addSubview:Phi2Label];
    Phi2Field = [[UITextField alloc] initWithFrame:CGRectMake(84, 52, 230, 20)];
    Phi2Field.text = @"0";
    Phi2Field.delegate = self;
    [bg1 addSubview:Phi2Field];
    
    UILabel *Phi3Label = [[UILabel alloc] initWithFrame:CGRectMake(10, 70, 80, 20)];
    Phi3Label.text = @"φ3(x) =";
    [bg1 addSubview:Phi3Label];
    Phi3Field = [[UITextField alloc] initWithFrame:CGRectMake(84, 72, 230, 20)];
    Phi3Field.text = @"exp(-x)*cos(x)";
    Phi3Field.delegate = self;
    [bg1 addSubview:Phi3Field];
    
    UILabel *Phi4Label = [[UILabel alloc] initWithFrame:CGRectMake(10, 90, 80, 20)];
    Phi4Label.text = @"φ4(x) =";
    [bg1 addSubview:Phi4Label];
    Phi4Field = [[UITextField alloc] initWithFrame:CGRectMake(84, 92, 230, 20)];
    Phi4Field.text = @"0";
    Phi4Field.delegate = self;
    [bg1 addSubview:Phi4Field];
    
    UILabel *ULabel = [[UILabel alloc] initWithFrame:CGRectMake(10, 160, 80, 20)];
    ULabel.text = @"U(x,t) =";
    [bg1 addSubview:ULabel];
    UField = [[UITextField alloc] initWithFrame:CGRectMake(84, 162, 230, 20)];
    UField.text = @"exp(-x-y)cos(x)cos(y)";
    UField.delegate = self;
    [bg1 addSubview:UField];
    
    //bg 2
    bg2 = [[UIView alloc] initWithFrame:CGRectMake(330, 160, 320, self.view.frame.size.height - 180)];
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
    order.text = @"Схема:";
    order.textAlignment = NSTextAlignmentCenter;
    
    UILabel *scheme = [[UILabel alloc] initWithFrame:CGRectMake(0, 160, 320, 15)];
    scheme.text = @"Зафиксировать:";
    scheme.textAlignment = NSTextAlignmentCenter;
    
    schemePicker = [[UiPickerViewSchemeElliptic alloc] initWithFrame:CGRectMake(0, 10, 320, 60)];
    variablePicker = [[UIPickerViewVariableElliptic alloc] initWithFrame:CGRectMake(0, 160, 320, 60)];
    
    [bg2 addSubview:order];
    [bg2 addSubview:scheme];
    [bg2 addSubview:variablePicker];
    [bg2 addSubview:schemePicker];
    
    UIButton *nextButton;
    
    if ([[UIScreen mainScreen] bounds].size.height < 500) {
        nextButton = [[UIButton alloc] initWithFrame:CGRectMake(230, 260, 60, 40)];
    }
    else {
        nextButton = [[UIButton alloc] initWithFrame:CGRectMake(230, 320, 60, 40)];
    }
    
    
    UILabel *blab2 = [[UILabel alloc] initWithFrame:CGRectMake(0, 0, 60, 20)];
    blab2.text = @"Далее";
    [nextButton addSubview:blab2];
    [nextButton addTarget:self action:@selector(proceedScheme:) forControlEvents:UIControlEventTouchUpInside];
    
    [bg2 addSubview:nextButton];
    
    
    UIButton *backButton;
    if ([[UIScreen mainScreen] bounds].size.height < 500) {
        backButton = [[UIButton alloc] initWithFrame:CGRectMake(40, 260, 60, 40)];
    }
    else {
        backButton = [[UIButton alloc] initWithFrame:CGRectMake(40, 320, 60, 40)];
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
