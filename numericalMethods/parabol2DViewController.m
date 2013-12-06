//
//  parabol2DViewController.m
//  numericalMethods
//
//  Created by Ilya on 05.12.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import "parabol2DViewController.h"

#import "CorePlotViewController.h"

#define ALPHABET3 @"qwertyuiopasdfghjklzxcvbnm1234567890-+*/()$"


double a,b,c,alpha1,alpha2,alpha3,alpha4,betta1,betta2,betta3,betta4,l,lx,ly;
long startApproximation;

@interface parabol2DViewController ()

@end

@implementation parabol2DViewController

@synthesize HUD;
@synthesize system;


double * processTridiagonalMatrixP(double *x, const size_t N, const double *a, const double *b, double *c) {
    
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
- (double)analyticFunction:(double)x y:(double)y t:(double) t{
//    return sin(x)*sin(y)*sin(t);
    NSArray *values = @[[NSNumber numberWithDouble:x], [NSNumber numberWithDouble:y], [NSNumber numberWithDouble:t], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"y", @"t", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionRealFunc evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) f:(double) x y:(double) y t:(double) t{
//    return sin(x)*sin(y)*(cos(t) + 2*sin(t));
    NSArray *values = @[[NSNumber numberWithDouble:x], [NSNumber numberWithDouble:y], [NSNumber numberWithDouble:t], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"y", @"t", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionF evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) psi:(double) x y:(double)y{
    NSArray *values = @[[NSNumber numberWithDouble:x], [NSNumber numberWithDouble:y], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"y", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPsi evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}


-(double) phi1:(double) y t:(double) t{
    //return 0;
    NSArray *values = @[[NSNumber numberWithDouble:y], [NSNumber numberWithDouble:t], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"y", @"t", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPhi1 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) phi2:(double) y t:(double) t{
//    return -sin(y)*sin(t);
    NSArray *values = @[[NSNumber numberWithDouble:y], [NSNumber numberWithDouble:t], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"y", @"t", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPhi2 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) phi3:(double) x t:(double) t{
//    return 0;
    NSArray *values = @[[NSNumber numberWithDouble:x], [NSNumber numberWithDouble:t], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"t", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPhi3 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

-(double) phi4:(double) x t:(double) t{
//    return -sin(x)*sin(t);
    NSArray *values = @[[NSNumber numberWithDouble:x], [NSNumber numberWithDouble:t], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"t", @"a", @"b", @"c"];
    
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    return [[expressionPhi4 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] doubleValue];
}

#pragma mark - fractional step method

-(double ***)fractional_step_method:(int) Nx Ny:(int) Ny K:(int)K  tau:(double)tau a:(double) a b:(double) b hx:(double) hx hy:(double) hy alpha1:(double) alpha1 betta1:(double) betta1 alpha2:(double) alpha2 betta2:(double) betta2 alpha3:(double) alpha3 betta3:(double) betta3 alpha4:(double) alpha4 betta4:(double) betta4{
    
    double ***U = (double ***)malloc((K +1) * sizeof(double **));
    for (int i = 0; i < (K +1); i++){
        U[i] = (double **)malloc((Nx+1) * sizeof(double *));
        for (int j = 0; j < (Nx +1); j++) {
            U[i][j] = (double *)malloc((Ny+1) * sizeof(double));
        }
    }

    double **tmp = (double **)malloc((Nx+1) * sizeof(double *));
    for (int j = 0; j < (Nx +1); j++) {
        tmp[j] = (double *)malloc((Ny+1) * sizeof(double));
    }

    
    double *lowerX = (double *)malloc((Nx+1) * sizeof(double));
    double *midX = (double *)malloc((Nx+1) * sizeof(double));
    double *upperX = (double *)malloc((Nx+1) * sizeof(double));
    double *answerX = (double *)malloc((Nx+1) * sizeof(double));
    
    memset(lowerX, 0, (Nx+1) * sizeof(double));
    memset(upperX, 0, (Nx+1) * sizeof(double));
    memset(midX, 0, (Nx+1) * sizeof(double));
    memset(answerX, 0, (Nx+1) * sizeof(double));

    double *lowerY = (double *)malloc((Ny+1) * sizeof(double));
    double *midY = (double *)malloc((Ny+1) * sizeof(double));
    double *upperY = (double *)malloc((Ny+1) * sizeof(double));
    double *answerY = (double *)malloc((Ny+1) * sizeof(double));
    
    memset(lowerY, 0, (Ny+1) * sizeof(double));
    memset(upperY, 0, (Ny+1) * sizeof(double));
    memset(midY, 0, (Ny+1) * sizeof(double));
    memset(answerY, 0, (Ny+1) * sizeof(double));
    
    
    for (int i = 0; i <= Nx; ++i) {
        for (int j = 0; j <= Ny; ++j) {
            U[0][i][j] = [self psi:i*hx y:j*hy];
        }
    }
    
    for (int k = 1; k <= K; ++k) {
        for (int j = 0; j <= Ny; ++j) {
            for (int i = 1; i < Nx; ++i) {
                lowerX[i] = -a/hx/hx;
                midX[i] = 1/tau + 2*a/hx/hx;
                upperX[i] = -a/hx/hx;
                answerX[i] = U[k-1][i][j]/tau + 0.5*[self f:i*hx y:j*hy t:(k+0.5)*tau];
            }
            lowerX[0] = 0;
            midX[0] = betta1 - alpha1/hx;
            upperX[0] = alpha1/hx;
            answerX[0] = [self phi1:j*hy t:(k+0.5)*tau];
            
            lowerX[Nx] = -alpha2/hx;
            midX[Nx] = betta2 + alpha2/hx;
            upperX[Nx] = 0;
            answerX[Nx] = [self phi2:j*hy t:(k+0.5)*tau];
            

            
            double * res = processTridiagonalMatrixP(answerX, Nx+1, lowerX, midX, upperX);
            for (int i = 0; i <= Nx; ++i) {
                tmp[i][j] = res[i];
            }
        }
        
        for (int i = 0; i <= Nx; ++i) {
            for (int j = 1; j < Ny; ++j) {
                lowerY[j] = -b/hy/hy;
                midY[j] = 1/tau + 2*b/hy/hy;
                upperY[j] = -b/hy/hy;
                answerY[j] = tmp[i][j]/tau + 0.5*[self f:i*hx y:j*hy t:(k+1)*tau];
            }
            lowerY[0] = 0;
            midY[0] = betta3 - alpha3/hy;
            upperY[0] = alpha3/hy;
            answerY[0] = [self phi3:i*hx t:(k+1)*tau];
            
            lowerY[Ny] = -alpha4/hy;
            midY[Ny] = betta4 + alpha4/hy;
            upperY[Ny] = 0;
            answerY[Ny] = [self phi4:i*hx t:(k+1)*tau];
            
            double * res = processTridiagonalMatrixP(answerY, Ny+1, lowerY, midY, upperY);
            for (int j = 0; j <= Ny; ++j) {
                U[k][i][j] = res[j];
            }
        }

    }
    
    return U;
}

#pragma mark - alternating direction method

-(double ***)alternating_direction_method:(int) Nx Ny:(int) Ny K:(int)K  tau:(double)tau a:(double) a b:(double) b hx:(double) hx hy:(double) hy alpha1:(double) alpha1 betta1:(double) betta1 alpha2:(double) alpha2 betta2:(double) betta2 alpha3:(double) alpha3 betta3:(double) betta3 alpha4:(double) alpha4 betta4:(double) betta4{
    
    double ***U = (double ***)malloc((K +1) * sizeof(double **));
    for (int i = 0; i < (K +1); i++){
        U[i] = (double **)malloc((Nx+1) * sizeof(double *));
        for (int j = 0; j < (Nx +1); j++) {
            U[i][j] = (double *)malloc((Ny+1) * sizeof(double));
        }
    }
    
    double **tmp = (double **)malloc((Nx+1) * sizeof(double *));
    for (int j = 0; j < (Nx +1); j++) {
        tmp[j] = (double *)malloc((Ny+1) * sizeof(double));
    }
    
    
    double *lowerX = (double *)malloc((Nx+1) * sizeof(double));
    double *midX = (double *)malloc((Nx+1) * sizeof(double));
    double *upperX = (double *)malloc((Nx+1) * sizeof(double));
    double *answerX = (double *)malloc((Nx+1) * sizeof(double));
    
    memset(lowerX, 0, (Nx+1) * sizeof(double));
    memset(upperX, 0, (Nx+1) * sizeof(double));
    memset(midX, 0, (Nx+1) * sizeof(double));
    memset(answerX, 0, (Nx+1) * sizeof(double));
    
    double *lowerY = (double *)malloc((Ny+1) * sizeof(double));
    double *midY = (double *)malloc((Ny+1) * sizeof(double));
    double *upperY = (double *)malloc((Ny+1) * sizeof(double));
    double *answerY = (double *)malloc((Ny+1) * sizeof(double));
    
    memset(lowerY, 0, (Ny+1) * sizeof(double));
    memset(upperY, 0, (Ny+1) * sizeof(double));
    memset(midY, 0, (Ny+1) * sizeof(double));
    memset(answerY, 0, (Ny+1) * sizeof(double));
    
    
    for (int i = 0; i <= Nx; ++i) {
        for (int j = 0; j <= Ny; ++j) {
            U[0][i][j] = [self psi:i*hx y:j*hy];
        }
    }
    
    for (int k = 1; k <= K; ++k) {
        for (int j = 1; j < Ny; ++j) {
            for (int i = 1; i < Nx; ++i) {
                lowerX[i] = -a/hx/hx;
                midX[i] = 2/tau + 2*a/hx/hx;
                upperX[i] = -a/hx/hx;
                answerX[i] = U[k-1][i][j]/tau*2 + b/hy/hy*(U[k-1][i][j+1] - 2*U[k-1][i][j] + U[k-1][i][j-1]) + [self f:i*hx y:j*hy t:(k+0.5)*tau];
            }
            lowerX[0] = 0;
            midX[0] = betta1 - alpha1/hx;
            upperX[0] = alpha1/hx;
            answerX[0] = [self phi1:j*hy t:(k+0.5)*tau];
            
            lowerX[Nx] = -alpha2/hx;
            midX[Nx] = betta2 + alpha2/hx;
            upperX[Nx] = 0;
            answerX[Nx] = [self phi2:j*hy t:(k+0.5)*tau];
            
            
            
            double * res = processTridiagonalMatrixP(answerX, Nx+1, lowerX, midX, upperX);
            for (int i = 0; i <= Nx; ++i) {
                tmp[i][j] = res[i];
            }
        }
        
        for (int i = 0; i <= Nx; ++i) {
            tmp[i][0] = U[k-1][i][0];
            tmp[i][Ny] = U[k-1][i][Ny];
        }
        
        for (int i = 1; i < Nx; ++i) {
            for (int j = 1; j < Ny; ++j) {
                lowerY[j] = -b/hy/hy;
                midY[j] = 2/tau + 2*b/hy/hy;
                upperY[j] = -b/hy/hy;
                answerY[j] = 2*tmp[i][j]/tau + a/hx/hx*(tmp[i+1][j] - 2*tmp[i][j] + tmp[i-1][j]) + [self f:i*hx y:j*hy t:(k+1)*tau];
            }
            lowerY[0] = 0;
            midY[0] = betta3 - alpha3/hy;
            upperY[0] = alpha3/hy;
            answerY[0] = [self phi3:i*hx t:(k+1)*tau];
            
            lowerY[Ny] = -alpha4/hy;
            midY[Ny] = betta4 + alpha4/hy;
            upperY[Ny] = 0;
            answerY[Ny] = [self phi4:i*hx t:(k+1)*tau];
            
            double * res = processTridiagonalMatrixP(answerY, Ny+1, lowerY, midY, upperY);
            for (int j = 0; j <= Ny; ++j) {
                U[k][i][j] = res[j];
            }
        }
        
        for (int j = 0; j <= Ny; ++j) {
            U[k][0][j] = tmp[0][j];
            U[k][Nx][j] = tmp[Nx][j];
        }
        
        double den = betta1 - alpha1/hx;
        
        U[k][0][0] = -U[k][1][0]*(alpha1 / hx / den) + [self phi1:0 t:(k+1.0)*tau] / den;
        U[k][0][Ny] = -U[k][1][Ny]*(alpha1 / hx / den) +  [self phi1:Ny*hy t:(k+1.0)*tau] / den;
        
        den = betta2 + alpha2/hx;
        
        U[k][Nx][0] = -U[k][Nx-1][0]*(-alpha2/hx/den) + [self phi2:0 t:(k+1.0)*tau] / den;
        U[k][Nx][Ny] = -U[k][Nx-1][Ny]*(-alpha2/hx/den) + [self phi2:Ny*hy t:(k+1.0)*tau] / den;
        
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
    
//    for (int i = 0; i <= [NxField text].intValue; ++i){
//        NSNumber *item = [[NSNumber alloc] initWithDouble:i * ([[lxField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue] / [NxField.text intValue])];
//        [_fixedXArray addObject:item];
//    }
//    for (int i = 0; i <= [NyField text].intValue; ++i){
//        NSNumber *item = [[NSNumber alloc] initWithDouble:i * ([[lyField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue] / [NyField.text intValue])];
//        [_fixedYArray addObject:item];
//    }
    for (int i = 0; i <= [KField text].intValue; ++i){
        NSNumber *item = [[NSNumber alloc] initWithDouble:i * ([TField text].doubleValue /  [KField text].doubleValue)];
        [_fixedTArray addObject:item];
    }
    
    fixedVariablePicker = [[UIPickerView alloc] initWithFrame:CGRectMake(100, 280, 100, 100)];
    fixedVariablePicker.showsSelectionIndicator = YES;
    fixedVariablePicker.delegate = self;
    [bg2 addSubview:fixedVariablePicker];

    
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
    } completion:^(BOOL finished) {
        [fixedVariablePicker removeFromSuperview];
//        [_fixedXArray removeAllObjects];
//        [_fixedYArray removeAllObjects];
        [_fixedTArray removeAllObjects];
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
    NSString* psi = [self proceedString:PsiField.text];

    
    expressionRealFunc = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:realFunc error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPhi1 = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:phi1 error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPhi2 = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:phi2 error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPhi3 = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:phi3 error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPhi4 = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:phi4 error:&error] error:&error] parsedExpressionWithError:&error];
    expressionF = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:f error:&error] error:&error] parsedExpressionWithError:&error];
    expressionPsi = [[DDParser parserWithTokenizer:[DDMathStringTokenizer tokenizerWithString:psi error:&error] error:&error] parsedExpressionWithError:&error];
    
    NSArray *values = @[[NSNumber numberWithDouble:1], [NSNumber numberWithDouble:1], [NSNumber numberWithDouble:1], [NSNumber numberWithDouble:a], [NSNumber numberWithDouble:b], [NSNumber numberWithDouble:c]];
    NSArray *keys = @[@"x", @"y", @"t", @"a", @"b", @"c"];
    NSDictionary *variableSubstitutions = [NSDictionary dictionaryWithObjects:values forKeys:keys];
    
    if ([expressionF evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция f(x,y,t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPhi1 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция φ1(y,t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPhi2 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция φ2(y,t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPhi3 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция φ3(x,t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    
    if ([expressionPhi4 evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция φ4(x,t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    if ([expressionPsi evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция ψ(x,y).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }

    if ([expressionRealFunc evaluateWithSubstitutions:variableSubstitutions evaluator:evaluator error:nil] == nil) {
        HUD.mode = MBProgressHUDModeText;
        HUD.labelText = @"Ошибка! Функция U(x,y,t).";
        HUD.detailsLabelText =  @"Невозможно распознать функцию.";
        [HUD show:YES];
        [HUD hide:YES afterDelay:3];
        return;
    }
    
    int Nx = [NxField.text intValue];
    int Ny = [NyField.text intValue];
    int K = [KField.text intValue];
    double T = [TField.text doubleValue];
    lx = [[lxField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    ly = [[lyField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    
    double hx = lx / Nx;
    double hy = ly / Ny;
    double tau = T / K;
    
    a = [[aField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    b = [[bField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    c = [[cField.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    alpha1 = [[alpha1Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    betta1 = [[betta1Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    alpha2 = [[alpha2Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    betta2 = [[betta2Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    alpha3 = [[alpha3Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    betta3 = [[betta3Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    alpha4 = [[alpha4Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    betta4 = [[betta4Field.text stringByReplacingOccurrencesOfString:@"," withString:@"."] doubleValue];
    long scheme = [schemePicker selectedRowInComponent:0]; // 0 -- дробных шагов, 1 -- переменных направлений
    long variable = [variablePicker selectedRowInComponent:0];   // 0 --  T, 1 -- X, 2 -- Y
    
    HUD.mode = MBProgressHUDAnimationFade;
    HUD.labelText = @"Идет расчет";
    [HUD show:YES];
    dispatch_async(dispatch_get_global_queue( DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^(void){
        double ***U;
        switch (scheme) {
            case 0:{
                U = [self fractional_step_method:Nx Ny:Ny K:K tau:tau a:a b:b hx:hx hy:hy alpha1:alpha1 betta1:betta1 alpha2:alpha2 betta2:betta2 alpha3:alpha3 betta3:betta3 alpha4:alpha4 betta4:betta4];
                break;
            }
            case 1:{
                U = [self alternating_direction_method:Nx Ny:Ny K:K tau:tau a:a b:b hx:hx hy:hy alpha1:alpha1 betta1:betta1 alpha2:alpha2 betta2:betta2 alpha3:alpha3 betta3:betta3 alpha4:alpha4 betta4:betta4];
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
                int TIME = (int)(_fixedT / tau);
                dataDict = [[NSMutableDictionary alloc] initWithCapacity:Nx];
                contentArrayErr = [[NSMutableArray alloc] init];
                dataDictAnalytic = [[NSMutableDictionary alloc] init];
                
                for (int i = 0; i <= Nx; ++i) {
                    NSMutableArray *contentArray = [[NSMutableArray alloc] init];
                    NSMutableArray *contentArrayAnalytic = [[NSMutableArray alloc] init];
                    
                    for (int j = 0; j <= Ny; ++j) {
                        id xAnalytic = [NSNumber numberWithDouble:(j)*hy];
                        id yAnalytic = [NSNumber numberWithDouble:[self analyticFunction:((i)*hx) y:j*hy t:TIME*tau]];
                        [contentArrayAnalytic addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:xAnalytic, @"x", yAnalytic, @"y", nil]];
                    }
                    id err = [NSNumber numberWithDouble:0.0];
                    for (int j = 0; j <= Ny; ++j) {
                        id x = [NSNumber numberWithDouble:j*hy];
                        id y = [NSNumber numberWithDouble:U[TIME][i][j]];
                        
                        err = [NSNumber numberWithDouble:([err doubleValue] > fabs([self analyticFunction:((i)*hx) y:j*hy t:TIME*tau] - U[TIME][i][j])) ? [err doubleValue] : fabs([self analyticFunction:((i)*hx) y:j*hy t:TIME*tau] - U[TIME][i][j])];
                        printf("%f	%f - %f\n", [x doubleValue], [y doubleValue], [self analyticFunction:((i)*hx) y:j*hy t:TIME*tau]);
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
                int TIME = (int)(_fixedT / tau);
                dataDict = [[NSMutableDictionary alloc] initWithCapacity:Ny];
                contentArrayErr = [[NSMutableArray alloc] init];
                dataDictAnalytic = [[NSMutableDictionary alloc] init];
                for (int i = 0; i <= Ny; ++i) {
                    NSMutableArray *contentArray = [[NSMutableArray alloc] init];
                    NSMutableArray *contentArrayAnalytic = [[NSMutableArray alloc] init];
                    
                    for (int j = 0; j <= Nx; ++j) {
                        id xAnalytic = [NSNumber numberWithDouble:j*hx];
                        id yAnalytic = [NSNumber numberWithDouble:[self analyticFunction:(j*hx) y:i*hy t:TIME*tau]];
                        [contentArrayAnalytic addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:xAnalytic, @"x", yAnalytic, @"y", nil]];
                    }
                    //epsiHlon[k][n] = max_i(u[k][i] - reshenie[k][i])
                    id err = [NSNumber numberWithDouble:0.0];
                    for (int j = 0; j <= Nx; ++j) {
                        id x = [NSNumber numberWithDouble:j*hx];
                        id y = [NSNumber numberWithDouble:U[TIME][j][i]];
                        
                        err = [NSNumber numberWithDouble:([err doubleValue] > fabs([self analyticFunction:(j*hx) y:i*hy t:TIME*tau] - U[TIME][j][i])) ? [err doubleValue] : fabs([self analyticFunction:(j*hx) y:i*hy t:TIME*tau] - U[TIME][j][i])];
                        printf("%f	%f\n", [x doubleValue], [y doubleValue]);
                        [contentArray addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:x, @"x", y, @"y", nil]];
                    }
                    id time = [NSNumber numberWithDouble:i*hy];
                    printf("K = %f\n", [time doubleValue]);
                    [dataDict setObject:contentArray forKey:time];
                    [dataDictAnalytic setObject:contentArrayAnalytic forKey:time];
                    [contentArrayErr addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:time, @"x", err, @"y", nil]];
                }
//            case 2:{
//                int X = (int)(_fixedX / (lx / [[NxField text] doubleValue]));
//                int Y = (int)(_fixedY / (ly / [[NyField text] doubleValue]));
//
//                dataDict = [[NSMutableDictionary alloc] initWithCapacity:K];
//                contentArrayErr = [[NSMutableArray alloc] init];
//                dataDictAnalytic = [[NSMutableDictionary alloc] init];
//                
//                for (int i = 0; i <= K; ++i) {
//                    NSMutableArray *contentArray = [[NSMutableArray alloc] init];
//                    NSMutableArray *contentArrayAnalytic = [[NSMutableArray alloc] init];
//                    
//                    for (int j = 0; j <= K; ++j) {
//                        id xAnalytic = [NSNumber numberWithDouble:_fixedX];
//                        id yAnalytic = [NSNumber numberWithDouble:[self analyticFunction:_fixedX y:_fixedY t:j*tau]];
//                        [contentArrayAnalytic addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:xAnalytic, @"x", yAnalytic, @"y", nil]];
//                    }
//                    id err = [NSNumber numberWithDouble:0.0];
//                    for (int j = 0; j <= K; ++j) {
//                        id x = [NSNumber numberWithDouble:_fixedX];
//                        id y = [NSNumber numberWithDouble:U[j][X][Y]];
//                        
//                        err = [NSNumber numberWithDouble:([err doubleValue] > fabs([self analyticFunction:_fixedX y:_fixedY t:j*tau] - U[j][X][Y])) ? [err doubleValue] : fabs([self analyticFunction:_fixedX y:_fixedY t:j*tau] - U[j][X][Y])];
//                        printf("%f	%f - %f\n", [x doubleValue], [y doubleValue], [self analyticFunction:_fixedX y:_fixedY t:j*tau]);
//                        [contentArray addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:x, @"x", y, @"y", nil]];
//                    }
//                    id time = [NSNumber numberWithDouble:i*tau];
//                    printf("K = %f\n", [time doubleValue]);
//                    [dataDict setObject:contentArray forKey:time];
//                    [dataDictAnalytic setObject:contentArrayAnalytic forKey:time];
//                    [contentArrayErr addObject:[NSMutableDictionary dictionaryWithObjectsAndKeys:time, @"x", err, @"y", nil]];
//                }
//                break;
//            }
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
                scrollBg.contentOffset = CGPointMake(0, 170);
            else
                scrollBg.contentOffset = CGPointMake(0, 70);
        }];
    }
    else if (bg2.frame.origin.x == 0){
        [UIView animateWithDuration:0.3 animations:^{
            if(self.view.frame.size.height < 500){
                CGRect frame = bg2.frame;
                frame.origin.y -= 190;
                bg2.frame = frame;
            }
            else{
                CGRect frame = bg2.frame;
                frame.origin.y -= 90;
                bg2.frame = frame;
            }
        }];
    }

}

-(void)keyboardHidden:(NSNotification *)note{
    if (bg1.frame.origin.x > 0){
        [UIView animateWithDuration:0.3 animations:^{
            scrollBg.contentOffset = CGPointMake(0, 0);
        }];
    }else if (bg2.frame.origin.x == 0){
        [UIView animateWithDuration:0.3 animations:^{
            CGRect frame = bg2.frame;
            frame.origin.y = 22;
            bg2.frame = frame;
        }];
    }
}

#pragma mark - viewController delegate
- (void)viewDidLoad
{
    [super viewDidLoad];
    
    _fixedXArray = [[NSMutableArray alloc] init];
    _fixedYArray = [[NSMutableArray alloc] init];
    _fixedTArray = [[NSMutableArray alloc] init];
    
    system = [[UIImageView alloc] initWithFrame:CGRectMake(63, 22, 194, 134)];
    system.image = [UIImage imageNamed:@"NM_sys4.png"];
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
    bxLabel.text = @"a =";
    [bg addSubview:bxLabel];
    aField = [[NMFTextField alloc] initWithFrame:CGRectMake(55, 12, 50, 20)];
    aField.text = @"1";
    [bg addSubview:aField];
    
    UILabel *byLabel = [[UILabel alloc] initWithFrame:CGRectMake(130, 10, 35, 20)];
    byLabel.text = @"b =";
    [bg addSubview:byLabel];
    bField = [[NMFTextField alloc] initWithFrame:CGRectMake(165, 12, 50, 20)];
    bField.text = @"1";
    [bg addSubview:bField];
    
    UILabel *cLabel = [[UILabel alloc] initWithFrame:CGRectMake(220, 10, 30, 20)];
    cLabel.text = @"c =";
    [bg addSubview:cLabel];
    cField = [[NMFTextField alloc] initWithFrame:CGRectMake(250, 12, 50, 20)];
    cField.text = @"0";
    [bg addSubview:cField];
    
    UILabel *lxLabel = [[UILabel alloc] initWithFrame:CGRectMake(130, 30, 30, 20)];
    lxLabel.text = @"lx =";
    [bg addSubview:lxLabel];
    lxField = [[NMFTextField alloc] initWithFrame:CGRectMake(160, 32, 50, 20)];
    lxField.text = @"3.14";
    [bg addSubview:lxField];
    
    UILabel *lyLabel = [[UILabel alloc] initWithFrame:CGRectMake(220, 30, 30, 20)];
    lyLabel.text = @"ly =";
    [bg addSubview:lyLabel];
    lyField = [[NMFTextField alloc] initWithFrame:CGRectMake(250, 32, 50, 20)];
    lyField.text = @"3.14";
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
    alpha2Field.text = @"1";
    [bg addSubview:alpha2Field];
    
    UILabel *betta2Label = [[UILabel alloc] initWithFrame:CGRectMake(150, 80, 60, 20)];
    betta2Label.text = @"β2 =";
    [bg addSubview:betta2Label];
    betta2Field = [[NMFTextField alloc] initWithFrame:CGRectMake(210, 82, 50, 20)];
    betta2Field.text = @"0";
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
    alpha4Field.text = @"1";
    [bg addSubview:alpha4Field];
    
    UILabel *betta4Label = [[UILabel alloc] initWithFrame:CGRectMake(150, 120, 60, 20)];
    betta4Label.text = @"β4 =";
    [bg addSubview:betta4Label];
    betta4Field = [[NMFTextField alloc] initWithFrame:CGRectMake(210, 122, 50, 20)];
    betta4Field.text = @"0";
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
    
    UILabel *KLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 190, 230, 20)];
    KLabel.text = @"Разбиений по времени:";
    [bg addSubview:KLabel];
    KField = [[NMFTextField alloc] initWithFrame:CGRectMake(250, 190, 50, 20)];
    KField.text = @"100";
    [bg addSubview:KField];
    
    UILabel *TLabel = [[UILabel alloc] initWithFrame:CGRectMake(20, 210, 230, 20)];
    TLabel.text = @"Финальное время T:";
    [bg addSubview:TLabel];
    TField = [[NMFTextField alloc] initWithFrame:CGRectMake(250, 210, 50, 20)];
    TField.text = @"1";
    [bg addSubview:TField];
    
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
    FLabel.text = @"f(x,y,t) =";
    [bg1 addSubview:FLabel];
    FField = [[UITextField alloc] initWithFrame:CGRectMake(84, 12, 230, 20)];
    FField.text = @"sin(x)*sin(y)*(cos(t) + 2*sin(t))";
    FField.delegate = self;
    [bg1 addSubview:FField];
    
    UILabel *Phi1Label = [[UILabel alloc] initWithFrame:CGRectMake(10, 30, 80, 20)];
    Phi1Label.text = @"φ1(y,t) =";
    [bg1 addSubview:Phi1Label];
    Phi1Field = [[UITextField alloc] initWithFrame:CGRectMake(84, 32, 230, 20)];
    Phi1Field.text = @"0";
    Phi1Field.delegate = self;
    [bg1 addSubview:Phi1Field];
    
    UILabel *Phi2Label = [[UILabel alloc] initWithFrame:CGRectMake(10, 50, 80, 20)];
    Phi2Label.text = @"φ2(y,t) =";
    [bg1 addSubview:Phi2Label];
    Phi2Field = [[UITextField alloc] initWithFrame:CGRectMake(84, 52, 230, 20)];
    Phi2Field.text = @"-sin(y)*sin(t)";
    Phi2Field.delegate = self;
    [bg1 addSubview:Phi2Field];
    
    UILabel *Phi3Label = [[UILabel alloc] initWithFrame:CGRectMake(10, 70, 80, 20)];
    Phi3Label.text = @"φ3(x,t) =";
    [bg1 addSubview:Phi3Label];
    Phi3Field = [[UITextField alloc] initWithFrame:CGRectMake(84, 72, 230, 20)];
    Phi3Field.text = @"0";
    Phi3Field.delegate = self;
    [bg1 addSubview:Phi3Field];
    
    UILabel *Phi4Label = [[UILabel alloc] initWithFrame:CGRectMake(10, 90, 80, 20)];
    Phi4Label.text = @"φ4(x,t) =";
    [bg1 addSubview:Phi4Label];
    Phi4Field = [[UITextField alloc] initWithFrame:CGRectMake(84, 92, 230, 20)];
    Phi4Field.text = @"-sin(x)*sin(t)";
    Phi4Field.delegate = self;
    [bg1 addSubview:Phi4Field];
    
    UILabel *PsiLabel = [[UILabel alloc] initWithFrame:CGRectMake(10, 110, 80, 20)];
    PsiLabel.text = @"ψ(x,y) =";
    [bg1 addSubview:PsiLabel];
    PsiField = [[UITextField alloc] initWithFrame:CGRectMake(84, 112, 230, 20)];
    PsiField.text = @"0";
    PsiField.delegate = self;
    [bg1 addSubview:PsiField];

    
    UILabel *ULabel = [[UILabel alloc] initWithFrame:CGRectMake(10, 160, 80, 20)];
    ULabel.text = @"U(x,y,t) =";
    [bg1 addSubview:ULabel];
    UField = [[UITextField alloc] initWithFrame:CGRectMake(84, 162, 230, 20)];
    UField.text = @"sin(x)*sin(y)*sin(t)";
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
    order.text = @"Метод:";
    order.textAlignment = NSTextAlignmentCenter;
    
    UILabel *scheme = [[UILabel alloc] initWithFrame:CGRectMake(0, 160, 320, 15)];
    scheme.text = @"Варировать:";
    scheme.textAlignment = NSTextAlignmentCenter;
    
    schemePicker = [[UIPickerViewSchemeParabolic2d alloc] initWithFrame:CGRectMake(0, 10, 320, 60)];
    variablePicker = [[UIPickerViewVariableParabolic2d alloc] initWithFrame:CGRectMake(0, 160, 320, 60)];
    
    [bg2 addSubview:order];
    [bg2 addSubview:scheme];
    [bg2 addSubview:variablePicker];
    [bg2 addSubview:schemePicker];
    
    UILabel *fixedTLabel = [[UILabel alloc] initWithFrame:CGRectMake(30, 350, 50, 20)];
    fixedTLabel.text = @"T:";
    fixedTLabel.textAlignment = NSTextAlignmentCenter;
    [bg2 addSubview:fixedTLabel];
    
    
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

- (BOOL)textField:(UITextField *)textField shouldChangeCharactersInRange:(NSRange)range replacementString:(NSString *)string{
    string = [string lowercaseString];
    NSCharacterSet *cs = [[NSCharacterSet characterSetWithCharactersInString:ALPHABET3] invertedSet];
    NSString *filtered = [[string componentsSeparatedByCharactersInSet:cs] componentsJoinedByString:@""];
    return (([string isEqualToString:filtered]));
}

#pragma mark - picker
- (void)pickerView:(UIPickerView *)pickerView didSelectRow: (NSInteger)row inComponent:(NSInteger)component {
    _fixedT = [[_fixedTArray objectAtIndex:row] doubleValue];

//    switch (component) {
//        case 0:{
//            _fixedX = [[_fixedXArray objectAtIndex:row] doubleValue];
//            break;
//        }
//        case 1:{
//            _fixedY = [[_fixedYArray objectAtIndex:row] doubleValue];
//            break;
//        }
//        case 2:{
//            _fixedT = [[_fixedTArray objectAtIndex:row] doubleValue];
//            break;
//        }
//        default:
//            break;
//    }
}


// tell the picker how many rows are available for a given component
- (NSInteger)pickerView:(UIPickerView *)pickerView numberOfRowsInComponent:(NSInteger)component {
    return [_fixedTArray count];
//
//    switch (component) {
//        case 0:{
//            return [_fixedXArray count];
//            break;
//        }
//        case 1:{
//            return [_fixedYArray count];
//            break;
//        }
//        case 2:{
//            return [_fixedTArray count];
//            break;
//        }
//        default:
//            break;
//    }
//    return 0;
}

// tell the picker how many components it will have
- (NSInteger)numberOfComponentsInPickerView:(UIPickerView *)pickerView {
    return 1;
}

// tell the picker the title for a given component
- (NSString *)pickerView:(UIPickerView *)pickerView titleForRow:(NSInteger)row forComponent:(NSInteger)component {
    return [NSString stringWithFormat:@"%.3f", [[_fixedTArray objectAtIndex:row] floatValue]];
//    switch (component) {
//        case 0:{
//            return [NSString stringWithFormat:@"%.3f", [[_fixedXArray objectAtIndex:row] floatValue]];
//            break;
//        }
//        case 1:{
//            return [NSString stringWithFormat:@"%.3f", [[_fixedYArray objectAtIndex:row] floatValue]];
//            break;
//        }
//        case 2:{
//            return [NSString stringWithFormat:@"%.3f", [[_fixedTArray objectAtIndex:row] floatValue]];
//            break;
//        }
//        default:
//            break;
//    }
//    return @"0";//[NSString stringWithFormat:@"%.3f", [[keys objectAtIndex:row] doubleValue]];
}

// tell the picker the width of each row for a given component
- (CGFloat)pickerView:(UIPickerView *)pickerView widthForComponent:(NSInteger)component {
    int sectionWidth = 75;
    
    return sectionWidth;
}


@end
