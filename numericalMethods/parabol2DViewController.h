//
//  parabol2DViewController.h
//  numericalMethods
//
//  Created by Ilya on 05.12.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import <UIKit/UIKit.h>
#import "UIPickerViewSchemeElliptic.h"
#import "UIPickerViewVariableElliptic.h"
#import "MBProgressHUD.h"
#import "DDMathParser.h"
#import "DDMathStringTokenizer.h"
#import "NMFTextField.h"

@interface parabol2DViewController : UIViewController<UITextFieldDelegate, UIPickerViewDelegate,MBProgressHUDDelegate>{
    IBOutlet NMFTextField *bxField;
    IBOutlet NMFTextField *byField;
    IBOutlet NMFTextField *cField;
    IBOutlet NMFTextField *lxField;
    IBOutlet NMFTextField *lyField;
    IBOutlet NMFTextField *alpha1Field;
    IBOutlet NMFTextField *betta1Field;
    IBOutlet NMFTextField *alpha2Field;
    IBOutlet NMFTextField *betta2Field;
    IBOutlet NMFTextField *alpha3Field;
    IBOutlet NMFTextField *betta3Field;
    IBOutlet NMFTextField *alpha4Field;
    IBOutlet NMFTextField *betta4Field;
    IBOutlet NMFTextField *lastXField;
    IBOutlet NMFTextField *lastYField;
    IBOutlet NMFTextField *NyField;
    IBOutlet NMFTextField *NxField;
    IBOutlet NMFTextField *epsilonField;
    
    IBOutlet UITextField *UField;
    IBOutlet UITextField *FField;
    IBOutlet UITextField *Phi1Field;
    IBOutlet UITextField *Phi2Field;
    IBOutlet UITextField *Phi3Field;
    IBOutlet UITextField *Phi4Field;
    UiPickerViewSchemeElliptic *schemePicker;
    UIPickerViewVariableElliptic *variablePicker;
    UIView *bg;
    UIView *bg1;
    UIView *bg2;
    UIScrollView *scrollBg;
    DDMathEvaluator *evaluator;
    DDMathStringTokenizer *tokenizer;
    DDExpression *expressionRealFunc;
    DDExpression *expressionPhi1;
    DDExpression *expressionPhi2;
    DDExpression *expressionPhi3;
    DDExpression *expressionPhi4;
    DDExpression *expressionF;
    
}

@property (nonatomic,retain) MBProgressHUD *HUD;
@property (nonatomic, strong) UIImageView *system;

@end
