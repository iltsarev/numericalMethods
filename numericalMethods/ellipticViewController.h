//
//  ellipticViewController.h
//  numericalMethods
//
//  Created by Ilya on 23.11.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import <UIKit/UIKit.h>
#import "UIPickerViewSchemeHyperbolic.h"
#import "UIPIckerViewOrder.h"
#import "UIPickerViewApproximation.h"
#import "MBProgressHUD.h"
#import "DDMathParser.h"
#import "DDMathStringTokenizer.h"
#import "NMFTextField.h"


@interface ellipticViewController : UIViewController<UITextFieldDelegate, UIPickerViewDelegate,MBProgressHUDDelegate>{
    IBOutlet NMFTextField *aField;
    IBOutlet NMFTextField *bField;
    IBOutlet NMFTextField *cField;
    IBOutlet NMFTextField *eField;
    IBOutlet NMFTextField *lField;
    IBOutlet NMFTextField *alphaField;
    IBOutlet NMFTextField *bettaField;
    IBOutlet NMFTextField *gammaField;
    IBOutlet NMFTextField *deltaField;
    IBOutlet NMFTextField *TField;
    IBOutlet NMFTextField *KField;
    IBOutlet NMFTextField *NField;
    
    IBOutlet UITextField *UField;
    IBOutlet UITextField *FField;
    IBOutlet UITextField *Phi0Field;
    IBOutlet UITextField *PhilField;
    IBOutlet UITextField *Psi1Field;
    IBOutlet UITextField *Psi2Field;
    IBOutlet UITextField *Psi1dField;
    IBOutlet UITextField *Psi1ddField;
    UIPIckerViewOrder *orderPicker;
    UIPickerViewSchemeHyperbolic *schemePicker;
    UIPickerViewApproximation *approxPicker;
    UIView *bg;
    UIView *bg1;
    UIView *bg2;
    UIScrollView *scrollBg;
    DDMathEvaluator *evaluator;
    DDMathStringTokenizer *tokenizer;
    DDExpression *expressionRealFunc;
    DDExpression *expressionPhi_0;
    DDExpression *expressionPhi_l;
    DDExpression *expressionPsi_1;
    DDExpression *expressionPsi_2;
    DDExpression *expressionPsi_1d;
    DDExpression *expressionPsi_1dd;
    DDExpression *expressionF;
    
}

@property (nonatomic,retain) MBProgressHUD *HUD;
@property (nonatomic, strong) UIImageView *system;

@end
