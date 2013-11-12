//
//  ViewController.h
//  numericalMethods
//
//  Created by Ilya on 02.10.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import <UIKit/UIKit.h>
#import "UIPickerViewScheme.h"
#import "UIPIckerViewOrder.h"
#import "MBProgressHUD.h"
#import "DDMathParser.h"
#import "DDMathStringTokenizer.h"
#import "NMFTextField.h"


@interface ViewController : UIViewController<UITextFieldDelegate, UIPickerViewDelegate,MBProgressHUDDelegate>{
    IBOutlet NMFTextField *aField;
    IBOutlet NMFTextField *bField;
    IBOutlet NMFTextField *cField;
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
    IBOutlet UITextField *PhiField;
    IBOutlet UITextField *Phi0Field;
    IBOutlet UITextField *PhilField;
    UIPIckerViewOrder *orderPicker;
    UIPickerViewScheme *schemePicker;
    UIView *bg;
    UIView *bg1;
    UIView *bg2;
    UIScrollView *scrollBg;
    DDMathEvaluator *evaluator;
    DDMathStringTokenizer *tokenizer;
    DDExpression *expressionRealFunc;
    DDExpression *expressionPhi;
    DDExpression *expressionPhi_0;
    DDExpression *expressionPhi_l;
    DDExpression *expressionF;

}

@property (nonatomic,retain) MBProgressHUD *HUD;
@property (nonatomic, strong) UIImageView *system;

@end
