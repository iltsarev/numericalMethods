//
//  hyperbolicViewController.h
//  numericalMethods
//
//  Created by Ilya on 03.11.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import <UIKit/UIKit.h>
#import "UIPickerViewSchemeHyperbolic.h"
#import "UIPIckerViewOrder.h"
#import "UIPickerViewApproximation.h"
#import "MBProgressHUD.h"

@interface hyperbolicViewController : UIViewController<UITextFieldDelegate, UIPickerViewDelegate,MBProgressHUDDelegate>{
    IBOutlet UITextField *aField;
    IBOutlet UITextField *bField;
    IBOutlet UITextField *cField;
    IBOutlet UITextField *eField;
    IBOutlet UITextField *lField;
    IBOutlet UITextField *alphaField;
    IBOutlet UITextField *bettaField;
    IBOutlet UITextField *gammaField;
    IBOutlet UITextField *deltaField;
    IBOutlet UITextField *TField;
    IBOutlet UITextField *KField;
    IBOutlet UITextField *NField;
    UIPIckerViewOrder *orderPicker;
    UIPickerViewSchemeHyperbolic *schemePicker;
    UIPickerViewApproximation *approxPicker;
    UIView *bg;
    UIView *bg2;
    UIScrollView *scrollBg;
}

@property (nonatomic,retain) MBProgressHUD *HUD;
@property (nonatomic, strong) UIImageView *system;

@end
