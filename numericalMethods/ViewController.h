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

@interface ViewController : UIViewController<UITextFieldDelegate, UIPickerViewDelegate,MBProgressHUDDelegate>{
    IBOutlet UITextField *aField;
    IBOutlet UITextField *bField;
    IBOutlet UITextField *cField;
    IBOutlet UITextField *lField;
    IBOutlet UITextField *alphaField;
    IBOutlet UITextField *bettaField;
    IBOutlet UITextField *gammaField;
    IBOutlet UITextField *deltaField;
    IBOutlet UITextField *TField;
    IBOutlet UITextField *KField;
    IBOutlet UITextField *NField;
    UIPIckerViewOrder *orderPicker;
    UIPickerViewScheme *schemePicker;
    UIView *bg;
    UIView *bg2;    
}

@property (nonatomic,retain)MBProgressHUD *HUD;

@end
