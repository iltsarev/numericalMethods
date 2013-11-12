//
//  NMFTextField.h
//  numericalMethods
//
//  Created by Ilya on 12.11.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import <UIKit/UIKit.h>

@interface NMFTextField : UITextField<UITextFieldDelegate>
-(BOOL) textFieldShouldReturn:(UITextField *)textField;
@end
