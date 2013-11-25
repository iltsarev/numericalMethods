//
//  UIPickerViewVariableElliptic.h
//  numericalMethods
//
//  Created by Ilya on 25.11.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import <UIKit/UIKit.h>

@interface UIPickerViewVariableElliptic : UIPickerView<UIPickerViewDelegate>{
    IBOutlet UIPickerView *schemePicker;
    IBOutlet NSArray *schemeArray;
}


@end
