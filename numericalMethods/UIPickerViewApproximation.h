//
//  UIPickerViewApproximation.h
//  numericalMethods
//
//  Created by Ilya on 03.11.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import <UIKit/UIKit.h>

@interface UIPickerViewApproximation : UIPickerView<UIPickerViewDelegate>{
    IBOutlet UIPickerView *aproximationPicker;
    IBOutlet NSArray *aproximationArray;
    
}

@end
