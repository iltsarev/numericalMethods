//
//  UIPickerViewApproximation.m
//  numericalMethods
//
//  Created by Ilya on 03.11.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import "UIPickerViewApproximation.h"

@implementation UIPickerViewApproximation

- (id)initWithFrame:(CGRect)frame
{
    self = [super initWithFrame:frame];
    if (self) {
        // Initialization code
        aproximationArray = [[NSArray alloc] init];
        aproximationArray = @[@"Первый порядок", @"Второй порядок"];
        
        self.showsSelectionIndicator = YES;
        self.delegate = self;
        
    }
    return self;
}

- (void)pickerView:(UIPickerView *)pickerView didSelectRow: (NSInteger)row inComponent:(NSInteger)component {
    // Handle the selection
}

// tell the picker how many rows are available for a given component
- (NSInteger)pickerView:(UIPickerView *)pickerView numberOfRowsInComponent:(NSInteger)component {
    return [aproximationArray count];
}

// tell the picker how many components it will have
- (NSInteger)numberOfComponentsInPickerView:(UIPickerView *)pickerView {
    return 1;
}

// tell the picker the title for a given component
- (NSString *)pickerView:(UIPickerView *)pickerView titleForRow:(NSInteger)row forComponent:(NSInteger)component {
    return [aproximationArray objectAtIndex:row];
}

// tell the picker the width of each row for a given component
- (CGFloat)pickerView:(UIPickerView *)pickerView widthForComponent:(NSInteger)component {
    int sectionWidth = 300;
    
    return sectionWidth;
}


@end
