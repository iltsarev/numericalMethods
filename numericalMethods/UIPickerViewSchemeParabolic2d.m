//
//  UIPickerViewSchemeParabolic2d.m
//  numericalMethods
//
//  Created by Ilya on 06.12.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import "UIPickerViewSchemeParabolic2d.h"

@implementation UIPickerViewSchemeParabolic2d

- (id)initWithFrame:(CGRect)frame
{
    self = [super initWithFrame:frame];
    if (self) {
        // Initialization code
        schemeArray = [[NSArray alloc] init];
        schemeArray = @[@"Дробных шагов", @"Переменных направлений"];
        
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
    return [schemeArray count];
}

// tell the picker how many components it will have
- (NSInteger)numberOfComponentsInPickerView:(UIPickerView *)pickerView {
    return 1;
}

// tell the picker the title for a given component
- (NSString *)pickerView:(UIPickerView *)pickerView titleForRow:(NSInteger)row forComponent:(NSInteger)component {
    return [schemeArray objectAtIndex:row];
}

// tell the picker the width of each row for a given component
- (CGFloat)pickerView:(UIPickerView *)pickerView widthForComponent:(NSInteger)component {
    int sectionWidth = 300;
    
    return sectionWidth;
}


@end
