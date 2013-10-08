//
//  UIPIckerViewOrder.m
//  numericalMethods
//
//  Created by Ilya on 08.10.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import "UIPIckerViewOrder.h"

@implementation UIPIckerViewOrder

- (id)initWithFrame:(CGRect)frame
{
    self = [super initWithFrame:frame];
    if (self) {
        // Initialization code
        orderArray = [[NSArray alloc] init];
        orderArray = @[@"2-ух точечная 1-го порядка", @"3-ех точечная 2-го порядка", @"2-ух точечная 2-го порядка"];
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
    return [orderArray count];
}

// tell the picker how many components it will have
- (NSInteger)numberOfComponentsInPickerView:(UIPickerView *)pickerView {
    return 1;
}

// tell the picker the title for a given component
- (NSString *)pickerView:(UIPickerView *)pickerView titleForRow:(NSInteger)row forComponent:(NSInteger)component {
    return [orderArray objectAtIndex:row];
}

// tell the picker the width of each row for a given component
- (CGFloat)pickerView:(UIPickerView *)pickerView widthForComponent:(NSInteger)component {
    int sectionWidth = 320;
    
    return sectionWidth;
}

/*
// Only override drawRect: if you perform custom drawing.
// An empty implementation adversely affects performance during animation.
- (void)drawRect:(CGRect)rect
{
    // Drawing code
}
*/

@end
