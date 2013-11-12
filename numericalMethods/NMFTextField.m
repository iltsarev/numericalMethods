//
//  NMFTextField.m
//  numericalMethods
//
//  Created by Ilya on 12.11.13.
//  Copyright (c) 2013 Ilya. All rights reserved.
//

#import "NMFTextField.h"

@implementation NMFTextField

- (id)initWithFrame:(CGRect)frame
{
    self = [super initWithFrame:frame];
    if (self) {
        UIToolbar *toolBar = [[UIToolbar alloc] initWithFrame:CGRectMake(0.0f,
                                                                         0.0f,
                                                                         44.0f,
                                                                         30.0f)];
        if (UI_USER_INTERFACE_IDIOM() == UIUserInterfaceIdiomPad)
        {
            toolBar.tintColor = [UIColor colorWithRed:0.6f
                                                green:0.6f
                                                 blue:0.64f
                                                alpha:1.0f];
        }
        else
        {
            toolBar.tintColor = [UIColor colorWithRed:0.56f
                                                green:0.59f
                                                 blue:0.63f
                                                alpha:1.0f];
        }
//        toolBar.backgroundColor = [UIColor grayColor];
//        toolBar.translucent = NO;
        toolBar.backgroundColor = [UIColor colorWithRed:0.15f
                                                 green:0.15f
                                                  blue:0.15f
                                                 alpha:0.1f];
        UIBarButtonItem *minusButton = [[UIBarButtonItem alloc] initWithTitle:@"-"
                                                                        style:UIBarButtonItemStyleBordered
                                                                       target:self
                                                                       action:@selector(barButtonAddText:)];
        [minusButton setTitleTextAttributes:[NSDictionary dictionaryWithObjectsAndKeys:[UIFont boldSystemFontOfSize:15], NSFontAttributeName,nil] forState:UIControlStateNormal];
        toolBar.items =   @[minusButton,
                             [[UIBarButtonItem alloc] initWithBarButtonSystemItem:UIBarButtonSystemItemFlexibleSpace
                                                                           target:nil
                                                                           action:nil],
                             [[UIBarButtonItem alloc] initWithBarButtonSystemItem:UIBarButtonSystemItemDone
                                                                           target:nil
                                                                           action:@selector(close:)],
                             // some more items could be added
                             ];
        self.inputAccessoryView = toolBar;
        [self setKeyboardType:UIKeyboardTypeDecimalPad];
        self.delegate = self;
    }
    return self;
}

-(IBAction)barButtonAddText:(UIBarButtonItem*)sender
{
    if (self.isFirstResponder)
    {
        [self insertText:sender.title];
    }
}

-(void)close:(id)sender{
    [self resignFirstResponder];
}

-(BOOL) textFieldShouldReturn:(UITextField *)textField{
    [textField resignFirstResponder];
    return YES;
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
