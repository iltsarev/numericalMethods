//
//  CorePlotExampleViewController.m
//  CorePlotExample
//
//  Created by Sylvain Marcotte on 11-02-16.
//  14 Oranges Software - http://www.14oranges.com
//
//  Original code obtained from http://www.switchonthecode.com/tutorials/using-core-plot-in-an-iphone-application
//

#import "CorePlotExampleViewController.h"

@implementation CorePlotExampleViewController



/*
// The designated initializer. Override to perform setup that is required before the view is loaded.
- (id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil {
    self = [super initWithNibName:nibNameOrNil bundle:nibBundleOrNil];
    if (self) {
        // Custom initialization
    }
    return self;
}
*/

/*
// Implement loadView to create a view hierarchy programmatically, without using a nib.
- (void)loadView {
}
*/



// Implement viewDidLoad to do additional setup after loading the view, typically from a nib.
- (void)viewDidLoad 
{
	[super viewDidLoad];
	
	defaultLayerHostingView.collapsesLayers = NO;
    [defaultLayerHostingView setAutoresizingMask:UIViewAutoresizingFlexibleWidth | UIViewAutoresizingFlexibleHeight];
#else
    [defaultLayerHostingView setAutoresizingMask:NSViewWidthSizable | NSViewHeightSizable];
#endif
    [defaultLayerHostingView setAutoresizesSubviews:YES];
    
    [hostingView addSubview:defaultLayerHostingView];
    [self generateData];
    [self renderInLayer:defaultLayerHostingView withTheme:theme animated:animated];}


-(NSUInteger)numberOfRecordsForPlot:(CPTPlot *)plot 
{
	return 51;
}

-(NSNumber *)numberForPlot:(CPTPlot *)plot field:(NSUInteger)fieldEnum  
			   recordIndex:(NSUInteger)index
{
    
	NSLog(@"Index = %d", index);
	double val = (index/5.0)-5;
	if(fieldEnum == CPTScatterPlotFieldX)
	{ return [NSNumber numberWithDouble:val]; }
	else
	{ 
		if(plot.identifier == @"X Squared Plot")
		{ return [NSNumber numberWithDouble:val*val]; }
		else
		{ return [NSNumber numberWithDouble:1/val]; }
	}
}




/*
// Override to allow orientations other than the default portrait orientation.
- (BOOL)shouldAutorotateToInterfaceOrientation:(UIInterfaceOrientation)interfaceOrientation {
    // Return YES for supported orientations
    return (interfaceOrientation == UIInterfaceOrientationPortrait);
}
*/

- (void)didReceiveMemoryWarning {
	// Releases the view if it doesn't have a superview.
    [super didReceiveMemoryWarning];
	
	// Release any cached data, images, etc that aren't in use.
}

- (void)viewDidUnload {
	// Release any retained subviews of the main view.
	// e.g. self.myOutlet = nil;
}




@end
