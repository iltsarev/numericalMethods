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
	
	graph = [[CPTXYGraph alloc] initWithFrame: self.view.bounds];
	
	CPTGraphHostingView *hostingView = (CPTGraphHostingView *)self.view;
	hostingView.hostedGraph = graph;
	graph.paddingLeft = 20.0;
	graph.paddingTop = 20.0;
	graph.paddingRight = 20.0;
	graph.paddingBottom = 20.0;
	
	CPTXYPlotSpace *plotSpace = (CPTXYPlotSpace *)graph.defaultPlotSpace;
	plotSpace.xRange = [CPTPlotRange plotRangeWithLocation:CPTDecimalFromFloat(-6 )
												   length:CPTDecimalFromFloat(12)];
	plotSpace.yRange = [CPTPlotRange plotRangeWithLocation:CPTDecimalFromFloat(-5)
												   length:CPTDecimalFromFloat(30)];
	
	CPTXYAxisSet *axisSet = (CPTXYAxisSet *)graph.axisSet;
	
	CPTLineStyle *lineStyle = [CPTLineStyle lineStyle];
	//lineStyle.lineColor = [CPTColor blackColor];
	//lineStyle.lineWidth = 2.0f;
	
	axisSet.xAxis.majorIntervalLength = [[NSDecimalNumber decimalNumberWithString:@"5"] decimalValue];
	axisSet.xAxis.minorTicksPerInterval = 4;
	axisSet.xAxis.majorTickLineStyle = lineStyle;
	axisSet.xAxis.minorTickLineStyle = lineStyle;
	axisSet.xAxis.axisLineStyle = lineStyle;
	axisSet.xAxis.minorTickLength = 5.0f;
	axisSet.xAxis.majorTickLength = 7.0f;
	axisSet.xAxis.labelOffset = 3.0f;
	
	axisSet.yAxis.majorIntervalLength = [[NSDecimalNumber decimalNumberWithString:@"5"] decimalValue];
	axisSet.yAxis.minorTicksPerInterval = 4;
	axisSet.yAxis.majorTickLineStyle = lineStyle;
	axisSet.yAxis.minorTickLineStyle = lineStyle;
	axisSet.yAxis.axisLineStyle = lineStyle;
	axisSet.yAxis.minorTickLength = 5.0f;
	axisSet.yAxis.majorTickLength = 7.0f;
	axisSet.yAxis.labelOffset = 3.0f;
	
	CPTScatterPlot *xSquaredPlot = [[CPTScatterPlot alloc]
									initWithFrame:self.view.bounds];
	xSquaredPlot.identifier = @"X Squared Plot";
	//xSquaredPlot.dataLineStyle.lineWidth = 1.0f;
	//xSquaredPlot.dataLineStyle.lineColor = [CPColor redColor];
	xSquaredPlot.dataSource = self;
	[graph addPlot:xSquaredPlot];
	
	CPTPlotSymbol *greenCirclePlotSymbol = [CPTPlotSymbol ellipsePlotSymbol];
	greenCirclePlotSymbol.fill = [CPTFill fillWithColor:[CPTColor greenColor]];
	greenCirclePlotSymbol.size = CGSizeMake(2.0, 2.0);
	xSquaredPlot.plotSymbol = greenCirclePlotSymbol;
	
	CPTScatterPlot *xInversePlot = [[CPTScatterPlot alloc]
									initWithFrame:self.view.bounds] ;
	xInversePlot.identifier = @"X Inverse Plot";
	//xInversePlot.dataLineStyle.lineWidth = 1.0f;
	//xInversePlot.dataLineStyle.lineColor = [CPColor blueColor];
	xInversePlot.dataSource = self;
	[graph addPlot:xInversePlot];
}

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
