#import "CorePlotViewController.h"
#import "CPViewControllerErrors.h"
@implementation CorePlotViewController

@synthesize dataForPlot;
@synthesize dictForPlot;
@synthesize dictForPlotErr;
@synthesize dictForPlotAnalytic;
@synthesize keys;
@synthesize size;
@synthesize time1;
@synthesize time2;
@synthesize time3;
@synthesize time4;

-(BOOL)shouldAutorotateToInterfaceOrientation:(UIInterfaceOrientation)toInterfaceOrientation
{
    return YES;
}

#pragma mark -
#pragma mark Initialization and teardown

-(void)scatterPlot:(CPTScatterPlot *)plot plotSymbolWasSelectedAtRecordIndex:(NSUInteger)index
{
    //CPTXYGraph *graph = [graphs objectAtIndex:0];
    
    if ( symbolTextAnnotation ) {
        [graph.plotAreaFrame.plotArea removeAnnotation:symbolTextAnnotation];
        symbolTextAnnotation = nil;
    }
    
    // Setup a style for the annotation
    CPTMutableTextStyle *hitAnnotationTextStyle = [CPTMutableTextStyle textStyle];
    hitAnnotationTextStyle.color    = [CPTColor whiteColor];
    hitAnnotationTextStyle.fontSize = 16.0f;
    hitAnnotationTextStyle.fontName = @"Helvetica-Bold";
    
    // Determine point of symbol in plot coordinates
    NSNumber *x          = [[dataForPlot objectAtIndex:index] valueForKey:@"x"];
    NSNumber *y          = [[dataForPlot objectAtIndex:index] valueForKey:@"y"];
    NSArray *anchorPoint = [NSArray arrayWithObjects:x, y, nil];
    
    // Add annotation
    // First make a string for the y value
    NSNumberFormatter *formatter = [[NSNumberFormatter alloc] init] ;
    [formatter setMaximumFractionDigits:2];
    NSString *yString = [formatter stringFromNumber:y];
    
    // Now add the annotation to the plot area
    CPTTextLayer *textLayer = [[CPTTextLayer alloc] initWithText:yString style:hitAnnotationTextStyle];
    symbolTextAnnotation              = [[CPTPlotSpaceAnnotation alloc] initWithPlotSpace:graph.defaultPlotSpace anchorPlotPoint:anchorPoint];
    symbolTextAnnotation.contentLayer = textLayer;
    symbolTextAnnotation.displacement = CGPointMake(0.0f, 20.0f);
    [graph.plotAreaFrame.plotArea addAnnotation:symbolTextAnnotation];
}

-(CPTPlotRange *)plotSpace:(CPTPlotSpace *)space willChangePlotRangeTo:(CPTPlotRange *)newRange forCoordinate:(CPTCoordinate)coordinate
{
    CPTXYAxisSet *axisSet = (CPTXYAxisSet *)space.graph.axisSet;
    
    CPTMutablePlotRange *changedRange = [newRange mutableCopy];
    
    switch ( coordinate ) {
        case CPTCoordinateX:
            [changedRange expandRangeByFactor:CPTDecimalFromDouble(1.025)];
            changedRange.location          = newRange.location;
            axisSet.xAxis.visibleAxisRange = changedRange;
            break;
            
        case CPTCoordinateY:
            [changedRange expandRangeByFactor:CPTDecimalFromDouble(1.05)];
            axisSet.yAxis.visibleAxisRange = changedRange;
            break;
            
        default:
            break;
    }
    
    return newRange;
}

-(void)viewDidLoad
{
    [super viewDidLoad];
    
    
    self.view.backgroundColor = [UIColor whiteColor];
    
    // Create graph from theme
    time1 = time2 = time3 = time4 = 0;
    
    CPTTheme *theme = [CPTTheme themeNamed:kCPTDarkGradientTheme];
    [graph applyTheme:theme];
    CPTGraphHostingView *hostingView = (CPTGraphHostingView *)self.view;
    CGRect bounds = hostingView.bounds;
    bounds.size.height +=20;
    graph = [[CPTXYGraph alloc] initWithFrame:bounds];
    hostingView.hostedGraph = graph;
    graph.title = self.title;
    CPTMutableTextStyle *textStyle = [CPTMutableTextStyle textStyle];
    textStyle.color                = [CPTColor grayColor];
    textStyle.fontName             = @"Helvetica-Bold";
    textStyle.fontSize             = round( bounds.size.height / CPTFloat(20.0) );
    graph.titleTextStyle           = textStyle;
    graph.titleDisplacement        = CPTPointMake( 0.0, textStyle.fontSize * CPTFloat(1.5) );
    graph.titlePlotAreaFrameAnchor = CPTRectAnchorTop;
    
    CGFloat boundsPadding = round( bounds.size.width / CPTFloat(20.0) ); // Ensure that padding falls on an integral pixel
    
    graph.paddingLeft = boundsPadding;
    
    if ( graph.titleDisplacement.y > 0.0 ) {
        graph.paddingTop = graph.titleTextStyle.fontSize * 2.0;
    }
    else {
        graph.paddingTop = boundsPadding;
    }
    
    graph.paddingRight  = boundsPadding;
    graph.paddingBottom = boundsPadding;
    
    // Add some initial data
   self.dataForPlot = self.dictForPlot[[NSNumber numberWithDouble:0.0]];
    
    graph.plotAreaFrame.paddingLeft   += 55.0;
//    graph.plotAreaFrame.paddingTop    += 40.0;
//    graph.plotAreaFrame.paddingRight  += 55.0;
    graph.plotAreaFrame.paddingBottom += 160.0;
   //graph.plotAreaFrame.masksToBorder  = NO;
    
    // Setup scatter plot space
    CPTXYPlotSpace *plotSpace = (CPTXYPlotSpace *)graph.defaultPlotSpace;
    plotSpace.allowsUserInteraction = YES;
    //plotSpace.delegate              = self;
    // Grid line styles
    CPTMutableLineStyle *majorGridLineStyle = [CPTMutableLineStyle lineStyle];
    majorGridLineStyle.lineWidth = 0.75;
    majorGridLineStyle.lineColor = [[CPTColor colorWithGenericGray:0.2] colorWithAlphaComponent:0.75];
    
    CPTMutableLineStyle *minorGridLineStyle = [CPTMutableLineStyle lineStyle];
    minorGridLineStyle.lineWidth = 0.25;
    minorGridLineStyle.lineColor = [[CPTColor whiteColor] colorWithAlphaComponent:0.1];
    
    CPTMutableLineStyle *redLineStyle = [CPTMutableLineStyle lineStyle];
    redLineStyle.lineWidth = 10.0;
    redLineStyle.lineColor = [[CPTColor redColor] colorWithAlphaComponent:0.5];
    
    CPTLineCap *lineCap = [CPTLineCap sweptArrowPlotLineCap];
    lineCap.size = CGSizeMake(15.0, 15.0);
    
    // Axes
    // Label x axis with a fixed interval policy
    CPTXYAxisSet *axisSet = (CPTXYAxisSet *)graph.axisSet;
    CPTXYAxis *x          = axisSet.xAxis;
    x.majorIntervalLength   = CPTDecimalFromDouble(1);
    x.minorTicksPerInterval = 4;
    x.majorGridLineStyle    = majorGridLineStyle;
    x.minorGridLineStyle    = minorGridLineStyle;
    x.axisConstraints       = [CPTConstraints constraintWithRelativeOffset:0.5];
    
    lineCap.lineStyle = x.axisLineStyle;
    lineCap.fill      = [CPTFill fillWithColor:lineCap.lineStyle.lineColor];
    x.axisLineCapMax  = lineCap;
    
    x.title       = @"x";
    x.titleOffset = 30.0;
    
    // Label y with an automatic label policy.
    CPTXYAxis *y = axisSet.yAxis;
    y.labelingPolicy              = CPTAxisLabelingPolicyAutomatic;
    y.minorTicksPerInterval       = 4;
    y.preferredNumberOfMajorTicks = 8;
    y.majorGridLineStyle          = majorGridLineStyle;
    y.minorGridLineStyle          = minorGridLineStyle;
    y.axisConstraints             = [CPTConstraints constraintWithLowerOffset:0.0];
    y.labelOffset                 = 10.0;
    
    lineCap.lineStyle = y.axisLineStyle;
    lineCap.fill      = [CPTFill fillWithColor:lineCap.lineStyle.lineColor];
    y.axisLineCapMax  = lineCap;
    y.axisLineCapMin  = lineCap;
    
    y.title       = @"U(x,t)";
    y.titleOffset = 32.0;
    
    // Set axes
    graph.axisSet.axes = [NSArray arrayWithObjects:x, y, nil];
    
    // Create a plot that uses the data source method
    CPTScatterPlot *dataSourceLinePlot = [[CPTScatterPlot alloc] init];
    dataSourceLinePlot.identifier = @"1st time";
    
    // Make the data source line use curved interpolation
    
    dataSourceLinePlot.interpolation = CPTScatterPlotInterpolationCurved;
    
    CPTMutableLineStyle *lineStyle = [dataSourceLinePlot.dataLineStyle mutableCopy] ;
    lineStyle.lineWidth              = 2.0;
    lineStyle.lineColor              = [CPTColor redColor];
    dataSourceLinePlot.dataLineStyle = lineStyle;
    
    dataSourceLinePlot.dataSource = self;
    [graph addPlot:dataSourceLinePlot];
    
    // First derivative
    CPTScatterPlot *firstPlot = [[CPTScatterPlot alloc] init] ;
    firstPlot.identifier    = @"2nd time";
    lineStyle.lineWidth     = 2.0;
    lineStyle.lineColor     = [CPTColor blueColor];
    firstPlot.dataLineStyle = lineStyle;
    firstPlot.dataSource    = self;
    
    [graph addPlot:firstPlot];
    
    
    // First derivative
    CPTScatterPlot *firstPlot2 = [[CPTScatterPlot alloc] init] ;
    firstPlot2.identifier    = @"3rd time";
    lineStyle.lineWidth     = 2.0;
    lineStyle.lineColor     = [CPTColor greenColor];
    firstPlot2.dataLineStyle = lineStyle;
    firstPlot2.dataSource    = self;
    
    [graph addPlot:firstPlot2];
    
    
    // Second derivative
    CPTScatterPlot *secondPlot = [[CPTScatterPlot alloc] init] ;
    secondPlot.identifier    = @"real";
    lineStyle.lineColor      = [CPTColor whiteColor];
    secondPlot.dataLineStyle = lineStyle;
    secondPlot.dataSource    = self;
    
    [graph addPlot:secondPlot];
    
    // Auto scale the plot space to fit the plot data
    [plotSpace scaleToFitPlots:[graph allPlots]];
    CPTMutablePlotRange *xRange = [plotSpace.xRange mutableCopy] ;
    CPTMutablePlotRange *yRange = [plotSpace.yRange mutableCopy] ;
    
    // Expand the ranges to put some space around the plot
    [xRange expandRangeByFactor:CPTDecimalFromDouble(1.2)];
    [yRange expandRangeByFactor:CPTDecimalFromDouble(1.2)];
    plotSpace.xRange = xRange;
    plotSpace.yRange = yRange;
    
    [xRange expandRangeByFactor:CPTDecimalFromDouble(1.025)];
    xRange.location = plotSpace.xRange.location;
    [yRange expandRangeByFactor:CPTDecimalFromDouble(1.05)];
    x.visibleAxisRange = xRange;
    y.visibleAxisRange = yRange;
    
    [xRange expandRangeByFactor:CPTDecimalFromDouble(3.0)];
    [yRange expandRangeByFactor:CPTDecimalFromDouble(3.0)];
    plotSpace.globalXRange = xRange;
    plotSpace.globalYRange = yRange;
    
    // Add plot symbols
    CPTMutableLineStyle *symbolLineStyle = [CPTMutableLineStyle lineStyle];
    symbolLineStyle.lineColor = [[CPTColor blackColor] colorWithAlphaComponent:0.5];
    CPTPlotSymbol *plotSymbol = [CPTPlotSymbol ellipsePlotSymbol];
    plotSymbol.fill               = [CPTFill fillWithColor:[[CPTColor blueColor] colorWithAlphaComponent:0.5]];
    plotSymbol.lineStyle          = symbolLineStyle;
    plotSymbol.size               = CGSizeMake(4.0, 4.0);
    
//    dataSourceLinePlot.plotSymbol = plotSymbol;
    firstPlot.plotSymbol = plotSymbol;
//    firstPlot2.plotSymbol = plotSymbol;
//    secondPlot.plotSymbol = plotSymbol;

    // Set plot delegate, to know when symbols have been touched
    // We will display an annotation when a symbol is touched
//    dataSourceLinePlot.delegate                        = self;
//    dataSourceLinePlot.plotSymbolMarginForHitDetection = 5.0f;
//    
//    firstPlot.delegate                        = self;
//    firstPlot.plotSymbolMarginForHitDetection = 5.0f;
//    
//    secondPlot.delegate                        = self;
//    secondPlot.plotSymbolMarginForHitDetection = 5.0f;
    
    // Add legend
    graph.legend                 = [CPTLegend legendWithGraph:graph];
    graph.legend.numberOfRows    = 1;
    graph.legend.textStyle       = x.titleTextStyle;
    graph.legend.fill            = [CPTFill fillWithColor:[CPTColor darkGrayColor]];
    graph.legend.borderLineStyle = x.axisLineStyle;
    graph.legend.cornerRadius    = 5.0;
    graph.legend.swatchSize      = CGSizeMake(25.0, 25.0);
    graph.legendAnchor           = CPTRectAnchorBottom;
    graph.legendDisplacement     = CGPointMake(0.0, 12.0);
    
#ifdef PERFORMANCE_TEST
    [NSTimer scheduledTimerWithTimeInterval:2.0 target:self selector:@selector(changePlotRange) userInfo:nil repeats:YES];
#endif
    
    UIButton *back;
    
    if ([[UIScreen mainScreen] bounds].size.height < 500) {
        back = [[UIButton alloc] initWithFrame:CGRectMake(5, self.view.frame.size.height - 50, 70, 30)];
        err = [[UIButton alloc] initWithFrame:CGRectMake(self.view.frame.size.width - 80, self.view.frame.size.height - 50, 70, 30)];
    }
    else {
        back = [[UIButton alloc] initWithFrame:CGRectMake(5, self.view.frame.size.height + 40, 70, 30)];
        err = [[UIButton alloc] initWithFrame:CGRectMake(self.view.frame.size.width - 80, self.view.frame.size.height + 40, 70, 30)];
    }
    
    [back setTitle:@"Назад" forState:UIControlStateNormal];
    back.titleLabel.textColor = [UIColor blackColor];
    //back.backgroundColor = [UIColor blackColor];
    [back  addTarget:self action:@selector(Back) forControlEvents:UIControlEventTouchUpInside];
    back.layer.transform = CATransform3DMakeRotation (M_PI, 1, 0, 0.f);
    [self.view addSubview:back];
    
    [err setTitle:@"Ошибки" forState:UIControlStateNormal];
    err.titleLabel.textColor = [UIColor blackColor];
    //back.backgroundColor = [UIColor blackColor];
    [err  addTarget:self action:@selector(Errors) forControlEvents:UIControlEventTouchUpInside];
    err.layer.transform = CATransform3DMakeRotation (M_PI, 1, 0, 0.f);
    [self.view addSubview:err];
    
    NSMutableArray *keys_m = [NSMutableArray arrayWithArray:[dictForPlot allKeys]];
    NSSortDescriptor *sortDescriptor;
    sortDescriptor = [[NSSortDescriptor alloc] initWithKey:@"self" ascending:YES];
    NSArray *sortDescriptors = [NSArray arrayWithObject:sortDescriptor];
    [keys_m sortUsingDescriptors:sortDescriptors];
    
    keys = [[NSArray alloc] initWithArray:keys_m];
    size = (int)keys.count;
    
    UIView *line = [[UIView alloc] initWithFrame:CGRectMake(0, 47, 320, 120)];
    line.backgroundColor = [UIColor grayColor];
    [self.view addSubview:line];
    
    UIPickerView *picker = [[UIPickerView alloc] initWithFrame:CGRectMake(0, 25, 320, 100)];
    picker.showsSelectionIndicator = YES;
    picker.delegate = self;
    picker.layer.transform = CATransform3DMakeRotation (M_PI, 1, 0, 0.f);
    [self.view addSubview:picker];

}

-(void)Back{
    [self dismissViewControllerAnimated:YES completion:^{}];
}


-(void)Errors{
    CorePlotViewControllerErrors *viewControllerToPresent = [[CorePlotViewControllerErrors alloc] initWithNibName:@"CPViewControllerErrors" bundle:nil];
    //viewControllerToPresent.dataForPlot = contentArray;
    viewControllerToPresent.dataForPlot = self.dictForPlotErr;
    NSLog(@"%@",self.dictForPlotErr);
    [self presentViewController:viewControllerToPresent animated:YES completion:^{}];
    viewControllerToPresent.view.backgroundColor = [UIColor grayColor];

}

-(void)changePlotRange
{
    // Setup plot space
    CPTXYPlotSpace *plotSpace = (CPTXYPlotSpace *)graph.defaultPlotSpace;
    
    plotSpace.xRange = [CPTPlotRange plotRangeWithLocation:CPTDecimalFromFloat(0.0) length:CPTDecimalFromFloat(3.0 + 2.0 * rand() / RAND_MAX)];
    plotSpace.yRange = [CPTPlotRange plotRangeWithLocation:CPTDecimalFromFloat(0.0) length:CPTDecimalFromFloat(3.0 + 2.0 * rand() / RAND_MAX)];
}

#pragma mark -
#pragma mark Plot Data Source Methods

-(NSUInteger)numberOfRecordsForPlot:(CPTPlot *)plot
{
    return [dataForPlot count];
}

-(NSNumber *)numberForPlot:(CPTPlot *)plot field:(NSUInteger)fieldEnum recordIndex:(NSUInteger)index
{
    NSArray *contentArray = nil;
    
    if ( [(NSString *)plot.identifier isEqualToString:@"1st time"] ) {
        contentArray = dictForPlot[[NSNumber numberWithDouble:time1]];
    }
    else if ( [(NSString *)plot.identifier isEqualToString:@"2nd time"] ) {
        contentArray = dictForPlot[[NSNumber numberWithDouble:time2]];
    }
    else if ( [(NSString *)plot.identifier isEqualToString:@"3rd time"] ) {
        contentArray = dictForPlot[[NSNumber numberWithDouble:time3]];
    }
    else if ( [(NSString *)plot.identifier isEqualToString:@"real"] ) {
        contentArray = dictForPlotAnalytic[[NSNumber numberWithDouble:time4]];
    }
    return [[contentArray objectAtIndex:index] valueForKey:(fieldEnum == CPTScatterPlotFieldX ? @"x" : @"y")];
}

#pragma mark -
#pragma mark Axis Delegate Methods

-(BOOL)axis:(CPTAxis *)axis shouldUpdateAxisLabelsAtLocations:(NSSet *)locations
{
    static CPTTextStyle *positiveStyle = nil;
    static CPTTextStyle *negativeStyle = nil;
    
    NSFormatter *formatter = axis.labelFormatter;
    CGFloat labelOffset    = axis.labelOffset;
    NSDecimalNumber *zero  = [NSDecimalNumber zero];
    
    NSMutableSet *newLabels = [NSMutableSet set];
    
    for ( NSDecimalNumber *tickLocation in locations ) {
        CPTTextStyle *theLabelTextStyle;
        
        if ( [tickLocation isGreaterThanOrEqualTo:zero] ) {
            if ( !positiveStyle ) {
                CPTMutableTextStyle *newStyle = [axis.labelTextStyle mutableCopy];
                newStyle.color = [CPTColor greenColor];
                positiveStyle  = newStyle;
            }
            theLabelTextStyle = positiveStyle;
        }
        else {
            if ( !negativeStyle ) {
                CPTMutableTextStyle *newStyle = [axis.labelTextStyle mutableCopy];
                newStyle.color = [CPTColor redColor];
                negativeStyle  = newStyle;
            }
            theLabelTextStyle = negativeStyle;
        }
        
        NSString *labelString       = [formatter stringForObjectValue:tickLocation];
        CPTTextLayer *newLabelLayer = [[CPTTextLayer alloc] initWithText:labelString style:theLabelTextStyle];
        
        CPTAxisLabel *newLabel = [[CPTAxisLabel alloc] initWithContentLayer:newLabelLayer];
        newLabel.tickLocation = tickLocation.decimalValue;
        newLabel.offset       = labelOffset;
        
        [newLabels addObject:newLabel];
    }
    
    axis.axisLabels = newLabels;
    
    return NO;
}

#pragma mark - picker
- (void)pickerView:(UIPickerView *)pickerView didSelectRow: (NSInteger)row inComponent:(NSInteger)component {
    switch (component) {
        case 0:{
            time1 = [[keys objectAtIndex:row] doubleValue];
            break;
        }
        case 1:{
            time2 = [[keys objectAtIndex:row] doubleValue];
            break;
        }
        case 2:{
            time3 = [[keys objectAtIndex:row] doubleValue];
            break;
        }
        case 3:{
            time4 = [[keys objectAtIndex:row] doubleValue];
            break;
        }
        default:
            break;
    }
    [graph reloadData];
}

// tell the picker how many rows are available for a given component
- (NSInteger)pickerView:(UIPickerView *)pickerView numberOfRowsInComponent:(NSInteger)component {
    return [keys count];
}

// tell the picker how many components it will have
- (NSInteger)numberOfComponentsInPickerView:(UIPickerView *)pickerView {
    return 4;
}

// tell the picker the title for a given component
- (NSString *)pickerView:(UIPickerView *)pickerView titleForRow:(NSInteger)row forComponent:(NSInteger)component {
    return [NSString stringWithFormat:@"%.3f", [[keys objectAtIndex:row] doubleValue]];
}

// tell the picker the width of each row for a given component
- (CGFloat)pickerView:(UIPickerView *)pickerView widthForComponent:(NSInteger)component {
    int sectionWidth = 75;
    
    return sectionWidth;
}


@end
