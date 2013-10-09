#import <UIKit/UIKit.h>
#import "CorePlot-CocoaTouch.h"
#import "CPTPlot.h"
#import "CPTAnnotationHostLayer.h"
#import "CPTDefinitions.h"
#import "CPTNumericDataType.h"

@interface CorePlotViewController : UIViewController <CPTPlotDataSource,CPTAxisDelegate,CPTScatterPlotDelegate>
{
    CPTXYGraph *graph;
    
    NSMutableArray *dataForPlot;
}

@property (readwrite, strong, nonatomic) NSMutableArray *dataForPlot;
@end

