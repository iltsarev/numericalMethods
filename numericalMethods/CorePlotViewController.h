#import <UIKit/UIKit.h>
#import "CorePlot-CocoaTouch.h"
#import "CPTPlot.h"
#import "CPTAnnotationHostLayer.h"
#import "CPTDefinitions.h"
#import "CPTNumericDataType.h"

@interface CorePlotViewController : UIViewController <CPTPlotDataSource,CPTAxisDelegate,CPTScatterPlotDelegate,CPTPlotSpaceDelegate>
{
    CPTXYGraph *graph;
    CPTPlotSpaceAnnotation *symbolTextAnnotation;
    CPTPlotSymbol *plotSymbolC;
    NSMutableArray *dataForPlot;
    NSMutableDictionary *dictForPlot;
}

@property (readwrite, strong, nonatomic) NSMutableArray *dataForPlot;
@property (readwrite, strong, nonatomic) NSMutableDictionary *dictForPlot;

@end

