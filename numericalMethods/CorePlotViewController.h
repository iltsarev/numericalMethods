#import <UIKit/UIKit.h>
#import "CorePlot-CocoaTouch.h"
#import "CPTPlot.h"
#import "CPTAnnotationHostLayer.h"
#import "CPTDefinitions.h"
#import "CPTNumericDataType.h"

@interface CorePlotViewController : UIViewController <CPTPlotDataSource,CPTAxisDelegate,CPTScatterPlotDelegate,CPTPlotSpaceDelegate, UIPickerViewDelegate>
{
    CPTXYGraph *graph;
    CPTPlotSpaceAnnotation *symbolTextAnnotation;
    CPTPlotSymbol *plotSymbolC;
    NSMutableArray *dataForPlot;
    NSMutableDictionary *dictForPlot;
    NSMutableDictionary *dictForPlotErr;
    NSMutableDictionary *dictForPlotAnalytic;

    NSArray *keys;
}
@property NSArray *keys;
@property (readwrite, strong, nonatomic) NSMutableArray *dataForPlot;
@property (readwrite, strong, nonatomic) NSMutableDictionary *dictForPlot;
@property (readwrite, strong, nonatomic) NSMutableDictionary *dictForPlotErr;
@property (readwrite, strong, nonatomic) NSMutableDictionary *dictForPlotAnalytic;

@property int size;

@property double time1;
@property double time2;
@property double time3;
@property double time4;

@end

