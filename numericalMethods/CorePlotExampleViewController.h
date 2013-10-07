//
//  CorePlotExampleViewController.h
//  CorePlotExample
//
//  Created by Sylvain Marcotte on 11-02-16.
//  14 Oranges Software - http://www.14oranges.com
//
//  Original code obtained from http://www.switchonthecode.com/tutorials/using-core-plot-in-an-iphone-application
//


#import <UIKit/UIKit.h>
#import "CorePlot-CocoaTouch.h"
#import "CPTPlot.h"
#import "CPTAnnotationHostLayer.h"
#import "CPTDefinitions.h"
#import "CPTNumericDataType.h"

@interface CorePlotExampleViewController : UIViewController <CPTPlotDataSource>
{
	CPTXYGraph *graph;
	
	int numberOfRecords;
}

@end

