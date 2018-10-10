#ifndef __VTK_BRAIN_EXTRACTION_FILTER_H__
#define __VTK_BRAIN_EXTRACTION_FILTER_H__
#pragma once
// me
#include "vtkbrainextractionfilter_export.h"
// vtk
#include <vtkPolyDataAlgorithm.h>
class vtkImageData;
class vtkPolyData;
class vtkBrainExtractionDecorator;
struct BET_Parameters;
class vtkBrainExtractionFilter : public vtkPolyDataAlgorithm
{
public:
	static vtkBrainExtractionFilter *New();
	vtkTypeMacro(vtkBrainExtractionFilter, vtkPolyDataAlgorithm);
	virtual void PrintSelf(ostream &os, vtkIndent indent) override;
	vtkImageData* GetOutputImage();
	vtkAlgorithmOutput* GetOutputPortImage();
	
	vtkGetMacro(Subdivision, int);
	vtkSetClampMacro(Subdivision, int, 0, 20);

	vtkGetMacro(NumOfIteration, int);
	vtkSetClampMacro(NumOfIteration, int, 0, VTK_INT_MAX);

	vtkGetMacro(IterationNumber, int);

	vtkGetMacro(SmoothArg, double);
	vtkSetMacro(SmoothArg, double);

	vtkGetVector3Macro(BrainCenter, double);
protected:
	vtkBrainExtractionFilter();
	virtual ~vtkBrainExtractionFilter() override;
	virtual int RequestData(vtkInformation *request,
		vtkInformationVector **inputVector,
		vtkInformationVector *outputVector) override; 
	int FillInputPortInformation(int port, vtkInformation *info) override;	
	int FillOutputPortInformation(int port, vtkInformation *info) override;
	void StepOfComputation(
		vtkImageData *data,
		vtkPolyData *polyData,
		const int pass,
		const double &smoothArg,
		const double &increaseSmoothing,
		const double &betMainParameter);
//////////////////////////////////////// Parameter ////////////////////////////////////////
	int Subdivision;
	int NumOfIteration;
	int IterationNumber;
	double SmoothArg;
	double BrainCenter[3];
	double InHomogeneityDirection[3];
private:
	vtkBrainExtractionFilter(const vtkBrainExtractionFilter&) VTK_DELETE_FUNCTION;
	vtkBrainExtractionFilter(vtkBrainExtractionFilter&&) VTK_DELETE_FUNCTION;
	void operator=(const vtkBrainExtractionFilter&) VTK_DELETE_FUNCTION;
	void operator=(vtkBrainExtractionFilter&&) VTK_DELETE_FUNCTION;
	vtkBrainExtractionDecorator *decorator;
	BET_Parameters *bp;
};

#endif // !__VTK_BRAIN_EXTRACTION_FILTER_H__
