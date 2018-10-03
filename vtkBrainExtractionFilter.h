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
class vtkBrainExtractionFilter : public vtkPolyDataAlgorithm
{
public:
	static vtkBrainExtractionFilter *New();
	vtkTypeMacro(vtkBrainExtractionFilter, vtkPolyDataAlgorithm);
	virtual void PrintSelf(ostream &os, vtkIndent indent) override;
	vtkImageData* GetOutputImage();
	vtkAlgorithmOutput* GetOutputPortImage();
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
		const double &betMainParameter,
		const double &t98,
		const double &t2,
		const double &t,
		const double &tm);
//////////////////////////////////////// Parameter ////////////////////////////////////////
	int Subdivision;
	int NumOfIteration;
	int IterationNumber;
	double BrainCenter[3];
	vtkBrainExtractionDecorator *decorator;
private:
	vtkBrainExtractionFilter(const vtkBrainExtractionFilter&) VTK_DELETE_FUNCTION;
	vtkBrainExtractionFilter(vtkBrainExtractionFilter&&) VTK_DELETE_FUNCTION;
	void operator=(const vtkBrainExtractionFilter&) VTK_DELETE_FUNCTION;
	void operator=(vtkBrainExtractionFilter&&) VTK_DELETE_FUNCTION;
};

#endif // !__VTK_BRAIN_EXTRACTION_FILTER_H__
