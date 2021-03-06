/**
 * @file	    vtkBrainExtractFilter.h
 * @language    C++
 * @author		WUZHUOBIN jiejin2022@163.com
 * @since 		Aug.26.2017
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *			This program is distributed in the hope that it will be useful, but	 *
 *			WITHOUT ANY WARRANTY; without even the implied warranty of			 * 
 *			MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.				 * 
 *			See the LICENSE for more detail.									 * 
 *			Copyright (c) WUZHUOBIN. All rights reserved.						 * 
 *			See COPYRIGHT for more detail.										 * 
 *			This software is distributed WITHOUT ANY WARRANTY; without even		 * 
 *			the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR	 * 
 *			PURPOSE.  See the above copyright notice for more information.		 *
 *			Internal usage only, without the permission of the author, please DO *
 *			NOT publish and distribute without the author's permission.  	     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */
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
/**
 * @class		vtkBrainExtractionFilter
 * @brief 		A vtk implementation of brain extraction algorithm. 	
 * @author		WUZHUOBIN 
 * @version		0.0
 * @date		Oct.13.2018
 * @since		Oct.13.2018
 * 
 * The vtkBrainExtractionFilter is a vtk implementation of BET(Brain Extraction Tool) of 
 * FSL. The original paper of this algorithm is reference in the README.MD. 
 * The algorithm is written in vtk and try to remove the dependency of FSL's library. 
*/
class VTKBRAINEXTRACTIONFILTER_EXPORT vtkBrainExtractionFilter : public vtkPolyDataAlgorithm
{
public:	
	static vtkBrainExtractionFilter *New();
	vtkTypeMacro(vtkBrainExtractionFilter, vtkPolyDataAlgorithm);
	virtual void PrintSelf(ostream &os, vtkIndent indent) override;
	/**
	 * @fn		vtkImageData* GetOutputImage();
	 * @brief	Get the brain segmentation image.
	 * @return	A vtkImageData of brain segmentation. 
	*/
	vtkImageData* GetOutputImage();
	/**
	 * @fn		vtkAlgorithmOutput* GetOutputPortImage();
	 * @brief	Get the brain segmentaiton outputport.
	 * @return 	A output port. 
	*/
	vtkAlgorithmOutput* GetOutputPortImage();
	/**
	 * @brief	Get Subdivision. 
	 * @return	Subdivision.
	 * @see 	Subdivision
	 * 
	 * The Subdivision is the for creating the sphere source. The sphere source is created 
	 * from an icosahedron by doing a few times of subdivision/tessellation.  
	*/
	vtkGetMacro(Subdivision, int);
	/**
	 * @brief		Set Subdivision.
	 * @param[in]	Subdivision Default is 4, it was clamp 0~20. 
	 * @see 		Subdivision
	 * 
	 * The Subdivision is the for creating the sphere source. The sphere source is created 
	 * from an icosahedron by doing a few times of subdivision/tessellation.  
	*/
	vtkSetClampMacro(Subdivision, int, 0, 20);
	/**
	 * @brief		Get the number of iteraction. 
	 * @return		NumOfIteraction
	 * @see			NumOfIteraction
	 * 
	 * The number of iteraction for executing #StepOfComputation(), which is 
	 * an iteraction number. 
	*/
	vtkGetMacro(NumOfIteration, int);
	/**
	 * @brief 		Set the number of iteraction. 
	 * @param[in]	NumberOfIteraction default is 1000, clamped to 0 - VTK_INT_MAX
	 * @see			NumOfIteraction
	 * 
	 * The number of iteraction for executing #StepOfComputation(), which is 
	 * an iteraction number. 
	*/
	vtkSetClampMacro(NumOfIteration, int, 0, VTK_INT_MAX);
	/**
	 * @brief		Get current iteraction. 
	 * @return		IteractionNumber
	 * @see			IteractionNumber
	 * 
	 * Get the current iteraction number while the function thread is blocking when 
	 * executing iteraction loop. So it only possible get the number when the 
	 * algorithm finish. 
	*/
	vtkGetMacro(IterationNumber, int);

	vtkGetMacro(SmoothArg, double);
	vtkSetMacro(SmoothArg, double);

	vtkGetVector3Macro(BrainCenter, double);

	vtkGetVector3Macro(InHomogeneityDirection, double);
	vtkSetVector3Macro(InHomogeneityDirection, double);

	vtkGetMacro(SelfIntersectionChecking, bool);
	vtkSetMacro(SelfIntersectionChecking, bool);
	vtkBooleanMacro(SelfIntersectionChecking, bool);

	vtkGetMacro(AlwaysRestartComputation, bool);
	vtkSetMacro(AlwaysRestartComputation, bool);
	vtkBooleanMacro(AlwaysRestartComputation, bool);
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
	static const int d1 = 7;				///< How far into the brain the minimum intensity is searched for. In parper, d1 = 20mm.
	static const int d2 = 3;				///< How far into the brain the minimum intensity is searched for. In paper, d2 = d1/2.
	int Subdivision;						///< Subdivision of the icosahedron to generate a sphere.
	int NumOfIteration;						///< Number of iteraction
	int IterationNumber;					///< Current iteraction number. 
	double SmoothArg;						///< Smooth arguement.
	double BrainCenter[3];					///< Center of brain.
	double InHomogeneityDirection[3];		///< Inhomogeneity direction of the MRI image.
	bool SelfIntersectionChecking;
	bool AlwaysRestartComputation;
private:
	vtkBrainExtractionFilter(const vtkBrainExtractionFilter&) VTK_DELETE_FUNCTION;
	vtkBrainExtractionFilter(vtkBrainExtractionFilter&&) VTK_DELETE_FUNCTION;
	void operator=(const vtkBrainExtractionFilter&) VTK_DELETE_FUNCTION;
	void operator=(vtkBrainExtractionFilter&&) VTK_DELETE_FUNCTION;
	vtkBrainExtractionDecorator *decorator;
	BET_Parameters *bp;
};

#endif // !__VTK_BRAIN_EXTRACTION_FILTER_H__
