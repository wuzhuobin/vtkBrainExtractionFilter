/**
 * @file    	vtkBrainExtractionDecorator.h
 * @language    C++
 * @author  	WUZHUOBIN jiejin2022@163.com
 * @since  		Oct.13.2018
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
 *  
 * 
 */
#ifndef __VTK_BRAIN_EXTRACTION_DECORATOR_H__
#define __VTK_BRAIN_EXTRACTION_DECORATOR_H__
#pragma once
// vtk
class vtkPlatonicSolidSource;
class vtkLinearSubdivisionFilter;
class vtkImageStencil;
class vtkPolyDataToImageStencil;
class vtkPolyDataNormalsCentroids;
class vtkPolyData;
class vtkImageData;
// std
#include <ostream>
/**
* @struct	BET_Parameters
* @brief	Brain Extraction Tool Parameters
* @author	WUZHUOBIN
*
* The Bet Parameters struct is using Plain Old Data(POD), following rules for standard-layout and
* trivial.
*/
struct BET_Parameters {
	double min;					///< The minimum intensity of the image.
	double max;					///< The maximum intensity of the image.
	double t98;					///< Caculated by looking at the intensity histogram, intensity above 98%.
	double t2;					///< Caculated by looking at the intensity histogram, intensity below 2%.
	double t;					///< Attemptes to distinguish between brain matter and background, 10% between #t2 and #t98.
	double tm;					///< The median intensity of all points with a sphere of the estimated #centerOfMass and #radius.
	double radius;				///< The estimated raduis of the sphere.
	struct {
		double x;
		double y;
		double z;
	} com;						 ///< The estimated center of the sphere. The memory layout just like an array. 
};
/**
* @fn				std::ostream& operator<<(std::ostream &os, BET_Parameters &bp)
* @brief			Serializing output of BET_Parameters
* @param[in]		os The std::ostream.
* @param[in]		bp The BET_Parameters need to be print.
* @return			The input std::ostream
*/
std::ostream& operator<<(std::ostream &os, BET_Parameters &bp);
/**
 * @class 	vtkBrainExtractionDecorator
 * @brief   A decorator class for vtkBrainExtractionFilter.
 * @author	WUZHUOBIN
 * @date	Oct.13.2018
 * @since 	Oct.13.2018 
*/
class vtkBrainExtractionDecorator
{
public:
	vtkBrainExtractionDecorator();
	~vtkBrainExtractionDecorator();
	void generateSphere(const int &subdivision, vtkPolyData *data);
	void generateLabelImage(vtkImageData *image, double label = 1.0);
	vtkImageData* polyDataToImage(vtkPolyData* polyData, vtkImageData *imageData);
	BET_Parameters initialParameters(vtkImageData *imageData, vtkPolyData *polyData, vtkPolyData *output);
	void mediumDistanceOfNeighbours(vtkPolyData *polyData) const;
	void normalsCentroids(vtkPolyData *input, vtkPolyData *output);
	const double selfIntersection(vtkPolyData *original, vtkPolyData *input) const;
private: 
	vtkPlatonicSolidSource *icosahedronSource;
	vtkLinearSubdivisionFilter *linearSubdivisionFilter;
	vtkImageStencil *imageStencil;
	vtkPolyDataToImageStencil *polyDataToImageStencil;
	vtkPolyDataNormalsCentroids *polyDataNormalsCentroids;

	vtkBrainExtractionDecorator(const vtkBrainExtractionDecorator &) = delete;
	vtkBrainExtractionDecorator(const vtkBrainExtractionDecorator &&) = delete;
	vtkBrainExtractionDecorator& operator=(const vtkBrainExtractionDecorator &) = delete;
	vtkBrainExtractionDecorator& operator=(vtkBrainExtractionDecorator&&) = delete;
};

#endif // !__VTK_BRAIN_EXTRACTION_DECORATOR_H__