#ifndef __VTK_BRAIN_EXTRACTION_DECORATOR_H__
#define __VTK_BRAIN_EXTRACTION_DECORATOR_H__
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
class vtkBrainExtractionDecorator
{
public:
	vtkBrainExtractionDecorator();
	~vtkBrainExtractionDecorator();
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
	void generateSphere(const int &subdivision, vtkPolyData *data);
	void generateLabelImage(vtkImageData *image, double label = 1.0);
	vtkImageData* polyDataToImage(vtkPolyData* polyData, vtkImageData *imageData);
	BET_Parameters initialParameters(vtkImageData *imageData, vtkPolyData *polyData, vtkPolyData *output);
	void normalsCentroidsNeighbourDistance(vtkPolyData *input, vtkPolyData *output);
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
/**
* @fn				std::ostream& operator<<(std::ostream &os, BET_Parameters &bp)
* @brief			Serializing output of BET_Parameters
* @param[in]		os The std::ostream.
* @param[in]		bp The BET_Parameters need to be print.
* @return			The input std::ostream
*/
std::ostream& operator<<(std::ostream &os, vtkBrainExtractionDecorator::BET_Parameters &bp);

//void binaryImageDataSource(
//	vtkImageData *image,
//	const double origin[3],
//	const double spacing[3],
//	const int extent[6],
//	int scalarType,
//	double label) {
//	image->SetOrigin(origin[0], origin[1], origin[2]);
//	image->SetSpacing(spacing[0], spacing[1], spacing[2]);
//	image->SetExtent(extent[0], extent[1], extent[2], extent[3], extent[4], extent[5]);
//	image->AllocateScalars(scalarType, 1);
//	typedef double VTK_TT;
//	for (vtkImageIterator<VTK_TT> it(image, image->GetExtent());
//		!it.IsAtEnd(); it.NextSpan()) {
//		for (VTK_TT *v = it.BeginSpan(); v != it.EndSpan(); ++v) {
//			*v = label;
//		}
//	}
//}

#endif // !__VTK_BRAIN_EXTRACTION_DECORATOR_H__