// me
#include "vtkBrainExtractionDecorator.h"
#include "vtkPolyDataNormalsCentroids.h"
// vtk 
#include <vtkNew.h>
#include <vtkPlatonicSolidSource.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkArrayIteratorTemplate.h>
#include <vtkImageStencil.h>
#include <vtkImageData.h>
#include <vtkImageHistogramStatistics.h>
#include <vtkImagePointIterator.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkImageIterator.h>
#include <vtkImagePointDataIterator.h>
// std
#include <ostream>
/**
 * @fn				template<typename vtkIterator> void generateSphereInternal(vtkIterator *it)
 * @brief			Internal template implementation for #generateSphere
 * @tparam			vtkInterator 
 * @param[in, out]	it
 */
template<typename vtkIterator>
static void generateSphereInternal(vtkIterator *it) {
	// template concept checking. 
	static_assert(!std::is_base_of<vtkIterator, vtkArrayIterator>::value, 
		"Error tempate type, type should be one of the followings: "
		"vtkArrayIteratorTemplate<double>, vtkArrayIteratorTemplate<float>. ");
	for (vtkIdType id = 0; id < it->GetNumberOfTuples(); ++id) {
		typename vtkIterator::ValueType *tuple = it->GetTuple(id);
		vtkMath::Normalize(tuple);
	}
}

std::ostream& operator<<(std::ostream &os, BET_Parameters &bp)
{
	os << "\t\t\t" << '\n';
	os << "\t\t\t" << "BET Parameters" << '\n';
	os << "\t\t\t" << "min " << bp.min << '\n';
	os << "\t\t\t" << "max " << bp.max << '\n';
	os << "\t\t\t" << "t98 " << bp.t98 << '\n';
	os << "\t\t\t" << "t2 " << bp.t2 << '\n';
	os << "\t\t\t" << "t " << bp.t << '\n';
	os << "\t\t\t" << "tm " << bp.tm << '\n';	
	os << "\t\t\t" << "radius " << bp.radius << '\n';
	os << "\t\t\t" << "center of mass " << bp.com.x << '\n';
	os << "\t\t\t" << "center of mass " << bp.com.y << '\n';
	os << "\t\t\t" << "center of mass " << bp.com.z << '\n';
	os << "\t\t\t" << '\n';
	return os;
}

vtkBrainExtractionDecorator::vtkBrainExtractionDecorator()
{
	// create Icosahedron, whose number of faces is twenty
	this->icosahedronSource = vtkPlatonicSolidSource::New();
	this->icosahedronSource->SetSolidTypeToIcosahedron();
	// re-tesselates.
	this->linearSubdivisionFilter = vtkLinearSubdivisionFilter::New();
	this->linearSubdivisionFilter->SetInputConnection(this->icosahedronSource->GetOutputPort());
	this->polyDataToImageStencil = vtkPolyDataToImageStencil::New();
	this->imageStencil = vtkImageStencil::New();
	this->imageStencil->SetStencilConnection(polyDataToImageStencil->GetOutputPort());
	this->imageStencil->ReverseStencilOff();
	this->imageStencil->SetBackgroundValue(0.0);
	this->polyDataNormalsCentroids = vtkPolyDataNormalsCentroids::New();
	this->polyDataNormalsCentroids->SetComputeCellCentroids(true);
	this->polyDataNormalsCentroids->SetComputeCellNormals(true);
	this->polyDataNormalsCentroids->SetComputePointCentroids(true);
	this->polyDataNormalsCentroids->SetComputePointNormals(true);
	this->polyDataNormalsCentroids->SetConsistency(true);
	this->polyDataNormalsCentroids->SetSplitting(false);
	this->polyDataNormalsCentroids->SetAutoOrientNormals(true);
	this->polyDataNormalsCentroids->SetFlipNormals(false);
}

vtkBrainExtractionDecorator::~vtkBrainExtractionDecorator()
{
	this->icosahedronSource->Delete();
	this->linearSubdivisionFilter->Delete();
	this->polyDataToImageStencil->Delete();
	this->imageStencil->Delete();
	this->polyDataNormalsCentroids->Delete();
}

void vtkBrainExtractionDecorator::generateSphere(const int & subdivision, vtkPolyData * data)
{
	this->linearSubdivisionFilter->SetNumberOfSubdivisions(subdivision);
	this->linearSubdivisionFilter->Update();
	data->ShallowCopy(this->linearSubdivisionFilter->GetOutput());
	data->GetCellData()->RemoveArray(0);
	vtkPoints *points = data->GetPoints();
	vtkArrayIterator *it = points->GetData()->NewIterator();
	switch (points->GetDataType())
	{
		vtkArrayIteratorTemplateMacroCase(VTK_DOUBLE, double, generateSphereInternal(static_cast<VTK_TT*>(it)));
		vtkArrayIteratorTemplateMacroCase(VTK_FLOAT, float, generateSphereInternal(static_cast<VTK_TT*>(it)));
	default:
	{
		throw std::invalid_argument(
		"Error tempate type, type should be one of the followings: "
		"vtkArrayIteratorTemplate<double>, vtkArrayIteratorTemplate<float>. ");

	}
	}
	it->Delete();
	data->Modified();
}

void vtkBrainExtractionDecorator::generateLabelImage(vtkImageData * image, double label)
{
	typedef double VTK_TT;
	switch (image->GetScalarType())
	{
		vtkTemplateMacro(
			for (vtkImageIterator<VTK_TT> it(image, image->GetExtent());
				!it.IsAtEnd(); it.NextSpan()) {
			for (VTK_TT *v = it.BeginSpan(); v != it.EndSpan(); ++v) {
				*v = label;
			}
		});
	default:
		break;
	}
}

vtkImageData * vtkBrainExtractionDecorator::polyDataToImage(vtkPolyData * polyData, vtkImageData * imageData)
{
	this->polyDataToImageStencil->SetInputData(polyData);
	this->polyDataToImageStencil->SetOutputOrigin(imageData->GetOrigin());
	this->polyDataToImageStencil->SetOutputSpacing(imageData->GetSpacing());
	this->polyDataToImageStencil->SetOutputWholeExtent(imageData->GetExtent());
	this->polyDataToImageStencil->Update();
	this->imageStencil->SetInputData(imageData);
	this->imageStencil->SetStencilConnection(this->polyDataToImageStencil->GetOutputPort());
	this->imageStencil->Update();
	return this->imageStencil->GetOutput();
}

BET_Parameters vtkBrainExtractionDecorator::initialParameters(vtkImageData * imageData, vtkPolyData * polyData, vtkPolyData * output)
{
	BET_Parameters bp;
	// calculate bp.max, bp.min, bp.t, bp.t2, bp.t98;
	vtkNew<vtkImageHistogramStatistics> imageHistogram;
	imageHistogram->SetInputData(imageData);
	imageHistogram->SetAutoRangePercentiles(2, 98);
	imageHistogram->SetAutoRangeExpansionFactors(0, 0);
	imageHistogram->SetGenerateHistogramImage(false);
	imageHistogram->Update();
	imageHistogram->GetAutoRange(bp.t2, bp.t98);
	bp.max = imageHistogram->GetMaximum();
	bp.min = imageHistogram->GetMinimum();
	bp.t = bp.t2 + .1*(bp.t98 - bp.t2);
	// calculate bp.radius
	vtkIdType numIsBrain = 0;
	const double *spacing = imageData->GetSpacing();
	// calculate bp.centerOfMass
	double centerOfMass[3]{ 0, 0, 0 };
	double counter = 0;
	vtkImagePointIterator it;
	for (it.Initialize(imageData); !it.IsAtEnd(); it.Next()) {
		double voxel;
		void *voidPointer = vtkImagePointIterator::GetVoidPointer(imageData, it.GetId());
		switch (imageData->GetScalarType())
		{
			vtkTemplateMacro(voxel = *static_cast<VTK_TT*>(voidPointer));
		}
		if (voxel > bp.t) {
			voxel = vtkMath::Min(voxel, bp.t98) - bp.t2;
			counter += voxel;
			centerOfMass[0] += voxel * it.GetPosition()[0];
			centerOfMass[1] += voxel * it.GetPosition()[1];
			centerOfMass[2] += voxel * it.GetPosition()[2];
			++numIsBrain;
		}
	}
	vtkMath::MultiplyScalar(centerOfMass, 1 / counter);
	//std::copy(centerOfMass, centerOfMass + 3, bp.centerOfMass.begin());
	memcpy(&bp.com, centerOfMass, sizeof(double) * 3);
	// using the brain's volume as the sphere's volume
	// basing on sphere volume equation, find radius.
	// sphere volume equation: 
	// volume = 4/3 * PI * radius^3
	double brainVolume = numIsBrain * spacing[0] * spacing[1] * spacing[2];
	bp.radius = pow(brainVolume * 0.75 / vtkMath::Pi(), 1.0 / 3.0);
	// transform
	vtkSmartPointer<vtkTransform> transform =
		vtkSmartPointer<vtkTransform>::New();
	transform->Identity();
	transform->Translate(centerOfMass);
	transform->Scale(bp.radius * 0.5, bp.radius * 0.5, bp.radius * 0.5);
	vtkSmartPointer<vtkTransformPolyDataFilter> transformPolyData =
		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformPolyData->SetInputData(polyData);
	transformPolyData->SetTransform(transform);
	transformPolyData->Update();
	output->ShallowCopy(transformPolyData->GetOutput());
	std::vector<double> voxels;
	double center[3];
	memcpy(center, static_cast<void*>(&bp.com), sizeof(bp.com));
	for (it.Initialize(imageData); !it.IsAtEnd(); it.Next()) {
		double voxel;
		void *voidPointer = vtkImagePointIterator::GetVoidPointer(imageData, it.GetId());
		if (vtkMath::Distance2BetweenPoints(it.GetPosition(), center) >= (bp.radius * bp.radius)) {
			continue;
		}
		switch (imageData->GetScalarType())
		{
			vtkTemplateMacro(voxel = *static_cast<VTK_TT*>(voidPointer));
		}
		if (voxel <= bp.t2 || bp.t98 <= voxel) {
			continue;
		}
		voxels.push_back(voxel);
	}
	// find medium;
	std::nth_element(voxels.begin(), voxels.begin() + floor(voxels.size() * 0.5) - 1, voxels.end());
	bp.tm = voxels[floor(voxels.size() * 0.5) - 1];
	return bp;
}

void vtkBrainExtractionDecorator::mediumDistanceOfNeighbours(vtkPolyData * polyData) const
{
	// Compute point neighbour distance, which are the means of a point's distance of neighbour points
	// The name of the array is "Mean distance"
	polyData->BuildCells();
	vtkNew<vtkFloatArray> meanDistance;
	meanDistance->Allocate(polyData->GetNumberOfPoints());
	meanDistance->SetNumberOfComponents(1);
	meanDistance->SetNumberOfTuples(polyData->GetNumberOfPoints());
	meanDistance->FillComponent(0, 0);
	meanDistance->SetName("Mean distance");
	polyData->GetPointData()->AddArray(meanDistance.Get());
	float *meanDistance_f = reinterpret_cast<float*>(meanDistance->GetVoidPointer(0));
	vtkIdType *counter = new vtkIdType[polyData->GetNumberOfPoints()]();
	for (vtkIdType cid = 0; cid < polyData->GetNumberOfCells(); ++cid) {
		vtkCell *cell = polyData->GetCell(cid);
		vtkIdType numPts;
		vtkIdType *points;
		polyData->GetCellPoints(cid, numPts, points);
		if (numPts < 2) {
			continue;
		}
		for (vtkIdType pid = 0; pid < numPts; ++pid) {
			// current point's neighbour
			vtkIdType before = (pid == 0 ? numPts - 1 : pid - 1);
			vtkIdType after = (pid == numPts - 1 ? 0 : pid + 1);
			counter[points[pid]] += 2;
			double beforePts[3];
			polyData->GetPoint(points[before], beforePts);
			double afterPts[3];
			polyData->GetPoint(points[after], afterPts);
			double currentPts[3];
			polyData->GetPoint(points[pid], currentPts);
			meanDistance_f[points[pid]] +=
				std::sqrt(vtkMath::Distance2BetweenPoints(currentPts, beforePts));
			meanDistance_f[points[pid]] +=
				std::sqrt(vtkMath::Distance2BetweenPoints(currentPts, afterPts));
		}
	}
	for (vtkIdType pid = 0; pid < polyData->GetNumberOfPoints(); ++pid) {
		if (counter[pid] != 0) {
			meanDistance_f[pid] /= counter[pid];
		}
	}
	delete counter;
}

void vtkBrainExtractionDecorator::normalsCentroids(vtkPolyData * input, vtkPolyData * output)
{
	// Compute point normals, which are the means of a points's all cell normals,
	// point inside to center of closed surface.
	// The name of the array is "Normal". 
	// Compute point centroids, which are the means of a point's neighbour points.
	// The name of the array is "Centroid". 
	this->polyDataNormalsCentroids->SetInputData(input);
	this->polyDataNormalsCentroids->Update();
	vtkPolyData *polyData = this->polyDataNormalsCentroids->GetOutput();
	output->ShallowCopy(polyData);
}

const double vtkBrainExtractionDecorator::selfIntersection(vtkPolyData * original, vtkPolyData * input) const
{
	double intersection = 0.0;
	vtkIdType numPoints_o = original->GetNumberOfPoints();
	vtkIdType numPoints = input->GetNumberOfPoints();
	if (numPoints_o != numPoints) {
		throw std::invalid_argument("The original polydata and the new polydata have different number of poitns. ");
	}
	// calculate the mean of medium distance of Neighbours of the original poly data. 
	this->mediumDistanceOfNeighbours(original);
	vtkDataArray *originalDistances = original->GetPointData()->GetArray("Mean distance");
	double mlo = 0.0;
	for (vtkIdType id = 0; id < numPoints_o; ++id) {
		mlo += originalDistances->GetTuple1(id);
	}
	mlo /= numPoints_o;
	// calculate the mean of medium distance of Neighbours of the input poly data. 
	this->mediumDistanceOfNeighbours(input);
	vtkDataArray *inputDistances = input->GetPointData()->GetArray("Mean distance");
	double ml = 0.0;
	for (vtkIdType id = 0; id < numPoints; ++id) {
		ml += originalDistances->GetTuple1(id);
	}
	ml /= numPoints;
	input->BuildLinks();
	input->BuildCells();
	vtkIdType io;
	vtkIdType i;
	double pio[3];
	double pi[3];
	for (io = 0, i = 0; io < numPoints_o && i < numPoints; ++io, ++i) {
		original->GetPoint(io, pio);
		input->GetPoint(i, pi);
		vtkIdType jo;
		vtkIdType j;
		double pjo[3];
		double pj[3];
		for (jo = 0, j = 0; jo < numPoints_o && j < numPoints; ++jo, ++j) {
			original->GetPoint(jo, pjo);
			input->GetPoint(io, pio);
			bool found = false;
			unsigned short ncells;
			vtkIdType* cells;
			input->GetPointCells(i, ncells, cells);
			for (vtkIdType cid = 0; cid < ncells && !found; ++cid) {
				vtkIdType npts; 
				vtkIdType* pts;
				input->GetCellPoints(cid, npts, pts);
				for (vtkIdType pid = 0; pid < npts && !found; ++pid) {
					double neighbourOrItself[3];
					const double TOLERANCE = 1e-8;
					input->GetPoint(pid, neighbourOrItself);
					if (fabs(pj[0] - neighbourOrItself[0]) >= TOLERANCE ||
						fabs(pj[1] - neighbourOrItself[1]) >= TOLERANCE ||
						fabs(pj[2] - neighbourOrItself[2]) >= TOLERANCE) {
						found = true;
					}
				}
			}
			if (!found) {
				double dist = vtkMath::Distance2BetweenPoints(pi, pj);
				if (dist < ml * ml) {
					dist = sqrt(dist) / ml;
					double disto = vtkMath::Distance2BetweenPoints(pio, pjo);
					disto = sqrt(disto) / mlo;
					intersection += (dist - disto) * (dist - disto);
				}
			}
		}
	}
	return intersection;
}
