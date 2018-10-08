// me 
#include "vtkBrainExtractionFilter.h"
#include "vtkBrainExtractionDecorator.h"
#include "vtkPolyDataNormalsCentroids.h"
// vtk
#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkImageData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
vtkStandardNewMacro(vtkBrainExtractionFilter);
void vtkBrainExtractionFilter::PrintSelf(ostream & os, vtkIndent indent)
{
	Superclass::PrintSelf(os, indent);
	os << indent << "Subdivision: " << this->Subdivision << "\n";
}

vtkImageData * vtkBrainExtractionFilter::GetOutputImage()
{
	return vtkImageData::SafeDownCast(this->GetOutputDataObject(1));
}

vtkAlgorithmOutput * vtkBrainExtractionFilter::GetOutputPortImage()
{
	return this->GetOutputPort(1);
}

vtkBrainExtractionFilter::vtkBrainExtractionFilter()
{
	this->SetNumberOfOutputPorts(2);
	this->Subdivision = 4;
	this->IterationNumber = 0;
	this->NumOfIteration = 1000;
	this->decorator = new vtkBrainExtractionDecorator;
}

vtkBrainExtractionFilter::~vtkBrainExtractionFilter()
{
	delete this->decorator;
}

int vtkBrainExtractionFilter::RequestData(vtkInformation * vtkNotUsed(request), vtkInformationVector ** inputVector, vtkInformationVector * outputVector)
{
	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
	vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
	// get the input and output
	vtkImageData *input = vtkImageData::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output0 = vtkPolyData::SafeDownCast(
		outInfo0->Get(vtkDataObject::DATA_OBJECT()));
	vtkImageData *output1 = vtkImageData::SafeDownCast(
		outInfo1->Get(vtkDataObject::DATA_OBJECT()));
	this->decorator->generateSphere(this->Subdivision, output0);
	vtkBrainExtractionDecorator::BET_Parameters bp = this->decorator->initialParameters(input, output0, output0);
	std::cerr << bp;
	//vtkBrainExtractionDecorator::normalsCentroidNeighbourDistance(output, output);
	double fraction_threshold = 0.5;
	vtkPolyData *originalPolyData = vtkPolyData::New();
	for (this->IterationNumber = 0; this->IterationNumber < this->NumOfIteration; ++this->IterationNumber) {
		vtkBrainExtractionFilter::StepOfComputation(input, output0, 0, 1.0, 0, pow(fraction_threshold, 0.275), bp.t98, bp.t2, bp.t, bp.tm);
		std::cerr << "Interation Number: " << this->IterationNumber << '\n';
	}
	
	originalPolyData->Delete();
	output1->DeepCopy(input);
	this->decorator->generateLabelImage(output1);
	output1->ShallowCopy(this->decorator->polyDataToImage(output0, output1));
	return 1;
}

int vtkBrainExtractionFilter::FillInputPortInformation(int port, vtkInformation * info)
{
	if (port == 0) {
		info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
		return 1;
	}
	else {
		return Superclass::FillInputPortInformation(port, info);
	}
}

int vtkBrainExtractionFilter::FillOutputPortInformation(int port, vtkInformation * info)
{
	if (port == 1) {
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
		return 1;
	}
	else {
		return Superclass::FillOutputPortInformation(port, info);
	}
}

void vtkBrainExtractionFilter::StepOfComputation(
	vtkImageData * data,
	vtkPolyData * polyData,
	const int pass,
	const double &smoothArg,
	const double & increaseSmoothing,
	const double & bt,
	const double & t98, 
	const double & t2, 
	const double & t, 
	const double & tm)
{	
	this->decorator->normalsCentroidsNeighbourDistance(polyData, polyData);
	vtkIdType numPoints = polyData->GetNumberOfPoints();
	vtkDataArray *points = polyData->GetPoints()->GetData();
	vtkDataArray *normals = polyData->GetPointData()->GetArray("Normals");
	vtkDataArray *centroids = polyData->GetPointData()->GetArray("Centroids");
	vtkDataArray *meanDistance = polyData->GetPointData()->GetArray("Mean distance");
	float *points_f = reinterpret_cast<float*>(points->GetVoidPointer(0));
	float *normals_f = reinterpret_cast<float*>(normals->GetVoidPointer(0));
	float *centroids_f = reinterpret_cast<float*>(centroids->GetVoidPointer(0));
	float *meanDistance_f = reinterpret_cast<float*>(meanDistance->GetVoidPointer(0));
	double l = 0;
//////////////////////////////////////// mean distance from vertex to neighboring ////////////////////////////////////////	
	for (vtkIdType id = 0; id < meanDistance->GetNumberOfValues(); ++id) {
		l += meanDistance->GetTuple1(id);
	}
	l /= meanDistance->GetNumberOfValues();
//////////////////////////////////////// s, st, sn ////////////////////////////////////////
	//std::cerr  << "s, st, sn\n";
	vtkSmartPointer<vtkFloatArray> s =
		vtkSmartPointer<vtkFloatArray>::New();
	s->SetNumberOfComponents(3);
	s->SetName("s");
	float *s_f = s->WritePointer(0, 3 * numPoints);
	polyData->GetPointData()->AddArray(s);
	vtkSmartPointer<vtkFloatArray> sn =
		vtkSmartPointer<vtkFloatArray>::New();
	sn->SetNumberOfComponents(3);
	sn->SetName("sn");
	float *sn_f = sn->WritePointer(0, 3 * numPoints);
	polyData->GetPointData()->AddArray(sn);
	memcpy(sn_f, normals_f, sizeof(float) * 3 * numPoints);
	vtkSmartPointer<vtkFloatArray> st =
		vtkSmartPointer<vtkFloatArray>::New();
	st->SetNumberOfComponents(3);
	st->SetName("st");
	float *st_f = st->WritePointer(0, 3 * numPoints);
	polyData->GetPointData()->AddArray(st);
	for (vtkIdType id = 0; id < numPoints; ++id) {
		vtkMath::Subtract(centroids_f + 3 * id, points_f + 3 * id, s_f + 3 * id);
		float tmp = vtkMath::Dot(s_f + 3 * id, normals_f + 3 * id);
		vtkMath::MultiplyScalar(sn_f + 3 * id, tmp);
		vtkMath::Subtract(s_f + 3 * id, sn_f + 3 * id, st_f + 3 * id);
	}
//////////////////////////////////////// u1 ////////////////////////////////////////
	//std::cerr  << "u1\n";
	constexpr float f1 = 0.5f;
	vtkSmartPointer<vtkFloatArray> u1 =
		vtkSmartPointer<vtkFloatArray>::New();
	u1->SetNumberOfComponents(3);
	u1->SetName("u1");
	float *u1_f = u1->WritePointer(0, 3 * numPoints);
	polyData->GetPointData()->AddArray(u1);
	memcpy(u1_f, st_f, 3 * numPoints * sizeof(float));
	for (vtkIdType id = 0; id < numPoints; ++id) {
		vtkMath::MultiplyScalar(u1_f + id * 3, f1);
	}	
//////////////////////////////////////// u2 ////////////////////////////////////////
		//std::cerr  << "u2\n";
		//const int pass = 0;
	const double rmin = 3.33 * smoothArg;
	const double rmax = 10 * smoothArg;
	const double E = (1 / rmin + 1 / rmax) / 2.;
	const double F = 6. / (1 / rmin - 1 / rmax);
	vtkSmartPointer<vtkFloatArray> u2 =
		vtkSmartPointer<vtkFloatArray>::New();
	u2->Allocate(3 * numPoints);
	u2->SetNumberOfComponents(3);
	u2->SetNumberOfTuples(numPoints);
	u2->SetName("u2");
	polyData->GetPointData()->AddArray(u2);
	float *u2_f = u2->WritePointer(0, 3 * numPoints);
	memcpy(u2_f, sn_f, 3 * numPoints * sizeof(float));
	for (vtkIdType id = 0; id < numPoints; ++id) {
		// @todo
		// calculation in source code
		//  rinv = (2 * fabs(sn|n))/(l*l);
		// means: 
		// r = l^2 / (2 * sqqrt(norm(sn)))
		// calculation in paper
		// r = l^2 / (2 * norm(sn))
		// paper's way is the following: 
		//float r_inv = (2 * vtkMath::Norm(sn_f + 3 * numPoints)) / (l * l);
		// source code is the following: 
		float r_inv = (2 * fabs(vtkMath::Dot(sn_f + id * 3, normals_f + id * 3))) / (l * l);
		float f2 = (1 + tanh(F*(r_inv - E)))*0.5;
		if (pass > 0) {
			if (vtkMath::Dot(s_f + 3 * id, normals_f + 3 * id) > 0) {
				f2 *= increaseSmoothing;
				f2 = vtkMath::Min(f2, 1.0f);
			}
		}
		vtkMath::MultiplyScalar(u2_f + 3 * id, f2);
	}
	//////////////////////////////////////// u3 ////////////////////////////////////////
	//std::cerr  << "u3\n";
	// main term of bet. 
	//if (bt != 0.0)
	//{
	//	bt = Min(1., Max(0., bet_main_parameter + local_th*((*i)->get_coord().Z - zcog) / radius));
	//}
	// @todo 
	// in paper d1 = 20mm, d2 = d1 / 2 ;
	// paper's way is the following: 
	//const int d1 = 20;
	//const int d2 = d1 / 2;
	// source code is the following: 
	const int d1 = 7;
	const int d2 = 3;
	const float normal_max_update_fraction = 0.5f;
	const float lambda_fit = 0.1;
	const double *bounds = data->GetBounds();
	const double *origin = data->GetOrigin();
	const double *spacing = data->GetSpacing();
	const int *extent = data->GetExtent();
	vtkSmartPointer<vtkFloatArray> u3 =
		vtkSmartPointer<vtkFloatArray>::New();
	u3->SetNumberOfComponents(3);
	u3->SetName("u3");
	float *u3_f = u3->WritePointer(0, 3 * numPoints);
	polyData->GetPointData()->AddArray(u3);
	memcpy(u3_f, normals_f, 3 * numPoints * sizeof(float));
	double dscale = vtkMath::Min(vtkMath::Min(vtkMath::Min(1.0, spacing[0]), spacing[1]), spacing[2]);
	for (vtkIdType id = 0; id < numPoints; ++id) {
		// @todo 
		// calculate in paper Imin = MAX(t2, MIN(tm, I(0), I(1), ..., I(d1))), Imax = MIN(tm, MAX(t, I(0), I(1), ...I(d2)))
		// papaer's code is the following: 
		float Imin = tm;
		float Imax = t;
		for (float step = 0; step < d1; ++step) {
			float p[3];
			float n[3];
			memcpy(n, normals_f + 3 * id, sizeof(float) * 3);
			vtkMath::MultiplyScalar(n, step);
			vtkMath::Subtract(points_f + 3 * id, n, p);
			int ijk[3];
			ijk[0] = vtkMath::Round((p[0] - origin[0]) / spacing[0]);
			ijk[1] = vtkMath::Round((p[1] - origin[1]) / spacing[1]);
			ijk[2] = vtkMath::Round((p[2] - origin[2]) / spacing[2]);
			if (extent[0] <= ijk[0] && ijk[0] <= extent[1] &&
				extent[2] <= ijk[1] && ijk[1] <= extent[3] &&
				extent[4] <= ijk[2] && ijk[2] <= extent[5]) {
				float im = data->GetScalarComponentAsFloat(ijk[0], ijk[1], ijk[2], 0);
				Imin = vtkMath::Min(Imin, im);
				if (step < d2) {
					Imax = vtkMath::Max(Imax, im);
				}
			}
		}
		Imin = vtkMath::Max((float)t2, Imin);
		Imax = vtkMath::Min((float)tm, Imax);
		const float tl = (Imax - t2) * bt + t2;
		float f3 = 2 * (Imin - tl) / (Imax - t2);
		f3 *= (normal_max_update_fraction * lambda_fit * l);
		vtkMath::MultiplyScalar(u3_f + 3 * id, f3);
		// source code is the following: 
		//float Imin = tm;
		//float Imax = t;
		//float f3 = 0;
		//float p[3];
		//vtkMath::Subtract(points_f + 3 * id, normals_f + 3 * id, p);
		//float iv = vtkMath::Round((p[0] - origin[0]) / spacing[0]);
		//float jv = vtkMath::Round((p[1] - origin[1]) / spacing[1]);
		//float kv = vtkMath::Round((p[2] - origin[2]) / spacing[2]);
		//if (extent[0] <= (int)iv && (int)iv <= extent[1] &&
		//	extent[2] <= (int)jv && (int)jv <= extent[3] &&
		//	extent[4] <= (int)kv && (int)kv <= extent[5]) {
		//	float im = data->GetScalarComponentAsFloat(iv, jv, kv, 0);
		//	Imin = vtkMath::Min(Imin, im);
		//	Imax = vtkMath::Max(Imax, im);
		//	float nxv = normals_f[id * 3 + 0] / spacing[0];
		//	float nyv = normals_f[id * 3 + 1] / spacing[1];
		//	float nzv = normals_f[id * 3 + 2] / spacing[2];
		//	int i2 = iv - (d1 - 1)*nxv;
		//	int j2 = jv - (d1 - 1)*nyv;
		//	int k2 = kv - (d1 - 1)*nzv;
		//	if (extent[0] <= i2 && i2 <= extent[1] &&
		//		extent[2] <= j2 && j2 <= extent[3] &&
		//		extent[4] <= k2 && k2 <= extent[5]) {
		//		im = data->GetScalarComponentAsFloat(i2, j2, k2, 0);
		//		Imin = vtkMath::Min(Imin, im);
		//		nxv *= dscale;
		//		nyv *= dscale;
		//		nzv *= dscale;
		//		for (double gi = 2.0; gi < d1; gi += dscale)
		//		{
		//			iv -= nxv; jv -= nyv; kv -= nzv;
		//			im = data->GetScalarComponentAsFloat(iv, jv, kv, 0);
		//			Imin = vtkMath::Min(Imin, im);

		//			if (gi < d2) {
		//				Imax = vtkMath::Max(Imax, im);
		//			}
		//		}

		//		Imin = vtkMath::Max((float)t2, Imin);
		//		Imax = vtkMath::Min((float)tm, Imax);

		//		const float tl = (Imax - t2) * bt + t2;
		//		if (Imax - t2 > 0) {
		//			f3 = 2 * (Imin - tl) / (Imax - t2);
		//		}
		//		else {
		//			f3 = (Imin - tl) * 2;
		//		}
		//	}
		//}
		//f3 *= (normal_max_update_fraction * lambda_fit * l);
		//vtkMath::MultiplyScalar(u3_f + 3 * id, f3);
	}
	//////////////////////////////////////// u ////////////////////////////////////////
	//std::cerr  << "u\n";
	vtkSmartPointer<vtkFloatArray> u =
		vtkSmartPointer<vtkFloatArray>::New();
	u->SetNumberOfComponents(3);
	u->SetName("u");
	float *u_f = u->WritePointer(0, 3 * numPoints);
	polyData->GetPointData()->AddArray(u);
	for (vtkIdType id = 0; id < numPoints; ++id) {
		vtkMath::Add(u1_f + 3 * id, u2_f + 3 * id, u_f + 3 * id);
		vtkMath::Add(u_f + 3 * id, u3_f + 3 * id, u_f + 3 * id);
		vtkMath::Add(points_f + 3 * id, u_f + 3 * id, points_f + 3 * id);
	}
}
