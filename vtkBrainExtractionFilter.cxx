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
	this->NumOfIteration = 1000;
	this->IterationNumber = 0;
	this->SmoothArg = 1.0;
	this->BrainCenter[0] = 0;
	this->BrainCenter[1] = 0;
	this->BrainCenter[2] = 0;
	this->InHomogeneityDirection[0] = 0;
	this->InHomogeneityDirection[1] = 0;
	this->InHomogeneityDirection[2] = 0;
	this->decorator = new vtkBrainExtractionDecorator;
	this->bp = new BET_Parameters;
}

vtkBrainExtractionFilter::~vtkBrainExtractionFilter()
{
	delete this->decorator;
	delete this->bp;
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
	(*this->bp) = this->decorator->initialParameters(input, output0, output0);
	cerr << *this->bp;
	this->UpdateProgress(0.1);
	double fraction_threshold = 0.5;
	vtkPolyData *originalPolyData = vtkPolyData::New();
	for (this->IterationNumber = 0; this->IterationNumber < this->NumOfIteration; ++this->IterationNumber) {
		vtkBrainExtractionFilter::StepOfComputation(
			input,
			output0,
			0,
			this->SmoothArg,
			0,
			pow(fraction_threshold, 0.275));
			//0.5);
		this->UpdateProgress(0.1 + 0.6 * this->IterationNumber / this->NumOfIteration);
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
	const double & bt)
{	
	this->decorator->normalsCentroids(polyData, polyData);
	this->decorator->mediumDistanceOfNeighbours(polyData);
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
	const double rmin = 3.33 * smoothArg;
	const double rmax = 10 * smoothArg;
	const double E = (1 / rmin + 1 / rmax) / 2.;
	const double F = 6. / (1 / rmin - 1 / rmax);
	vtkSmartPointer<vtkFloatArray> u2 =
		vtkSmartPointer<vtkFloatArray>::New();
	u2->SetNumberOfComponents(3);
	u2->SetName("u2");
	float *u2_f = u2->WritePointer(0, 3 * numPoints);
	polyData->GetPointData()->AddArray(u2);
	memcpy(u2_f, sn_f, 3 * numPoints * sizeof(float));
	for (vtkIdType id = 0; id < numPoints; ++id) {
		// @todo
		// calculation in source code
		//  rinv = (2 * fabs(sn|n))/(l*l);
		// means: 
		// r = l^2 / (2 * sqrt(norm(sn)))
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
		//float Imin = tm;
		//float Imax = t;
		//for (float step = 0; step < d1; ++step) {
		//	float p[3];
		//	float n[3];
		//	memcpy(n, normals_f + 3 * id, sizeof(float) * 3);
		//	vtkMath::MultiplyScalar(n, step);
		//	vtkMath::Subtract(points_f + 3 * id, n, p);
		//	int ijk[3];
		//	ijk[0] = vtkMath::Round((p[0] - origin[0]) / spacing[0]);
		//	ijk[1] = vtkMath::Round((p[1] - origin[1]) / spacing[1]);
		//	ijk[2] = vtkMath::Round((p[2] - origin[2]) / spacing[2]);
		//	if (extent[0] <= ijk[0] && ijk[0] <= extent[1] &&
		//		extent[2] <= ijk[1] && ijk[1] <= extent[3] &&
		//		extent[4] <= ijk[2] && ijk[2] <= extent[5]) {
		//		float im = data->GetScalarComponentAsFloat(ijk[0], ijk[1], ijk[2], 0);
		//		Imin = vtkMath::Min(Imin, im);
		//		if (step < d2) {
		//			Imax = vtkMath::Max(Imax, im);
		//		}
		//	}
		//}
		//Imin = vtkMath::Max((float)t2, Imin);
		//Imax = vtkMath::Min((float)tm, Imax);
		//const float tl = (Imax - t2) * bt + t2;
		//float f3 = 2 * (Imin - tl) / (Imax - t2);
		//f3 *= (normal_max_update_fraction * lambda_fit * l);
		//vtkMath::MultiplyScalar(u3_f + 3 * id, f3);
		// source code is the following: 
		float Imin = this->bp->tm;
		float Imax = this->bp->t;
		float f3 = 0;
		float p[3];
		vtkMath::Subtract(points_f + 3 * id, normals_f + 3 * id, p);
		float iv = (p[0] - origin[0]) / spacing[0];
		float jv = (p[1] - origin[1]) / spacing[1];
		float kv = (p[2] - origin[2]) / spacing[2];
		if (extent[0] <= vtkMath::Round(iv) && vtkMath::Round(iv) <= extent[1] &&
			extent[2] <= vtkMath::Round(jv) && vtkMath::Round(jv) <= extent[3] &&
			extent[4] <= vtkMath::Round(kv) && vtkMath::Round(kv) <= extent[5]) {
			float im = data->GetScalarComponentAsFloat(iv, jv, kv, 0);
			Imin = vtkMath::Min(Imin, im);
			Imax = vtkMath::Max(Imax, im);
			float nxv = normals_f[id * 3 + 0] / spacing[0];
			float nyv = normals_f[id * 3 + 1] / spacing[1];
			float nzv = normals_f[id * 3 + 2] / spacing[2];
			int i2 = vtkMath::Round(iv - (d1 - 1)*nxv);
			int j2 = vtkMath::Round(jv - (d1 - 1)*nyv);
			int k2 = vtkMath::Round(kv - (d1 - 1)*nzv);
			if (extent[0] <= i2 && i2 <= extent[1] &&
				extent[2] <= j2 && j2 <= extent[3] &&
				extent[4] <= k2 && k2 <= extent[5]) {
				im = data->GetScalarComponentAsFloat(i2, j2, k2, 0);
				Imin = vtkMath::Min(Imin, im);
				nxv *= dscale;
				nyv *= dscale;
				nzv *= dscale;
				for (double gi = 2.0; gi < d1; gi += dscale)
				{
					iv -= nxv; jv -= nyv; kv -= nzv;
					im = data->GetScalarComponentAsFloat(vtkMath::Round(iv), vtkMath::Round(jv), vtkMath::Round(kv), 0);
					Imin = vtkMath::Min(Imin, im);
					if (gi < d2) {
						Imax = vtkMath::Max(Imax, im);
					}
				}
				Imin = vtkMath::Max((float)(this->bp->t2), Imin);
				Imax = vtkMath::Min((float)(this->bp->tm), Imax);
				// main term of bet. 
				double _bt = bt;
				if (vtkMath::Norm(this->InHomogeneityDirection) != 0.0) {
					//bt = Min(1., Max(0., bet_main_parameter + local_th*((*i)->get_coord().Z - zcog) / radius));
					float offsetFromCenter[3];
					float center[3]{ this->bp->com.x, this->bp->com.y, this->bp->com.z };
					vtkMath::Subtract(points_f + 3 * id, center, offsetFromCenter);
					_bt = vtkMath::Min(1.0, vtkMath::Max(0.0, bt + vtkMath::Dot(offsetFromCenter, center) / this->bp->radius));
				}
				const float tl = (Imax - this->bp->t2) * _bt + this->bp->t2;
				if (Imax - this->bp->t2 > 0) {
					f3 = 2 * (Imin - tl) / (Imax - this->bp->t2);
				}
				else {
					f3 = (Imin - tl) * 2;
				}
			}
		}
		f3 *= (normal_max_update_fraction * lambda_fit * l);
		vtkMath::MultiplyScalar(u3_f + 3 * id, f3);
	}
	//////////////////////////////////////// u ////////////////////////////////////////
	vtkSmartPointer<vtkFloatArray> u =
		vtkSmartPointer<vtkFloatArray>::New();
	u->SetNumberOfComponents(3);
	u->SetName("u");
	float *u_f = u->WritePointer(0, 3 * numPoints);
	polyData->GetPointData()->AddArray(u);
	double meanMovement = 0;
	for (vtkIdType id = 0; id < numPoints; ++id) {
		vtkMath::Add(u1_f + 3 * id, u2_f + 3 * id, u_f + 3 * id);
		vtkMath::Add(u_f + 3 * id, u3_f + 3 * id, u_f + 3 * id);
		vtkMath::Add(points_f + 3 * id, u_f + 3 * id, points_f + 3 * id);
		meanMovement += vtkMath::Norm(u_f + 3 * id);
	}
	meanMovement /= numPoints;
	cerr << "meanMovement" << meanMovement << '\n';
}
