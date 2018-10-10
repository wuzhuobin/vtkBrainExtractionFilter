#include <vtkBrainExtractionFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkNIFTIImageReader.h>
#include <vtkNIFTIImageWriter.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkCallbackCommand.h>
int main(int argc, char **argv)
{
	vtkSmartPointer<vtkPolyData> polyData =
		vtkSmartPointer<vtkPolyData>::New();

	auto reader =
		vtkSmartPointer<vtkNIFTIImageReader>::New();
	reader->SetFileName("T2.nii");
	reader->Update();

	auto progress =
		vtkSmartPointer<vtkCallbackCommand>::New();
	progress->SetCallback([](vtkObject *caller, unsigned long eid,
                                     void *clientdata, void *calldata) {
		vtkBrainExtractionFilter *bef = reinterpret_cast<vtkBrainExtractionFilter*>(caller);
		vtkWarningWithObjectMacro(bef, << "Progress: " << bef->GetProgress() << '.');
	});

	vtkSmartPointer<vtkBrainExtractionFilter> bef =
		vtkSmartPointer<vtkBrainExtractionFilter>::New();
	bef->SetInputData(reader->GetOutput());
	bef->AddObserver(vtkCommand::ProgressEvent, progress);
	bef->SetNumOfIteration(1000);
	bef->SetSmoothArg(1);
	bef->Update();
	polyData = bef->GetOutput();

	vtkSmartPointer<vtkPolyDataWriter> polyDataWriter =
		vtkSmartPointer<vtkPolyDataWriter>::New();
	polyDataWriter->SetInputData(polyData);
	polyDataWriter->SetFileName("sphere.vtk");
	polyDataWriter->Write();

	vtkSmartPointer<vtkNIFTIImageWriter> imageWriter =
		vtkSmartPointer<vtkNIFTIImageWriter>::New();
	imageWriter->SetInputConnection(bef->GetOutputPortImage());
	imageWriter->SetFileName("sphere.nii");
	imageWriter->Write();

	cin.get();

	return 0;
}