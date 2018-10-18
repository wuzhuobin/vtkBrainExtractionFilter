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
	argc = 3;
	argv[1] = "T2.nii";
	argv[2] = "output.vtk";
	if (argc != 3) {
		cerr << "Number of input parameters should be 2. ";
		return 1;
	}
	cerr << "Since the vtkNIFTIImageReader does not care about the orientation "
		<< "and origin information of the original image. \n"
		<< "Please make sure they are zeros. \n"
		<< "OR THE RESULT MAY BE INCORRECT. \n";
	auto reader =
		vtkSmartPointer<vtkNIFTIImageReader>::New();
	reader->SetFileName("T2.nii");
	reader->Update();

	auto progress =
		vtkSmartPointer<vtkCallbackCommand>::New();
	progress->SetCallback([](vtkObject *caller, unsigned long eid,
                                     void *clientdata, void *calldata) {
		vtkBrainExtractionFilter *bef = reinterpret_cast<vtkBrainExtractionFilter*>(caller);
		//vtkWarningWithObjectMacro(bef, << "Progress: " << bef->GetProgress() << '.');
		cout << "Progress: " << bef->GetProgress() << '.' << '\n';
	});

	vtkSmartPointer<vtkBrainExtractionFilter> bef =
		vtkSmartPointer<vtkBrainExtractionFilter>::New();
	bef->SetInputData(reader->GetOutput());
	bef->AddObserver(vtkCommand::ProgressEvent, progress);
	bef->SetNumOfIteration(1000);
	bef->SetSmoothArg(1);
	bef->Update();

	vtkSmartPointer<vtkPolyDataWriter> polyDataWriter =
		vtkSmartPointer<vtkPolyDataWriter>::New();
	polyDataWriter->SetInputData(bef->GetOutput());
	polyDataWriter->SetFileName(argv[2]);
	polyDataWriter->Write();

	vtkSmartPointer<vtkNIFTIImageWriter> imageWriter =
		vtkSmartPointer<vtkNIFTIImageWriter>::New();
	imageWriter->SetInputConnection(bef->GetOutputPortImage());
	imageWriter->SetFileName((std::string(argv[2]) + ".nii").c_str());
	imageWriter->Write();

	cin.get();

	return 0;
}