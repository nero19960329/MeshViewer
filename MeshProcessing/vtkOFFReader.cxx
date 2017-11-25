/*=========================================================================

Program:   Visualization Toolkit
Module:    $RCSfile: vtkOFFReader.cxx,v $

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOFFReader.h"

#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTriangle.h"

#include <sstream>
#include <fstream>

vtkStandardNewMacro(vtkOFFReader);

// Description:
// Instantiate object with NULL filename.
vtkOFFReader::vtkOFFReader()
{
	this->FileName = NULL;

	this->SetNumberOfInputPorts(0);
}

vtkOFFReader::~vtkOFFReader()
{
	if (this->FileName)
	{
		delete[] this->FileName;
		this->FileName = NULL;
	}
}

/*--------------------------------------------------------

This is only partial support for the OFF format, which is
quite complicated. To find a full specification,
search the net for "OFF format", eg.:

http://people.sc.fsu.edu/~burkardt/data/off/off.html

We support only vertices and faces composed of 3 vertices.

---------------------------------------------------------*/

// a replacement for isspace()
int is_whitespace(char c)
{
	if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f')
		return 1;
	else
		return 0;
}

int vtkOFFReader::RequestData(
	vtkInformation *vtkNotUsed(request),
	vtkInformationVector **vtkNotUsed(inputVector),
	vtkInformationVector *outputVector)
{
	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the ouptut
	vtkPolyData *output = vtkPolyData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	if (!this->FileName)
	{
		vtkErrorMacro(<< "A FileName must be specified.");
		return 0;
	}

	std::ifstream in(this->FileName);

	if (!in)
	{
		vtkErrorMacro(<< "File " << this->FileName << " not found");
		return 0;
	}

	vtkDebugMacro(<< "Reading file");

	// intialise some structures to store the file contents in
	vtkPoints *points = vtkPoints::New();
	vtkCellArray *polys = vtkCellArray::New();

	// -- work through the file line by line, assigning into the above structures as appropriate --

	float xyz[3];

	std::string line;

	getline(in, line); //throw the 1st line away, it always says "OFF"

					   //Read until we get to the first line that is not a comment or a blank line.
					   //This line second states the number of vertices, faces, and edges. 
	while (getline(in, line))
	{
		//std::cout << line << std::endl;
		if (line.size() == 0)
			continue;
		if (line[0] == '#')
			continue;

		//if we get to here, this is the info line
		break;
	}
	//std::cout << "Found info line." << std::endl;
	//At this point, the line were are interested in is stored in 'line'
	//We are only interested in vertices and faces.
	std::stringstream ss;
	ss << line;
	unsigned int NumberOfVertices, NumberOfFaces;
	ss >> NumberOfVertices >> NumberOfFaces;
	//std::cout << "Vertices: " << NumberOfVertices << " Faces: " << NumberOfFaces << std::endl;
	//Now we can start to read the vertices
	unsigned int VertexCounter = 0;

	while (getline(in, line) && VertexCounter < NumberOfVertices)
	{
		if (line.size() == 0) //skip blank lines (they should only occur before the vertices start)
			continue;

		std::stringstream ssVertex;
		ssVertex << line;
		float x, y, z;
		ssVertex >> x >> y >> z;
		points->InsertNextPoint(x, y, z);
		//std::cout << "adding vertex: " << x << " " << y << " " << z << std::endl;
		VertexCounter++;
	} // (end of vertex while loop)

	unsigned int FaceCounter = 0;

	//read faces
	do
	{
		std::stringstream ssFace;
		ssFace << line;
		unsigned int NumFaceVerts, Vert0, Vert1, Vert2;
		ssFace >> NumFaceVerts >> Vert0 >> Vert1 >> Vert2;
		if (NumFaceVerts != 3)
		{
			vtkErrorMacro(<< "File " << this->FileName << " contains a face with >3 (" << NumFaceVerts << ") vertices.");
			return 0;
		}

		//std::cout << "adding triangle: " << Vert0 << " " << Vert1 << " " << Vert2 << std::endl;
		vtkTriangle* triangle = vtkTriangle::New();

		triangle->GetPointIds()->SetId(0, Vert0);
		triangle->GetPointIds()->SetId(1, Vert1);
		triangle->GetPointIds()->SetId(2, Vert2);

		polys->InsertNextCell(triangle);

		triangle->Delete();

		FaceCounter++;

	} while (getline(in, line) && FaceCounter < NumberOfFaces);

	// we have finished with the file
	in.close();

	output->SetPoints(points);
	output->SetPolys(polys);

	points->Delete();
	polys->Delete();

	return 1;
}

void vtkOFFReader::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);

	os << indent << "File Name: "
		<< (this->FileName ? this->FileName : "(none)") << "\n";
}
