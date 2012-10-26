// -*- C++ -*-
/**
		@file		hotOcean.cpp
		@since		2012-10-26

		@author		Nico Rehberg, Imre Tuske, Szabolcs Horvatth

		@brief		Maya deformer to displace a surface using Houdini Ocean Toolkit (implementation).

		Implementation: Nico Rehberg <mail@nico-rehberg.de>

		The Houdini Ocean Toolkit is copyrighted by Drew Whitehouse
		see http://odforce.net/wiki/index.php/HoudiniOceanToolkit

		This program is free software; you can redistribute it and/or
		modify it under the terms of the GNU General Public License
		as published by the Free Software Foundation; either version 2
		of the License, or (at your option) any later version.

		This program is distributed in the hope that it will be useful,
		but WITHOUT ANY WARRANTY; without even the implied warranty of
		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
		GNU General Public License for more details.
*/

#include "hotOcean.h"
#include <maya/MGlobal.h>
#include <stdexcept>

#define M_PI 3.14159265358979323846

MTypeId		hotOceanDeformer::id( 0x0007443a );

MObject		hotOceanDeformer::globalScale;
MObject		hotOceanDeformer::resolution;
MObject		hotOceanDeformer::size;
MObject		hotOceanDeformer::windSpeed;
MObject		hotOceanDeformer::waveHeight;
MObject		hotOceanDeformer::shortestWave;
MObject		hotOceanDeformer::choppiness;
MObject		hotOceanDeformer::windDirection;
MObject		hotOceanDeformer::dampReflections;
MObject		hotOceanDeformer::windAlign;
MObject		hotOceanDeformer::oceanDepth;
MObject		hotOceanDeformer::time;
MObject		hotOceanDeformer::seed;
MObject		hotOceanDeformer::interpolation;
MObject		hotOceanDeformer::deformSpace;
MObject		hotOceanDeformer::vertexColor;
MObject		hotOceanDeformer::doJMinus;
MObject		hotOceanDeformer::doJPlus;
MObject		hotOceanDeformer::doEMinus;
MObject		hotOceanDeformer::doEPlus;



hotOceanDeformer::hotOceanDeformer()
: _ocean(0)
, _ocean_context(0)
, _ocean_scale(1.0)
, _ocean_needs_rebuild(true)
, _mesh_changed(true)
, _initTangentSpace(true)
{
}



hotOceanDeformer::~hotOceanDeformer()
{
	if (_ocean) {
		delete _ocean;
	}

	if (_ocean_context) {
		delete _ocean_context;
	}
}



/**	Attribute initialization.
*/
MStatus hotOceanDeformer::initialize()
{
	MFnNumericAttribute nAttr;
	globalScale=nAttr.create( "globalScale", "scale", MFnNumericData::kDouble  );
	nAttr.setDefault(1.0);
	nAttr.setKeyable(false);
	nAttr.setMin(0.00001);
	nAttr.setSoftMax(10);
	nAttr.setChannelBox(true);
	addAttribute( globalScale );

	resolution=nAttr.create( "resolution", "res", MFnNumericData::kInt  );
	nAttr.setDefault(6);
	nAttr.setKeyable(false);
	nAttr.setMin(4);
	nAttr.setMax(11);
	nAttr.setChannelBox(true);
	addAttribute( resolution );

	size=nAttr.create( "size", "size", MFnNumericData::kDouble  );
	nAttr.setDefault(200.0);
	nAttr.setKeyable(false);
	nAttr.setMin(0.01);
	nAttr.setSoftMax(1000);
	nAttr.setChannelBox(true);
	addAttribute( size );

	windSpeed = nAttr.create( "windSpeed", "ws", MFnNumericData::kDouble );
	nAttr.setDefault(15.0);
	nAttr.setKeyable(true);
	nAttr.setMin(0.01);
	nAttr.setSoftMax(100);
	addAttribute( windSpeed );

	waveHeight = nAttr.create( "waveHeight", "wh", MFnNumericData::kDouble );
	nAttr.setDefault(2.0);
	nAttr.setKeyable(true);
	nAttr.setMin(0.01);
	nAttr.setSoftMax(20);
	addAttribute( waveHeight );

	shortestWave = nAttr.create( "shortestWave", "sw", MFnNumericData::kDouble );
	nAttr.setDefault(0.01);
	nAttr.setKeyable(true);
	nAttr.setMin(0.00001);
	nAttr.setSoftMax(20);
	addAttribute( shortestWave );

	choppiness = nAttr.create( "choppiness", "chop", MFnNumericData::kDouble );
	nAttr.setDefault(0.7);
	nAttr.setKeyable(true);
	nAttr.setMin(0.0);
	nAttr.setSoftMax(4.0);
	addAttribute( choppiness );

	windDirection = nAttr.create( "windDirection", "dir", MFnNumericData::kDouble );
	nAttr.setDefault(0);
	nAttr.setKeyable(true);
	nAttr.setMin(0.0);
	nAttr.setSoftMax(360.0);
	addAttribute( windDirection );

	dampReflections = nAttr.create( "dampReflections", "damp", MFnNumericData::kDouble );
	nAttr.setDefault(0.5);
	nAttr.setKeyable(true);
	nAttr.setMin(0.0);
	nAttr.setMax(1.0);
	addAttribute( dampReflections );

	windAlign = nAttr.create( "windAlign", "wa", MFnNumericData::kDouble );
	nAttr.setDefault(2.0);
	nAttr.setKeyable(true);
	nAttr.setMin(0.01);
	nAttr.setSoftMax(10.0);
	addAttribute( windAlign );

	oceanDepth = nAttr.create( "oceanDepth", "dpt", MFnNumericData::kDouble );
	nAttr.setDefault(200.0);
	nAttr.setKeyable(true);
	nAttr.setMin(0.01);
	nAttr.setSoftMax(1000.0);
	addAttribute( oceanDepth );

	time = nAttr.create( "time", "t", MFnNumericData::kDouble );
	nAttr.setDefault(0.0);
	nAttr.setKeyable(true);
	nAttr.setSoftMin(0.0);
	nAttr.setSoftMax(1000.0);
	addAttribute( time );

	seed = nAttr.create( "seed", "s", MFnNumericData::kInt );
	nAttr.setDefault(0);
	nAttr.setKeyable(true);
	nAttr.setMin(0);
	nAttr.setSoftMax(9999);
	addAttribute( seed );

	interpolation = nAttr.create( "interpolation", "int", MFnNumericData::kBoolean );
	nAttr.setDefault(false);
	nAttr.setKeyable(false);
	nAttr.setChannelBox(true);
	addAttribute( interpolation );

	MFnEnumAttribute eAttr;
	deformSpace = eAttr.create( "space", "space", 0 );  // can set default 0, 1, 2, etc..
	eAttr.setStorable(true);
	eAttr.setKeyable(false);
	eAttr.setChannelBox(true);
	eAttr.addField("world", 0);
	eAttr.addField("object", 1);
	eAttr.addField("tangent", 2);
	addAttribute(deformSpace);

	vertexColor = nAttr.create("vertexColor", "vertexColor", MFnNumericData::kBoolean, false);
	nAttr.setKeyable(false); nAttr.setChannelBox(true);
	addAttribute(vertexColor);

	doJMinus = nAttr.create("doJMinus", "dojmn", MFnNumericData::kBoolean, false);
	nAttr.setKeyable(false); nAttr.setChannelBox(true);
	addAttribute(doJMinus);

	doJPlus = nAttr.create("doJPlus", "dojps", MFnNumericData::kBoolean, false);
	nAttr.setKeyable(false); nAttr.setChannelBox(true);
	addAttribute(doJPlus);

	doEMinus = nAttr.create("doEMinus", "doemn", MFnNumericData::kBoolean, false);
	nAttr.setKeyable(false); nAttr.setChannelBox(true);
	addAttribute(doEMinus);

	doEPlus = nAttr.create("doEPlus", "doeps", MFnNumericData::kBoolean, false);
	nAttr.setKeyable(false); nAttr.setChannelBox(true);
	addAttribute(doEPlus);

	attributeAffects( globalScale, outputGeom );
	attributeAffects( resolution, outputGeom );
	attributeAffects( size, outputGeom );
	attributeAffects( windSpeed, outputGeom );
	attributeAffects( waveHeight, outputGeom );
	attributeAffects( shortestWave, outputGeom );
	attributeAffects( choppiness, outputGeom );
	attributeAffects( windDirection, outputGeom );
	attributeAffects( dampReflections, outputGeom );
	attributeAffects( windAlign, outputGeom );
	attributeAffects( oceanDepth, outputGeom );
	attributeAffects( time, outputGeom );
	attributeAffects( seed, outputGeom );
	attributeAffects( interpolation, outputGeom );
	attributeAffects( deformSpace, outputGeom );
	attributeAffects( vertexColor, outputGeom );
	attributeAffects( doJMinus, outputGeom );
	attributeAffects( doJPlus, outputGeom );
	attributeAffects( doEMinus, outputGeom );
	attributeAffects( doEPlus, outputGeom );

	return MS::kSuccess;
}



MStatus hotOceanDeformer::setDependentsDirty(
	const MPlug &plugBeingDirtied,
	MPlugArray &affectedPlugs )
{
	if (!(( plugBeingDirtied.partialName() == "t" ) || plugBeingDirtied.partialName() == "en")) {
		//cout << "######### HotOcean need rebuild because of = " << plugBeingDirtied.partialName() << std::endl;
		_ocean_needs_rebuild = true;
	}

	// in geo changed ip[0].ig

	if (plugBeingDirtied.partialName() == "ip[0].ig") {
		_mesh_changed = true;
		_initTangentSpace = true;
	}

	if ( plugBeingDirtied.partialName() == "vertexColor" ) {
		_mesh_changed = true;
	}

	if ( plugBeingDirtied.partialName() == "space" ) {
		_initTangentSpace = true;
	}

	return MS::kSuccess;
}



/*
MStatus hotOceanDeformer::deform( MDataBlock& block,
	MItGeometry& iter,
	const MMatrix& worldSpace,
	unsigned int multiIndex )
{
}
*/



MStatus hotOceanDeformer::compute( const MPlug & plug, MDataBlock & block )
{
#define CHK(msg) if ( MS::kSuccess!=status ) { throw(msg); }

	MStatus			status;

	MFnDependencyNode	this_dnode(thisMObject(), &status);
	MString			dnode_name = this_dnode.name(),
				FN_dnode = dnode_name+": ";

	try
	{
		if (plug.attribute() == outputGeom)
		{
			// get all attributes
			//
			MDataHandle globalScaleData = block.inputValue(globalScale,&status);
			CHK("Error getting globalScale data handle");
			double globalScale = globalScaleData.asDouble();

			MDataHandle resolutionData = block.inputValue(resolution,&status);
			CHK("Error getting resolution data handle");
			int resolution = resolutionData.asInt();
			resolution = (int) pow(2.0,resolution);

			MDataHandle sizeData = block.inputValue(size,&status);
			CHK("Error getting size data handle");
			double size = sizeData.asDouble();

			MDataHandle windSpeedData = block.inputValue(windSpeed,&status);
			CHK("Error getting windSpeed data handle");
			double windSpeed = windSpeedData.asDouble();

			MDataHandle waveHeightData = block.inputValue(waveHeight, &status);
			CHK("Error getting waveHeight data handle");
			double waveHeight = waveHeightData.asDouble();

			MDataHandle shortestWaveData = block.inputValue(shortestWave,&status);
			CHK("Error getting shortestWave data handle");
			double shortestWave = shortestWaveData.asDouble();

			MDataHandle choppinessData = block.inputValue(choppiness,&status);
			CHK("Error getting choppiness data handle");
			double choppiness = choppinessData.asDouble();

			MDataHandle windDirectionData = block.inputValue(windDirection,&status);
			CHK("Error getting windDirection data handle");
			double windDirection = windDirectionData.asDouble();

			MDataHandle dampReflectionsData = block.inputValue(dampReflections,&status);
			CHK("Error getting dampReflection data handle");
			double dampReflections = dampReflectionsData.asDouble();

			MDataHandle windAlignData = block.inputValue(windAlign,&status);
			CHK("Error getting windAlign data handle");
			double windAlign = windAlignData.asDouble();

			MDataHandle oceanDepthData = block.inputValue(oceanDepth,&status);
			CHK("Error getting oceanDepth data handle");
			double oceanDepth = oceanDepthData.asDouble();

			MDataHandle timeData = block.inputValue(time,&status);
			CHK("Error getting time data handle");
			double time = timeData.asDouble();

			MDataHandle seedData = block.inputValue(seed,&status);
			CHK("Error getting seed data handle");
			int seed = seedData.asInt();

			MDataHandle interpolationData = block.inputValue(interpolation,&status);
			CHK("Error getting interpolation data handle");
			bool interpolation = interpolationData.asBool();

			MDataHandle deformSpaceData = block.inputValue(deformSpace, &status );
			CHK("Error getting deformation space data handle");
			int deformSpace = deformSpaceData.asShort();	// Now get it as an SHORT

			MDataHandle vertexColorsData = block.inputValue(vertexColor,&status);
			CHK("Error getting do Vertex Colors data handle");
			bool doVertexColors = vertexColorsData.asBool();


			bool do_jacobian = true || doVertexColors; // for vertex colors


			// determine the envelope (this is a global scale factor)
			//
			MDataHandle envData = block.inputValue(envelope,&status);
			CHK("Error getting envelope data handle");
			float env = envData.asFloat();


			// if we need to (re)initialize the ocean, do this
			//
			if ( !_ocean || _ocean_needs_rebuild )
			{
				if (_ocean)
				{
					delete _ocean;
				}

				if (_ocean_context)
				{
					delete _ocean_context;
				}

				_ocean = new drw::Ocean(resolution, resolution,
					size/float(resolution), size/float(resolution),
					windSpeed, shortestWave, 0.00001f, windDirection/180.0f * M_PI,
					1.0f-dampReflections, windAlign, oceanDepth, seed);

				_ocean_scale   = _ocean->get_height_normalize_factor();
				_ocean_context = _ocean->new_context(true, choppiness>0, false, do_jacobian);

				_ocean_needs_rebuild = false;
			}


			// sum up the waves at this timestep
			//
			_ocean->update( time, *_ocean_context, true, (choppiness>0), false, do_jacobian,
				_ocean_scale * waveHeight, choppiness);

			unsigned int mIndex = plug.logicalIndex();
			MObject thisNode = this->thisMObject();
			MPlug inPlug(thisNode,input);
			inPlug.selectAncestorLogicalIndex(mIndex,input);

			MDataHandle hInput = block.inputValue(inPlug);
			MDataHandle inputGeomDataH = hInput.child(inputGeom);
			MDataHandle hOutput = block.outputValue(plug);
			hOutput.copy(inputGeomDataH);

			MFnMesh inputMesh( inputGeomDataH.asMesh() );

			MMatrix worldSpace = inputGeomDataH.geometryTransformMatrix();
			
			MPointArray verts;
			unsigned int nPoints = inputMesh.numVertices();
			inputMesh.getPoints(verts);


			// some statics to speed things up
			//
			const float envGlobalScale = env * globalScale;
			const float oneOverGlobalScale = 1.0/globalScale;

			const MColor black(0.0, 0.0, 0.0);

			MColorArray
				jMinus = MColorArray(nPoints, black),
				jPlus = MColorArray(nPoints, black),
				eMinus = MColorArray(nPoints, black),
				ePlus = MColorArray(nPoints, black);

			if (doVertexColors)
			{
				//Create one set per mel
				char buffer[512];
				//MString meshName = inputMesh.name();
				//cout << inputMesh.name().asChar() << "  " << std::endl;
				//sprintf_s( buffer, "polyColorSet -create -colorSet \"jMinus1\" %s",inputMesh.name().asChar());
				//status = MGlobal::executeCommand(buffer);
				//if (status != MS::kSuccess) MGlobal::displayError("Error creating colorSet");
				//sprintf_s( buffer, "polyColorPerVertex -rgb 0.5 0.0 0.0 %s",inputMesh.name().asChar());
				//status = MGlobal::executeCommand("polyColorPerVertex -rgb 0.5 0.0 0.0 " + inputMesh.name());
				//if (status != MS::kSuccess) MGlobal::displayError("Error coloring colorSet");

				MString cset_name;
				cset_name = "jMinus";
				status = inputMesh.createColorSetDataMesh(cset_name);
				CHK("Error creating jMinus colorset.");

				//sprintf_s( buffer, "polyColorPerVertex -rgb 0.5 0.0 0.0 %s",this->name().asChar());
				//status = MGlobal::executeCommand(buffer);

				//if (status != MS::kSuccess) MGlobal::displayError("Error coloring colorSet");
				//sprintf_s( buffer, "setAttr \"%s.difs\" 0",inputMesh.name().asChar());
				//cout << inputMesh.name() << std::endl;
				//status = MGlobal::executeCommand("setAttr \"" + meshName + ".difs\" 0");


				cset_name = "jPlus";
				status = inputMesh.createColorSetDataMesh(cset_name);
				CHK("Error creating jPlus colorset.");

				cset_name = "eMinus";
				status = inputMesh.createColorSetDataMesh(cset_name);
				CHK("Error creating eMinus colorset.");

				cset_name = "ePlus";
				status = inputMesh.createColorSetDataMesh(cset_name);
				CHK("Error creating ePlus colorset.");
			}


			if ((_initTangentSpace) && (deformSpace == 2))
			{
				//get the tangents, normas, uvs
				_tangents.setLength(nPoints);
				_normals.setLength(nPoints);
				_binormals.setLength(nPoints);
				_uList.setLength(nPoints);
				_vList.setLength(nPoints);

				MVector vec;
				MIntArray vertexList;

				for (int i=0; i<inputMesh.numPolygons(); i++)
				{
					vertexList.clear();
					inputMesh.getPolygonVertices(i,vertexList);
					for (unsigned int j=0; j<vertexList.length(); j++) {
						inputMesh.getFaceVertexTangent (i, vertexList[j], vec);
						_tangents[vertexList[j]] = vec.normal();
						inputMesh.getFaceVertexNormal (i, vertexList[j], vec);
						_normals[vertexList[j]] = vec.normal();
						vec = _tangents[vertexList[j]]^_normals[vertexList[j]];
						_binormals[vertexList[j]] = vec.normal();
						inputMesh.getPolygonUV(i,j,_uList[vertexList[j]],_vList[vertexList[j]]);
					}
				}

				_initTangentSpace = false;

			} //if (_initTangentSpace)) && (deformSpace == 2)



			if (_mesh_changed && doVertexColors)
			{
				_vertexNumberList.setLength(inputMesh.numFaceVertices());

				MIntArray vertexList;
				int colorIndex;

				//create a list where [vertex per face number] = vertexNumber
				for (int i=0; i<inputMesh.numPolygons(); i++)
				{

					// GET THE VERTEX INDICES FOR CURRENT FACE:
					vertexList.clear();
					inputMesh.getPolygonVertices(i,vertexList);

					// ITERATE THROUGH EACH FACE VERTEX TO SEE IF EACH ONE BELONGS TO THE ORIGINAL VERTEX SELECTION:
					for (unsigned j=0; j<vertexList.length(); j++)
					{
						inputMesh.getFaceVertexColorIndex(i,j,colorIndex);
						_vertexNumberList[colorIndex] = vertexList[j];
						//cout << "vertexId: " << vertexList[j] << " per Face Id " << colorIndex << std::endl;
					}
				}

			} // if (doVertexColors)


			_mesh_changed = false;

			if (deformSpace == 2)
			{
#pragma omp parallel for
				for(int i=0; i<nPoints; ++i)
				{
					drw::EvalData evaldata;

					// do the waves
					//
					if (interpolation)
						_ocean_context->eval2_xz(oneOverGlobalScale *_uList[i], oneOverGlobalScale *_vList[i], evaldata);
					else
						_ocean_context->eval_xz(oneOverGlobalScale *_uList[i], oneOverGlobalScale *_vList[i], evaldata);

					verts[i] += evaldata.disp[0] * envGlobalScale * _tangents[i];
					verts[i] += evaldata.disp[1] * envGlobalScale * _normals[i];
					verts[i] += evaldata.disp[2] * envGlobalScale * _binormals[i];

					if (doVertexColors)
					{
						jMinus.set(i,evaldata.Jminus,evaldata.Jminus,evaldata.Jminus);
						jPlus.set(i,evaldata.Jplus,evaldata.Jplus,evaldata.Jplus);
						eMinus.set(i,evaldata.Eminus[0],evaldata.Eminus[1],evaldata.Eminus[2]);
						ePlus.set(i,evaldata.Eplus[0],evaldata.Eplus[1],evaldata.Eplus[2]);
					}


				} // for

				//for (int i = 0; i < vertexNumberList.length(); i++)
				//cout << "i: " << i << " colorIndex " << vertexNumberList[i] << std::endl;



				//if (setExists == 0)
				//{

				//}

				/*inputMesh.getColorSetNames(existingColorSets);
				for (int i=0; i<existingColorSets.length(); i++)
				{
					cout << "Found Set: " << existingColorSets[i].asChar() << std::endl;
				}*/
				//status = inputMesh.setCurrentColorSetName(setName);
				//CHK("Error switching to colorset.");

				//status = inputMesh.setFaceColors(jMinus, vertexNumberList);
				//cout << "ColorNum " << jMinus.length() << " IndexNum "<< vertexNumberList.length() <<  std::endl;
				//cout << "Current Set " << inputMesh.currentColorSetName().asChar() <<  std::endl;
				//status = inputMesh.setFaceColors(jMinus, vertexNumberList);



			}
			else // object or worldspace
			{
#pragma omp parallel for
				for (int i=0; i<nPoints; ++i)
				{
					drw::EvalData evaldata;
					MPoint pt = verts[i];

					if (deformSpace == 0)
						pt *= worldSpace;

					// do the waves
					//
					if (interpolation)
						_ocean_context->eval2_xz(oneOverGlobalScale *pt.x,oneOverGlobalScale *pt.z, evaldata);
					else
						_ocean_context->eval_xz(oneOverGlobalScale *pt.x,oneOverGlobalScale *pt.z, evaldata);

					pt.x += evaldata.disp[0] * envGlobalScale;
					pt.y += evaldata.disp[1] * envGlobalScale;
					pt.z += evaldata.disp[2] * envGlobalScale;

					if (doVertexColors)
					{
						jMinus.set(i,evaldata.Jminus,evaldata.Jminus,evaldata.Jminus);
						jPlus.set(i,evaldata.Jplus,evaldata.Jplus,evaldata.Jplus);
						eMinus.set(i,evaldata.Eminus[0],evaldata.Eminus[1],evaldata.Eminus[2]);
						ePlus.set(i,evaldata.Eplus[0],evaldata.Eplus[1],evaldata.Eplus[2]);
						//cout << "Eminus " << evaldata.Eminus[0] << std::endl;
					}


					if (deformSpace == 0)
						pt *= worldSpace.inverse();

					verts[i] = pt;

				}
			}

			// write values back onto output using fast set method on iterator
			inputMesh.setPoints(verts);
			/*for (int i = 0; i < _vertexNumberList.length(); i++)
						cout << "i: " << i << " colorIndex " << _vertexNumberList[i] << std::endl;*/

			if (doVertexColors)
			{
				//status = inputMesh.createColorSetDataMesh(setName);
				//CHK("Error creating colorset.");
				//MStringArray existingColorSets;
				//inputMesh.getColorSetNames(existingColorSets);
				////setExists = 0;
				//for (int i=0; i<existingColorSets.length(); i++)
				//{
				//	cout << "Found Set: " << existingColorSets[i].asChar() << std::endl;

				//}
				MString setName;

				setName = MString("jMinus");
				status = inputMesh.clearColors(&setName);
				CHK("Error clearing colors.");
				status = inputMesh.setColors(jMinus,&setName);
				CHK("Error setting colors.");
				status = inputMesh.assignColors(_vertexNumberList,&setName);
				CHK("Error assigning colors.");

				setName = MString("jPlus");
				status = inputMesh.clearColors(&setName);
				CHK("Error clearing colors.");
				status = inputMesh.setColors(jPlus,&setName);
				CHK("Error setting colors.");
				status = inputMesh.assignColors(_vertexNumberList,&setName);
				CHK("Error assigning colors.");

				setName = MString("eMinus");
				status = inputMesh.clearColors(&setName);
				CHK("Error clearing colors.");
				status = inputMesh.setColors(eMinus,&setName);
				CHK("Error setting colors.");
				status = inputMesh.assignColors(_vertexNumberList,&setName);
				CHK("Error assigning colors.");

				setName = MString("ePlus");
				status = inputMesh.clearColors(&setName);
				CHK("Error clearing colors.");
				status = inputMesh.setColors(ePlus,&setName);
				CHK("Error setting colors.");
				status = inputMesh.assignColors(_vertexNumberList,&setName);
				CHK("Error assigning colors.");

				//char buffer[512];
				//sprintf_s( buffer, "polyOptions -colorShadedDisplay true %s",inputMesh.name().asChar());
				//status = MGlobal::executeCommand(buffer);
				//if (status != MS::kSuccess) MGlobal::displayError("Error displaying colors");

				//inputMesh.cleanupEdgeSmoothing();
				//status = inputMesh.updateSurface();
				//CHK("Error redrawing.");
				//status = inputMesh.setDisplayColors(true);
				//CHK("Error displaying colors.");
				//MColor temp = MColor(1.0f,0.0f,0.0f);
				//inputMesh.setFaceColor(temp,0);
				//status = inputMesh.setFaceColors(ePlus,_vertexNumberList);
				//CHK("test.");
			}
		}
	}

	// error handling
	//
	catch( char const *e )
	{
		MGlobal::displayWarning(FN_dnode+MString(e));
	}

	catch( std::runtime_error & e )
	{
		MGlobal::displayWarning(FN_dnode+MString(e.what()));
	}

	catch(...)
	{
		MGlobal::displayError(FN_dnode+"unknown error occurred");
		throw;
	}

	return status;
#undef CHK
}



/**	Registration of deformer node.
*/
MStatus initializePlugin(MObject obj)
{
	MStatus s;
	MFnPlugin plugin(obj, "Nico Rehberg" , "1.1", "Any");

	s = plugin.registerNode("hotOceanDeformer",
		hotOceanDeformer::id, hotOceanDeformer::creator,
		hotOceanDeformer::initialize, MPxNode::kDeformerNode);

	return s;
}



/**	De-registration of deformer node.
*/
MStatus uninitializePlugin(MObject obj)
{
	MStatus s;
	MFnPlugin plugin(obj);
	s = plugin.deregisterNode(hotOceanDeformer::id);
	return s;
}

