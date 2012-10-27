// -*- C++ -*-
/**
		@file		hotOcean.cpp
		@since		2012-10-26

		@author		Nico Rehberg, Imre Tuske

		@brief		Maya deformer to displace a surface using Houdini Ocean Toolkit (implementation).

		Implementation: Nico Rehberg <mail@nico-rehberg.de>

		Thanks to Szabolcs Horvatth and John Patrick for tips and
		hints on colorsets-in-deformers hackery.

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

#include <maya/MFnNumericAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MItMeshPolygon.h>


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
	MFnNumericAttribute	nAttr;
	MFnEnumAttribute	eAttr;
	MFnUnitAttribute	uAttr;

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

	/*
	time = uAttr.create( "time", "t", MFnUnitAttribute::kTime );
	uAttr.setKeyable(true);
	*/
	time = nAttr.create( "time", "t", MFnNumericData::kDouble );
	nAttr.setKeyable(true);
	addAttribute(time);

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



void hotOceanDeformer::postConstructor()
{
	setDeformationDetails(kDeformsColors);
}



MStatus hotOceanDeformer::accessoryNodeSetup( MDagModifier & cmd )
{
	MStatus			s;
	MFnDependencyNode	this_dnode(thisMObject(), &s);
	
	s=cmd.commandToExecute("expression -s \""+this_dnode.name()+".time=time\";");
	return s;
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

			MDataHandle d;
			d = block.inputValue(doJMinus, &status); bool do_jminus = d.asBool() && doVertexColors;
			d = block.inputValue(doJPlus,  &status); bool do_jplus  = d.asBool() && doVertexColors;
			d = block.inputValue(doEMinus, &status); bool do_eminus = d.asBool() && doVertexColors;
			d = block.inputValue(doEPlus,  &status); bool do_eplus  = d.asBool() && doVertexColors;

			doVertexColors = doVertexColors && (do_jminus || do_jplus || do_eminus || do_eplus);

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
			_ocean->update( time, *_ocean_context, true, choppiness>0, false, do_jacobian,
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

			// jacobian arrays
			//
			MFloatArray jMinus, jPlus;
			MColorArray eMinus, ePlus;

			if (do_jminus) jMinus.setLength(nPoints);
			if (do_jplus)  jPlus.setLength(nPoints);
			if (do_eminus) eMinus.setLength(nPoints);
			if (do_eplus)  ePlus.setLength(nPoints);


			_mesh_changed = false;

			if (deformSpace == 2)
			{
				if ( _initTangentSpace )
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
				}


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

					if (doVertexColors) {
						if (do_jminus) jMinus.set(evaldata.Jminus, i);
						if (do_jplus)  jPlus.set(evaldata.Jplus, i);
						if (do_eminus) eMinus.set(i, evaldata.Eminus[0], evaldata.Eminus[1], evaldata.Eminus[2]);
						if (do_eplus)  ePlus.set( i, evaldata.Eplus[0],  evaldata.Eplus[1],  evaldata.Eplus[2]);
					}
				}
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

					if (doVertexColors) {
						if (do_jminus) jMinus.set(evaldata.Jminus, i);
						if (do_jplus)  jPlus.set(evaldata.Jplus, i);
						if (do_eminus) eMinus.set(i, evaldata.Eminus[0], evaldata.Eminus[1], evaldata.Eminus[2]);
						if (do_eplus)  ePlus.set( i, evaldata.Eplus[0],  evaldata.Eplus[1],  evaldata.Eplus[2]);
					}

					if (deformSpace == 0)
						pt *= worldSpace.inverse();

					verts[i] = pt;

				}
			}

			// write values back onto output using fast set method on iterator
			//
			inputMesh.setPoints(verts);


			// write vertex colors (per-face vertex method)
			//
			if (doVertexColors)
			{
				MColorArray	jm, jp, em, ep;
				MIntArray	indices;
				unsigned	nth=0;

				MItMeshPolygon pit(inputGeomDataH.asMesh(), &status);
				CHK("couldn't get mesh poly iterator");

				for( ;  !pit.isDone();  pit.next() )
				{
					for( unsigned v=0, vc=pit.polygonVertexCount();  v<vc;  ++v, ++nth )
					{
						unsigned int p = pit.vertexIndex(v);

						if (do_jminus) jm.append(jMinus[p], jMinus[p], jMinus[p]);
						if (do_jplus)  jp.append(jPlus[p],  jPlus[p],  jPlus[p]);
						if (do_eminus) em.append(eMinus[p]);
						if (do_eplus)  ep.append(ePlus[p]);

						indices.append(nth);
					}
				}

				MString cset;

				if (do_jminus) {
					cset = "jMinus";
					status = inputMesh.createColorSetDataMesh(cset);
					inputMesh.setColors(jm, &cset);
					inputMesh.assignColors(indices, &cset);
				}
				
				if (do_jplus) {
					cset = "jPlus";
					status = inputMesh.createColorSetDataMesh(cset);
					inputMesh.setColors(jp, &cset);
					inputMesh.assignColors(indices, &cset);
				}
				
				if (do_eminus) {
					cset = "eMinus";
					status = inputMesh.createColorSetDataMesh(cset);
					inputMesh.setColors(em, &cset);
					inputMesh.assignColors(indices, &cset);
				}
				
				if (do_eplus) {
					cset = "ePlus";
					status = inputMesh.createColorSetDataMesh(cset);
					inputMesh.setColors(ep, &cset);
					inputMesh.assignColors(indices, &cset);
				}
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

