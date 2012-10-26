// -*- C++ -*-
/**
		@file		hotOcean.h
		@since		2012-10-26

		@author		Nico Rehberg, Imre Tuske, Szabolcs Horvatth

		@brief		Maya deformer to displace a surface using Houdini Ocean Toolkit (headers).

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

#include "Ocean.h"

#include <string.h>
#include <math.h>
#include <omp.h>

#include <maya/MIOStream.h>

#include <maya/MPxDeformerNode.h>
#include <maya/MItGeometry.h>

#include <maya/MTypeId.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>

#include <maya/MFnNumericAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnPlugin.h>
#include <maya/MFnDependencyNode.h>

#include <maya/MFnMesh.h>

#include <maya/MPoint.h>
#include <maya/MMatrix.h>

#include <maya/MPointArray.h>
#include <maya/MVectorArray.h>
#include <maya/MFloatArray.h>
#include <maya/MIntArray.h>


class hotOceanDeformer : public MPxDeformerNode
{
public:

	hotOceanDeformer();
	virtual			~hotOceanDeformer();

	static void *		creator() { return new hotOceanDeformer(); }
	static MStatus		initialize();

	virtual MStatus		setDependentsDirty( const MPlug &, MPlugArray & );
	virtual MStatus		compute( const MPlug &, MDataBlock & );

/*
	virtual MStatus		deform( MDataBlock &, MItGeometry &, const MMatrix &, unsigned int );
*/

public:
	static MObject globalScale;
	static MObject resolution;
	static MObject size;
	static MObject windSpeed;
	static MObject waveHeight;
	static MObject shortestWave;
	static MObject choppiness;
	static MObject windDirection;
	static MObject dampReflections;
	static MObject windAlign;
	static MObject oceanDepth;
	static MObject time;
	static MObject seed;
	static MObject interpolation;
	static MObject deformSpace;
	static MObject vertexColor;

	static MTypeId id;

protected:
	// This is where all the wave action takes place
	drw::Ocean *		_ocean;
	drw::OceanContext *	_ocean_context;
	float              	_ocean_scale;

	// If this is true we will create a new instance of drw::Ocean
	// next time it runs.
	bool			_ocean_needs_rebuild;
	bool			_mesh_changed;
	bool			_initTangentSpace;

	MVectorArray		_tangents;
	MVectorArray		_normals;
	MVectorArray		_binormals;
	MFloatArray		_uList, _vList;
	MIntArray		_vertexNumberList;
};


