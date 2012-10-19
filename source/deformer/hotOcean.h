// -*- C++ -*-
/***************************************************************************
 * hotOcean.cpp      a maya deformer to displace a surface using
 * the Houdini Ocean Toolkit code
 *
 * implementation Nico Rehberg <mail@nico-rehberg.de>
 *
 *
 * The Houdini Ocean Toolkit is copyrighted by Drew Whitehouse
 * see http://odforce.net/wiki/index.php/HoudiniOceanToolkit
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 *
 ***************************************************************************/


#include <string.h>
#include <maya/MIOStream.h>
#include <math.h>

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

#include "Ocean.h"

#include <omp.h>
#include <maya/MPointArray.h>
#include <maya/MVectorArray.h>
#include <maya/MFloatArray.h>
#include <maya/MIntArray.h>


class hotOceanDeformer : public MPxDeformerNode
{
public:

	hotOceanDeformer();
	virtual			~hotOceanDeformer();

	static void *		creator();
	static MStatus		initialize();

	virtual MStatus		setDependentsDirty( const MPlug &, MPlugArray & );
	virtual MStatus		compute( const MPlug &, MDataBlock & );

/*
	virtual MStatus		deform( MDataBlock &, MItGeometry &, const MMatrix &, unsigned int );
*/

public:
	// attributes
	//
	static  MObject globalScale;
	static  MObject resolution;		//grid aufloesung
	static  MObject size;			// The grid mentiond above is computed for and applied to the input geometry in tiles of this size.
	static  MObject windSpeed;		// Wind Speed - Affects the shape of the waves, "Windspeed (m/s)"
	static  MObject waveHeigth;
	static  MObject shortestWave;		// Shortest Wavelength(m)
	static  MObject choppiness;
	static  MObject windDirection;		// Wind direction in degrees
	static  MObject dampReflections;	// Damp reflections - In a 'fully developed' ocean you will have waves travelling in both the forward and backwards directions. This parameter damps out the negative direcion waves.
	static  MObject windAlign;		// Wind Alignment - Controls how closely the waves travel in the direction of the wind.
	static  MObject oceanDepth;		// Ocean Depth - Affects the spectrum of waves generated. Visually in doesnt seem to have that great an influence.
	static  MObject time;
	static  MObject seed;			// Seed - Seeds the random number generator.
	static  MObject interpolation;		// interpolation: linear or smooth
	static  MObject deformSpace;
	static  MObject vertexColor;
	static  MTypeId id;

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


