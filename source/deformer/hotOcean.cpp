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
 

 
#include "hotOcean.h" 

MTypeId     hotOceanDeformer::id( 0x0007443a );

// attributes  //
MObject		hotOceanDeformer::globalScale;
MObject     hotOceanDeformer::resolution;             //grid aufloesung
MObject     hotOceanDeformer::size;             // The grid mentiond above is computed for and applied to the input geometry in tiles of this size.
MObject		hotOceanDeformer::windSpeed;     // Wind Speed - Affects the shape of the waves, "Windspeed (m/s)"
MObject		hotOceanDeformer::waveHeigth;
MObject		hotOceanDeformer::shortestWave;    // Shortest Wavelength(m)
MObject		hotOceanDeformer::choppiness;
MObject		hotOceanDeformer::windDirection;    // Wind direction in degrees
MObject		hotOceanDeformer::dampReflections;    // Damp reflections - In a ìfully developedÅEocean you will have waves travelling in both the forward and backwards directions. This parameter damps out the negative direcion waves.
MObject		hotOceanDeformer::windAlign;    // Wind Alignment - Controls how closely the waves travel in the direction of the wind.
MObject		hotOceanDeformer::oceanDepth;  // Ocean Depth - Affects the spectrum of waves generated. Visually in doesnít seem to have that great an influence.
MObject		hotOceanDeformer::time;
MObject		hotOceanDeformer::seed;   //Seed - Seeds the random number generator. 
MObject		hotOceanDeformer::interpolation;   
MObject		hotOceanDeformer::deformSpace;
MObject		hotOceanDeformer::vertexColor;

hotOceanDeformer::hotOceanDeformer() 
//
//      Description:
//              constructor
//
{ 
	_ocean = 0;
    _ocean_context = 0;
    _ocean_scale = 1.0f;

    _ocean_needs_rebuild = true;
	_mesh_changed = true;
	_initTangentSpace = true;
}

hotOceanDeformer::~hotOceanDeformer()
//
//      Description:
//              destructor
//
{
    if (_ocean)
    {
        delete _ocean;
    }
    if (_ocean_context)
    {
        delete _ocean_context;
    }
}

void* hotOceanDeformer::creator()
//
//      Description:
//              create the yTwist
//
{
        return new hotOceanDeformer();
}

MStatus hotOceanDeformer::initialize()
//
//      Description:
//              initialize the attributes
//
{
        // local attribute initialization
        //
	    MFnNumericAttribute nAttr;
        globalScale=nAttr.create( "globalScale", "scale", MFnNumericData::kDouble  );
            nAttr.setDefault(1.0); nAttr.setKeyable(false); addAttribute( globalScale ); 
			nAttr.setMin(0.00001); nAttr.setSoftMax(10); nAttr.setChannelBox(true);

        resolution=nAttr.create( "resolution", "res", MFnNumericData::kInt  );
            nAttr.setDefault(6); nAttr.setKeyable(false); addAttribute( resolution ); 
			nAttr.setMin(4); nAttr.setMax(11); nAttr.setChannelBox(true);

		size=nAttr.create( "size", "size", MFnNumericData::kDouble  );
            nAttr.setDefault(200.0); nAttr.setKeyable(false); addAttribute( size ); 
			nAttr.setMin(0.01);  nAttr.setSoftMax(1000); nAttr.setChannelBox(true);

        windSpeed = nAttr.create( "windSpeed", "ws", MFnNumericData::kDouble );
            nAttr.setDefault(15.0); nAttr.setKeyable(true); addAttribute( windSpeed );
			nAttr.setMin(0.01); nAttr.setSoftMax(100); 

		waveHeigth = nAttr.create( "waveHeigth", "wh", MFnNumericData::kDouble );
            nAttr.setDefault(2.0); nAttr.setKeyable(true); addAttribute( waveHeigth );
			nAttr.setMin(0.01); nAttr.setSoftMax(20); 

		shortestWave = nAttr.create( "shortestWave", "sw", MFnNumericData::kDouble );
            nAttr.setDefault(0.01); nAttr.setKeyable(true); addAttribute( shortestWave );
			nAttr.setMin(0.00001);  nAttr.setSoftMax(20); 

		choppiness = nAttr.create( "choppiness", "chop", MFnNumericData::kDouble );
            nAttr.setDefault(0.7); nAttr.setKeyable(true); addAttribute( choppiness );
			nAttr.setMin(0.0); nAttr.setMax(1.0);

		windDirection = nAttr.create( "windDirection", "dir", MFnNumericData::kDouble );
            nAttr.setDefault(0); nAttr.setKeyable(true); addAttribute( windDirection );
			nAttr.setMin(0.0); nAttr.setSoftMax(360.0); 

		dampReflections = nAttr.create( "dampReflections", "damp", MFnNumericData::kDouble );
            nAttr.setDefault(0.5); nAttr.setKeyable(true); addAttribute( dampReflections );
			nAttr.setMin(0.0); nAttr.setMax(1.0);
 
		windAlign = nAttr.create( "windAlign", "wa", MFnNumericData::kDouble );
            nAttr.setDefault(2.0); nAttr.setKeyable(true); addAttribute( windAlign );
			nAttr.setMin(0.01); nAttr.setSoftMax(10.0);
  
		oceanDepth = nAttr.create( "oceanDepth", "dpt", MFnNumericData::kDouble );
            nAttr.setDefault(200.0); nAttr.setKeyable(true); addAttribute( oceanDepth );
			nAttr.setMin(0.01); nAttr.setSoftMax(1000.0);

		time = nAttr.create( "time", "t", MFnNumericData::kDouble );
            nAttr.setDefault(0.0); nAttr.setKeyable(true); addAttribute( time );
			nAttr.setSoftMin(0.0); nAttr.setSoftMax(1000.0);

		seed = nAttr.create( "seed", "s", MFnNumericData::kInt );
            nAttr.setDefault(0); nAttr.setKeyable(true); addAttribute( seed );
			nAttr.setMin(0); nAttr.setSoftMax(9999);
		
		interpolation = nAttr.create( "interpolation", "int", MFnNumericData::kBoolean );
            nAttr.setDefault(false); nAttr.setKeyable(false); addAttribute( interpolation ); nAttr.setChannelBox(true);
		
		MFnEnumAttribute eAttr;		
		deformSpace = eAttr.create( "space", "space", 0 );  // can set default 0, 1, 2, etc..
		eAttr.setStorable(true); eAttr.setKeyable(false); addAttribute(deformSpace); eAttr.setChannelBox(true);
		eAttr.addField("world", 0); eAttr.addField("object", 1); eAttr.addField("tangent", 2) ;

		vertexColor = nAttr.create( "doVertexColor", "vertexColor", MFnNumericData::kBoolean ); 
			nAttr.setDefault(false); nAttr.setKeyable(false); addAttribute( vertexColor ); nAttr.setChannelBox(true);


	

        // affects
        //
		attributeAffects( hotOceanDeformer::globalScale, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::resolution, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::size, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::windSpeed, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::waveHeigth, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::shortestWave, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::choppiness, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::windDirection, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::dampReflections, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::windAlign, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::oceanDepth, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::time, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::seed, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::interpolation, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::deformSpace, hotOceanDeformer::outputGeom );
		attributeAffects( hotOceanDeformer::vertexColor, hotOceanDeformer::outputGeom );

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
		if ( plugBeingDirtied.partialName() == "vertexColor" ){
			_mesh_changed = true;
		}
		if ( plugBeingDirtied.partialName() == "space" ){
			_initTangentSpace = true;
		}
        return( MS::kSuccess );
    }


//MStatus
//hotOceanDeformer::deform( MDataBlock& block,
//                            MItGeometry& iter,
//                            const MMatrix& worldSpace,
//                            unsigned int multiIndex)

MStatus hotOceanDeformer::compute( const MPlug& plug, MDataBlock& block )

//
// Method: deform
//
// 
//
// Arguments:
//   block              : the datablock of the node
//       iter           : an iterator for the geometry to be deformed
//   worldSpace                  : matrix to transform the point into world space
//       multiIndex : the index of the geometry that we are deforming
//
//
{
		
	MStatus status = MS::kSuccess;
	if (plug.attribute() == outputGeom) {
        
        // get all attributes
        //
		MDataHandle globalScaleData = block.inputValue(globalScale,&status);
        McheckErr(status, "Error getting globalScale data handle\n");
        double globalScale = globalScaleData.asDouble();

        MDataHandle resolutionData = block.inputValue(resolution,&status);
        McheckErr(status, "Error getting resolution data handle\n");
        int resolution = resolutionData.asInt();
		resolution = (int) pow(2.0,resolution);

		MDataHandle sizeData = block.inputValue(size,&status);
        McheckErr(status, "Error getting size data handle\n");
        double size = sizeData.asDouble();

		MDataHandle windSpeedData = block.inputValue(windSpeed,&status);
        McheckErr(status, "Error getting windSpeed data handle\n");
        double windSpeed = windSpeedData.asDouble();

		MDataHandle waveHeigthData = block.inputValue(waveHeigth,&status);
        McheckErr(status, "Error getting waveHeigth data handle\n");
        double waveHeigth = waveHeigthData.asDouble();

		MDataHandle shortestWaveData = block.inputValue(shortestWave,&status);
        McheckErr(status, "Error getting shortestWave data handle\n");
        double shortestWave = shortestWaveData.asDouble();

		MDataHandle choppinessData = block.inputValue(choppiness,&status);
        McheckErr(status, "Error getting choppiness data handle\n");
        double choppiness = choppinessData.asDouble();

		MDataHandle windDirectionData = block.inputValue(windDirection,&status);
        McheckErr(status, "Error getting windDirection data handle\n");
        double windDirection = windDirectionData.asDouble();

		MDataHandle dampReflectionsData = block.inputValue(dampReflections,&status);
        McheckErr(status, "Error getting dampReflection data handle\n");
        double dampReflections = dampReflectionsData.asDouble();

		MDataHandle windAlignData = block.inputValue(windAlign,&status);
        McheckErr(status, "Error getting windAlign data handle\n");
        double windAlign = windAlignData.asDouble();

		MDataHandle oceanDepthData = block.inputValue(oceanDepth,&status);
        McheckErr(status, "Error getting oceanDepth data handle\n");
        double oceanDepth = oceanDepthData.asDouble();

		MDataHandle timeData = block.inputValue(time,&status);
        McheckErr(status, "Error getting time data handle\n");
        double time = timeData.asDouble();

		MDataHandle seedData = block.inputValue(seed,&status);
        McheckErr(status, "Error getting seed data handle\n");
        int seed = seedData.asInt();

		MDataHandle interpolationData = block.inputValue(interpolation,&status);
        McheckErr(status, "Error getting interpolation data handle\n");
        bool interpolation = interpolationData.asBool();

		MDataHandle deformSpaceData = block.inputValue(deformSpace, &status );
		McheckErr(status, "Error getting deformation space data handle\n");
		int deformSpace = deformSpaceData.asShort();	// Now get it as an SHORT

		MDataHandle vertexColorsData = block.inputValue(vertexColor,&status);
        McheckErr(status, "Error getting do Vertex Colors data handle\n");
        bool doVertexColors = false; //vertexColorsData.asBool();


        // determine the envelope (this is a global scale factor)
        //
        MDataHandle envData = block.inputValue(envelope,&status);
        McheckErr(status, "Error getting envelope data handle\n");      
        float env = envData.asFloat();  
	
        // if we need to (re)initialize the ocean, do this
		if (!_ocean || _ocean_needs_rebuild)
        {
            if (_ocean)
            {
                delete _ocean;
            }

            if (_ocean_context)
            {
                delete _ocean_context;
            }

            _ocean = new drw::Ocean(resolution,resolution,size/float(resolution),size/float(resolution),
                          windSpeed,shortestWave,0.00001f, windDirection/180.0f * M_PI,
                          1.0f-dampReflections,windAlign,oceanDepth,seed);
            _ocean_scale   = _ocean->get_height_normalize_factor();

            _ocean_context = _ocean->new_context(true,(choppiness>0),false,true);

            _ocean_needs_rebuild = false;
           // cout << "######### HotOcean, rebuilt ocean, norm_factor = " << _ocean_scale   
			//<<	"resolution = " << resolution
			//	 <<	"windSpeed = " << windSpeed
			//<< std::endl;
        }
		
		// sum up the waves at this timestep
        _ocean->update( time,*_ocean_context,true,(choppiness>0),false,true,
                       _ocean_scale * waveHeigth,choppiness);
		
		unsigned int mIndex = plug.logicalIndex();
		MObject thisNode = this->thisMObject();
		MPlug inPlug(thisNode,input);
		inPlug.selectAncestorLogicalIndex(mIndex,input);
		MDataHandle hInput = block.inputValue(inPlug);
		MDataHandle inputGeomDataH = hInput.child(inputGeom);
		MDataHandle hOutput = block.outputValue(plug);
		hOutput.copy(inputGeomDataH);
		MFnMesh inputMesh = inputGeomDataH.asMesh();
		MMatrix worldSpace;
		worldSpace = inputGeomDataH.geometryTransformMatrix();
		MPointArray verts;
		unsigned int nPoints = inputMesh.numVertices();
		inputMesh.getPoints(verts);
		
		
		// some statics to speed things up

		const float envGlobalScale = env * globalScale;
		const float oneOverGlobalScale = 1.0/globalScale;
		const MColor black = MColor(0.0, 0.0, 0.0);
		MColorArray jMinus = MColorArray(nPoints, black);
		MColorArray jPlus = MColorArray(nPoints, black);
		MColorArray eMinus = MColorArray(nPoints, black);
		MColorArray ePlus = MColorArray(nPoints, black);

		//cout << "thisNode " << this->name() << std::endl;
	

		//bool setExists = 0;
		//MStringArray existingColorSets;
		//inputMesh.getColorSetNames(existingColorSets);
		//for (int i=0; i<existingColorSets.length(); i++)
		//{
		//	cout << "Found Set: " << existingColorSets[i].asChar() << std::endl;
		//	setExists = 1;
		//}


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

			cout << "Creating Sets"  << std::endl;
			MString tmp;
			tmp = "jMinus"; status = inputMesh.createColorSetDataMesh(tmp);
			McheckErr(status, "Error creating colorset.\n");
			
			//sprintf_s( buffer, "polyColorPerVertex -rgb 0.5 0.0 0.0 %s",this->name().asChar());
			//status = MGlobal::executeCommand(buffer);
			
		    //if (status != MS::kSuccess) MGlobal::displayError("Error coloring colorSet");
			//sprintf_s( buffer, "setAttr \"%s.difs\" 0",inputMesh.name().asChar());
			cout << inputMesh.name() << std::endl;
			//status = MGlobal::executeCommand("setAttr \"" + meshName + ".difs\" 0");
			

			tmp = "jPlus"; status = inputMesh.createColorSetDataMesh(tmp);
			McheckErr(status, "Error creating colorset.\n");

			tmp = "eMinus"; status = inputMesh.createColorSetDataMesh(tmp);
			McheckErr(status, "Error creating colorset.\n");

			tmp = "ePlus"; status = inputMesh.createColorSetDataMesh(tmp);
			McheckErr(status, "Error creating colorset.\n");
		}
			


		if ((_initTangentSpace) && (deformSpace == 2)){

				//get the tangents, normas, uvs
				_tangents.setLength(nPoints);
				_normals.setLength(nPoints);
				_binormals.setLength(nPoints);
				_uList.setLength(nPoints);
				_vList.setLength(nPoints);

				MVector vec;
				MIntArray vertexList;

		
				for(int i=0;i<inputMesh.numPolygons();i++)
				{
						vertexList.clear();					
						inputMesh.getPolygonVertices(i,vertexList);
						for(unsigned int j=0; j<vertexList.length();j++) {
							inputMesh.getFaceVertexTangent (i, vertexList[j], vec);
							_tangents[vertexList[j]] = vec.normal();
							inputMesh.getFaceVertexNormal (i, vertexList[j], vec);
							_normals[vertexList[j]] = vec.normal();
							vec = _tangents[vertexList[j]]^_normals[vertexList[j]];
							_binormals[vertexList[j]] = vec.normal();
							inputMesh.getPolygonUV(i,j,_uList[vertexList[j]],_vList[vertexList[j]]);
						}
				}

				_initTangentSpace == false;

		} //if (_initTangentSpace)) && (deformSpace == 2)

		if ((_mesh_changed) && (doVertexColors))	
		{
				
			_vertexNumberList.setLength(inputMesh.numFaceVertices());

			cout << "Listlegth "  << inputMesh.numFaceVertices() <<std::endl;
			
			MIntArray vertexList;
			int colorIndex;

			//create a list where [vertex per face number] = vertexNumber
			for(int i=0;i<inputMesh.numPolygons();i++)
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

		}//if (doVertexColors)
		_mesh_changed = false;
		
		if (deformSpace == 2)
		{
			#pragma omp parallel for
			for(int i=0; i<nPoints; i++) 
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
				

			}//for
	
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
				//McheckErr(status, "Error switching to colorset.\n");

				//status = inputMesh.setFaceColors(jMinus, vertexNumberList);
				//cout << "ColorNum " << jMinus.length() << " IndexNum "<< vertexNumberList.length() <<  std::endl;
				//cout << "Current Set " << inputMesh.currentColorSetName().asChar() <<  std::endl; 				
				//status = inputMesh.setFaceColors(jMinus, vertexNumberList);


		
		}
		else //object or worldspace
		{
			#pragma omp parallel for
 			for(int i=0; i<nPoints; i++) {
	                
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
			//McheckErr(status, "Error creating colorset.\n");
				//MStringArray existingColorSets;
				//inputMesh.getColorSetNames(existingColorSets);
				////setExists = 0;
				//for (int i=0; i<existingColorSets.length(); i++)
				//{
				//	cout << "Found Set: " << existingColorSets[i].asChar() << std::endl;

				//}
			MString setName = MString("jMinus");
			status = inputMesh.clearColors(&setName);
			McheckErr(status, "Error clearing colors.\n");
			status = inputMesh.setColors(jMinus,&setName);
			McheckErr(status, "Error setting colors.\n");
			status = inputMesh.assignColors(_vertexNumberList,&setName);
			McheckErr(status, "Error assigning colors.\n");
			setName = MString("jPlus");
			status = inputMesh.clearColors(&setName);
			McheckErr(status, "Error clearing colors.\n");
			status = inputMesh.setColors(jPlus,&setName);
			McheckErr(status, "Error setting colors.\n");
			status = inputMesh.assignColors(_vertexNumberList,&setName);
			McheckErr(status, "Error assigning colors.\n");
			setName = MString("eMinus");
			status = inputMesh.clearColors(&setName);
			McheckErr(status, "Error clearing colors.\n");
			status = inputMesh.setColors(eMinus,&setName);
			McheckErr(status, "Error setting colors.\n");
			status = inputMesh.assignColors(_vertexNumberList,&setName);
			McheckErr(status, "Error assigning colors.\n");
			setName = MString("ePlus");
			status = inputMesh.clearColors(&setName);
			McheckErr(status, "Error clearing colors.\n");
			status = inputMesh.setColors(ePlus,&setName);
			McheckErr(status, "Error setting colors.\n");
			status = inputMesh.assignColors(_vertexNumberList,&setName);
			McheckErr(status, "Error assigning colors.\n");

			//char buffer[512];
			//sprintf_s( buffer, "polyOptions -colorShadedDisplay true %s",inputMesh.name().asChar());
			//status = MGlobal::executeCommand(buffer);
			//if (status != MS::kSuccess) MGlobal::displayError("Error displaying colors");


    
			//inputMesh.cleanupEdgeSmoothing();
			//status = inputMesh.updateSurface();
			//McheckErr(status, "Error redrawing.\n");
			//status = inputMesh.setDisplayColors(true);
			//McheckErr(status, "Error displaying colors.\n");
			//MColor temp = MColor(1.0f,0.0f,0.0f);
			//inputMesh.setFaceColor(temp,0);
			//status = inputMesh.setFaceColors(ePlus,_vertexNumberList);
			//McheckErr(status, "test.\n");
		}
	}

		
    return status;
}

// standard initialization procedures
//

MStatus initializePlugin( MObject obj )
{
        MStatus result;
        MFnPlugin plugin( obj, "Nico Rehberg" , "1.1", "Any");
        result = plugin.registerNode( "hotOceanDeformer", hotOceanDeformer::id, hotOceanDeformer::creator, 
                                                                  hotOceanDeformer::initialize, MPxNode::kDeformerNode );
        return result;
}

MStatus uninitializePlugin( MObject obj)
{
        MStatus result;
        MFnPlugin plugin( obj );
        result = plugin.deregisterNode( hotOceanDeformer::id );
        return result;
}
 
