// -*- C++ -*-
/***************************************************************************
 * hotOcean.cpp      a mental ray surface/displacement shader to use the
 * Houdini Ocean Toolkit code inside mental ray
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

   
  
#include "shader.h"
#include "math.h" 
#include "stdlib.h"
#include "mayaapi.h"

// this is where all the magic happens
#include "Ocean.h"


extern "C" DLLEXPORT int hotOcean_version(void) {return 1;}

struct hotOcean{
	miScalar globalScale;			
	// todo int textureSpace; 	// 0= worldspace, 1 = uv Space
	int resolution;			//grid aufloesung
	miScalar size;			// The grid mentiond above is computed for and applied to the input geometry in tiles of this size.
	miScalar windSpeed;     // Wind Speed - Affects the shape of the waves, "Windspeed (m/s)"
	miScalar waveHeigth;
	miScalar shortestWave;    // Shortest Wavelength(m)
	miScalar choppiness;
	miScalar windDirection;    // Wind direction in degrees
	miScalar dampReflections;    // Damp reflections - In a “fully developedEocean you will have waves travelling in both the forward and backwards directions. This parameter damps out the negative direcion waves.
	miScalar windAlign;    // Wind Alignment - Controls how closely the waves travel in the direction of the wind.
	miScalar oceanDepth;  // Ocean Depth - Affects the spectrum of waves generated. Visually in doesn’t seem to have that great an influence.
	miScalar time;
	int   seed;   //Seed - Seeds the random number generator.
	int   textureSpace;  
	miVector inDisplacement; //an input for layering multiple displacements
	
};  

struct hotOceanOutputs{
	    miScalar		displacement;
		miVector		displacementVector;
		miVector		normal;
		miScalar		jminus;
		miScalar		jplus;
		miVector		eminus;
		miVector		eplus;

};

struct pointSample
{
	double x;
	double z;
	int samples;
};

struct OceanHolder  
{
	float              normalize_factor;
    drw::Ocean        *ocean;
    drw::OceanContext *context;
	bool				empty;
	bool				doBuffer;
	double				stepSize;
	pointSample			sampleSequence[];

}; 


//////

void* miaux_user_memory_pointer(miState *state, int allocation_size) 
{ 
    void **user_pointer; 
    mi_query(miQ_FUNC_USERPTR, state, 0, &user_pointer); 
    if (allocation_size > 0) { 
		*user_pointer = mi_mem_allocate(allocation_size); 
    } 	
    return *user_pointer; 
} 

miBoolean miaux_release_user_memory(char* shader_name, miState *state, void *params) 
{ 
    if (params != NULL) {  /* Shader instance exit */ 
		
		void **user_pointer;                            
		if (!mi_query(miQ_FUNC_USERPTR, state, 0, &user_pointer)) 
			mi_fatal("Could not get user pointer in shader exit function %s_exit",  
				 shader_name); 	
		
		//OceanHolder *oh = reinterpret_cast<OceanHolder*>(*user_pointer);
		//delete (oh->context); 
		//delete (oh->ocean); 
		mi_mem_release(*user_pointer);
    } 
    return miTRUE; 
} 

extern "C" DLLEXPORT miBoolean hotOcean_init(
					miState *state, struct hotOcean *params, miBoolean *instance_init_required)
{
	if (params == NULL) //Haupt Shader Init
		*instance_init_required = miTRUE;
	else // Init fuer die aktuelle Shader Instanz
	{
			int res		= *mi_eval_integer(&params->resolution);
			res = (int)pow (2.0, res);
			float size	= *mi_eval_scalar(&params->size);
			float V		= *mi_eval_scalar(&params->windSpeed);
			float l		= *mi_eval_scalar(&params->shortestWave);
			float  height_scale   = *mi_eval_scalar(&params->waveHeigth);
			float  chop_amount    =  *mi_eval_scalar(&params->choppiness);
			float w		= *mi_eval_scalar(&params->windDirection);
			float damp	= *mi_eval_scalar(&params->dampReflections);
			float align	= *mi_eval_scalar(&params->windAlign);
			float depth	= *mi_eval_scalar(&params->oceanDepth);
			int   seed	= *mi_eval_integer(&params->seed);  
			//int textureSpace		= *mi_eval_integer(&params->textureSpace);

			bool    do_chop        =  false;			
			if (chop_amount > 0.0001 )
				do_chop    =  true;			
			bool    do_jacobian    = true;    
			bool    do_normal      = true;
			
			int do_buffer = 0;
			if (do_chop)
				do_buffer = 1;
			
			void *uData = miaux_user_memory_pointer(state,sizeof(OceanHolder) + do_buffer * res * res * sizeof(pointSample));
			OceanHolder *oh = reinterpret_cast<OceanHolder*>(uData);

			miScalar  now            = *mi_eval_scalar(&params->time);
			mi_info("Creating Ocean with resolution %i for time %f",res,now);
			if (do_buffer)
				mi_info("Using point buffer size %i kB",do_buffer * res * res * sizeof(pointSample)/1024);

			oh->ocean = new drw::Ocean(res,res,size/float(res),size/float(res),
                          V,l,0.000001f, w/180.0f * M_PI,  
                          1.0f-damp,align,depth,seed);   
			oh->normalize_factor = oh->ocean->get_height_normalize_factor();
			oh->context = oh->ocean->new_context(true,do_chop,do_normal,do_jacobian);     
			oh->ocean->update (now,*oh->context, true,do_chop,do_normal,do_jacobian,
							   height_scale * oh->normalize_factor,chop_amount);
			oh->empty = false; 
			oh->doBuffer = (bool) do_buffer;
			miScalar globalScale = *mi_eval_scalar(&params->globalScale);
			oh->stepSize = size * globalScale/res; 
	
			if ((state->type == miRAY_DISPLACE) && (oh->doBuffer))
			{
				// init the Buffer as empty
				pointSample *dispBuffer = (pointSample*) oh->sampleSequence;
				for (int x = 0; x < res*res; x++)
					dispBuffer[x].samples = 0; 
			}

			/*
			mi_info(" Time %f", now);
			mi_info(" res %f", res);
			mi_info(" size/float(res) %f", size/float(res));
			mi_info(" V %f", V);
			mi_info(" l %f", l);
			mi_info(" w/180.0f * M_PI %f", w/180.0f * M_PI);
			mi_info(" 1.0f-damp %f", 1.0f-damp);
			mi_info(" align %f", align);
			mi_info(" depth %f", depth);
			mi_info(" seed %f", seed); 

			oh->context->eval2_xz(10.0f,10.0f);
			mi_info(" Test 10 10 %f %f %f",oh->context->disp[0], oh->context->disp[1], oh->context->disp[2]);

			oh->context->eval2_xz(-10.0f,10.0f);
			mi_info(" Test -10 10 %f %f %f",oh->context->disp[0], oh->context->disp[1], oh->context->disp[2]);
			*/

		}

	return miTRUE;
}//init
 
extern "C" DLLEXPORT miBoolean hotOcean_exit(
					miState *state, struct hotOcean *params)
{


	if (params != NULL)	{  /* Shader instance exit */ 
		
		void **user_pointer;                            
		if (!mi_query(miQ_FUNC_USERPTR, state, 0, &user_pointer)) 
			mi_fatal("Could not get user pointer in shader exit function hotOcean_exit"); 	
		
		OceanHolder *oh = reinterpret_cast<OceanHolder*>(*user_pointer);
		delete (oh->context); 
		delete (oh->ocean); 

		mi_mem_release(*user_pointer);
    } 
    return miTRUE;	

	//void *uData = miaux_user_memory_pointer(state, 0);
	//OceanHolder *oh = reinterpret_cast<OceanHolder*>(uData);	
	//return miaux_release_user_memory("hotOcean", state, params);
 
}//exit

pointSample * getPointBuffer(miState *state, int bufferX, int bufferZ, pointSample * dispBuffer, int res)
{	
	//mi_info("%i %i",bufferX, bufferZ);
	if (bufferX < 0)
		bufferX += res;
	if (bufferZ < 0)
		bufferZ += res;
	bufferX = bufferX % res;
	bufferZ = bufferZ % res;
	
	pointSample * result = &dispBuffer[bufferX * res + bufferZ];
	// TODO: if this sample is empty interpolate it from neighbours
	if (result->samples == 0)
		{
			miLock *lock;
			mi_query(miQ_FUNC_LOCK, state, miNULLTAG, &lock);
			mi_lock(*lock);

			//mi_info("No Samples for %i %i", bufferX, bufferZ);
			result->x = 0;
			result->z = 0;

			mi_unlock(*lock);

		} 

	return result;
}

extern "C" DLLEXPORT miBoolean hotOcean(
					struct hotOceanOutputs *result, miState *state, struct hotOcean *params)
{	
	 
//mi_info("tangent1 %f %f", state->derivs[0],state->derivs[1]);
		/* check for illegal calls */
        if (state->type == miRAY_SHADOW) {
			return(miFALSE); }
		int textureSpace		= *mi_eval_integer(&params->textureSpace);
		miScalar globalScale = *mi_eval_scalar(&params->globalScale);
		miVector inDisplacement = *mi_eval_vector(&params->inDisplacement);
		void *uData = miaux_user_memory_pointer(state, 0);		
		OceanHolder *oh = reinterpret_cast<OceanHolder*>(uData);
		float  x,z;
		
		int res		= *mi_eval_integer(&params->resolution);
		res = (int)pow (2.0, res);
		float size	= *mi_eval_scalar(&params->size);		

		pointSample *dispBuffer = (pointSample*) oh->sampleSequence; 
	

		miVector pointLocal;
		miVector    tangent, binormal;

		if (textureSpace == 0) {
			mi_point_to_world(state,&pointLocal,&state->point);
		}
		else if (textureSpace == 1){
			mi_point_to_object(state,&pointLocal,&state->point);
		}
		x              = pointLocal.x;  // Point to sample
		z              = pointLocal.z;  // Point to sample

		//if (textureSpace == 2) {
		//	x = state->tex_list[0].x;
		//	z = state->tex_list[0].y;
		//
		//	// get normal, tangent and binormal
		//	

		//	
 


		//	miUint  num = 0;
		//	/* we like to deal with triangles only, like in mi_tri_vectors() */
		//	if ( mi_query(miQ_PRI_NUM_VERTICES, state, 0, &num)  &&  num == 3) {
		//		mi_info("verts %i", num);
		//		miVector    tangents[3], binormals[3];				
		//		miVector *q[] = { &tangents[0], &tangents[1],  &tangents[2]};
		//		miVector *q1[] = { &binormals[0], &binormals[1],  &binormals[2]};


		//		if ( mi_query(miQ_PRI_VERTICES_BUMPS_U, state, 0, &q, 2) &&  mi_query(miQ_PRI_VERTICES_BUMPS_V, state, 0, &q1)) {
		//			mi_info("tangent1 %f %f %f", tangents[0].x,tangents[0].y,tangents[0].z );
		//			/* average between 3 for now*/
		//			mi_vector_add(&tangent,&tangents[0],&tangents[1]);
		//			mi_vector_add(&tangent,&tangent,&tangents[2]);
		//			mi_vector_normalize(&tangent);
	
		//			mi_vector_add(&binormal,&binormals[0],&binormals[1]);
		//			mi_vector_add(&binormal,&binormal,&binormals[2]);
		//			mi_vector_normalize(&binormal);


		//		}
		//		else {
		//			mi_fatal("Error getting derivs.");
		//		}
	
		//	}
		//	else {
		//			mi_info("Error nums.");
		//	}
		//		
		//	//}
		//	//else {
		//	//	mi_fatal("This is not a triangle. Number of tris %i", num);
		//	//}
		//}

		/* the displacement part */
        if (state->type == miRAY_DISPLACE) {	
			miVector displacement;
			//miVector	u, v, n;
			double newX, newZ;   		 


			drw::EvalData evaldata; 
			oh->context->eval2_xz((1.0/globalScale)* x,(1.0/globalScale)* z, evaldata);

			displacement.x = evaldata.disp[0] * globalScale;
			displacement.y = evaldata.disp[1] * globalScale;
			displacement.z = evaldata.disp[2] * globalScale;

			
			if (textureSpace == 1){
				mi_vector_to_world(state, &displacement, &displacement);	
			}
			//if (textureSpace == 2){
			//	miVector tempDisplacement;
			//	// disp = x * tangent + y * normal + z * binormal
			//	tempDisplacement.x = displacement.x * tangent.x + displacement.y * state->normal.x + displacement.z * binormal.x;
			//	tempDisplacement.y = displacement.x * tangent.y + displacement.y * state->normal.y + displacement.z * binormal.y;
			//	tempDisplacement.z = displacement.x * tangent.z + displacement.y * state->normal.z + displacement.z * binormal.z;
			//	displacement = tempDisplacement;

			//}
	
			displacement.x +=  inDisplacement.x;
			displacement.y +=  inDisplacement.y;
			displacement.z +=  inDisplacement.z; 			
			  
			//doesn't compile for maya 8.5 but seems to work nonetheless, so just disable next 4 lines for 8.5 
			miBoolean is_root;
			mi_query(miQ_FUNC_IS_ROOT,state,0,&is_root); 
			//this is strange, if I don't check for root the assignment crashes Maya 32bit but renders fine with Maya 64
			if(! is_root) 
				result->displacementVector = displacement;

			//mi_info("dX %f dZ %f R %f G %f B %f",x, z, result->displacementVector.x,result->displacementVector.y,result->displacementVector.z);
			mi_vector_from_world(state, &displacement, &displacement);							
	
			// get where this point ends up after displaceing for buffering
			miVector newPoint = state->point;
			miVector worldDisplace;
			if (textureSpace == 0) {
				mi_point_to_world (state,&newPoint,&newPoint);
				mi_vector_to_world (state,&worldDisplace,&displacement);
			}
			else if (textureSpace == 1){
				mi_point_to_object(state,&newPoint,&newPoint);
				mi_vector_to_object (state,&worldDisplace,&displacement);
			}
			 
			
			mi_vector_add(&newPoint, &newPoint, &worldDisplace);
			//mi_info("%f moved by %f to %f", x, worldDisplace.x, newPoint.x);
			
			// change normal and give back displacement length for the actual displacement
			miScalar dispLength = mi_vector_norm(&displacement);
			mi_vector_normalize(&displacement);			
			state->normal.x = displacement.x;
			state->normal.y = displacement.y; 
			state->normal.z = displacement.z;

			result->displacement = dispLength; 
 
			//save where the displaced point came from, so we can later get the correct texture placement
			if (oh->doBuffer)
			{	
  
				int bufferX, bufferZ;
				if (newPoint.x > 0)	bufferX = (int (newPoint.x /oh->stepSize))%res;
				else bufferX = res - 1 - ((int (-newPoint.x /oh->stepSize))%res);				
				if (newPoint.z > 0)	bufferZ = (int (newPoint.z /oh->stepSize))%res;
				else bufferZ = res - 1 - ((int (-newPoint.z /oh->stepSize))%res);

				pointSample *thisPointSample = &dispBuffer[bufferX *res + bufferZ];
				//mi_info("Buffering %i %i stepsize %f res %i",bufferX,bufferZ,oh->stepSize,res );	
				miLock *lock;
				mi_query(miQ_FUNC_LOCK, state, miNULLTAG, &lock);
				mi_lock(*lock);
				
				if (thisPointSample->samples == 0)
				{
					thisPointSample->x = newPoint.x - x;
					thisPointSample->z = newPoint.z - z;										
				}
				else //interpolate
				{
					thisPointSample->x = ((thisPointSample->samples * thisPointSample->x) + ( newPoint.x - x))/(thisPointSample->samples + 1);
					thisPointSample->z = ((thisPointSample->samples * thisPointSample->z) + ( newPoint.z - z))/(thisPointSample->samples + 1);					
				}
				thisPointSample->samples ++;
				//mi_info("Buffering %i %i with %f %f",bufferX,bufferZ,thisPointSample->x,thisPointSample->z);				
				mi_unlock(*lock);
			}

		}
		else
		{//mi_info("here");

			// the color part
			////////////////////

			// read where the point was before being displaced
			if (oh->doBuffer)
			{
				int bufferX, bufferZ;
				if (x > 0)	bufferX = (int (x /oh->stepSize))%res;
				else bufferX = res - 1 - ((int (-x /oh->stepSize))%res);				
				if (z > 0)	bufferZ = (int (z /oh->stepSize))%res;
				else bufferZ = res - 1 - ((int (-z /oh->stepSize))%res);
				
				pointSample * pointDif = getPointBuffer(state, bufferX, bufferZ, dispBuffer, res);
				//mi_info("Point %f %f with %i %i came from %f %f", x,z, bufferX,bufferZ, x - pointDif->x, z - pointDif->z);
				
				//smoothly interpolate from sample to sample
				float deltaX = ((x - ((int (x /oh->stepSize)) * oh->stepSize)) / oh->stepSize) - 0.5;
				if (x < 0) deltaX += 1;
				float deltaZ = ((z - ((int (z /oh->stepSize)) * oh->stepSize)) / oh->stepSize) - 0.5;
				if (z < 0) deltaZ += 1; 
				//mi_info("stepsize %f delta %f %f " , oh->stepSize, deltaX, deltaZ);				
				pointSample * verticalNb, * horizontalNb, * diagonalNb;

				int nbBufferX, nbBufferZ;
				if (deltaX > 0)
					nbBufferX = bufferX + 1;	 				
				else
				{
					nbBufferX = bufferX - 1;
					deltaX *= -1;
				} 				
				if (deltaZ > 0)
					nbBufferZ = bufferZ + 1;		
				else
				{ 
					nbBufferZ = bufferZ - 1;		
					deltaZ *= -1;
				}  
				
				horizontalNb = getPointBuffer(state, nbBufferX, bufferZ, dispBuffer, res);
				verticalNb = getPointBuffer(state, bufferX, nbBufferZ, dispBuffer, res);
				diagonalNb = getPointBuffer(state, nbBufferX, nbBufferZ, dispBuffer, res);
				//filter bilinear
				//first the two horizontal directions
				float horizX1, horizX2, horizZ1, horizZ2;
				horizX1 = (1-deltaX) * pointDif->x + deltaX * horizontalNb->x;
				horizZ1 = (1-deltaX) * pointDif->z + deltaX * horizontalNb->z;
				horizX2 = (1-deltaX) * verticalNb->x + deltaX * diagonalNb->x;
				horizZ2 = (1-deltaX) * verticalNb->z + deltaX * diagonalNb->z;
				// now interpolate vertical, result is the displacement, subract this from our point
				x -= (1 - deltaZ) * horizX1 + deltaZ * horizX2;
				z -= (1 - deltaZ) * horizZ1 + deltaZ * horizZ2;


			} 
			
			//evaluate and get the results from the ocean

			drw::EvalData evaldata; 
			oh->context->eval2_xz((1.0/globalScale)* x,(1.0/globalScale)* z, evaldata);

			result->displacementVector.x = evaldata.disp[0] * globalScale;
			result->displacementVector.y = evaldata.disp[1] * globalScale;
			result->displacementVector.z = evaldata.disp[2] * globalScale;

			result->normal.x = evaldata.normal[0]; 
			result->normal.y = evaldata.normal[1];
			result->normal.z = evaldata.normal[2];

			result->jminus  = evaldata.Jminus;
			result->jplus   = evaldata.Jplus;
			result->eminus.x = evaldata.Eminus[0];
			result->eminus.y = evaldata.Eminus[1];
			result->eminus.z = evaldata.Eminus[2];

			result->eplus.x = evaldata.Eplus[0];
			result->eplus.y = evaldata.Eplus[1];
			result->eplus.z = evaldata.Eplus[2];
			

			// find out how much the calculated normal is rotated and apply that rotation to state->normal to get our resulting normal
			mi_vector_normalize(&result->normal);
			miVector  rotAxis, yUp;
			yUp.x = yUp.z = 0; yUp.y = 1;
			mi_vector_prod(&rotAxis,&yUp,&result->normal);
			mi_vector_normalize(&rotAxis);
			miScalar angle = acos(mi_vector_dot(&yUp,&result->normal));
			miMatrix rotMatrix;
			mi_matrix_rotate_axis(rotMatrix,&rotAxis,angle);
			mi_vector_transform(&result->normal, &state->normal, rotMatrix);
			mi_vector_normalize(&result->normal);
		
		}//the color part
		return miTRUE;
}




/***************************************************************************
// A vector displacement shader for test purposes. Not yet fully functional. Use at own risk
/***************************************************************************/

extern "C" DLLEXPORT int nr_vectorDisplace_version(void) {return 1;}


struct nr_vectorDisplace{ 
	//int textureSpace; 	// 0= worldspace, 1 = uv Space
	// integer "textureSpace" default 0,		#:	enum "worldspace=0:texturespace=1"
	miScalar factor;
	miColor displacement;
	 
};

extern "C" DLLEXPORT miBoolean nr_vectorDisplace(
					miScalar *result, miState *state, struct nr_vectorDisplace *params)
{	
	//int textureSpace = *mi_eval_integer(&params->textureSpace);
	//mi_point_to_world(state,&state->point, &state->point);
	miColor inColor = *mi_eval_color(&params->displacement);
	miScalar factor = *mi_eval_scalar(&params->factor);	
	//mi_info("X %f Z %f R %f G %f B %f",state->point.x, state->point.z, inColor.r,inColor.g,inColor.b);

	miVector displacement; 
	miVector	u, v, n; 

	//if (textureSpace == 0)	
	//{
		displacement.x = inColor.r;
		displacement.y = inColor.g;
		displacement.z = inColor.b; 
		mi_vector_from_world(state, &displacement, &displacement);
	/*}
	/else
	{


		u = state->derivs[0];
		//v = state->derivs[1];

		mi_vector_normalize(&u);
		// Compute v to be perpendicular to u
		// (in the tangent plane) 
		mi_vector_prod(&v, &state->normal, &u);
		mi_vector_normalize(&v);
		// rot is u, gruen is normal, blau ist v
		//disp = r*u + g*normal + b*v
		mi_vector_mul(&u, inColor.r);
		mi_vector_mul(&v, inColor.b);
		n = state->normal;
		mi_vector_mul(&n, inColor.g);
		mi_vector_add(&displacement, &u, &n); 
		mi_vector_add(&displacement, &displacement, &v);
	}*/
	miScalar dispLength = mi_vector_norm(&displacement);
	mi_vector_normalize(&displacement);
	//if (abs(displacement.x) < 0.001)
	//{
		//mi_info("X %f Y %f Z %f Length %f", displacement.x, displacement.y, displacement.z, dispLength);
		//displacement.y = 2;
	//}
	state->normal = displacement;	
	*result = dispLength*factor;
	return miTRUE;
}