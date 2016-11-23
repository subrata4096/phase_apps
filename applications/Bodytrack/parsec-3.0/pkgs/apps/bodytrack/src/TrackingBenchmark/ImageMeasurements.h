//-------------------------------------------------------------
//      ____                        _      _
//     / ___|____ _   _ ____   ____| |__  | |
//    | |   / ___| | | |  _  \/ ___|  _  \| |
//    | |___| |  | |_| | | | | |___| | | ||_|
//     \____|_|  \_____|_| |_|\____|_| |_|(_) Media benchmarks
//                           
//	 © 2006, Intel Corporation, licensed under Apache 2.0 
//
//  file : ImageMeasurements.h
//  author :	Scott Ettinger	- scott.m.ettinger@intel.com
//				Jean-Yves Bouguet - jean-yves.bouguet@intel.com
//  description : Image Measurements (silhouette and edges)
//  modified : 
//--------------------------------------------------------------


#ifndef IMAGEMEASUREMENTS_H
#define IMAGEMEASUREMENTS_H

#if defined(HAVE_CONFIG_H)
# include "config.h"
#endif

#include <FlexLib.h>
#include "BinaryImage.h"
#include "ImageProjection.h"

#define FlexImage8u FlexImage<Im8u,1>
#define FlexImage32f FlexImage<Im32f,1>

#include <vector>

//Default granularity of sampling of the projected cylinder edges:
#define STEP_DEFAULT 5.00f //0.50f //5f

//Default horizontal and vertical granularities of sampling of the inside of the projected cylinder:
#define H_STEP_DEFAULT 10.00f // 1.00f //10f
#define V_STEP_DEFAULT 10.00f //1.00f //10f
#define CHECKPHASE() ((phase_manish != 0) && ((phase_manish == -1) || ((whileLoopCounter >= (total_iterations_manish * (phase_manish-1) / total_phases_manish)) && (whileLoopCounter < (total_iterations_manish * (phase_manish) / total_phases_manish))))) 
#define CHECKPHASE_MULTI( phase ) ((phase != 0) && ((phase == -1) || ((whileLoopCounter >= (total_iterations_manish * (phase-1) / total_phases_manish)) && (whileLoopCounter < (total_iterations_manish * (phase) / total_phases_manish)))))
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>
// extern 

extern double time_taken_edgeError_1      ; 
extern double time_taken_edgeError_2      ; 
extern double time_taken_insideError_inner; 
extern double time_taken_insideError_outer; 
extern double time_taken_ieEdge_inner     ; 
extern double time_taken_ieEdge_outer     ; 
extern double time_taken_ieInside_outer   ; 
extern double time_taken_ieInside_inner   ; 

extern int count_time_taken_edgeError_1 ;
extern int count_time_taken_edgeError_2 ;
extern int count_time_taken_insideError_inner ;
extern int count_time_taken_insideError_outer ;
extern int count_time_taken_ieEdge_inner ;
extern int count_time_taken_ieEdge_outer ;
extern int count_time_taken_ieInside_outer ;
extern int count_time_taken_ieInside_inner ;

extern int phase_manish ;
extern int whileLoopCounter ;
extern int total_phases_manish;
extern int total_iterations_manish;
extern int init_num_particles;

long
perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
               int cpu, int group_fd, unsigned long flags);

long long perfProfileRead(int* fdptr);
void perfProfileInit(struct perf_event_attr* peptr, int* fdptr);

void perfProfileStart(int* fdptr);

void perfProfileEnd(int* fdptr);


class ImageMeasurements
{
private:
	std::vector<Point> mSamples;											//pixel samples (for each body part, edge and inside)

	float mStep;															//Sampling resolution of the edges in pixels (default: STEP_DEFAULT)
	float mHstep,mVstep;													//Horizontal and vertical sampling resolutions of inside the limbs (defaults: H_STEP_DEFAULT,V_STEP_DEFAULT)
	
	//compute error at the edges of a projected cylinder (body part)
	void EdgeError(const ProjectedCylinder &ProjCyl, const FlexImage8u &EdgeMap, float &error, int &samplePoints);	

	//compute error inside of a projected cylinder (body part)
	void InsideError(const ProjectedCylinder &ProjCyl, const BinaryImage &FGmap, int &error, int &samplePoints);

public:
	ImageMeasurements(){SetSamplingResolutions(STEP_DEFAULT, H_STEP_DEFAULT, V_STEP_DEFAULT); };
	~ImageMeasurements(){};

	//Set sampling densities
	void SetSamplingResolutions(float edge, float inside_h, float inside_v){mStep=edge; mHstep=inside_h; mVstep=inside_v;};
	
	//Edge error of a complete body on all camera images
	float ImageErrorEdge(std::vector<FlexImage8u> &ImageMaps, MultiCameraProjectedBody &ProjBodies);

	//Silhouette error of a complete body on all camera images
	float ImageErrorInside(std::vector<BinaryImage> &ImageMaps, MultiCameraProjectedBody &ProjBodies);
};

#endif
