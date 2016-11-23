//-------------------------------------------------------------
//      ____                        _      _
//     / ___|____ _   _ ____   ____| |__  | |
//    | |   / ___| | | |  _  \/ ___|  _  \| |
//    | |___| |  | |_| | | | | |___| | | ||_|
//     \____|_|  \_____|_| |_|\____|_| |_|(_) Media benchmarks
//                         
//	 © 2006, Intel Corporation, licensed under Apache 2.0 
//
//  file : ImageMeasurements.cpp
//  author : Scott Ettinger - scott.m.ettinger@intel.com
//			 Jean-Yves Bouguet - jean-yves.bouguet@intel.com
//			 
//  description : Image Measurements (silhouette and edges)
//  modified : 
//--------------------------------------------------------------

#if defined(HAVE_CONFIG_H)
# include "config.h"
#endif

#include "ImageMeasurements.h"
#include <math.h>
#include <time.h>

using namespace std;

double time_taken_edgeError_1 =       0;
double time_taken_edgeError_2 =       0;
double time_taken_insideError_inner = 0;
double time_taken_insideError_outer = 0;
double time_taken_ieEdge_inner =      0;
double time_taken_ieEdge_outer =      0;
double time_taken_ieInside_outer =    0;
double time_taken_ieInside_inner =    0;


int count_time_taken_edgeError_1 =       0;
int count_time_taken_edgeError_2 =       0;
int count_time_taken_insideError_inner = 0;
int count_time_taken_insideError_outer = 0;
int count_time_taken_ieEdge_inner =      0;
int count_time_taken_ieEdge_outer =      0;
int count_time_taken_ieInside_outer =    0;
int count_time_taken_ieInside_inner =    0;

int phase_manish =    0;
int whileLoopCounter =    0;
int total_phases_manish = 0;
int total_iterations_manish = 0;
int init_num_particles = 0;

long
perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
               int cpu, int group_fd, unsigned long flags)
{
   int ret;

   ret = syscall(__NR_perf_event_open, hw_event, pid, cpu,
                  group_fd, flags);
   return ret;
}

long long perfProfileRead(int* fdptr)
{
   int fd = *fdptr;
   long long count;
   //printf("%d\n", fd);   
   read(fd, &count, sizeof(long long));

   printf("Number of HW instructions Used=%lld\n", count);
   close(fd);

   return count;
}

void perfProfileInit(struct perf_event_attr* peptr, int* fdptr)
{
  int fd = *fdptr;
   //printf("%d\n", fd);   

   peptr->type = PERF_TYPE_HARDWARE;
   peptr->size = sizeof(struct perf_event_attr);
   peptr->config = PERF_COUNT_HW_INSTRUCTIONS;
   peptr->disabled = 1;
   peptr->exclude_kernel = 1;
   peptr->exclude_hv = 1;
   peptr->exclude_idle = 1;

   fd = perf_event_open(peptr, 0, -1, -1, 0);
   if (fd == -1) {
      fprintf(stderr, "Error opening leader %llx\n", peptr->config);
      exit(EXIT_FAILURE);
   }
   *fdptr = fd;
}

void perfProfileStart(int* fdptr)
{
   int fd = *fdptr;
   //printf("%d\n", fd);   
   ioctl(fd, PERF_EVENT_IOC_RESET, 0);
   ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
}

void perfProfileEnd(int* fdptr)
{
   int fd = *fdptr;
   //printf("%d\n", fd);   
   ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
}


//2D vector magnitude
inline double mag(Point &p)
{	return sqrt((double)(p.x * p.x + p.y * p.y));
}

//accumulate error at a given edge sample point (round to nearest integral point)
inline void SampleEdgePoint(float xf, float yf, const FlexImage8u &EdgeMap, int &error, int &samplePoints)
{
	int x = int(xf + 0.5f), y = int(yf + 0.5f);
	if((x >= 0) && (x < EdgeMap.Width()) && (y >= 0) && (y < EdgeMap.Height())) //check image bounds
	{	int e = 255 - EdgeMap(x,y);												//get value from image map and compute difference
		// printf("SEP x = %d, y = %d, call = %d\n", x, y, EdgeMap(x,y));
		error += (e * e);														//sum squared error values
		samplePoints++;															//count points sampled
	}
}

//Generate Samples for points along the non-joint edges of the cylinder
void ImageMeasurements::EdgeError(const ProjectedCylinder &ProjCyl, const FlexImage8u &EdgeMap, float &error, int &samplePoints)
{
	int ErrorSSD = 0;
	const Point &p1 = ProjCyl.mPts[0];
	Point s1;																//get direction vector of side 1 of the 2D cylinder projection
	s1.Set(ProjCyl.mPts[1].x - p1.x, ProjCyl.mPts[1].y - p1.y);
	int n1 = max((int)(mag(s1) / mStep + 0.5), 4);							//compute number of points sampled (sample at least 4)
	float d1 = 1.0f / (float)n1++;											//get fraction of side length per sample

	const Point &p2 = ProjCyl.mPts[2];
	Point s2;																//repeat for side 2 of cylinder
	s2.Set(ProjCyl.mPts[3].x - p2.x, ProjCyl.mPts[3].y - p2.y);
	int n2 = max((int)(mag(s2) / mStep + 0.5), 4);					
	float d2 = 1.0f / (float)n2++;

	float delta = 0;
    //Manish
    int loopchange1 = atoi(getenv("LOOPCHANGE_EDGEERROR_1"));
    int loopskip1 = atoi(getenv("LOOPITERSKIP_EDGEERROR_1"));

    int loopchange2 = atoi(getenv("LOOPCHANGE_EDGEERROR_2"));
    int loopskip2 = atoi(getenv("LOOPITERSKIP_EDGEERROR_2"));

    int all_phase_change = (getenv("CHANGE_ALL_PHASE")) ? atoi(getenv("CHANGE_ALL_PHASE")) : 0;

    loopchange1 = (CHECKPHASE()) ? loopchange1 : 1;
    loopskip1 = CHECKPHASE() ? loopskip1 : 1;
    loopchange2 = CHECKPHASE() ? loopchange2 : 1;
    loopskip2 = CHECKPHASE() ? loopskip2 : 1;

    if(all_phase_change != 0)
    {
        if(CHECKPHASE_MULTI(1) == 1)
        {
    	    loopchange1 = (CHECKPHASE()) ? loopchange1 : 1;
    	    loopskip1 = CHECKPHASE() ? loopskip1 : 5;
    	    loopchange2 = CHECKPHASE() ? loopchange2 : 1;
   	    loopskip2 = CHECKPHASE() ? loopskip2 : 5;
        }
        if(CHECKPHASE_MULTI(2) == 1)
        {
    	    loopchange1 = (CHECKPHASE()) ? loopchange1 : 1;
    	    loopskip1 = CHECKPHASE() ? loopskip1 : 3;
    	    loopchange2 = CHECKPHASE() ? loopchange2 : 1;
   	    loopskip2 = CHECKPHASE() ? loopskip2 : 3;
        }
        if(CHECKPHASE_MULTI(3) == 1)
        {
    	    loopchange1 = (CHECKPHASE()) ? loopchange1 : 2;
    	    loopskip1 = CHECKPHASE() ? loopskip1 : 2;
    	    loopchange2 = CHECKPHASE() ? loopchange2 : 3;
   	    loopskip2 = CHECKPHASE() ? loopskip2 : 2;
        }
        if(CHECKPHASE_MULTI(4) == 1)
        {
    	    loopchange1 = (CHECKPHASE()) ? loopchange1 : 4;
    	    loopskip1 = CHECKPHASE() ? loopskip1 : 2;
    	    loopchange2 = CHECKPHASE() ? loopchange2 : 4;
   	    loopskip2 = CHECKPHASE() ? loopskip2 : 2;
        }
    }


    clock_t manish_Clock;
    manish_Clock = clock();

    // printf("n1 = %d, n2 = %d, p1.x  = %f, p1.y = %f, d1 = %f, s1.x=%f, s1.y = %f\n", n1, n2, p1.x, p1.y, d1, s1.x, s1.y);
	printf("ApproxLog: AB: EDGEERROR_1\n");
	for(int i = 0; i < n1; i+=loopchange1)												//generate sample points along each side of cylinder projection
	{
		if((loopskip1 != 1) && ((i % loopskip1) == 0))
			continue;
		float x = p1.x + delta * s1.x;
		float y = p1.y + delta * s1.y;
	    // printf("i = %d, delta = %f, x = %f, y = %f, p1.x  = %f, p1.y = %f, d1 = %f, s1.x=%f, s1.y = %f\n",i, delta, x, y, p1.x, p1.y, d1, s1.x, s1.y);
		SampleEdgePoint(x, y, EdgeMap, ErrorSSD, samplePoints);				//accumulate error at computed edge points on side 1
  		delta += d1 * loopchange1;
	}

    manish_Clock = clock() - manish_Clock;
    time_taken_edgeError_1 += ((double)manish_Clock)/CLOCKS_PER_SEC; // in seconds
    count_time_taken_edgeError_1 += 1;
    // printf("Time taken by edge error loop 1 = %f\n",
    //            time_taken);

	
    manish_Clock = clock();
	delta = 0;
	
	printf("ApproxLog: AB: EDGEERROR_2\n");
	for(int i = 0; i < n2; i+=loopchange2)
	{
		if((loopskip2 != 1) && ((i % loopskip2) == 0))
			continue;
		float x = p2.x + delta * s2.x;
		float y = p2.y + delta * s2.y;
		SampleEdgePoint(x, y, EdgeMap, ErrorSSD, samplePoints);				//accumulate error at comptued edge points on side 2
		delta += d2 * loopchange2;
	}
	error += (float)ErrorSSD / (255.0f * 255.0f);


    manish_Clock = clock() - manish_Clock;
    time_taken_edgeError_2 += ((double)manish_Clock)/CLOCKS_PER_SEC; // in seconds
    count_time_taken_edgeError_2 += 1;;
    // printf("Time taken by loop edge error loop 2 = %f\n",
    //            time_taken);
}

//check bounds of nearest integral point and accumulate error value
inline void SampleInsidePoint(float xf, float yf, const BinaryImage &FGmap, int &error, int &samplePoints)
{
	int x = int(xf + 0.5f), y = int(yf + 0.5f);
	if((x >= 0) && (x < FGmap.Width()) && (y >= 0) && (y < FGmap.Height())) //check image bounds
	{	int e = 1 - FGmap(x,y);												//get value from image map and compute difference
		error += e;															//sum squared error values (since err = {1,0} same as sum of errors)
		samplePoints++;														//count points sampled
	}
}

//Sample points inside the projected cylinder
void ImageMeasurements::InsideError(const ProjectedCylinder &ProjCyl, const BinaryImage &FGmap, int &error, int &samplePoints)
{
	const Point &p1 = ProjCyl.mPts[0], &p2 = ProjCyl.mPts[3];
	Point s1, s2;
	s1.Set(ProjCyl.mPts[1].x - p1.x, ProjCyl.mPts[1].y - p1.y);				//get vectors along sides
	s2.Set(ProjCyl.mPts[2].x - p2.x, ProjCyl.mPts[2].y - p2.y);
	Point m(p1.x + s1.x / 2.0f - (p2.x + s2.x / 2.0f), p1.y + s1.y / 2.0f - (p2.y + s2.y / 2.0f));
	int n1 = max((int)(mag(s1) / mVstep + 0.5), 4);							//compute number of points sampled on each side (sample at least 4)
	int n2 = max((int)(mag(m) / mHstep + 0.5f), 4);							//compute number of interior samples using the midpoint length
	float d1 = 1.0f / n1++;													//get fraction of side lengths per sample
	float d2 = 1.0f / n2;
	float delta1 = 0;
	Point e1, e2;

	    //Manish
    int loopchange1 = atoi(getenv("LOOPCHANGE_INSIDEERROR_OUTER"));
    int loopskip1 = atoi(getenv("LOOPITERSKIP_INSIDEERROR_OUTER"));

    int loopchange2 = atoi(getenv("LOOPCHANGE_INSIDEERROR_INNER"));
    int loopskip2 = atoi(getenv("LOOPITERSKIP_INSIDEERROR_INNER"));
    
    int all_phase_change = (getenv("CHANGE_ALL_PHASE")) ? atoi(getenv("CHANGE_ALL_PHASE")) : 0;

    loopchange1 = CHECKPHASE() ? loopchange1 : 1;
    loopskip1 = CHECKPHASE() ? loopskip1 : 1;
    loopchange2 = CHECKPHASE() ? loopchange2 : 1;
    loopskip2 = CHECKPHASE() ? loopskip2 : 1;
    
    if(all_phase_change != 0)
    {
        if(CHECKPHASE_MULTI(1) == 1)
        {
    	    loopchange1 = (CHECKPHASE()) ? loopchange1 : 1;
    	    loopskip1 = CHECKPHASE() ? loopskip1 : 6;
    	    loopchange2 = CHECKPHASE() ? loopchange2 : 1;
   	    loopskip2 = CHECKPHASE() ? loopskip2 : 6;
        }
        if(CHECKPHASE_MULTI(2) == 1)
        {
    	    loopchange1 = (CHECKPHASE()) ? loopchange1 : 1;
    	    loopskip1 = CHECKPHASE() ? loopskip1 : 5;
    	    loopchange2 = CHECKPHASE() ? loopchange2 : 1;
   	    loopskip2 = CHECKPHASE() ? loopskip2 : 5;
        }
        if(CHECKPHASE_MULTI(3) == 1)
        {
    	    loopchange1 = (CHECKPHASE()) ? loopchange1 : 2;
    	    loopskip1 = CHECKPHASE() ? loopskip1 : 4;
    	    loopchange2 = CHECKPHASE() ? loopchange2 : 3;
   	    loopskip2 = CHECKPHASE() ? loopskip2 : 4;
        }
        if(CHECKPHASE_MULTI(4) == 1)
        {
    	    loopchange1 = (CHECKPHASE()) ? loopchange1 : 4;
    	    loopskip1 = CHECKPHASE() ? loopskip1 : 4;
    	    loopchange2 = CHECKPHASE() ? loopchange2 : 4;
   	    loopskip2 = CHECKPHASE() ? loopskip2 : 4;
        }
    }


    clock_t manish_Clock1;
    clock_t manish_Clock2;
    // time_taken2 = 0;
	    manish_Clock1 = clock();

	printf("ApproxLog: AB: InsideError\n");
	for(int i = 0; i < n1; i+=loopchange1)					//(atoi(getenv("LOOPCHANGE")))							//generate sample points along each side of cylinder projection
	{ 		
		if((loopskip1 != 1) && ((i % loopskip1) == 0))
			continue;


	    e1.Set(p1.x + delta1 * s1.x, p1.y + delta1 * s1.y);
		e2.Set(p2.x + delta1 * s2.x, p2.y + delta1 * s2.y);
		m.Set(e2.x - e1.x, e2.y - e1.y);									//get vector between edge points
		delta1 += d1 * loopchange1;
		float delta2 = 0;
		
		// clock_t manish_Clock2;
	    manish_Clock2 = clock();
	  	
		for(int j = 0; j < n2; j+=loopchange2)	//(atoi(getenv("LOOPCHANGE")))										//generate interior samples
		{
			if((loopskip2 != 1) && ((i % loopskip2) == 0))
				continue;
			SampleInsidePoint(e1.x + delta2 * m.x, e1.y + delta2 * m.y, FGmap, error, samplePoints);
			delta2 += d2 * loopchange2;
		}
		manish_Clock2 = clock() - manish_Clock2;
	    time_taken_insideError_inner += ((double)manish_Clock2)/CLOCKS_PER_SEC; // in seconds
	    count_time_taken_insideError_inner += 1;
	}


    // printf("Time taken by loop inner inside error  = %f\n",
    //            time_taken2);
    manish_Clock1 = clock() - manish_Clock1;
    time_taken_insideError_outer += ((double)manish_Clock1)/CLOCKS_PER_SEC; // in seconds
    count_time_taken_insideError_outer += 1;
    // printf("Time taken by loop pf outer inside error = %f\n",
    //            time_taken);

}

//compute edge map error term for all cameras given the set of 2D body geometry projections
float ImageMeasurements::ImageErrorEdge(std::vector<FlexImage8u> &ImageMaps, MultiCameraProjectedBody &ProjBodies) 
{
	int samples = 0;
	float error = 0;
    //Manish
    int loopchange1 = atoi(getenv("LOOPCHANGE_IEEDGE_OUTER"));
    int loopskip1 = atoi(getenv("LOOPITERSKIP_IEEDGE_OUTER"));

    int loopchange2 = atoi(getenv("LOOPCHANGE_IEEDGE_INNER"));
    int loopskip2 = atoi(getenv("LOOPITERSKIP_IEEDGE_INNER"));

    loopchange1 = CHECKPHASE() ? loopchange1 : 1;
    loopskip1 = CHECKPHASE() ? loopskip1 : 1;
    loopchange2 = CHECKPHASE() ? loopchange2 : 1;
    loopskip2 = CHECKPHASE() ? loopskip2 : 1;

	clock_t manish_Clock1;
	clock_t manish_Clock2;
	// double time_taken2 = 0;
    manish_Clock1 = clock();


    // printf("Im size = %d\n", (int)ImageMaps.size());
	for(int i = 0; i < (int)ImageMaps.size(); i+=loopchange1)							//for each camera, compute the edge map error term
	{			
		if((loopskip1 != 1) && ((i % loopskip1) == 0))
			continue;
		int nParts = ProjBodies(i).Size();
		// printf("nParts = %d\n", nParts);
	    manish_Clock2 = clock();
		for(int j = 0; j < nParts; j+=loopchange2)	
		{		
			if((loopskip2 != 1) && ((j % loopskip2) == 0))
				continue;
									//accumulate edge error for each body part, counting samples
			EdgeError(ProjBodies(i)(j), ImageMaps[i], error, samples);
		}

		manish_Clock2 = clock() - manish_Clock2;
	    time_taken_ieEdge_inner += ((double)manish_Clock2)/CLOCKS_PER_SEC; // in seconds
	    count_time_taken_ieEdge_inner += 1;
	}


    // printf("Time taken by loop inner error edge = %f\n",
    //            time_taken2);
    manish_Clock1 = clock() - manish_Clock1;
    time_taken_ieEdge_outer += ((double)manish_Clock1)/CLOCKS_PER_SEC; // in seconds
    count_time_taken_ieEdge_outer += 1;
    // printf("Time taken by loop outer error edge = %f\n",
	    //            time_taken);

	//cout << "Samples = " << samples << endl;
	return (float)error / samples;											//normalize to number of samples
}

//compute silhouette error term for all cameras given the set of 2D body geometry projections
float ImageMeasurements::ImageErrorInside(std::vector<BinaryImage> &ImageMaps, MultiCameraProjectedBody &ProjBodies)
{
    //Manish
    int loopchange1 = atoi(getenv("LOOPCHANGE_IEINSIDE_OUTER"));
    int loopskip1 = atoi(getenv("LOOPITERSKIP_IEINSIDE_OUTER"));

    int loopchange2 = atoi(getenv("LOOPCHANGE_IEINSIDE_INNER"));
    int loopskip2 = atoi(getenv("LOOPITERSKIP_IEINSIDE_INNER"));

    loopchange1 = CHECKPHASE() ? loopchange1 : 1;
    loopskip1 = CHECKPHASE() ? loopskip1 : 1;
    loopchange2 = CHECKPHASE() ? loopchange2 : 1;
    loopskip2 = CHECKPHASE() ? loopskip2 : 1;

	int samples = 0;
	int error = 0;

	// double time_taken2 = 0;
	clock_t manish_Clock1;
	clock_t manish_Clock2;
    manish_Clock1 = clock();
	for(int i = 0; i < (int)ImageMaps.size(); i+=loopchange1)							//for each camera, compute the edge map error term
	{	
		if((loopskip1 != 1) && ((i % loopskip1) == 0))
			continue;
		int nParts = ProjBodies(i).Size();

		// clock_t manish_Clock2;
	    manish_Clock2 = clock();
		for(int j = 0; j < nParts; j+=loopchange2)	
		{			
			if((loopskip2 != 1) && ((j % loopskip2) == 0))
				continue;
											//accumulate edge error for each body part, counting samples
			InsideError(ProjBodies(i)(j), ImageMaps[i], error, samples);
		}

		manish_Clock2 = clock() - manish_Clock2;
	    time_taken_ieInside_inner += ((double)manish_Clock1)/CLOCKS_PER_SEC; // in seconds
	    count_time_taken_ieInside_inner += 1;
	}
    // printf("Time taken by loop inner error inside = %f\n",
    //            time_taken2);
    manish_Clock1 = clock() - manish_Clock1;
	    time_taken_ieInside_outer += ((double)manish_Clock2)/CLOCKS_PER_SEC; // in seconds
	    count_time_taken_ieInside_outer += 1;
    // printf("Time taken by loop outer error inside = %f\n",
    //            time_taken);

	//cout << "Samples = " << samples << endl;
	return (float)error / samples;
}
