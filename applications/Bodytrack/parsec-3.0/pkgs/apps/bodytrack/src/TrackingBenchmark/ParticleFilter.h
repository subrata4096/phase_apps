//------------------------------------------------------------------------
//      ____                        _      _
//     / ___|____ _   _ ____   ____| |__  | |
//    | |   / ___| | | |  _  \/ ___|  _  \| |
//    | |___| |  | |_| | | | | |___| | | ||_|
//     \____|_|  \_____|_| |_|\____|_| |_|(_) Media benchmarks
//                           
//	 © 2006, Intel Corporation, licensed under Apache 2.0 
//
//  file :	 ParticleFilter.h
//  author : Scott Ettinger - scott.m.ettinger@intel.com
//
//  description : Generic Particle Filter class.  Supports annealed 
//				  particle filtering.  Templated on a model object
//				  which evaluates the particle likelihoods. Vector width
//				  is determined by the model initial state.  Number of 
//				  annealing steps is determined by model StdDevs() function.
//
//				  Model object must support member functions : 
//					std::vector<fpType> InitialState(); 
//					fpType LogLikelihood(std::vector<fpType> &v, bool &valid);
//					std::vector<std::vector<fpType> > StdDevs();
//					void GetObservation(fpType timeval);
//		
//  modified : 
//--------------------------------------------------------------------------

#ifndef PARTICLEFILTER_H
#define PARTICLEFILTER_H

#if defined(HAVE_CONFIG_H)
# include "config.h"
#endif

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sys/types.h>
#include "RandomGenerator.h"
#include "AnnealingFactor.h"

//Manish
#include "ImageMeasurements.h"

//Manish START
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//Manish END

#ifndef uint
#define uint unsigned int
#endif

#undef min

//Generic particle filter class templated on model object
template<class T> 
class ParticleFilter{	

public:
//Types
	typedef float fpType;
	typedef std::vector<fpType> Vectorf;

protected:
//variables
	T *mModel;												//templated model object evaluates particle likelihoods
	std::vector<Vectorf > mParticles, mNewParticles;		//lists of particles
	Vectorf	mWeights, mCdf;									//particle weights, cumulative distribution
	Vectorf mBins, mSamples, mLikelihoods;					//resampling bins, new samples, particle likelihoods

	int mNParticles;										//number of particles used
	int old_mNParticles_mParts;										//number of particles used
	int mBestParticle;										//index of particle with highest likelihood
	int mMinParticles;										//minimum number of particles allowed
	std::vector<RandomGenerator> mRnd;						//random number generators - should be replaced with a single parallel leapfrog generator for better quality random numbers.
	bool mInitialized;										//particle initialization flag

//functions
	void CalcCDF(const Vectorf &weights, Vectorf &dst);										//calculate the cumulative distribution function from particle weights
	void AddGaussianNoise(Vectorf &p, const Vectorf &stdDevs, RandomGenerator &rnd) const;	//distribute particle randomly according to given standard deviations
	void Resample(Vectorf &cdf, Vectorf &bins, Vectorf &samples, int n);					//Monte Carlo resampling given a CDF
	virtual void CalcWeights(std::vector<Vectorf > &particles);								//calculate particle weights based on model likelihood
	virtual void GenerateNewParticles(int k);												//generate new particles distributed by model annealing level std dev



public:
	//Constructors
	ParticleFilter()							{mMinParticles = 5; mInitialized = false;};			
	ParticleFilter(T &model)					{mModel = &model; mMinParticles = 5; mInitialized = false;}; 
	virtual ~ParticleFilter() {};									

	//Get Functions
	T &Model()									{return *mModel; };
	int NumParticles() const					{return mNParticles; };
	int GoodParticles()	const					{return (int)mParticles.size(); };
	std::vector<Vectorf > &StdDevs() const		{return mModel->StdDevs(); };
	
	//Set Functions
	void SetModel(T &model)						{mModel = &model; };
	void SetMinimumParticles(int n)				{mMinParticles = n;};

	//Set number of particles to n and generate initial values
	void InitializeParticles(int n);						

	//Update filter to a new set of particles - returns false if model observation fails
	//calls model to get observation at given time, and to get likelihoods of each particle
	virtual bool Update(fpType timeVal);
	
	//State estimate (weighted mean of all particles)
	void Estimate(Vectorf &estimate);
	
	//Return particle with highest likelihood
	void BestParticle(Vectorf &p) {p = mParticles[mBestParticle]; };
	
};

//------------------------------------------ Implementation -----------------------------------------------------

//vector multiply by scalar
template<class T1, class T2>
inline void operator*=(std::vector<T1> &v, const T2 c)
{	for(uint i = 0; i < v.size(); i++)
		v[i] *= c;
}

//vector constant subtraction
template<class T1, class T2>
inline void operator-=(std::vector<T1> &v, const T2 c)
{	for(uint i = 0; i < v.size(); i++)
		v[i] -= c;
}

//distribute particle randomly according to given standard deviations
template<class T>
inline void ParticleFilter<T>::AddGaussianNoise(Vectorf &p, const Vectorf &stdDevs, RandomGenerator &rnd) const
{	for(uint i = 0; i < stdDevs.size(); i++)
		p[i] += (fpType)rnd.RandN() * stdDevs[i];				
}

//calculate the CDF from particle weights
template<class T>
inline void ParticleFilter<T>::CalcCDF(const Vectorf &weights, Vectorf &dst)
{	dst.resize(weights.size());
	dst[0] = weights[0];
	for(uint i = 1; i < weights.size(); i++)							//compute cumulative sum
		dst[i] = dst[i - 1] + weights[i];
	dst *= fpType(1.0) / dst[dst.size() - 1];							//normalize cdf
}

//Monte Carlo resampling given a cdf.  Uses incremental eponential random values to
//create a sorted uniform random sample for fast inverse of the cdf for all samples in 1 pass
template<class T>
inline void ParticleFilter<T>::Resample(Vectorf &cdf, Vectorf &bins, Vectorf &samples, int n)
{
	RandomGenerator r;
	samples[0] = (fpType)r.RandExp();
	for(int i = 1; i < n; i++)											//generate a new set of sorted random samples
		samples[i] = samples[i - 1] + (fpType)r.RandExp();
	samples *= (fpType)1.0 / (samples[samples.size() - 1]);				//normalize 
	cdf[cdf.size() - 1] += 1;											//prevent overrun due to numerical error
	int p = 0;
	bins.assign(cdf.size(), 0);
	for(uint i = 0; i < samples.size(); i++)							//populate sample bins based on probability distribution
	{	while(cdf[p] < samples[i])
			p++;
		bins[p]++;
	}
}

//calculate particle weights (mWeights) and find highest likelihood particle. 
//computes an optimal annealing factor and scales the likelihoods. 
//Also removes any particles reported as invalid by the model.
template<class T>
void ParticleFilter<T>::CalcWeights(std::vector<Vectorf > &particles)
{	
	bool valid;
	mBestParticle = 0;
	fpType total = 0, best = 0, minWeight = 1e30f, annealingFactor = 1;
	mWeights.resize(particles.size());
	uint i = 0;
	int manishCounter =0;
	while(i < particles.size())											//compute likelihood weights for each particle
	{	
		manishCounter++;
		mWeights[i] = mModel->LogLikelihood(particles[i], valid);
		if(!valid)														//if not valid(model prior), remove the particle from the list
		{	particles[i] = particles[particles.size() - 1];
			particles.pop_back(); mWeights.pop_back();
		}
		else
			minWeight = std::min(mWeights[i++], minWeight);				//find minimum weight
	}
	// printf("CALCWEIGHT OUTER, iter while = %d, mWeights = %d, particles = %d, mMinParticles = %d, mstdevs = %d\n", manishCounter,  (int)mWeights.size(), (int)particles.size(), mMinParticles, (int)mModel->StdDevs().size());
	if((int)particles.size() < mMinParticles) return;					//bail out if not enough valid particles
	mWeights -= minWeight;												//shift weights to zero for numerical stability
	if(mModel->StdDevs().size() > 1) 
		annealingFactor = BetaAnnealingFactor(mWeights, 0.5f);			//calculate annealing factor if more than 1 step
	// printf("AnnealingFactor = %f\n", annealingFactor);
	for(uint i = 0; i < mWeights.size(); i++)
	{	double wa = annealingFactor * mWeights[i];
		mWeights[i] = (float)exp(wa);
		total += mWeights[i];											//save sum of all weights
		if(i == 0 || mWeights[i] > best)								//find highest likelihood particle
		{	best = mWeights[i];
			mBestParticle = i;
		}
	}
	mWeights *= fpType(1.0) / total;									//normalize weights
}

//Generate a set of inital particles
template<class T>
void ParticleFilter<T>::InitializeParticles(int n)
{	
	// printf("\n\n\n\nINITIALIZZING %d\n", n)
	mRnd.resize(n);														//initialize random number generators
	for(int i = 0; i < n; i++)
		mRnd[i].Seed(i * 2);
	mNParticles = n;
	Vectorf initVector = mModel->InitialState();						//get inital state vector 
	if(initVector.size() != mModel->StdDevs()[0].size() )
		std::cout << "Warning!! PF::Model initial Vector and stdDev vector are not equal lengths!" << std::endl;
	bool minValid = false;
	while(!minValid)
	{	mParticles.resize(n);
		for(int i = 0; i < n; i++)
		{	Vectorf &p = mParticles[i];									//create new particle with randomized initial value
			p = initVector;
			AddGaussianNoise(p, mModel->StdDevs()[0], mRnd[i]);			//distribute particles randomly about the initial vector
		}
		CalcWeights(mParticles);										//calculate initial weights and remove any particles that are invalid by model prior
		minValid = (int)mParticles.size() >= mMinParticles;				//repeat until minimum number of valid particles is met
		if(!minValid)							
			std::cout << "Warning : initial particle set does not meet minimum number of particles. Resampling.." << std::endl;
	}
	mCdf.resize(n);														//allocate space 
	mSamples.resize(n);
	mBins.resize(n);
	mInitialized = true;
}

//generate new particles distributed with std deviation given by the model annealing parameter
template<class T> 
void ParticleFilter<T>::GenerateNewParticles(int k)
{	int p = 0;
	mNewParticles.resize(mNParticles);
	for(uint i = 0; i < mBins.size(); i++)									//distribute new particles randomly according to model stdDevs
	{
		for(uint j = 0; j < mBins[i]; j++)
		{	
			// printf("Particle mNew p = %d, mP i = %d, j = %d, mBins.size() = %d, mBins[%d] = %f\n", p, i, j, mBins.size(), i, mBins[i]);
			// if(p >= mNParticles) 
			// 	return;
			mNewParticles[p] = mParticles[i];								//add new particle for each entry in each bin distributed randomly about duplicated particle
			AddGaussianNoise(mNewParticles[p], mModel->StdDevs()[k], mRnd[p]);
			p++;
			// printf("EXITParticle mNew %d, mP %d\n", p, i);
		}
	}
	// printf("EXITING\n\n\n");
}

//Particle filter update (model and observation updates must be called first)  
template<class T>
bool ParticleFilter<T>::Update(fpType timeval)								//weights have already been computed from previous step or initialization
{						
	if(!mInitialized)														//check for proper initialization
	{	std::cout << "Update Error : Particles not initialized" << std::endl; 
		return false;
	}	
	if(!mModel->GetObservation(timeval))
	{	std::cout << "Update Error : Model observation failed for time : " << timeval << std::endl;
		return false;
	}

    //Manish

    int loopchange = atoi(getenv("LOOPCHANGE_PFUPDATE"));
    int loopskip = atoi(getenv("LOOPITERSKIP_PFUPDATE"));
    int particleApproxKnob = atoi(getenv("PARTICLE_FACTOR_CHANGE"));
    int particleApproxKnob_final = atoi(getenv("PARTICLE_FACTOR_CHANGE"));

    clock_t manish_Clock;
    manish_Clock = clock();

    int manishforLoopCounter = 0;
	int manishCounter = 0;


	// printf("mMinParticles = %d, mNParticles = %d\n", mMinParticles, mNParticles);

                        int all_phase_change = (getenv("CHANGE_ALL_PHASE")) ? atoi(getenv("CHANGE_ALL_PHASE")) : 0;
  
                        if(all_phase_change != 0)
   			{
        		   if(CHECKPHASE_MULTI(1) == 1)
                           {
  				loopchange = 1;
				loopskip = 1;
                           }
                           if(CHECKPHASE_MULTI(2) == 1)
                           {
  				loopchange = 1;
				loopskip = 1;
 				
                           }
                           if(CHECKPHASE_MULTI(3) == 1)
                           {
  				loopchange = 1;
				loopskip = 1;
                           }
                           if(CHECKPHASE_MULTI(4) == 1)
                           {
  				loopchange = 1;
				loopskip = 1;
                           }
                        } 

    printf("ApproxLog: AB: ParticleFilter-Update\n");	

    for(int k = (int)mModel->StdDevs().size() - 1; k >= 0 ; k-=loopchange)//atoi(getenv("LOOPCHANGE")))			//loop over all annealing steps starting with highest
	{
		if((loopskip != 1) && ((k % loopskip) == 0))
			continue;

		// if(particleApproxKnob_final != -1)
		// {
			// mNParticles = temp_mMinParticles / particleApproxKnob_final;
		// }
		// else
		// {
			// mNParticles = temp_mMinParticles;
		// }



	    manishforLoopCounter += 1;
        CalcCDF(mWeights, mCdf);											//Monte Carlo re-sampling 
		Resample(mCdf, mBins, mSamples, mNParticles);		
		bool minValid = false;
		int old_mNParticles = mNParticles;
		old_mNParticles_mParts = mNParticles;
		while(!minValid)
		{
			
	// int mNParticles;										//number of particles used
	// int mBestParticle;										//index of particle with highest likelihood
	// int mMinParticles;										//minimum number of particles allowed
			particleApproxKnob_final = (CHECKPHASE()) ? particleApproxKnob : 1;
		
                        int all_phase_change = (getenv("CHANGE_ALL_PHASE")) ? atoi(getenv("CHANGE_ALL_PHASE")) : 0;
  
                        if(all_phase_change != 0)
   			{
        		   if(CHECKPHASE_MULTI(1) == 1)
                           {
  				particleApproxKnob_final = 1;
                           }
                           if(CHECKPHASE_MULTI(2) == 1)
                           {
  				particleApproxKnob_final = 1;
 				
                           }
                           if(CHECKPHASE_MULTI(3) == 1)
                           {
  				particleApproxKnob_final = 1;
                           }
                           if(CHECKPHASE_MULTI(4) == 1)
                           {
  				particleApproxKnob_final = 1;
                           }
                        } 
                        
			mNParticles = 	init_num_particles / particleApproxKnob_final;

			if(old_mNParticles != mNParticles)
			{
				printf("Changing particles from %d to %d\n", old_mNParticles, mNParticles);
				InitializeParticles(mNParticles);
	        	CalcCDF(mWeights, mCdf);											//Monte Carlo re-sampling 
				Resample(mCdf, mBins, mSamples, mNParticles);		

					// mCdf.resize(mNParticles);														//allocate space 
					// mSamples.resize(mNParticles);
					// mBins.resize(mNParticles);
					// mRnd.resize(mNParticles);
					// mParticles.resize(mNParticles);
			}
			old_mNParticles = mNParticles;
			// printf("whileLoopCounter = %d, checkphase = %d, mNParticles = %d, particleApproxKnob_final = %d\n", whileLoopCounter, CHECKPHASE(), mNParticles, particleApproxKnob_final);
			manishCounter++;
			whileLoopCounter++;
			// mNParticles = temp_mMinParticles;
			GenerateNewParticles(k);
			//MANISH
			// particleApproxKnob_final = (CHECKPHASE()) ? particleApproxKnob : 1;
			// mNParticles = temp_mMinParticles / particleApproxKnob_final;
			// mNewParticles.resize(mNParticles);
			//ENDMANISH
			CalcWeights(mNewParticles);										//calculate particle weights and remove any invalid particles
			minValid = (int)mNewParticles.size() >= mMinParticles;			//repeat if not enough valid particles
			if(!minValid) 
				std::cout << "Not enough valid particles - Resampling!!!" << std::endl;
		}
		mParticles = mNewParticles;											//save new particle set
	}

	printf("while loop minvalid, iter = %d\n", manishCounter);
	printf("for loop annealing layers, iter = %d\n", manishforLoopCounter);


    manish_Clock = clock() - manish_Clock;
    double time_taken = ((double)manish_Clock)/CLOCKS_PER_SEC; // in seconds

    printf("Time taken by loop pf update = %f\n",
               time_taken);

    // Global prints
    printf("Time taken by edge error loop 1 = %f\n",
               (time_taken_edgeError_1*100000) / count_time_taken_edgeError_1);
    // printf("init time edge1 = %f, count_time_taken_edgeError_1 = %d\n", time_taken_edgeError_1, count_time_taken_edgeError_1);

    printf("Time taken by loop edge error loop 2 = %f\n",
               (time_taken_edgeError_2*100000) / count_time_taken_edgeError_2);
    // printf("init time edge2 = %f, count_time_taken_edgeError_2 = %d\n", time_taken_edgeError_2, count_time_taken_edgeError_2);

    printf("Time taken by loop inner inside error  = %f\n",
               (time_taken_insideError_inner*100000) / count_time_taken_insideError_inner);


    printf("Time taken by loop pf outer inside error = %f\n",
               (time_taken_insideError_outer*100000) / count_time_taken_insideError_outer);

    printf("Time taken by loop inner error edge = %f\n",
               (time_taken_ieEdge_inner*100000) / count_time_taken_ieEdge_inner);


    printf("Time taken by loop outer error edge = %f\n",
               (time_taken_ieEdge_outer*100000) / count_time_taken_ieEdge_outer);


    printf("Time taken by loop inner error inside = %f\n",
               (time_taken_ieInside_inner*100000) / count_time_taken_ieInside_inner);


    printf("Time taken by loop outer error inside = %f\n",
               (time_taken_ieInside_outer*100000) / count_time_taken_ieInside_outer);

	time_taken_edgeError_1 =       0;
	time_taken_edgeError_2 =       0;
	time_taken_insideError_inner = 0;
	time_taken_insideError_outer = 0;
	time_taken_ieEdge_inner =      0;
	time_taken_ieEdge_outer =      0;
	time_taken_ieInside_outer =    0;
	time_taken_ieInside_inner =    0;

	return true;
}

//calculate expected value of the particle set
template<class T> 
void ParticleFilter<T>::Estimate(Vectorf &estimate)
{
	estimate.assign(mParticles[0].size(), 0);								//clear estimate values
	for(uint i = 0; i < mParticles.size(); i++)								//calculate expected value 
		for(uint j = 0; j < estimate.size(); j++)
			estimate[j] += mParticles[i][j] * mWeights[i];
}

#endif
