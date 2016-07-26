//#include "randomKernel.h"
//#include "Macro.h"



typedef struct __attribute__ ((aligned)) ParticleStatus{
	float3 pos, dir;
	float energy, maxSigma, mass, charge;
	int ifPrimary;
}PS;


__kernel void initParticles(__global PS * particle, float T, float2 width, float3 sourceCenter, float m, float c, int randSeed){
	size_t gid = get_global_id(0);
	
//	if(gid == 1){
//		int size = sizeof(PS);
//		//printf("size of PS: %d\n", size);
//	}

	particle[gid].pos.z = sourceCenter.z;
	int iseed[2];
	iseed[0] = randSeed;
	iseed[1] = gid*gid;
	MTrng(iseed);
	particle[gid].pos.x = (MTrng(iseed) - 0.5f) * width.s0;
	particle[gid].pos.y = (MTrng(iseed) - 0.5f) * width.s1;
	
	particle[gid].dir = normalize((float3)(0.0f, 0.0f, 0.0f) - sourceCenter);
	
//	particle[gid].dir.x = -1.0f;
//	particle[gid].dir.y = 0.0f;
//	particle[gid].dir.z = 0.0f;

	particle[gid].energy = T;
	particle[gid].maxSigma = 0.0f;
	particle[gid].mass = m;
	particle[gid].charge = c;
	particle[gid].ifPrimary = 1;
}


bool ifInsidePhantom(int3 voxIndex, int3 phantomSize){
	if(voxIndex.x < 0 || voxIndex.x >= phantomSize.x
		|| voxIndex.y < 0 || voxIndex.y >= phantomSize.y
		|| voxIndex.z < 0 || voxIndex.z >= phantomSize.z)
		return false;
	return true;
}


float step2VoxBoundary(float3 pos, float3 dir, float3 voxSize, int * cb, int3 phantomSize, int3 voxIndex) {
	float stepX, stepY, stepZ;

	float3 phantomBoundary = convert_float3(phantomSize) * voxSize;

//	//printf("floor(0.0) = %f, ceil(0.0) = %f\n", floor(0.0f), ceil(0.0f));
	if(fabs(dir.x) < EPSILON)
		stepX = INF;
	else if(dir.x > 0)
		stepX = ((voxIndex.x + 1)*voxSize.x - phantomBoundary.x * 0.5f - pos.x)/dir.x;
	else
		stepX = (voxIndex.x * voxSize.x - phantomBoundary.x * 0.5f  - pos.x)/dir.x;

	if(fabs(dir.y) < EPSILON)
		stepY = INF;
	else if(dir.y > 0)
		stepY = ((voxIndex.y + 1)*voxSize.y - phantomBoundary.y * 0.5f - pos.y)/dir.y;
	else
		stepY = (voxIndex.y * voxSize.y - phantomBoundary.y * 0.5f - pos.y)/dir.y;

	if(fabs(dir.z) < EPSILON)
		stepZ = INF;
	else if(dir.z > 0)
		stepZ = ((voxIndex.z + 1)*voxSize.z - phantomBoundary.z * 0.5f - pos.z)/dir.z;
	else
		stepZ = (voxIndex.z * voxSize.z - phantomBoundary.z * 0.5f - pos.z)/dir.z;

	float minStep;
	if(stepX < stepY){
		minStep = stepX;
		if(minStep < stepZ)
			*cb = 1;
		else{
			minStep = stepZ;
			*cb = 3;
		}
	}
	else{
		minStep = stepY;
		if(minStep < stepZ)
			*cb = 2;
		else{
			minStep = stepZ;
			*cb = 3;
		}
	}
	
//	if(isnan(stepX) || isnan(stepY) || isnan(stepZ))
		//printf("stepX = %f, stepY = %f, stepZ = %f\n", stepX, stepY, stepZ);

	return fabs(minStep);
}

inline float totalLinearSigma(float4 vox, float4 eleWeight, read_only image2d_t MCS, read_only image2d_t ICS, float e) {
	sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;
	float ics = read_imagef(ICS, sampler, (float2)(e - 1.0f + 0.5f, 0.5f)).s0;
	float4 mcs = read_imagef(MCS, sampler, (float2)(e - 1.0f + 0.5f, 0.5f));
	#if __ONLYEM__ == 1
		return ics*vox.s3;
	#endif
	return ics*vox.s3 + dot(mcs, eleWeight)*vox.s2;
}

inline float gamma(PS * particle) {
	return (particle->energy + particle->mass) / particle->mass;
}
inline float beta(PS * particle) {
	return sqrt(1 - particle->mass*particle->mass/((particle->energy + particle->mass)*(particle->energy + particle->mass)));
}

inline float maxDeltaElectronEnergy(PS * particle) {
	return (2 * ME*(gamma(particle)*gamma(particle) - 1.0f)) / (1 + 2 * gamma(particle)*ME / particle->mass + ME*ME / (particle->mass*particle->mass));
}

inline float maxOxygenEnergy(PS * particle){
	return (2 * MO*beta(particle)*beta(particle)*gamma(particle)*gamma(particle)) / (1 + 2 * gamma(particle)*MO / particle->mass + MO*MO / (particle->mass*particle->mass));
}

inline float momentumSquare(PS * particle) {
	return particle->energy*particle->energy + 2 * particle->energy*particle->mass;
}

float radiationLength(float4 vox, int ifWater)  
//	calculate the radiation length ratio \rho_wX_w/(\rhoX_0(\rho)
{
	if(ifWater)
		return XW;
	float ratio;
	if (vox.s2 >= 0.9)
	{
		ratio = 1.19f + 0.44f*log(vox.s2 - 0.44f);
	}
	else if (vox.s2 >= 0.26)
	{
		ratio = 1.0446f - 0.2180f*vox.s2;
	}
	else
	{
		ratio = 0.9857f + 0.0085f*vox.s2;
	}
	return WATERDENSITY*XW / (ratio*vox.s2);
}



float energyInOneStep(float4 vox, PS * particle, read_only image2d_t RSPW, read_only image2d_t MSPR, float stepLength, int ifWater) {
	//calculate equivalent step in water
	sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;
	float stepInWater;
	float4 mspr = read_imagef(MSPR, sampler, (float2)(particle->energy - 1.0f + 0.5f, vox.s1 + 0.5f));
	stepInWater = mspr.s0*stepLength*vox.s2/WATERDENSITY;

	if(ifWater)
		stepInWater = stepLength;

	//calculate energy transfer
	float4 rspw = read_imagef(RSPW, sampler, (float2)(particle->energy - 1.0f + 0.5f, 0.5f));
//	printf("energy = %f, rspw = %v4f\n", particle->energy, rspw);
	float de1 = stepInWater*rspw.s0;
	float b = rspw.s1;

	float temp = particle->energy/particle->mass;
	float eps = de1/particle->energy;
	float t = de1*(1.0f + eps/(1.0f+temp)/(2.0f+temp) + eps*eps*(2.0f+2.0f*temp+temp*temp)/(1.0f+temp)/(1.0f+temp)/(2.0f+temp)/(2.0f+temp)
			- b*eps*(0.5f+2.0f*eps/3.0f/(1.0f+temp)/(2.0f+temp) + (1.0f-b)*eps/6.0f) );
//	if(isnan(t))
		//printf("de1 = %f, b = %f, temp = %f, eps = %f\n", de1, b, temp, eps);
	return t;
}


void transform(float3 * dir, float theta, float phi){
	// if original direction is along z-axis
//	//printf("direction before transform %f %f %f\n", (*dir).x, (*dir).y, (*dir).z);
	float temp = 1.0 - ZERO;
	if ((*dir).z*(*dir).z >= temp){
		if ((*dir).z > 0){
			(*dir).x = sin(theta)*cos(phi);
			(*dir).y = sin(theta)*sin(phi);
			(*dir).z = cos(theta);
		}
		else{
			(*dir).x = -sin(theta)*cos(phi);
			(*dir).y = -sin(theta)*sin(phi);
			(*dir).z = -cos(theta);
		}
	}
	else{
		float u, v, w;
		u = (*dir).x*cos(theta) + sin(theta)*((*dir).x*(*dir).z*cos(phi) - (*dir).y*sin(phi)) / sqrt(1.0 - (*dir).z*(*dir).z);
		v = (*dir).y*cos(theta) + sin(theta)*((*dir).y*(*dir).z*cos(phi) + (*dir).x*sin(phi)) / sqrt(1.0 - (*dir).z*(*dir).z);
		w = (*dir).z*cos(theta) - sqrt(1.0f - (*dir).z*(*dir).z)*sin(theta)*cos(phi);

		(*dir).x = u;
		(*dir).y = v;
		(*dir).z = w;
	}
	
//	if(any(isnan(*dir))){
		//printf("transform result: %v3f\n", *dir);
		//printf("theta %f, phi %f\n", theta, phi);
//	}
//	//printf("direction after transform %f %f %f\n", (*dir).x, (*dir).y, (*dir).z);

	*dir = normalize(*dir);
//	//printf("direction after normalization %f %f %f\n", (*dir).x, (*dir).y, (*dir).z);

}

inline void update(PS * thisOne, float stepLength, float energyTransfer, float theta, float phi, int crossBound, int3 * voxIndex){	
//	//printf("before update dir %f %f %f, pos %f %f %f\n", thisOne->dir.x, thisOne->dir.y, thisOne->dir.z, thisOne->pos.x, thisOne->pos.y, thisOne->pos.z);
	
	thisOne->pos += thisOne->dir*stepLength;

	switch(crossBound){
	case(1):
		if(thisOne->dir.x > 0)
			(*voxIndex).x += 1;
		else
			(*voxIndex).x -= 1;
		break;

	case(2):
		if(thisOne->dir.y > 0)
			(*voxIndex).y += 1;
		else
			(*voxIndex).y -= 1;
		break;

	case(3):
		if(thisOne->dir.z > 0)
			(*voxIndex).z += 1;
		else
			(*voxIndex).z -= 1;
		break;

	case(0):
		break;
	}

	transform(&(thisOne->dir), theta, phi);
	thisOne->energy -= energyTransfer;
//	//printf("after update dir %f %f %f, pos %f %f %f\n", thisOne->dir.x, thisOne->dir.y, thisOne->dir.z, thisOne->pos.x, thisOne->pos.y, thisOne->pos.z);

}

inline void atomicAdd(volatile global float * source, const float operand) {
	union {
		unsigned int intVal;
		float floatVal;
	} newVal;
	union {
		unsigned int intVal;
		float floatVal;
	} prevVal;

	do {
		prevVal.floatVal = *source;
		newVal.floatVal = prevVal.floatVal + operand;
	} while (atomic_cmpxchg((volatile global unsigned int *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);

}

void getMutex(global int * mutex){
	int occupied = atomic_xchg(mutex, 1);
	while(occupied > 0){
		occupied = atomic_xchg(mutex, 1);
	}
}

void releaseMutex(global int * mutex){
	int prevVal = atomic_xchg(mutex, 0);
}


void score(global float* doseCounter, int scoringType, int absIndex, int nVoxels, float energyTransfer, float stepLength, PS * thisOne, int * iseed, read_only image2d_t MSPR, float materialID){
	// choose a dose counter
	int doseCounterId = convert_int_rtn(MTrng(iseed)*NDOSECOUNTERS);
	volatile global float * counter = &doseCounter[absIndex + doseCounterId * nVoxels];
	sampler_t dataSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

	// score: 0 is DOSE2MEDIUM, 1 is DOSE2WATTER, 2 is FLUENCE, 3 is LET
	switch(scoringType){
	case(0):
		atomicAdd(counter, energyTransfer);
		break;
	case(1):
		energyTransfer = energyTransfer/read_imagef(MSPR, dataSampler, (float2)(thisOne->energy - 0.5f, materialID + 0.5f)).s0;
		atomicAdd(counter, energyTransfer);
		break;
	case(3):
		atomicAdd(counter, energyTransfer*energyTransfer/stepLength);
		break;
	}
}

void scoreFluence(global float* doseCounter, int scoringType, int absIndex, int nVoxels, int * iseed){
	if(scoringType != 2)
		return;
	int doseCounterId = convert_int_rtn(MTrng(iseed)*NDOSECOUNTERS);
	volatile global float * counter = &doseCounter[absIndex + doseCounterId * nVoxels];
	atomicAdd(counter, 1.0f);
}


void store(PS * newOne, __global PS * secondary, volatile __global uint * nSecondary, global int * mutex2Secondary){
	if(*nSecondary == 0){
//		printf("\n secondary particle overflow!!!\n");
		return;	
	}
	
//	getMutex(mutex2Secondary);
	uint ind = atomic_dec(nSecondary);
	secondary[ind-1] = *newOne;
//	releaseMutex(mutex2Secondary);

//	//printf("store to # %d\n", *nSecondary);
//	//printf("stored secondary status: energy %f, pos %v3f, dir %v3f, ifPrimary %d, mass %f, charge %f, maxSigma %f\n", newOne->energy, newOne->pos, newOne->dir, newOne->ifPrimary, newOne->mass, newOne->charge, newOne->maxSigma);

}


void ionization(PS * thisOne, int * iseed, float * energyTransfer){
	

	float E = thisOne->energy + thisOne->mass;
	float Te;
	float rand;
	while(true){
		rand = MTrng(iseed);
		Te = MINELECTRONENERGY*maxDeltaElectronEnergy(thisOne)
			/((1-rand)*maxDeltaElectronEnergy(thisOne)+rand*MINELECTRONENERGY);
		if(MTrng(iseed) < 1.0-beta(thisOne)*beta(thisOne)*Te/maxDeltaElectronEnergy(thisOne)+Te*Te/(2.0*E*E))
			break;
	}
	*energyTransfer += Te;
}

void inelastic(PS * thisOne, read_only image2d_t yields, read_only image3d_t products, global PS * secondary, volatile __global uint * nSecondary, int * iseed,  global int * mutex2Secondary){
	uint nParticles = get_image_height(yields);
	PS newOne;
	float4 thisYield, thisProduct;
	sampler_t linearSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;
	sampler_t nearSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;
//	//printf("thisOne energy = %f, nParticles = %d\n", thisOne->energy, nParticles);
	for(uint i = 0; i < nParticles; i++){
		thisYield = read_imagef(yields, nearSampler, (float2)(thisOne->energy/12.0/10.0 - 1.0 + 0.5, (float)(i) + 0.5f));
		if(fabs(thisYield.s2) < 0.1f)
			continue;
//		//printf("yield = %f, %f, %f\n", thisYield.s0, thisYield.s1, thisYield.s2);
		while(MTrng(iseed) < thisYield.s0){
			float t = MTrng(iseed);
			thisProduct = read_imagef(products, linearSampler, (float4)(thisOne->energy/12.0/10.0 - 1.0 + 0.5, t*256 - 1.0 + 0.5, (float)(i) + 0.5f, 0));
//			//printf("particle index = %d, nParticles = %d, energy = %f, random = %f\n product = %v4f\n", i, nParticles, thisOne->energy, t, thisProduct);
			newOne.energy = thisProduct.s0;
			newOne.ifPrimary = 0;
			newOne.mass = thisYield.s1;
			newOne.charge = thisYield.s2;
			newOne.maxSigma = 0.0f;
			newOne.pos = thisOne->pos;
			newOne.dir = thisOne->dir;
			thisProduct = read_imagef(products, linearSampler, (float4)(thisOne->energy/12.0/10.0 - 1.0 + 0.5, MTrng(iseed)*256 - 1.0 + 0.5, (float)(i) + 0.5f, 0));
			transform(&(newOne.dir), thisProduct.s1, MTrng(iseed)*2.0f*PI);
			
			store(&newOne, secondary, nSecondary, mutex2Secondary);
			thisYield.s0 -= 1.0f; 
		}	
	
	}
}



void hardEvent(PS * thisOne, float * energyTransfer, float4 vox, float4 eleWeight, read_only image2d_t MCS, read_only image2d_t ICS, global PS * secondary, volatile __global uint * nSecondary, 
				int * iseed, global int * mutex2Secondary, __read_only image2d_t CHYield, __read_only image2d_t COYield, __read_only image2d_t CCYield,  __read_only image2d_t CCaYield,
				__read_only image3d_t CHProducts, __read_only image3d_t COProducts, __read_only image3d_t CCProducts, __read_only image3d_t CCaProducts){
	sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;
	float ics = read_imagef(ICS, sampler, (float2)(thisOne->energy - 1.0f + 0.5f, 0.5f)).s0;
	float4 mcs = read_imagef(MCS, sampler, (float2)(thisOne->energy - 1.0f + 0.5f, 0.5f));
	float sigIon = ics*vox.s3;
	float sigCHI = mcs.s0*eleWeight.s0*vox.s2;
	float sigCOI = mcs.s1*eleWeight.s1*vox.s2;
	float sigCCI = mcs.s2*eleWeight.s2*vox.s2;
	float sigCCaI = mcs.s3*eleWeight.s3*vox.s2;
	float sigma = sigIon + sigCHI + sigCOI + sigCCI + sigCCaI;
	#if __ONLYEM__ == 1
		sigma = sigIon;
	#endif
//	sigma = sigIon;
//	//printf("energy = %f, mcs = %v4f\n", thisOne->energy, mcs);
	float rand = MTrng(iseed)*(thisOne->maxSigma > sigma ? thisOne->maxSigma : sigma);
//	//printf("rand = %f, sigIon = %f, sigCHI = %f, sigCOI = %f, thisOne->maxSigma = %f\n", rand, sigIon, sigCHI, sigCOI, thisOne->maxSigma);
	#if __ONLYEM__ == 1
	if(rand < sigIon){
		ionization(thisOne, iseed, energyTransfer);
		return;
	}
	#else
	if(rand < sigIon){
		ionization(thisOne, iseed, energyTransfer);
		return;
	}
	else if(rand < sigIon + sigCHI){
		inelastic(thisOne, CHYield, CHProducts, secondary, nSecondary, iseed, mutex2Secondary);
		thisOne->energy = 0.0f;
		return;
	}
	else if(rand < sigIon + sigCHI + sigCOI){
		inelastic(thisOne, COYield, COProducts, secondary, nSecondary, iseed, mutex2Secondary);
		thisOne->energy = 0.0f;
		return;
	}
	else if(rand < sigIon + sigCHI + sigCOI + sigCCI){
		inelastic(thisOne, CCYield, CCProducts, secondary, nSecondary, iseed, mutex2Secondary);
		thisOne->energy = 0.0f;
		return;
	}
	else if(rand < sigIon + sigCHI + sigCOI + sigCCI + sigCCaI){
		inelastic(thisOne, CCaYield, CCaProducts, secondary, nSecondary, iseed, mutex2Secondary);
		thisOne->energy = 0.0f;
		return;
	}
	else
		return;
	#endif
}

void rayTrace(PS * particle, int3 phantomSize, float3 voxSize){
	float3 phantomBoundary1, phantomBoundary2;
	phantomBoundary1 = -convert_float3(phantomSize)*voxSize/2.0f;
	phantomBoundary2 = convert_float3(phantomSize)*voxSize/2.0f;

	float3 delta1, delta2, delta;
	delta1 = (phantomBoundary1 - particle->pos)/particle->dir;
	delta2 = (phantomBoundary2 - particle->pos)/particle->dir;
	delta =	fmin(delta1, delta2);

	float translation = fmax(fmax(delta.x, delta.y), delta.z);
//	//printf("particle pos = %f, %f, %f, dir = %f, %f, %f, delta = %v3f, translation = %f\n", 
//		particle->pos.x, particle->pos.y, particle->pos.z, particle->dir.x, particle->dir.y, particle->dir.z, delta, translation);	
	particle->pos += (translation + 1e-5) * particle->dir;
//	//printf("particle pos = %f, %f, %f, dir = %f, %f, %f, delta = %v3f, translation = %f\n", 
//		particle->pos.x, particle->pos.y, particle->pos.z, particle->dir.x, particle->dir.y, particle->dir.z, delta, translation);
}

__kernel void propagate(__global PS * particle, __global float * doseCounter, int scoringType, __read_only image3d_t voxels, __read_only image3d_t elementWeight, float3 voxSize, __read_only image2d_t MCS, 
		__read_only image2d_t ICS, __read_only image2d_t RSPW, __read_only image2d_t MSPR, __read_only image2d_t CHYield, __read_only image2d_t COYield, __read_only image2d_t CCYield, 
		__read_only image2d_t CCaYield, __read_only image3d_t CHProducts, __read_only image3d_t COProducts, __read_only image3d_t CCProducts, __read_only image3d_t CCaProducts, 
		__global PS * secondary, volatile __global uint * nSecondary, int randSeed, int ifWater, __global int * mutex){

	size_t gid = get_global_id(0);
	PS thisOne = particle[gid];
	

	int3 phantomSize = (int3)(get_image_width(voxels), get_image_height(voxels), get_image_depth(voxels));
	int nVoxels = phantomSize.x * phantomSize.y * phantomSize.z;

	float stepLength, stepInWater, thisMaxStep, step2bound, energyTransfer, sigma1, sigma2, sigma, sampledStep, variance, theta0, theta, phi, energyStrag, es, residualRange;
	int3 voxIndex;
	int absIndex, absIndex2, crossBound;
	int iseed[2];
	iseed[0] = gid*gid;
	iseed[1] = randSeed;
	MTrng(iseed);
	sampler_t voxSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;
	sampler_t dataSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;
	
	float4 vox, eleWeight;
	bool ifHard;
	absIndex = -1;
	absIndex2 = -1;

	voxIndex = convert_int3_rtn(thisOne.pos / voxSize + convert_float3(phantomSize) / 2.0f);

	if(!ifInsidePhantom(voxIndex, phantomSize)){
		rayTrace(&thisOne, phantomSize, voxSize);
		voxIndex = convert_int3_rtn(thisOne.pos / voxSize + convert_float3(phantomSize) / 2.0f);
	}

	while (ifInsidePhantom(voxIndex, phantomSize)){
		absIndex = voxIndex.x + voxIndex.y*phantomSize.x + voxIndex.z*phantomSize.x*phantomSize.y;
		vox = read_imagef(voxels, voxSampler, (float4)(convert_float3(voxIndex), 0.0f));
		eleWeight = read_imagef(elementWeight, voxSampler, (float4)(convert_float3(voxIndex), 0.0f));

		if(ifWater){
			eleWeight.s0 = 2.0f*HYDROGEN/(2.0f*HYDROGEN + OXYGEN);
			eleWeight.s1 = OXYGEN/(2.0f*HYDROGEN + OXYGEN);
			eleWeight.s2 = 0.0f;
			eleWeight.s3 = 0.0f;
		}


		if(absIndex != absIndex2){
			absIndex2 = absIndex;
			scoreFluence(doseCounter, scoringType, absIndex, nVoxels, iseed);
		}

//		if(gid == 0)
//			printf("position %v3f, energy %f, voxel %v3d, abs voxel %d\n", thisOne.pos, thisOne.energy, voxIndex, absIndex);
		
		step2bound = step2VoxBoundary(thisOne.pos, thisOne.dir, voxSize, &crossBound, phantomSize, voxIndex);
		residualRange = read_imagef(RSPW, dataSampler, (float2)(thisOne.energy - 1.0f + 0.5f, 0.5f)).s2 / vox.s2;

		if (residualRange < step2bound || thisOne.energy < 50.0f){ // 100MeV correspon to range of 0.01 cm
//			if(gid == 0)
//				printf("residualRange = %v4f, step2bound = %f\n", read_imagef(RSPW, dataSampler, (float2)(thisOne.energy - 1.0f + 0.5f, 0.5f)), step2bound);
			stepInWater = thisOne.energy / read_imagef(RSPW, dataSampler, (float2)(thisOne.energy - 1.0 + 0.5f, 0.5f)).s0;
			stepLength = stepInWater*WATERDENSITY / vox.s2 / read_imagef(MSPR, dataSampler, (float2)(thisOne.energy - 0.5f, vox.s1 + 0.5f)).s0;
			if(ifWater)
				stepLength = stepInWater;
			
			score(doseCounter, scoringType, absIndex, nVoxels, thisOne.energy, stepLength, &thisOne, iseed, MSPR, vox.s1);
			return;
		}

		// rescale maxStep to let energy transferred in one step < MAXENERGYRATIO
		thisMaxStep = MAXSTEP;
		energyTransfer = energyInOneStep(vox, &thisOne, RSPW, MSPR, thisMaxStep, ifWater);

		if (energyTransfer > MAXENERGYRATIO*thisOne.energy){
			stepInWater = MAXENERGYRATIO*thisOne.energy / read_imagef(RSPW, dataSampler, (float2)((1 - 0.5f*MAXENERGYRATIO)*thisOne.energy - 1.0f + 0.5f, 0.5f)).s0;
			thisMaxStep = stepInWater*WATERDENSITY / vox.s2 / read_imagef(MSPR, dataSampler, (float2)((1 - 0.5f*MAXENERGYRATIO)*thisOne.energy - 0.5f, vox.s1 + 0.5f)).s0;
			if(ifWater)
				thisMaxStep = stepInWater;
		}		
		

		// get linear macro cross section to sample a step
		energyTransfer = energyInOneStep(vox, &thisOne, RSPW, MSPR, thisMaxStep, ifWater);
		sigma1 = totalLinearSigma(vox, eleWeight, MCS, ICS, thisOne.energy);
		sigma2 = totalLinearSigma(vox, eleWeight, MCS, ICS, thisOne.energy - energyTransfer);
		sigma = sigma1 > sigma2 ? sigma1 : sigma2;


		// sample one step
		if(sigma < ZERO)
			sampledStep = INF;
		else
			sampledStep = -log(MTrng(iseed)) / sigma;
		if (sampledStep > step2bound || sampledStep > thisMaxStep){
			ifHard = false;
			if(step2bound < thisMaxStep)
				sampledStep = step2bound;
			else{
				sampledStep = thisMaxStep;
				crossBound = 0;
			}
		}
		else{
			ifHard = true;
			crossBound = 0;
		}
		stepLength = sampledStep;
		stepLength = stepLength > ZERO ? stepLength : ZERO;


//		//printf("sample step: if crossBound %d, steplength %f\n", crossBound, stepLength);

		// get energy transferred (plus energy straggling) in this sampled step
		energyTransfer = energyInOneStep(vox, &thisOne, RSPW, MSPR, stepLength, ifWater);
		variance = TWOPIRE2MENEW*vox.s3*stepLength*thisOne.charge*thisOne.charge*fmin(MINELECTRONENERGY, maxDeltaElectronEnergy(&thisOne))*(1.0f / (beta(&thisOne)*beta(&thisOne)) - 0.5f);
		do{
			energyStrag = MTGaussian(iseed) * sqrt(variance);
		}while(energyTransfer + energyStrag < 0);
		energyTransfer += energyStrag;
		energyTransfer = energyTransfer < thisOne.energy ? energyTransfer : thisOne.energy;

		es = 19.8 + 0.0023*thisOne.energy;
		theta0 = es*sqrt(stepLength/radiationLength(vox, ifWater))*thisOne.charge/beta(&thisOne)/sqrt(momentumSquare(&thisOne))
				*(1.0f + 0.038f*log(stepLength/radiationLength(vox, ifWater)));
		theta = MTGaussian(iseed) * theta0;
		phi = 2.0f*PI*MTrng(iseed);
		thisOne.maxSigma = sigma;

		//hard event
		if(ifHard)
			hardEvent(&thisOne, &energyTransfer, vox, eleWeight, MCS, ICS, secondary, nSecondary, iseed, mutex, CHYield, COYield, CCYield, 
				CCaYield, CHProducts, COProducts, CCProducts, CCaProducts);
		
		update(&thisOne, stepLength, energyTransfer, theta, phi, crossBound, &voxIndex);
		score(doseCounter, scoringType, absIndex, nVoxels, energyTransfer, stepLength, &thisOne, iseed, MSPR, vox.s1);
	}
}