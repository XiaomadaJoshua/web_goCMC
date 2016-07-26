//#include "randomKernel.h"
//#include "Macro.h"

typedef struct __attribute__ ((aligned)) ParticleStatus{
	float3 pos, dir;
	float energy, maxSigma, mass, charge;
	int ifPrimary;
}PS;


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

inline float gamma(PS * particle) {
	return (particle->energy + particle->mass) / particle->mass;
}
inline float beta(PS * particle) {
	return sqrt(1 - particle->mass*particle->mass/((particle->energy + particle->mass)*(particle->energy + particle->mass)));
}

inline float maxDeltaElectronEnergy(PS * particle) {
	return (2 * ME*beta(particle)*beta(particle)*gamma(particle)*gamma(particle)) / (1 + 2 * gamma(particle)*ME / particle->mass + ME*ME / (particle->mass*particle->mass));
}

inline float momentumSquare(PS * particle) {
	return particle->energy*particle->energy + 2 * particle->energy*particle->mass;
}


float radiationLength(float4 vox)  
//	calculate the radiation length ratio \rho_wX_w/(\rhoX_0(\rho)
{
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


inline float RSP(PS * particle, float4 vox){
	float Tup = min(MINELECTRONENERGY, maxDeltaElectronEnergy(particle));
	float beta2 = beta(particle)*beta(particle);
	float gamma2 = gamma(particle)*gamma(particle);
	return vox.s3*TWOPIRE2MENEW*particle->charge*particle->charge/beta2*(log(2.0*ME*beta2*gamma2*Tup/vox.s0/vox.s0) 
															- beta2*(1.0 + Tup/maxDeltaElectronEnergy(particle)));
}

float MSPR(PS * particle, float4 vox){
	float Temax = maxDeltaElectronEnergy(particle);
	float Tup = min(MINELECTRONENERGY, Temax);
	float beta2 = beta(particle)*beta(particle);
	float gamma2 = gamma(particle)*gamma(particle);
	return (log(2.0*ME*beta2*gamma2*Tup/vox.s0/vox.s0) - beta2*(1.0 + Tup/Temax))/(log(2.0*ME*beta2*gamma2*Tup/EXCITATION/EXCITATION) - beta2*(1.0 + Tup/Temax));
}



void score(global float * doseCounter, int scoringType, int absIndex, int nVoxels, float energyTransfer, float stepLength, PS * thisOne, int * iseed, float4 vox){
	// choose a dose counter
	int doseCounterId = convert_int_rtn(MTrng(iseed)*NDOSECOUNTERS);
	volatile global float * counter = &doseCounter[absIndex + doseCounterId * nVoxels];
	// find scoring type: 0 is DOSE2MEDIUM, 1 is DOSE2WATTER, 2 is FLUENCE, 3 is LET
	switch(scoringType){
	case(0):
		atomicAdd(counter, energyTransfer);
		break;
	case(1):
		energyTransfer = energyTransfer/MSPR(thisOne, vox);
		atomicAdd(counter, energyTransfer);
		break;
	}
}



void ionization(PS * thisOne, int * iseed, float * energyTransfer){
	

	float E = thisOne->energy + thisOne->mass;
	float Te;
	float rand;
	while(true){
		rand = MTrng(iseed);
		Te = MINELECTRONENERGY*maxDeltaElectronEnergy(thisOne)
			/((1-rand)*maxDeltaElectronEnergy(thisOne)+rand*MINELECTRONENERGY);
		if(MTrng(iseed) < 1-beta(thisOne)*beta(thisOne)*Te/maxDeltaElectronEnergy(thisOne)+Te*Te/(2*E*E))
			break;
	}

	*energyTransfer += Te;

}


inline float sigmaIon(PS * particle, float eDensity){
	float Temax = maxDeltaElectronEnergy(particle);

	if(Temax <= MINELECTRONENERGY)
		return 0.0;
	else{
		float beta2 = beta(particle)*beta(particle);
		float E = particle->energy + particle->mass;
		return eDensity*TWOPIRE2MENEW*particle->charge*particle->charge/beta2*((-1.0/Temax - beta2/Temax*log(Temax) + Temax/2.0/(E*E)) 
                                    - (-1.0/MINELECTRONENERGY - beta2/Temax*log(MINELECTRONENERGY) + MINELECTRONENERGY/2/(E*E)) );
	}
}


inline float energyInOneStep(float4 vox, PS * particle, float stepLength) {
	return stepLength*RSP(particle, vox);
}



void hardEvent(PS * thisOne, float *energyTransfer, int * iseed, float4 vox){

	float sigma = sigmaIon(thisOne, vox.s3);
	float rand = MTrng(iseed);
	rand *= thisOne->maxSigma > sigma ? thisOne->maxSigma : sigma;
	if(rand < sigma){
		ionization(thisOne, iseed, energyTransfer);
		return;
	}
	else 
		return;
}


__kernel void propagate(__global PS * particle, __global float * doseCounter, int scoringType, __read_only image3d_t voxels, float3 voxSize,int randSeed, __read_only image2d_t RSPW){

	size_t gid = get_global_id(0);
	PS thisOne = particle[gid];
//	if(isnan(thisOne.energy) || isnan(thisOne.mass)){
//		//printf("simulate secondary particle\n");
//		//printf("simulated secondary status: energy %f, pos %v3f, dir %v3f, ifPrimary %d, mass %f, charge %f, maxSigma %f\n", thisOne.energy, thisOne.pos, thisOne.dir, thisOne.ifPrimary, thisOne.mass, thisOne.charge, thisOne.maxSigma);
//	}

	int3 phantomSize = (int3)(get_image_width(voxels), get_image_height(voxels), get_image_depth(voxels));
	int nVoxels = phantomSize.x * phantomSize.y * phantomSize.z;

	float stepLength, thisMaxStep, step2bound, energyTransfer, sigma1, sigma2, sigma, sampledStep, variance, theta0, theta, phi, residualRange;
	int3 voxIndex;
	int absIndex, absIndex2, crossBound;
	int iseed[2];
	iseed[0] = gid;
	iseed[1] = randSeed;
	MTrng(iseed);
	sampler_t dataSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;
	sampler_t voxSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;
	float4 vox;
	bool ifHard;
	int step = 0;
	absIndex = -1;
	absIndex2 = -1;

	voxIndex = convert_int3_rtn(thisOne.pos / voxSize + convert_float3(phantomSize) / 2.0f);

	while (ifInsidePhantom(voxIndex, phantomSize)){
		absIndex = voxIndex.x + voxIndex.y*phantomSize.x + voxIndex.z*phantomSize.x*phantomSize.y;
		vox = read_imagef(voxels, voxSampler, (float4)(convert_float3(voxIndex), 0.0f));		
		residualRange = read_imagef(RSPW, dataSampler, (float2)(MC*thisOne.energy/thisOne.mass - 1.0f + 0.5f, 0.5f)).s2 / vox.s2;
		residualRange *= thisOne.mass*CC*CC/(MC*thisOne.charge*thisOne.charge);
		step2bound = step2VoxBoundary(thisOne.pos, thisOne.dir, voxSize, &crossBound, phantomSize, voxIndex);
		if (residualRange <= step2bound || thisOne.energy <= MINSECONDENERGY){
			stepLength = thisOne.energy / RSP(&thisOne, vox);
			score(doseCounter, scoringType, absIndex, nVoxels, thisOne.energy, stepLength, &thisOne, iseed, vox);
			return;
		}

		// rescale maxStep to let energy transferred in one step < MAXENERGYRATIO
		thisMaxStep = MAXSTEP;
		energyTransfer = energyInOneStep(vox, &thisOne, thisMaxStep);
		if (energyTransfer > MAXENERGYRATIO*thisOne.energy)
			thisMaxStep = MAXENERGYRATIO*thisOne.energy / RSP(&thisOne, vox);


		// get linear macro cross section to sample a step
		energyTransfer = energyInOneStep(vox, &thisOne, thisMaxStep);
		thisOne.energy -= energyTransfer;
		sigma1 = sigmaIon(&thisOne, vox.s3);
		thisOne.energy += energyTransfer;
		sigma2 = sigmaIon(&thisOne, vox.s3);
		sigma = sigma1 > sigma2 ? sigma1 : sigma2;


		// sample one step
		if(sigma < ZERO)
			sampledStep = INF;
		else
			sampledStep = -log(MTrng(iseed)) / sigma;
		stepLength = sampledStep < thisMaxStep ? sampledStep : thisMaxStep;
		if (stepLength >= step2bound){
			ifHard = false;
			stepLength = step2bound;
		}
		else{
			ifHard = true;
			crossBound = 0;
		}

//		printf("sample step %f, if crossBound %d, steplength %f\n", sampledStep, crossBound, stepLength);

		// get energy transferred (plus energy straggling) in this sampled step
		energyTransfer = energyInOneStep(vox, &thisOne, stepLength);
		variance = TWOPIRE2MENEW*vox.s3*stepLength*thisOne.charge*thisOne.charge*fmin(MINELECTRONENERGY, maxDeltaElectronEnergy(&thisOne))*(1.0f / (beta(&thisOne)*beta(&thisOne)) - 0.5f);
		energyTransfer += MTGaussian(iseed) * sqrt(variance);
		energyTransfer = energyTransfer > 0 ? energyTransfer : -energyTransfer;
		energyTransfer = energyTransfer < thisOne.energy ? energyTransfer : thisOne.energy;


		// deflection
		theta0 = ES*thisOne.charge*sqrt(stepLength/radiationLength(vox))/beta(&thisOne)/sqrt(momentumSquare(&thisOne));
//		if(isnan(theta0))
//			printf("stepLength = %f, beta = %f, momentumSquare = %f, energy = %f, mass = %f\n", stepLength, beta(&thisOne), momentumSquare(&thisOne), thisOne.energy, thisOne.mass);
		theta = MTGaussian(iseed) * theta0;
		phi = 2.0f*PI*MTrng(iseed);
		thisOne.maxSigma = sigma;

		//hard event
		if(ifHard)
			hardEvent(&thisOne, &energyTransfer, iseed, vox);

		update(&thisOne, stepLength, energyTransfer, theta, phi, crossBound, &voxIndex);
		score(doseCounter, scoringType, absIndex, nVoxels, energyTransfer, stepLength, &thisOne, iseed, vox);

	}
	
}