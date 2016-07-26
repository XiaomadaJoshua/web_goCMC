//#include "Macro.h"

__kernel void initialize(__global float * doseCounter){
	size_t gid = get_global_id(0);
	doseCounter[gid] = 0;
}


__kernel void finalize(__global float * doseBuff, __global float * meanCounter, __global float * uncertaintyCounter, int scoringType, read_only image3d_t voxels, float3 voxSize, ulong nPaths){
	size_t idx = get_global_id(0);
	size_t idy = get_global_id(1);
	size_t idz = get_global_id(2);
	size_t absId = idx + idy*get_global_size(0) + idz*get_global_size(0)*get_global_size(1);
	int nVoxels = get_global_size(0)*get_global_size(1)*get_global_size(2);

	sampler_t voxSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;
	float4 vox = read_imagef(voxels, voxSampler, (float4)(idx, idy, idz, 0.0f));
	float mean = 0.0f, var = 0.0f, std = 0.0f;
	float volume = voxSize.x*voxSize.y*voxSize.z;
	float mass = vox.s2*volume;


	for(int i = 0; i < NDOSECOUNTERS; i++){
		mean += doseBuff[absId + i*nVoxels];
		var += doseBuff[absId + i*nVoxels]*doseBuff[absId + i*nVoxels];
	}
	std = sqrt(NDOSECOUNTERS*var - mean*mean);

	mean = mean/nPaths;
	std = std/nPaths;
	if(scoringType == 0 || scoringType == 1){
		mean = mean*1e3*MEV2JOULES/mass; // 1e3 to convert cm^3 to mm^3
		std = std*1e3*MEV2JOULES/mass;
	}

	meanCounter[absId] = mean;
	uncertaintyCounter[absId] = std;
}

__kernel void tempStore(__global float * doseCounter, __global float * doseBuff){
	size_t idx = get_global_id(0);
	size_t idy = get_global_id(1);
	size_t idz = get_global_id(2);
	size_t absId = idx + idy*get_global_size(0) + idz*get_global_size(0)*get_global_size(1);
	int nVoxels = get_global_size(0)*get_global_size(1)*get_global_size(2);
	for(int i  = 0; i < NDOSECOUNTERS; i++){
		doseBuff[absId + i*nVoxels] += doseCounter[absId + i*nVoxels];
		doseCounter[absId + i*nVoxels] = 0;
	}
}