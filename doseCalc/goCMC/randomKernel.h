#pragma once
//#include "Random123\threefry.h"

#define M_RAN_INVM32 2.32830643653869628906e-010

float inverseCumulativeNormal(float u){
	float a[4] = { 2.50662823884f, -18.61500062529f, 41.39119773534f, -25.44106049637f };
	float b[4] = { -8.47351093090f, 23.08336743743f, -21.06224101826f, 3.13082909833f };
	float c[9] = { 0.3374754822726147f, 0.9761690190917186f, 0.1607979714918209f, 0.0276438810333863f, 0.0038405729373609f, 0.0003951896511919f,
		0.0000321767881768f, 0.0000002888167364f, 0.0000003960315187f };
	float x = u - 0.5f;
	float r;
	if (x < 0.42f && x > -0.42f)  // Beasley-Springer
	{
		float y = x*x;
		r = x*(((a[3] * y + a[2])*y + a[1])*y + a[0]) / ((((b[3] * y + b[2])*y + b[1])*y + b[0])*y + 1.0f);
	}
	else // Moro
	{
		r = u;
		if (x>0.0f)
			r = 1.0f - u;
		r = log(-log(r));
		r = c[0] + r*(c[1] + r*(c[2] + r*(c[3] + r*(c[4] + r*(c[5] +
			r*(c[6] + r*(c[7] + r*c[8])))))));
		if (x<0.0f)
			r = -r;
	}
	return r;
}

float MTrng(int * iseed){
	int I1 = iseed[0] / 53668;
	iseed[0] = 40014 * (iseed[0] - I1 * 53668) - I1 * 12211;
	if (iseed[0] < 0) iseed[0] = iseed[0] + 2147483563;

	int I2 = iseed[1] / 52774;
	iseed[1] = 40692 * (iseed[1] - I2 * 52774) - I2 * 3791;
	if (iseed[1] < 0) iseed[1] = iseed[1] + 2147483399;

	int IZ = iseed[0] - iseed[1];
	if (IZ < 1) IZ = IZ + 2147483562;
	return (float)(IZ*4.656612873077392578125e-10);
}

float MTGaussian(int * iseed){
	float u1 = MTrng(iseed);
	float u2 = MTrng(iseed);
	return sqrt(-2.0f*log(u1))*cos(2.0f*PI*u2);
}

/*
float MTGaussian(int * iseed){
return inverseCumulativeNormal(MTrng(iseed));
}
*/

float MTExp(int * iseed, float lambda){

	return -log(1 - MTrng(iseed))*lambda;
}













/*
void getUniform(threefry2x32_ctr_t * c, float * random){
	threefry2x32_key_t k = { {0} };

	*c = threefry2x32(*c, k);
	
	random[0] = (int)c->v[0] * M_RAN_INVM32 + 0.5f;
	random[1] = (int)c->v[1] * M_RAN_INVM32 + 0.5f;
}


void getGaussian(threefry2x32_ctr_t * c, float * variates){
	getUniform(c, variates);
	for (unsigned int i = 0; i<2; i++){
		float x = *(variates + i);
		*(variates + i) = inverseCumulativeNormal(x);
	}
}

void getExp(threefry2x32_ctr_t * c, float * variates, float p){
	getUniform(c, variates);
	for (unsigned int i = 0; i<2; i++){
		*(variates + i) = -log(1 - *(variates + i)) / p;
	}
}*/