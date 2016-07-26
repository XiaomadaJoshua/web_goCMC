#define INF 1.0e20
#define WATERDENSITY 1.000 // g/cm^3
#define ME 0.510998928  //electron mass, in MeV
#define MO 14903.3460795634 //oxygen mass in MeV
#define MINELECTRONENERGY 0.1 // MeV
#define TWOPIRE2MENEW 0.08515495201157892 //2pi*r_e^2*m_e*n_{ew}, where r_e in cm, m_e in eV, n_ew = 3.34e23/cm^3
#define XW 36.514 	//radiation length of water, in cm (used in Geant4)
#define PI 3.1415926535897932384626433
#define ES 25 //MeV, ES parameter for secondary particles
#define MAXSTEP 0.2 //in cm
#define MAXENERGYRATIO 0.25 //Max energy decay ratio of initial energy in a step
#define ZERO 1e-6
#define EPSILON 1e-20
#define SECONDARYNUMBERRATIO 6 // ratio of secondary stack size over the number of a batch primary
#define NDOSECOUNTERS 10 // number of dosecounters
#define MC 11174.86339 // carbon mass in MeV	
#define CC 6 // carbon charge number
#define MINCARBONENERGY 1.0 //Min carbon energy to transport in MeV
#define MEV2JOULES 1.6021773e-13
#define EXCITATION 70e-6  //mean exitation enery of water, in MeV (used in secondary stopping power calculation)
#define NE 3.342795e23  // edensity in water, in 1/cm^3
#define RE 2.8179403267e-13 // electron classical radius, in cm
#define MINSECONDENERGY 1.0 // Min secondary energy to transport in MeV
#define OXYGEN 15.9994 // oxygen atomic mass
#define HYDROGEN 1.00794 //hydrogen atomic mass
#define CARBON 12.0107 // carbon atomic mass

#define MIN(a,b) (a > b ? b : a)
#define MIN3(a,b,c) (a > b ? b : a) > c ? c : (a > b ? b : a)
#define ABS(a) a > 0 ? a : -a


#define __ONLYEM__ 0
//	if only EM interaction is turned on, 1--on, 0--off
#define __GPU__ 1
//	if GPU simulation is turned on, 1--on, 0--off
#define __LINUX__ 0
//  if linux is used
#define __DOWNSAMPLE__ 0
//	if downsample CT