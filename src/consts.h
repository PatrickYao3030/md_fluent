#define MAXCELLNUM 9999
#define EPS 5.0e-4
#define PVDF 0
#define PTFE 1
#define PP 2
#define PES 3
#define SOLID 0
#define GAS 1

struct PorousMaterials{
	real thickness;
	real porosity;
	real tortuosity;
	real conductivity;
	real MDcoeff;
};
struct Components{
	real water;
	real salts;
};
struct CellInfos{
	int seq;
	int index;
	real centroid[ND_ND];
	real temperature;
	struct Components massfraction;
};

