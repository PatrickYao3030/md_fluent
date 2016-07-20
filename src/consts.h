struct PorousMaterials{
	real thickness;
	real porosity;
	real tortuosity;
	real conductivity;
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