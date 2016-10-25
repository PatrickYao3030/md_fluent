#define MAXRECLINE 99999
#define MAXCELLNUM 9999
#define EPS 5.0e-4
#define PVDF 0
#define PTFE 1
#define PP 2
#define PES 3
#define SOLID 0
#define GAS 1
// check the zones' id in FLUENT
#define id_domain 1
#define id_FeedFluid 5
#define id_PermFluid 4
#define id_FeedInterface 13
#define id_PermInterface 14
/* the customized messages displayed in command-line window of FLUENT have 3 levels
   Lvl.0 stands for serious problems needed care
	 Lvl.1            warnings
	 Lvl.2            minor messages
	 Lvl.3            debug messages
	 these messages are controled to be display by the switch of id_message
*/
#define id_message 2

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
struct Flux{
	real mass;
	real heat;
};
struct CellInfos{
	int seq;
	int index;
	real centroid[ND_ND];
	real temperature;
	real volume;
	real area;
	struct Flux flux;
	struct Components massfraction;
};
struct MessageInfos{
	int flag;
	char content[79];
};
