
template<class TETRA_MESH_TYPE>
class TetraStats
{
typedef typename TETRA_MESH_TYPE::TetraType TetraType;

TETRA_MESH_TYPE * Tetra;
TetraType* T;

public:
double volume;
double ratio;

TetraStats(){T=0;}

~TetraStats(){}

static double ComputeVolume(TETRA_MESH_TYPE *tm)
{

	double total=0.f;
	TETRA_MESH_TYPE::TetraIterator ti;
	for (ti=tm->tetra.begin();ti<tm->tetra.end();ti++)
	{
		if (!ti->IsD())
			total+=ti->ComputeVolume();
	}
	return total;
}

static double ComputeRatioMedia(TETRA_MESH_TYPE *tm)
{
	double total=0.f;
	TETRA_MESH_TYPE::TetraIterator ti;
	int i=0;
	for (ti=tm->tetra.begin();ti<tm->tetra.end();ti++)
	{
		if (!ti->IsD())
		{
			total+=ti->AspectRatio();
			i++;
		}
	}
	return (total/i);
}

void SetTetraMesh(TETRA_MESH_TYPE* T)
{
	Tetra=T;
}

void Update()
{
	ratio=ComputeRatioMedia(Tetra);
	volume=ComputeVolume(Tetra);
}

void SetTetraInfo(TetraType *Te)
{
	if (T!=0)
		T->ClearS();
	T=Te;
}

void ClearTetraInfo(TetraType *Te)
{
	T=0;
}

TetraType * TCurrent()
{
	return T;
}
};
