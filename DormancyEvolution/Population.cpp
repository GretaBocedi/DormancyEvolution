#include "Population.h"

std::random_device rd3;
std::mt19937 gen3(rd3());

Population::Population() {
	N = 0;
	Nseeds = 0;
	NseedsB = 0;
	mean_dorm = 0.0;
	std_dorm = 0.0;
}
//---------------------------------------------------------------------------
Population::~Population()
{
}
//---------------------------------------------------------------------------
void Population::reproduce(int fec, int xx, int yy, bool evol, std::bernoulli_distribution mut, std::uniform_real_distribution<> mutation) {
	double genotype;

	Individual* seed;

	for (int i = 0; i < N; i++) {
		genotype = inds[i].dorm;

		for (int j = 0; j < fec; j++) {
			seed = new Individual(genotype, xx, yy);
			if (evol && mut(gen3)) seed->mutate(mutation);
			seeds.push_back(*seed);
			delete seed;
			Nseeds++;			
		}
	}
}
//---------------------------------------------------------------------------
void Population::outPop(bool evol, int r, int g, int x, int y, std::ofstream *pop) {

	*pop << r << "\t" << g << "\t" << x << "\t" << y << "\t" << N << "\t" << NseedsB;

	if (evol) *pop << "\t" << mean_dorm << "\t" << std_dorm;

	*pop << endl;

}