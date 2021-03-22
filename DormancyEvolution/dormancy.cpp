#include "dormancy.h"


#if CLUSTER
int main(int argc, char* argv[])
{
	// Get the current directory.
	char* buffer = getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "/"; //Current directory path
	dirOut = dir + "Outputs/"; //Outpus folder path

	SimNr = std::atoi(argv[1]);
	replicates = std::atoi(argv[2]);
	years = std::atoi(argv[3]);
	expansion_start = std::atoi(argv[4]);
	init_dorm = std::atof(argv[5]);
	dorm_evol = std::atoi(argv[6]);
	Allee = std::atoi(argv[7]);
	stochasticity = std::atoi(argv[8]);
	env_std = std::atof(argv[9]);
	env_ac = std::atof(argv[10]);


	RunModel();

	cout << "Simulation completed" << endl;

	return 0;
}
#else
int _tmain(int argc, _TCHAR* argv[])
{

	// Get the current directory.
	char* buffer = _getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "\\"; //Current directory path
	dirOut = dir + "Outputs\\"; //Outpus folder path

	RunModel();

	std::cout << "Simulation completed" << endl;

	return 0;
}
#endif

void RunModel(void) {

	int p;

	outPop_header(d);

	for (int r = 0; r < replicates; r++) {
		cout << "------------------------------" << endl;
		cout << "REP = " << r << endl;

		Initialise();
		//cout << "init_OK" << endl;

		for (int g = 0; g < years; g++) {
			if (g % 50 == 0) cout << "gen = " << g << endl;

			if (g < expansion_start) x_m = x_max0;
			else x_m = x_max;

			if (stochasticity) env_stoch();

			reproduction();
			//cout << "rep_OK" << endl;

			dispersal();
			//cout << "disp_OK" << endl;

			survival(r, g);
			//cout << "surv_OK" << endl;
			//cout << "Ntot = " << Ntot << "  NseedsB = " << NseedsB_tot << endl;

			if (Ntot < 1) break;
		}
		delete_landscape();
	}
	pop.close();

}

//---------------------------------------------------------------------------
const string Int2Str(const int x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
//---------------------------------------------------------------------------
const string Float2Str(const double x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
//---------------------------------------------------------------------------
//Initialise
void Initialise(void) {

	eps = 0.0;

	for (int i = 0; i < (int)K; i++) {
		land[init_x][init_y].inds.push_back(Individual(init_dorm, init_x, init_y));
		land[init_x][init_y].N++;
	}
}
//---------------------------------------------------------------------------
void env_stoch(void)
{
	eps = eps * env_ac + normEnv(rdgen) * sqrt(1.0 - env_ac*env_ac);
	fec = s * (1.0 + eps);
	if (fec < 0.0) fec = 0.0;
	if (fec > gamma) fec = gamma;
}
//---------------------------------------------------------------------------
//Reproduction
void reproduction(void) {

	int minN = 0;

	if (Allee) minN = 1;

	for (int x = 0; x < x_m; x++) {
		for (int y = 0; y < y_max; y++) {
			if (land[x][y].N > minN) {
				land[x][y].reproduce((int)fec, x, y, dorm_evol, mutate, mutation);
			}
			land[x][y].N = 0; //adult plants die
			land[x][y].inds.clear();
		}
	}
}
//---------------------------------------------------------------------------
void dispersal(void) {
	int new_x, new_y;

	for (int x = 0; x < x_m; x++) {
		for (int y = 0; y < y_max; y++) {
			if (land[x][y].Nseeds > 0) {

				for (int i = 0; i < land[x][y].Nseeds; i++) {
					if (disperse(rdgen)) {
						land[x][y].seeds[i].disp(x_m, y_max);
						if (land[x][y].seeds[i].alive) {
							new_x = land[x][y].seeds[i].x;
							new_y = land[x][y].seeds[i].y;
							land[new_x][new_y].seedsB.push_back(land[x][y].seeds[i]);
							land[new_x][new_y].NseedsB++;						
						}						
					}
					else {
						land[x][y].seedsB.push_back(land[x][y].seeds[i]);
						land[x][y].NseedsB++;						
					}
				}
				land[x][y].seeds.clear();
				land[x][y].Nseeds = 0;
			}
		}
	}
}
//---------------------------------------------------------------------------
//germination and survival
void survival(int rr, int gg) {

	int Nseed = 0;
	int Nseedlings = 0;
	double ps = 0.0; //survival probability

	vector<Individual>::iterator iter;

	Ntot = 0;
	NseedsB_tot = 0;

	for (int x = 0; x < x_m; x++) {
		for (int y = 0; y < y_max; y++) {
			if (land[x][y].NseedsB > 0) {

				Nseed = 0;
				Nseedlings = 0;
				for (int i = 0; i < land[x][y].NseedsB; i++) {
					if (land[x][y].seedsB[i].germinate()) {
						land[x][y].seeds.push_back(land[x][y].seedsB[i]); //seeds vector here used for seedlings
						Nseedlings++;
					}
					else {
						if (dorm_surv(rdgen)) {
							land[x][y].seedsB_tmp.push_back(land[x][y].seedsB[i]); //seeds that do not germoinate and survive
							Nseed++;
						}
					}
				}

				//density-dependent seedling survival
				if (Nseedlings > 0) {
					ps = std::fmin(K / (double)Nseedlings, 1.0);
					std::bernoulli_distribution survive(ps);

					for (iter = land[x][y].seeds.begin(); iter != land[x][y].seeds.end(); iter++) {
						if (survive(rdgen)) {
							land[x][y].inds.push_back(*iter);
							land[x][y].N++;
						}
					}
				}

				land[x][y].seeds.clear();
				land[x][y].seedsB.clear();
				land[x][y].seedsB = land[x][y].seedsB_tmp;
				land[x][y].seedsB_tmp.clear();
				land[x][y].NseedsB = Nseed;

				Ntot += land[x][y].N;
				NseedsB_tot += land[x][y].NseedsB;

				if (dorm_evol) {
					land[x][y].mean_dorm = 0.0;
					land[x][y].std_dorm = 0.0;
					for (int i = 0; i < land[x][y].N; i++) {
						land[x][y].mean_dorm += land[x][y].inds[i].dorm; //sum of genotypes
						land[x][y].std_dorm += sqrt(land[x][y].inds[i].dorm); //sum of squares
					}
					land[x][y].std_dorm = land[x][y].std_dorm - sqrt(land[x][y].mean_dorm) / land[x][y].N; //standard deviation
					land[x][y].mean_dorm /= land[x][y].N; //mean genotype
				}

				if (gg % out_interval == 0) land[x][y].outPop(dorm_evol, rr, gg, x, y, &pop);
			}
		}
	}

}

//---------------------------------------------------------------------------
void outPop_header(double d) {
	string name;

	name = "Sim" + Int2Str(SimNr) + "_Pop.txt";

	pop.open(name.c_str());

	if (dorm_evol) { //dormancy evolution
		pop << "rep\tgen\tx\ty\tN\tNseedsB\tmean_dorm\tstd_dorm" << endl;
	}
	else { //fix dormancy
		pop << "rep\tgen\tx\ty\tN\tNseedsB" << endl;
	}
}
//---------------------------------------------------------------------------
void delete_landscape(void) {
	for (int x = 0; x < x_max; x++) {
		for (int y = 0; y < y_max; y++) {
			land[x][y].N = 0;
			land[x][y].Nseeds = 0;
			land[x][y].NseedsB = 0;
			land[x][y].mean_dorm = 0.0;
			land[x][y].std_dorm = 0.0;

			land[x][y].inds.clear();
			land[x][y].seeds.clear();
			land[x][y].seedsB.clear();
			land[x][y].seedsB_tmp.clear();
		}
	}
}