#include "dormancy.h"

int main(void) {

	int p;

	if (experiment == 1) {
		for (double d = min_fix_dorm; d < 1.0; d += 0.05) {

			outPop_header(d);

			init_dorm = d;

			for (int r = 0; r < replicates; r++) {
				Initialise();

				//cout << "init_OK" << endl;

				for (int g = 0; g < gen; g++) {
					cout << "gen = " << g << endl;

					if (fixed_s == false) s = sample_s(rdgen);

					reproduction();
					//cout << "rep_OK" << endl;

					dispersal();
					//cout << "disp_OK" << endl;

					survival(r, g);
					//cout << "surv_OK" << endl;

					cout << "Ntot = " << Ntot << "  NseedsB = " << NseedsB_tot << endl;
				}

				delete_landscape();
			}
			pop.close();
		}
	}
	else {

		outPop_header(0.0);

		for (int r = 0; r < replicates; r++) {
			
			InitialiseB(0);

			cout << "init_OK" << endl;

			for (int g = 0; g < genB; g++) {
				cout << "gen = " << g << endl;

				if (fixed_s == false) s = sample_s(rdgen);

				reproductionB();
				cout << "rep_OK" << endl;

				dispersalB();
				cout << "disp_OK" << endl;

				survivalB(r, g);
				cout << "surv_OK" << endl;

				cout << "Ntot = " << Ntot << "  NseedsB = " << NseedsB_tot << endl;
			}

			InitialiseB(1);

			for (int g = 0; g < gen; g++) {
				cout << "gen = " << g << endl;

				if (fixed_s == false) s = sample_s(rdgen);

				reproduction();
				//cout << "rep_OK" << endl;

				dispersal();
				//cout << "disp_OK" << endl;

				survival(r, g);
				//cout << "surv_OK" << endl;

				cout << "Ntot = " << Ntot << "  NseedsB = " << NseedsB_tot << endl;
			}
			delete_landscape();
		}
		pop.close();
	}

	std::cin >> p;

	return 0;
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

//Initialise
void Initialise(void) {
		for (int i = 0; i < K; i++) {
			land[init_x][init_y].inds.push_back(Individual(init_dorm, init_x, init_y));
			land[init_x][init_y].N++;
		}
}

void InitialiseB(int type) {

	if (type == 0) { //initialise every cell at k
		for (int x = 0; x < 50; x++) {
			for (int y = 0; y < 50; y++) {
				land0[x][y].inds.push_back(Individual(init_dorm, x, y));
				land0[x][y].N++;
			}
		}
	}
	else { //initialise one cell at K with individuals sampled from land0 for invasion experiment
		std::uniform_int_distribution<> coord(0, 49);

		int x, y, z;

		for (int i = 0; i < K; i++) {
			do {
				x = coord(rdgen);
				y = coord(rdgen);
			} while (land0[x][y].N == 0);

			std::uniform_int_distribution<> sample(0, land0[x][y].N - 1);
			z = sample(rdgen);

			land[init_x][init_y].inds.push_back(land0[x][y].inds[z]);
			land[init_x][init_y].N++;
		}
	}
}

//Reproduction
void Individual::mutate(void) {
	dorm += mutation(rdgen);
	if (dorm < 0.0) dorm = 0.0;
	if (dorm > 1.0) dorm = 1.0;
}

void Population::reproduce(int xx, int yy) {
	double genotype;

	for (int i = 0; i < N; i++) {
		genotype = inds[i].dorm;

		for (int j = 0; j < s; j++) {
			seeds.push_back(Individual(genotype, xx, yy));
			Nseeds++;
			if (mutate(rdgen)) {
				seeds[Nseeds - 1].mutate();
			}
		}
	}
}

void reproduction(void) {

	for (int x = 0; x < x_max; x++) {
		for (int y = 0; y < y_max; y++) {
			if (land[x][y].N > 0) {
				land[x][y].reproduce(x, y);
			}
			land[x][y].N = 0;
			land[x][y].inds.clear();
		}
	}
}

void reproductionB(void) { //reproduction on land0

	for (int x = 0; x < 50; x++) {
		for (int y = 0; y < 50; y++) {
			if (land0[x][y].N > 0) {
				land0[x][y].reproduce(x, y);
			}
			land0[x][y].N = 0;
			land0[x][y].inds.clear();
		}
	}
}


void Individual::disp(int x_m, int y_m) {
	int new_x, new_y;
	do {
		new_x = x + step(rdgen);
		new_y = y + step(rdgen);
	} while (new_x == x && new_y == y);
	x = new_x;
	y = new_y;

	if (new_x < 0 || new_x > x_m - 1 || new_y < 0 || new_y > y_m - 1) alive = false;
}

void dispersal(void) {
	int new_x, new_y;

	for (int x = 0; x < x_max; x++) {
		for (int y = 0; y < y_max; y++) {
			if (land[x][y].Nseeds > 0) {
				for (int i = 0; i < land[x][y].Nseeds; i++) {
					if (disperse(rdgen)) {
						land[x][y].seeds[i].disp(x_max, y_max);
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

void dispersalB(void) {
	int new_x, new_y;

	for (int x = 0; x < 50; x++) {
		for (int y = 0; y < 50; y++) {
			if (land0[x][y].Nseeds > 0) {
				for (int i = 0; i < land0[x][y].Nseeds; i++) {
					if (disperse(rdgen)) {
						land0[x][y].seeds[i].disp(50,50);
						if (land0[x][y].seeds[i].alive) {
							new_x = land0[x][y].seeds[i].x;
							new_y = land0[x][y].seeds[i].y;
							land0[new_x][new_y].seedsB.push_back(land0[x][y].seeds[i]);
							land0[new_x][new_y].NseedsB++;

						}
					}
					else {
						land0[x][y].seedsB.push_back(land0[x][y].seeds[i]);
						land0[x][y].NseedsB++;

					}
				}
				land0[x][y].seeds.clear();
				land0[x][y].Nseeds = 0;

			}
		}
	}
}

//germination
bool Individual::germinate(void) {

	std::bernoulli_distribution germ(1.0 - dorm);

	return germ(rdgen);
}

//germination and survival
void survival(int rr, int gg) {

	int Ns;
	int Nseedlings;

	Ntot = 0;
	NseedsB_tot = 0;

	for (int x = 0; x < x_max; x++) {
		for (int y = 0; y < y_max; y++) {
			if (land[x][y].NseedsB > 0) {
				//cout << land[x][y].NseedsB << endl;
				//int pp = land[x][y].seedsB.size();
				//cout << pp << endl;

				Ns = 0;
				Nseedlings = 0;
				for (int i = 0; i < land[x][y].NseedsB; i++) {
					if (land[x][y].seedsB[i].germinate()) {
						land[x][y].seeds.push_back(land[x][y].seedsB[i]); //seeds vector here used for seedlings
						Nseedlings++;
					}
					else {
						if (dorm_surv(rdgen)) {
							land[x][y].seedsB_tmp.push_back(land[x][y].seedsB[i]);
							Ns++;
						}
					}
				}

				if (Nseedlings > 0) {
					if (Nseedlings > K) {
						random_shuffle(land[x][y].seeds.begin(), land[x][y].seeds.end());
						for (int i = 0; i < K; i++) {
							land[x][y].inds.push_back(land[x][y].seeds[i]);
							land[x][y].N++;
						}
					}
					else {
						for (int i = 0; i < Nseedlings; i++) {
							land[x][y].inds.push_back(land[x][y].seeds[i]);
							land[x][y].N++;
						}
					}
				}

				land[x][y].seeds.clear();
				land[x][y].seedsB.clear();
				land[x][y].seedsB = land[x][y].seedsB_tmp;
				land[x][y].seedsB_tmp.clear();
				land[x][y].NseedsB = Ns;

				Ntot += land[x][y].N;
				NseedsB_tot += land[x][y].NseedsB;

				if (beta > 0.0) {
					land[x][y].mean_dorm = 0.0;
					land[x][y].std_dorm = 0.0;
					for (int i = 0; i < land[x][y].N; i++) {
						land[x][y].mean_dorm += land[x][y].inds[i].dorm; //sum of genotypes
						land[x][y].std_dorm += sqrt(land[x][y].inds[i].dorm); //sum of squares
					}
					land[x][y].std_dorm = land[x][y].std_dorm - sqrt(land[x][y].mean_dorm) / land[x][y].N; //standard deviation
					land[x][y].mean_dorm /= land[x][y].N; //mean genotype
				}

				if (gg % out_interval == 0) land[x][y].outPop(rr, gg, x, y);

			}
		}
	}

}

//germination and survival
void survivalB(int rr, int gg) {

	int Ns;
	int Nseedlings;

	Ntot = 0;
	NseedsB_tot = 0;

	for (int x = 0; x < 50; x++) {
		for (int y = 0; y < 50; y++) {
			if (land0[x][y].NseedsB > 0) {
				//cout << land[x][y].NseedsB << endl;
				//int pp = land[x][y].seedsB.size();
				//cout << pp << endl;

				Ns = 0;
				Nseedlings = 0;
				for (int i = 0; i < land0[x][y].NseedsB; i++) {
					if (land0[x][y].seedsB[i].germinate()) {
						land0[x][y].seeds.push_back(land0[x][y].seedsB[i]); //seeds vector here used for seedlings
						Nseedlings++;
					}
					else {
						if (dorm_surv(rdgen)) {
							land0[x][y].seedsB_tmp.push_back(land0[x][y].seedsB[i]);
							Ns++;
						}
					}
				}

				if (Nseedlings > 0) {
					if (Nseedlings > K) {
						random_shuffle(land0[x][y].seeds.begin(), land0[x][y].seeds.end());
						for (int i = 0; i < K; i++) {
							land0[x][y].inds.push_back(land0[x][y].seeds[i]);
							land0[x][y].N++;
						}
					}
					else {
						for (int i = 0; i < Nseedlings; i++) {
							land0[x][y].inds.push_back(land0[x][y].seeds[i]);
							land0[x][y].N++;
						}
					}
				}

				land0[x][y].seeds.clear();
				land0[x][y].seedsB.clear();
				land0[x][y].seedsB = land0[x][y].seedsB_tmp;
				land0[x][y].seedsB_tmp.clear();
				land0[x][y].NseedsB = Ns;

				Ntot += land0[x][y].N;
				NseedsB_tot += land0[x][y].NseedsB;

				if (beta > 0.0) {
					land0[x][y].mean_dorm = 0.0;
					land0[x][y].std_dorm = 0.0;
					for (int i = 0; i < land0[x][y].N; i++) {
						land0[x][y].mean_dorm += land0[x][y].inds[i].dorm; //sum of genotypes
						land0[x][y].std_dorm += sqrt(land0[x][y].inds[i].dorm); //sum of squares
					}
					land0[x][y].std_dorm = land0[x][y].std_dorm - sqrt(land0[x][y].mean_dorm) / land0[x][y].N; //standard deviation
					land0[x][y].mean_dorm /= land0[x][y].N; //mean genotype
				}
				if (gg % out_interval == 0) land0[x][y].outPop(rr, gg, x, y);
			}

		}
	}

}

void outPop_header(double d) {
	string name;

	if (experiment == 1)  name = "Sim" + Int2Str(SimNr) + "_Pop_d" + Float2Str(d) + ".txt";
	else name = "Sim" + Int2Str(SimNr) + "_Pop.txt";

	pop.open(name.c_str());

	if (beta > 0.0) { //dormancy evolution
		pop << "rep\tgen\tx\ty\tN\tNseedsB\tmean_dorm\tstd_dorm" << endl;
	}
	else { //fix dormancy
		pop << "rep\tgen\tx\ty\tN\tNseedsB" << endl;
	}

}

void Population::outPop(int r, int g, int x, int y) {

	pop << r << "\t" << g << "\t" << x << "\t" << y << "\t" << N << "\t" << NseedsB;

	if (beta > 0.0) pop << "\t" << mean_dorm << "\t" << std_dorm;

	pop << endl;
	
}

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