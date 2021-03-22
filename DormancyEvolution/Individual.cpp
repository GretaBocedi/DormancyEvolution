#include "Individual.h"

std::random_device rd2;
std::mt19937 gen(rd2());

std::uniform_int_distribution<> step(-1, 1);

Individual::Individual(double d, int xx, int yy) {
	alive = true;
	dorm = d;
	x = xx;
	y = yy;
}

//---------------------------------------------------------------------------
Individual::~Individual()
{
}
//---------------------------------------------------------------------------
void Individual::mutate(std::uniform_real_distribution<> mut) {
	dorm += mut(gen);
	if (dorm < 0.0) dorm = 0.0;
	if (dorm > 1.0) dorm = 1.0;
}
//---------------------------------------------------------------------------
//Nearest neighbour dispersal
void Individual::disp(int x_m, int y_m) {
	int new_x = 0;
	int new_y = 0;

	do {
		new_x = x + step(gen);
		new_y = y + step(gen);
	} while (new_x == x && new_y == y);

	//wrap y-axis
	if (new_x > -1 && new_x < x_m) {
		if (new_y > y_m - 1 || new_y < 0) new_y = y_m - std::abs(new_y);
	}
	
	if (new_x < 0 || new_x > x_m - 1) alive = false;

	x = new_x;
	y = new_y;
}
//---------------------------------------------------------------------------
//germination
bool Individual::germinate(void) {

	std::bernoulli_distribution germ(1.0 - dorm);

	return germ(gen);
}