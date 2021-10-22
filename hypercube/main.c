#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <windows.h>

// width and height of screen
#define ww 100
#define wh 50

void clr(CHAR_INFO* d)
{
	for (int i = 0; i < ww * wh; ++i)
	{
		d[i].Attributes = FOREGROUND_GREEN;
		d[i].Char.UnicodeChar = ' ';
	}
}

void set(CHAR_INFO* d, COORD pt, char c)
{
	d[pt.Y * ww + pt.X].Char.UnicodeChar = c;
}

char getp(CHAR_INFO* d, COORD* pts, double err)
{
	if (abs(pts[0].Y - pts[2].Y) < 2)
	{
		if (err > 0.5)
		{
			return '-';
		}
		return '_';
	}
	
	if (abs(pts[0].X - pts[2].X) < 2 &&
		(pts[0].X >= pts[2].X || pts[1].X != pts[2].X) &&
		(pts[0].X <= pts[2].X || pts[1].X != pts[0].X))
	{
		return '|';
	}

	int mX = pts[0].Y < pts[2].Y ? pts[0].X : pts[2].X;
	return mX < pts[1].X ? '\\' : '/';\
}

void ln(CHAR_INFO* d, COORD a, COORD b)
{
	set(d, a, '@');
	set(d, b, '@');

	int dx = abs(b.X - a.X), sx = a.X < b.X ? 1 : -1;
	int dy = abs(b.Y - a.Y), sy = a.Y < b.Y ? 1 : -1;
	int err = (dx > dy ? dx : -dy) / 2, e2;

	COORD pts[3];
	double ers[3];

	for (int i = 0; i < 3; ++i)
	{
		pts[i] = a;
		ers[i] = ((double)err - dx) / ((double)dy - dx);
		ers[i] = sy == 1 ? 1.0f - ers[i] : ers[i];

		if (a.X == b.X && a.Y == b.Y) {
			return;
		}

		e2 = err;
		if (e2 > -dx) { err -= dy; a.X += sx; }
		if (e2 < dy) { err += dx; a.Y += sy; }
	}

	for (;;)
	{
		set(d, pts[1], getp(d, pts, ers[1]));

		pts[0] = pts[1];
		pts[1] = pts[2];
		pts[2] = a;

		ers[0] = ers[1];
		ers[1] = ers[2];
		ers[2] = ((double)err - dx) / ((double)dy - dx);
		ers[2] = sy == 1 ? 1.0f - ers[2] : ers[2];
		
		if (a.X == b.X && a.Y == b.Y) {
			break;
		}

		e2 = err;
		if (e2 > -dx) { err -= dy; a.X += sx; }
		if (e2 < dy) { err += dx; a.Y += sy; }
	}

	// add the final point
	set(d, pts[1], getp(d, pts, ers[1]));
}

// hypercube vertices in 4D
double V4[16][4] =
{
	{-1, -1, -1, -1},
	{ 1, -1, -1, -1},
	{-1,  1, -1, -1},
	{ 1,  1, -1, -1},
	{-1, -1,  1, -1},
	{ 1, -1,  1, -1},
	{-1,  1,  1, -1},
	{ 1,  1,  1, -1},
	{-1, -1, -1,  1},
	{ 1, -1, -1,  1},
	{-1,  1, -1,  1},
	{ 1,  1, -1,  1},
	{-1, -1,  1,  1},
	{ 1, -1,  1,  1},
	{-1,  1,  1,  1},
	{ 1,  1,  1,  1},
};

// store the vertices once they have been projected to 3D
double V3[16][3];

// final 2D projection
double V2[16][2];

// the indices for each line
int indices[32][2] =
{
	// cube #1
	{0, 1},
	{0, 2},
	{0, 4},
	{1, 3},
	{1, 5},
	{2, 3},
	{2, 6},
	{3, 7},
	{4, 5},
	{4, 6},
	{5, 7},
	{6, 7},

	// in-between lines
	{0,	8},
	{1,	9},
	{2,	10},
	{3,	11},
	{4,	12},
	{5,	13},
	{6,	14},
	{7,	15},

	// cube #2
	{8, 9},
	{8, 10},
	{8, 12},
	{9, 11},
	{9, 13},
	{10, 11},
	{10, 14},
	{11, 15},
	{12, 13},
	{12, 14},
	{13, 15},
	{14, 15},
};

double dot4(const double* V, const double* U)
{
	return (V[0] * U[0]) + (V[1] * U[1]) + (V[2] * U[2]) + (V[3] * U[3]);
}

double norm4(const double* V)
{
	return sqrt(dot4(V, V));
}

// cross4 computes the four-dimensional cross product of the three vectors
// U, V and W, in that order.
// returns the resulting four-vector.
void cross4(double* result, const double* U, const double* V, const double* W)
{
	// intermediate values
	double A, B, C, D, E, F;       

	// calculate intermediate values
	A = (V[0] * W[1]) - (V[1] * W[0]);
	B = (V[0] * W[2]) - (V[2] * W[0]);
	C = (V[0] * W[3]) - (V[3] * W[0]);
	D = (V[1] * W[2]) - (V[2] * W[1]);
	E = (V[1] * W[3]) - (V[3] * W[1]);
	F = (V[2] * W[3]) - (V[3] * W[2]);

	// calculate the result-vector components
	result[0] = (U[1] * F) - (U[2] * E) + (U[3] * D);
	result[1] = -(U[0] * F) + (U[2] * C) - (U[3] * B);
	result[2] = (U[0] * E) - (U[1] * C) + (U[3] * A);
	result[3] = -(U[0] * D) + (U[1] * B) - (U[2] * A);
}

void vecSub4(double* result, const double* a, const double* b)
{
	result[0] = a[0] - b[0];
	result[1] = a[1] - b[1];
	result[2] = a[2] - b[2];
	result[3] = a[3] - b[3];
}

void vecScale4(double* vec, double m)
{
	vec[0] *= m;
	vec[1] *= m;
	vec[2] *= m;
	vec[3] *= m;
}

void matVecMul4(double* result, const double* mat, const double* vec)
{
	for (int row = 0; row < 4; ++row)
	{
		result[row] = 0;
		for (int col = 0; col < 4; ++col)
		{
			result[row] += mat[col * 4 + row] * vec[col];
		}
	}
}

// creates a rotation matrix for the XW plane
// T is the angle in radians
void rotXW4(double* result, double T)
{
	// column vectors
	double* Wa = result + 4 * 0;
	double* Wb = result + 4 * 1;
	double* Wc = result + 4 * 2;
	double* Wd = result + 4 * 3;

	Wa[0] = cos(T);
	Wa[1] = 0;
	Wa[2] = 0;
	Wa[3] = -sin(T);

	Wb[0] = 0;
	Wb[1] = 1;
	Wb[2] = 0;
	Wb[3] = 0;

	Wc[0] = 0;
	Wc[1] = 0;
	Wc[2] = 1;
	Wc[3] = 0;

	Wd[0] = sin(T);
	Wd[1] = 0;
	Wd[2] = 0;
	Wd[3] = cos(T);
}

double from4[4] = { 5, 0, 0, 0 };
double to4[4] = { 0, 0, 0, 0 };
double up4[4] = { 0, 1, 0, 0 };
double over4[4] = { 0, 0, 1, 0 };

// generate a 4D view matrix
void view4(double* result)
{
	// column vectors
	double* Wa = result + 4 * 0;
	double* Wb = result + 4 * 1;
	double* Wc = result + 4 * 2;
	double* Wd = result + 4 * 3;

	// vector norm
	double norm;    

	// get the normalized Wd column-vector.
	vecSub4(Wd, to4, from4);
	norm = norm4(Wd);
	vecScale4(Wd, 1 / norm);

	// calculate the normalized Wa column-vector.
	cross4(Wa, up4, over4, Wd);
	norm = norm4(Wa);
	vecScale4(Wa, 1 / norm);

	// calculate the normalized Wb column-vector.
	cross4(Wb, over4, Wd, Wa);
	norm = norm4(Wb);
	vecScale4(Wb, 1 / norm);

	// calculate the Wc column-vector.
	cross4(Wc, Wd, Wa, Wb);
}

void projectTo3D(double vAngle, const double* matView, const double* matRotation)
{
	// column vectors
	const double* Wa = matView + 4 * 0;
	const double* Wb = matView + 4 * 1;
	const double* Wc = matView + 4 * 2;
	const double* Wd = matView + 4 * 3;

	// divisor Values
	double S, T;    

	T = 1 / tan(vAngle / 2);

	for (int i = 0; i < 16; ++i)
	{
		double V[4];
		matVecMul4(V, matRotation, V4[i]);

		double Vf[4];
		vecSub4(Vf, V, from4);

		S = T / dot4(Vf, Wd);

		V3[i][0] = S * dot4(Vf, Wa);
		V3[i][1] = S * dot4(Vf, Wb);
		V3[i][2] = S * dot4(Vf, Wc);
	}
}

double dot3(const double* V, const double* U)
{
	return (V[0] * U[0]) + (V[1] * U[1]) + (V[2] * U[2]);
}

double norm3(const double* V)
{
	return sqrt(dot3(V, V));
}

void cross3(double* result, const double* U, const double* V)
{
	result[0] = (U[1] * V[2]) - (U[2] * V[1]);
	result[1] = (U[2] * V[0]) - (U[0] * V[2]);
	result[2] = (U[0] * V[1]) - (U[1] * V[0]);
}
void vecSub3(double* result, const double* a, const double* b)
{
	result[0] = a[0] - b[0];
	result[1] = a[1] - b[1];
	result[2] = a[2] - b[2];
}

void vecScale3(double* vec, double m)
{
	vec[0] *= m;
	vec[1] *= m;
	vec[2] *= m;
}

void matVecMul3(double* result, const double* mat, const double* vec)
{
	for (int row = 0; row < 3; ++row)
	{
		result[row] = 0;
		for (int col = 0; col < 3; ++col)
		{
			result[row] += mat[col * 3 + row] * vec[col];
		}
	}
}

void rotXZ3(double* result, double T)
{
	// column vectors
	double* Va = result + 3 * 0;
	double* Vb = result + 3 * 1;
	double* Vc = result + 3 * 2;

	Va[0] = cos(T);
	Va[1] = 0;
	Va[2] = -sin(T);

	Vb[0] = 0;
	Vb[1] = 1;
	Vb[2] = 0;

	Vc[0] = sin(T);
	Vc[1] = 0;
	Vc[2] = cos(T);
}

double from3[3] = { 3.00f, 0.99f, 1.82f };
double to3[3] = { 0, 0, 0 };
double up3[3] = { 0, -1, 0 };

// generate a 3D view matrix
void view3(double* result)
{
	double* Va = result + 3 * 0;
	double* Vb = result + 3 * 1;
	double* Vc = result + 3 * 2;

	double norm;

	// Get the normalized Vc column-vector.
	vecSub3(Vc, to3, from3);
	norm = norm3(Vc);
	vecScale3(Vc, 1 / norm);

	// Calculate the normalized Va column-vector.
	cross3(Va, Vc, up3);
	norm = norm3(Va);
	vecScale3(Va, 1 / norm);

	// Calculate the Vb column-vector.
	cross3(Vb, Va, Vc);
}

void projectTo2D(double vAngle, const double* matView, const double* matRotation)
{
	// column vectors
	const double* Va = matView + 3 * 0;
	const double* Vb = matView + 3 * 1;
	const double* Vc = matView + 3 * 2;

	// divisor values
	double  S, T;

	T = 1 / tan(vAngle / 2);

	for (int i = 0; i < 16; ++i)
	{
		double V[3];
		matVecMul3(V, matRotation, V3[i]);

		double Vf[3];
		vecSub3(Vf, V, from3);

		S = T / dot3(Vf, Vc);

		V2[i][0] = (ww / 2) + (ww * S * dot3(Vf, Va));
		V2[i][1] = (wh / 2) + (wh * S * dot3(Vf, Vb));
	}
}

int main(int argc, const char* argv[])
{
	// get the console handle
	HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);

	// set console dimensions
	COORD s = { ww, wh };
	SMALL_RECT r = { 0, 0, ww, wh };
	COORD z = { 0, 0 };
	SetConsoleScreenBufferSize(h, s);
	SetConsoleWindowInfo(h, TRUE, &r);

	CHAR_INFO d[wh * ww];

	double viewMat4[4 * 4];
	view4(viewMat4);
	double rot4[4 * 4];

	double viewMat3[3 * 3];
	view3(viewMat3);
	double rot3[3 * 3];

	double rotation = 0;

	for (;;)
	{
		rotation += 0.01f;

		rotXW4(rot4, rotation);
		projectTo3D(M_PI / 3, viewMat4, rot4);

		rotXZ3(rot3, rotation * 0.3);
		projectTo2D(M_PI / 4, viewMat3, rot3);

		clr(d);

		for (int i = 0; i < 32; ++i)
		{
			int a = indices[i][0];
			int b = indices[i][1];
			COORD c1 = { (SHORT)V2[a][0], (SHORT)V2[a][1] };
			COORD c2 = { (SHORT)V2[b][0], (SHORT)V2[b][1] };
			ln(d, c1, c2);
		}

		WriteConsoleOutput(h, d, s, z, &r);

		Sleep(1);
	}
}