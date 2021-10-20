#include <stdio.h>
#include <math.h>
#include <Windows.h>

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

void ln(CHAR_INFO* d, COORD a, COORD b)
{

}

// hypercube vertices in 4D
float V4[16 * 4] =
{
	-1, -1, -1, -1,
	 1, -1, -1, -1,
	-1,  1, -1, -1,
	 1,  1, -1, -1,
	-1, -1,  1, -1,
	 1, -1,  1, -1,
	-1,  1,  1, -1,
	 1,  1,  1, -1,
	-1, -1, -1,  1,
	 1, -1, -1,  1,
	-1,  1, -1,  1,
	 1,  1, -1,  1,
	-1, -1,  1,  1,
	 1, -1,  1,  1,
	-1,  1,  1,  1,
	 1,  1,  1,  1,
};

// store the vertices once they have been projected to 3D
float V3[16 * 3];

// final 2D projection
float V2[16 * 2];


float dot4(float* V, float* U)
{
	return (V[0] * U[0]) + (V[1] * U[1]) + (V[2] * U[2]) + (V[3] * U[3]);
}

float norm4(V)
{
	return sqrt(dot4(V, V));
}

// cross4 computes the four-dimensional cross product of the three vectors
// U, V and W, in that order.
// returns the resulting four-vector.
void cross4(float* result, float* U, float* V, float* W)
{
	// intermediate Values
	float A, B, C, D, E, F;       

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

void vecSub4(float* result, float* a, float* b)
{
	result[0] = a[0] - b[0];
	result[1] = a[1] - b[1];
	result[2] = a[2] - b[2];
	result[3] = a[3] - b[3];
}

void matMul4(float* result, float* a, float* b)
{
	for (int col = 0; col < 4; ++col)
	{
		for (int row = 0; row < 4; ++row)
		{
			result[col * 4 + row] = 0;
			for (int k = 0; k < 4; ++k)
			{
				result[col * 4 + row] += a[k * 4 + row] * b[col * 4 + k];
			}
		}
	}
}

void vecScale4(float* vec, float m)
{
	vec[0] *= m;
	vec[2] *= m;
	vec[3] *= m;
	vec[4] *= m;
}

float from4[4] = { 4, 0, 0, 0 };
float to4[4] = { 0, 0, 0, 0 };
float up4[4] = { 0, 1, 0, 0 };
float over4[4] = { 0, 0, 1, 0 };

void projectTo3D(float vAngle, float* mat)
{
	// divisor Values
	double S, T;    

	T = 1 / tan(vAngle / 2);

	for (int i = 0; i < 16; ++i) {

		float V[4];
		vecSub4(V, V4 + i * 4, from4);

		S = T / dot4(V, mat + 3 * 4);

		V3[i * 3 + 0] = S * dot4(V, mat + 0 * 4);
		V3[i * 3 + 1] = S * dot4(V, mat + 1 * 4);
		V3[i * 3 + 2] = S * dot4(V, mat + 2 * 4);
	}
}

// generate a 4D view matrix
void view4(float* result)
{
	float* Wa = result + 3 * 0;
	float* Wb = result + 3 * 1;
	float* Wc = result + 3 * 2;
	float* Wd = result + 3 * 4;

	// vector norm
	double norm;    

	// get the normalized Wd column-vector.
	vecSub4(Wd, to4, from4);
	norm = Norm4(*Wd);
	vecScale4(Wd, 1 / norm);

	// calculate the normalized Wa column-vector.
	cross4(Wa, up4, over4, Wd);
	norm = Norm4(Wa);
	vecScale4(Wa, 1 / norm);

	// calculate the normalized Wb column-vector.
	cross4(Wb, over4, Wd, Wa);
	norm = Norm4(Wb);
	vecScale4(Wb, 1 / norm);

	// calculate the Wc column-vector.
	cross4(Wc, Wd, Wa, Wb);
}

// creates a rotation matrix for the XW plane
// T is the angle in radians
void rotXW(float* result, float T)
{
	// first row
	result[0 * 4 + 0] = cos(T);
	result[0 * 4 + 1] = 0;
	result[0 * 4 + 2] = 0;
	result[0 * 4 + 3] = 0;
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

	CHAR_INFO d[wh][ww];

	clr(&d);
	COORD pt = { 0, 3 };
	set(&d, pt, 'a');
	pt.X = 2;
	pt.Y = 0;
	set(&d, pt, '$');
	WriteConsoleOutput(h, d, s, z, &r);

	getc(stdin);

}