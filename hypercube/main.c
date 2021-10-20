#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <Windows.h>

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

float V42[16 * 4];

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
	// intermediate values
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

void vecScale4(float* vec, float m)
{
	vec[0] *= m;
	vec[1] *= m;
	vec[2] *= m;
	vec[3] *= m;
}

void matMul4(float* result, const float* a, const float* b)
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

void matVecMul4(float* result, const float* mat, const float* vec)
{
	for (int row = 0; row < 4; ++row)
	{
		result[row] = 0;
		for (int col = 0; col < 4; ++col)
		{
			result[row] += mat[col * 4 + row] * vec[col];
		}
	}
	return result;
}

// creates a rotation matrix for the XW plane
// T is the angle in radians
void rotXW(float* result, float T)
{
	// column vectors
	float* Wa = result + 4 * 0;
	float* Wb = result + 4 * 1;
	float* Wc = result + 4 * 2;
	float* Wd = result + 4 * 3;

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

// creates a rotation matrix for the ZW plane
// T is the angle in radians
void rotZW(float* result, float T)
{
	// column vectors
	float* Wa = result + 4 * 0;
	float* Wb = result + 4 * 1;
	float* Wc = result + 4 * 2;
	float* Wd = result + 4 * 3;

	Wa[0] = 1;
	Wa[1] = 0;
	Wa[2] = 0;
	Wa[3] = 0;

	Wb[0] = 0;
	Wb[1] = 1;
	Wb[2] = 0;
	Wb[3] = 0;

	Wc[0] = 0;
	Wc[1] = 0;
	Wc[2] = cos(T);
	Wc[3] = sin(T);

	Wd[0] = 0;
	Wd[1] = 0;
	Wd[2] = -sin(T);
	Wd[3] = cos(T);
}

float from4[4] = { 8, 0, 0, 0 };
float to4[4] = { 0, 0, 0, 0 };
float up4[4] = { 0, 1, 0, 0 };
float over4[4] = { 0, 0, 1, 0 };

// generate a 4D view matrix
void view4(float* result)
{
	// column vectors
	float* Wa = result + 4 * 0;
	float* Wb = result + 4 * 1;
	float* Wc = result + 4 * 2;
	float* Wd = result + 4 * 3;

	// vector norm
	double norm;    

	// get the normalized Wd column-vector.
	vecSub4(Wd, &to4, &from4);
	norm = norm4(Wd);
	vecScale4(Wd, 1 / norm);

	// calculate the normalized Wa column-vector.
	cross4(Wa, &up4, &over4, Wd);
	norm = norm4(Wa);
	vecScale4(Wa, 1 / norm);

	// calculate the normalized Wb column-vector.
	cross4(Wb, &over4, Wd, Wa);
	norm = norm4(Wb);
	vecScale4(Wb, 1 / norm);

	// calculate the Wc column-vector.
	cross4(Wc, Wd, Wa, Wb);
}

void projectTo3D(float vAngle, const float* mat)
{
	// column vectors
	const float* Wa = mat + 4 * 0;
	const float* Wb = mat + 4 * 1;
	const float* Wc = mat + 4 * 2;
	const float* Wd = mat + 4 * 3;

	// divisor Values
	double S, T;    

	T = 1 / tan(vAngle / 2);

	for (int i = 0; i < 16; ++i)
	{
		float V[4];
		vecSub4(&V, V42 + i * 4, from4);

		S = T / dot4(&V, Wd);

		V3[i * 3 + 0] = S * dot4(&V, Wa);
		V3[i * 3 + 1] = S * dot4(&V, Wb);
		V3[i * 3 + 2] = S * dot4(&V, Wc);
	}
}

float dot3(float* V, float* U)
{
	return (V[0] * U[0]) + (V[1] * U[1]) + (V[2] * U[2]);
}

float norm3(V)
{
	return sqrt(dot3(V, V));
}

void cross3(float* result, float* U, float* V)
{
	result[0] = (U[1] * V[2]) - (U[2] * V[1]);
	result[1] = (U[2] * V[0]) - (U[0] * V[2]);
	result[2] = (U[0] * V[1]) - (U[1] * V[0]);
}
void vecSub3(float* result, float* a, float* b)
{
	result[0] = a[0] - b[0];
	result[1] = a[1] - b[1];
	result[2] = a[2] - b[2];
}

void vecScale3(float* vec, float m)
{
	vec[0] *= m;
	vec[1] *= m;
	vec[2] *= m;
}

float from3[3] = { 3.00, 0.99, 1.82 };
float to3[3] = { 0, 0, 0 };
float up3[3] = { 0, -1, 0 };

// generate a 3D view matrix
void view3(float* result)
{
	float* Va = result + 3 * 0;
	float* Vb = result + 3 * 1;
	float* Vc = result + 3 * 2;

	double norm;

	// Get the normalized Vc column-vector.
	vecSub3(Vc, &to3, &from3);
	norm = norm3(Vc);
	vecScale3(Vc, 1 / norm);

	// Calculate the normalized Va column-vector.
	cross3(Va, Vc, &up3);
	norm = norm3(Va);
	vecScale3(Va, 1 / norm);

	// Calculate the Vb column-vector.
	cross3(Vb, Va, Vc);
}

void projectTo2D(float vAngle, const float* mat)
{
	// column vectors
	const float* Va = mat + 3 * 0;
	const float* Vb = mat + 3 * 1;
	const float* Vc = mat + 3 * 2;

	// divisor values
	double  S, T;

	T = 1 / tan(vAngle / 2);

	for (int i = 0; i < 16; ++i)
	{
		float V[3];
		vecSub3(&V, V3 + i * 3, from3);

		S = T / dot3(&V, Vc);

		V2[i * 2 + 0] = (ww / 2) + (ww * S * dot3(&V, Va));
		V2[i * 2 + 1] = (wh / 2) + (wh * S * dot3(&V, Vb));
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

	CHAR_INFO d[wh][ww];

	float viewMat4[4][4];
	float modelMat4[4][4];
	float MV4[4][4];

	float viewMat3[3][3];

	//float modelMat3[3][3];
	//float MV3[3][3];

	float rotation = 0;

	for (;;)
	{
		rotation += 0.001;

		float rot4[4][4];
		rotXW(&rot4, rotation);
		for (int i = 0; i < 16; ++i)
		{
			matVecMul4(&V42[i * 4], &rot4, &V4[i * 4]);
		}

		view4(&viewMat4);
		//rotZW(&modelMat4, rotation);
		//matMul4(&MV4, &viewMat4, &modelMat4);

		projectTo3D(M_PI / 4, &viewMat4);

		view3(&viewMat3);
		projectTo2D(M_PI / 4, &viewMat3);

		clr(&d);

		for (int i = 0; i < 16; ++i)
		{
			COORD pt = { V2[i * 2 + 0], V2[i * 2 + 1] };

			if (pt.X >= 0 && pt.X < ww && pt.Y >= 0 && pt.Y <= wh)
			{
				set(&d, pt, '@');
			}

		}

		WriteConsoleOutput(h, d, s, z, &r);
	}


	//getc(stdin);

}