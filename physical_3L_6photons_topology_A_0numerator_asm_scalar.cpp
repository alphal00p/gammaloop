#include <iostream>
#include <complex>
#include <cmath>

extern "C" {
	void physical_3L_6photons_topology_A_0numerator_asm_scalar_complex(std::complex<double>* params, std::complex<double>* out) {
	std::complex<double> Z[114];
	Z[0] = params[0];
	Z[1] = params[1];
	Z[2] = params[2];
	Z[3] = params[3];
	Z[4] = params[4];
	Z[5] = params[5];
	Z[6] = params[6];
	Z[7] = params[7];
	Z[8] = params[8];
	Z[9] = params[9];
	Z[10] = params[10];
	Z[11] = params[11];
	Z[12] = params[12];
	Z[13] = params[13];
	Z[14] = params[14];
	Z[15] = params[15];
	Z[16] = params[16];
	Z[17] = params[17];
	Z[18] = params[18];
	Z[19] = params[19];
	Z[20] = params[20];
	Z[21] = params[21];
	Z[22] = params[22];
	Z[23] = params[23];
	Z[24] = params[24];
	Z[25] = params[25];
	Z[26] = params[26];
	Z[27] = params[27];
	Z[28] = params[28];
	Z[29] = params[29];
	Z[30] = params[30];
	Z[31] = params[31];
	Z[32] = params[32];
	Z[33] = params[33];
	Z[34] = params[34];
	Z[35] = params[35];
	Z[36] = params[36];
	Z[37] = params[37];
	Z[38] = params[38];
	Z[39] = params[39];
	Z[40] = params[40];
	Z[41] = params[41];
	Z[42] = params[42];
	Z[43] = params[43];
	Z[44] = params[44];
	Z[45] = params[45];
	Z[46] = params[46];
	Z[47] = params[47];
	Z[48] = params[48];
	Z[49] = params[49];
	Z[50] = params[50];
	Z[51] = params[51];
	Z[52] = params[52];
	Z[53] = params[53];
	Z[54] = params[54];
	Z[55] = params[55];
	Z[56] = params[56];
	Z[57] = params[57];
	Z[58] = params[58];
	Z[59] = params[59];
	Z[60] = params[60];
	Z[61] = params[61];
	Z[62] = params[62];
	Z[63] = params[63];
	Z[64] = params[64];
	Z[65] = params[65];
	Z[66] = params[66];
	Z[67] = params[67];
	Z[68] = params[68];
	Z[69] = params[69];
	Z[70] = params[70];
	Z[71] = params[71];
	Z[72] = params[72];
	Z[73] = params[73];
	Z[74] = params[74];
	Z[75] = params[75];
	Z[76] = params[76];
	Z[77] = params[77];
	Z[78] = params[78];
	Z[79] = params[79];
	Z[80] = params[80];
	Z[81] = params[81];
	Z[82] = params[82];
	Z[83] = params[83];
	Z[84] = params[84];
	Z[85] = params[85];
	Z[86] = params[86];
	Z[87] = params[87];
	Z[88] = params[88];
	Z[89] = params[89];
	Z[90] = params[90];
	Z[91] = params[91];
	Z[92] = params[92];
	Z[93] = params[93];
	Z[94] = params[94];
	Z[95] = params[95];
	Z[96] = params[96];
	Z[97] = 0.5;
	Z[98] = 89.89849108367626;
	Z[99] = -1;
	Z[100] = 3;
	Z[101] = 1;
	Z[102] = 0.0075467711139788835;
	Z[103] = 0.118;
	Z[104] = 173;
	Z[105] = 3.141592653589793;
	__asm__(
		"xorpd xmm0, xmm0\n\t"
		"addpd xmm0, XMMWORD PTR [%0+1584]\n\t"
		"addpd xmm0, XMMWORD PTR [%0+1600]\n\t"
		"movapd XMMWORD PTR [%0+1696], xmm0\n\t"
		"xorpd xmm0, xmm0\n\t"
		"addpd xmm0, XMMWORD PTR [%0+1616]\n\t"
		"addpd xmm0, XMMWORD PTR [%0+1600]\n\t"
		"movapd XMMWORD PTR [%0+1712], xmm0\n\t"
		"movapd xmm1, XMMWORD PTR [%0+1584]\n\t"
		"movapd xmm2, XMMWORD PTR [%0+1712]\n\t"
		"movapd xmm0, xmm1\n\t"
		"unpckhpd xmm0, xmm0\n\t"
		"unpcklpd xmm1, xmm1\n\t"
		"mulpd xmm0, xmm2\n\t"
		"mulpd xmm1, xmm2\n\t"
		"shufpd xmm0, xmm0, 1\n\t"
		"addsubpd xmm1, xmm0\n\t"
		"movapd XMMWORD PTR [%0+1712], xmm1\n\t"
		"xorpd xmm0, xmm0\n\t"
		"addpd xmm0, XMMWORD PTR [%0+1616]\n\t"
		"addpd xmm0, XMMWORD PTR [%0+1664]\n\t"
		"movapd XMMWORD PTR [%0+1728], xmm0\n\t"
		"movapd xmm1, XMMWORD PTR [%0+1584]\n\t"
		"movapd xmm2, XMMWORD PTR [%0+1728]\n\t"
		"movapd xmm0, xmm1\n\t"
		"unpckhpd xmm0, xmm0\n\t"
		"unpcklpd xmm1, xmm1\n\t"
		"mulpd xmm0, xmm2\n\t"
		"mulpd xmm1, xmm2\n\t"
		"shufpd xmm0, xmm0, 1\n\t"
		"addsubpd xmm1, xmm0\n\t"
		"movapd XMMWORD PTR [%0+1744], xmm1\n\t"
:
                        : "r"(Z)
                        : "memory");
                        	Z[110] = pow(Z[100], -1);
	Z[111] = sqrt(Z[102]);
	Z[112] = sqrt(Z[103]);
	Z[113] = sqrt(Z[105]);
	Z[106] = Z[97]*Z[97]*Z[98]*Z[106]*Z[107]*Z[109]*Z[110]*Z[111]*Z[111]*Z[111]*Z[111]*Z[111]*Z[111]*Z[112]*Z[112]*Z[112]*Z[112]*Z[108]*Z[108]*Z[108]*Z[108]*Z[108]*Z[108]*Z[108]*Z[108]*Z[108]*Z[113]*Z[113]*Z[113]*Z[113]*Z[113]*Z[113]*Z[113]*Z[113]*Z[113]*Z[113];
	out[0] = Z[106];
	return;
}

	void physical_3L_6photons_topology_A_0numerator_asm_scalar_double(double* params, double* out) {
	double Z[114];
	Z[0] = params[0];
	Z[1] = params[1];
	Z[2] = params[2];
	Z[3] = params[3];
	Z[4] = params[4];
	Z[5] = params[5];
	Z[6] = params[6];
	Z[7] = params[7];
	Z[8] = params[8];
	Z[9] = params[9];
	Z[10] = params[10];
	Z[11] = params[11];
	Z[12] = params[12];
	Z[13] = params[13];
	Z[14] = params[14];
	Z[15] = params[15];
	Z[16] = params[16];
	Z[17] = params[17];
	Z[18] = params[18];
	Z[19] = params[19];
	Z[20] = params[20];
	Z[21] = params[21];
	Z[22] = params[22];
	Z[23] = params[23];
	Z[24] = params[24];
	Z[25] = params[25];
	Z[26] = params[26];
	Z[27] = params[27];
	Z[28] = params[28];
	Z[29] = params[29];
	Z[30] = params[30];
	Z[31] = params[31];
	Z[32] = params[32];
	Z[33] = params[33];
	Z[34] = params[34];
	Z[35] = params[35];
	Z[36] = params[36];
	Z[37] = params[37];
	Z[38] = params[38];
	Z[39] = params[39];
	Z[40] = params[40];
	Z[41] = params[41];
	Z[42] = params[42];
	Z[43] = params[43];
	Z[44] = params[44];
	Z[45] = params[45];
	Z[46] = params[46];
	Z[47] = params[47];
	Z[48] = params[48];
	Z[49] = params[49];
	Z[50] = params[50];
	Z[51] = params[51];
	Z[52] = params[52];
	Z[53] = params[53];
	Z[54] = params[54];
	Z[55] = params[55];
	Z[56] = params[56];
	Z[57] = params[57];
	Z[58] = params[58];
	Z[59] = params[59];
	Z[60] = params[60];
	Z[61] = params[61];
	Z[62] = params[62];
	Z[63] = params[63];
	Z[64] = params[64];
	Z[65] = params[65];
	Z[66] = params[66];
	Z[67] = params[67];
	Z[68] = params[68];
	Z[69] = params[69];
	Z[70] = params[70];
	Z[71] = params[71];
	Z[72] = params[72];
	Z[73] = params[73];
	Z[74] = params[74];
	Z[75] = params[75];
	Z[76] = params[76];
	Z[77] = params[77];
	Z[78] = params[78];
	Z[79] = params[79];
	Z[80] = params[80];
	Z[81] = params[81];
	Z[82] = params[82];
	Z[83] = params[83];
	Z[84] = params[84];
	Z[85] = params[85];
	Z[86] = params[86];
	Z[87] = params[87];
	Z[88] = params[88];
	Z[89] = params[89];
	Z[90] = params[90];
	Z[91] = params[91];
	Z[92] = params[92];
	Z[93] = params[93];
	Z[94] = params[94];
	Z[95] = params[95];
	Z[96] = params[96];
	Z[97] = 0.5;
	Z[98] = 89.89849108367626;
	Z[99] = -1;
	Z[100] = 3;
	Z[101] = 1;
	Z[102] = 0.0075467711139788835;
	Z[103] = 0.118;
	Z[104] = 173;
	Z[105] = 3.141592653589793;
	__asm__(
		"movsd xmm0, QWORD PTR [%0+792]\n\t"
		"addsd xmm0, QWORD PTR [%0+800]\n\t"
		"movsd QWORD PTR [%0+848], xmm0\n\t"
		"movsd xmm0, QWORD PTR [%0+808]\n\t"
		"addsd xmm0, QWORD PTR [%0+800]\n\t"
		"movsd QWORD PTR [%0+856], xmm0\n\t"
		"movsd xmm0, QWORD PTR [%0+792]\n\t"
		"mulsd xmm0, QWORD PTR [%0+856]\n\t"
		"movsd QWORD PTR [%0+856], xmm0\n\t"
		"movsd xmm0, QWORD PTR [%0+808]\n\t"
		"addsd xmm0, QWORD PTR [%0+832]\n\t"
		"movsd QWORD PTR [%0+864], xmm0\n\t"
		"movsd xmm0, QWORD PTR [%0+792]\n\t"
		"mulsd xmm0, QWORD PTR [%0+864]\n\t"
		"movsd QWORD PTR [%0+872], xmm0\n\t"
:
                        : "r"(Z)
                        : "memory");
                        	Z[110] = pow(Z[100], -1);
	Z[111] = sqrt(Z[102]);
	Z[112] = sqrt(Z[103]);
	Z[113] = sqrt(Z[105]);
	__asm__(
		"movsd xmm0, QWORD PTR [%0+776]\n\t"
		"mulsd xmm0, QWORD PTR [%0+776]\n\t"
		"mulsd xmm0, QWORD PTR [%0+784]\n\t"
		"mulsd xmm0, QWORD PTR [%0+848]\n\t"
		"mulsd xmm0, QWORD PTR [%0+856]\n\t"
		"mulsd xmm0, QWORD PTR [%0+872]\n\t"
		"mulsd xmm0, QWORD PTR [%0+880]\n\t"
		"mulsd xmm0, QWORD PTR [%0+888]\n\t"
		"mulsd xmm0, QWORD PTR [%0+888]\n\t"
		"mulsd xmm0, QWORD PTR [%0+888]\n\t"
		"mulsd xmm0, QWORD PTR [%0+888]\n\t"
		"mulsd xmm0, QWORD PTR [%0+888]\n\t"
		"mulsd xmm0, QWORD PTR [%0+888]\n\t"
		"mulsd xmm0, QWORD PTR [%0+896]\n\t"
		"mulsd xmm0, QWORD PTR [%0+896]\n\t"
		"mulsd xmm0, QWORD PTR [%0+896]\n\t"
		"mulsd xmm0, QWORD PTR [%0+896]\n\t"
		"mulsd xmm0, QWORD PTR [%0+864]\n\t"
		"mulsd xmm0, QWORD PTR [%0+864]\n\t"
		"mulsd xmm0, QWORD PTR [%0+864]\n\t"
		"mulsd xmm0, QWORD PTR [%0+864]\n\t"
		"mulsd xmm0, QWORD PTR [%0+864]\n\t"
		"mulsd xmm0, QWORD PTR [%0+864]\n\t"
		"mulsd xmm0, QWORD PTR [%0+864]\n\t"
		"mulsd xmm0, QWORD PTR [%0+864]\n\t"
		"mulsd xmm0, QWORD PTR [%0+864]\n\t"
		"mulsd xmm0, QWORD PTR [%0+904]\n\t"
		"mulsd xmm0, QWORD PTR [%0+904]\n\t"
		"mulsd xmm0, QWORD PTR [%0+904]\n\t"
		"mulsd xmm0, QWORD PTR [%0+904]\n\t"
		"mulsd xmm0, QWORD PTR [%0+904]\n\t"
		"mulsd xmm0, QWORD PTR [%0+904]\n\t"
		"mulsd xmm0, QWORD PTR [%0+904]\n\t"
		"mulsd xmm0, QWORD PTR [%0+904]\n\t"
		"mulsd xmm0, QWORD PTR [%0+904]\n\t"
		"mulsd xmm0, QWORD PTR [%0+904]\n\t"
		"movsd QWORD PTR [%0+848], xmm0\n\t"
:
            : "r"(Z)
            : "memory");
            	out[0] = Z[106];
	return;
}
}
