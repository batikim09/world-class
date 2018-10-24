//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise,
// Copyright 2018 Yukara Ikemiya
//
// test.exe input.wav outout.wav f0 spec
// input.wav  : Input file
//
// output.wav : Output file
// f0         : F0 scaling (a positive number)
// spec       : Formant scaling (a positive number)
//
// Note: This version output three speech synthesized by different algorithms.
//       When the filename is "output.wav", "01output.wav", "02output.wav" and
//       "03output.wav" are generated. They are almost all the same.
//-----------------------------------------------------------------------------
#include <chrono>
#include "world.hpp"
#include "audioio.hpp"

using chrono_tp = chrono::system_clock::time_point;

// time measurement methods
inline chrono_tp get_time_now(){
	return chrono::system_clock::now();
}

inline double get_elapsed_msec(chrono_tp st, chrono_tp end) {
	return chrono::duration_cast<chrono::nanoseconds>(end-st).count() / 1.0e6;
}

//-----------------------------------------------------------------------------
// Test program.
// test.exe input.wav outout f0 spec flag
// input.wav  : argv[1] Input file path
// output	  : argv[2] Output file name
// f0         : argv[3] F0 scaling (a positive number)
// spec       : argv[4] Formant shift (a positive number)
//-----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
	if (argc != 2 && argc != 3 && argc != 4 && argc != 5) {
		cout << "Usage : ./test [input wav path] [output wav name]" << endl;
		cout << "Example: ./test ./aiueo.wav output" << endl;
		return -2;
	}

	// Memory allocation is carried out in advanse.
	// This is for compatibility with C language.
	int x_length = GetAudioLength(argv[1]);
	if (x_length <= 0) {
		if (x_length == 0) { cerr << "error: File not found.\n" << endl; }
		else { cerr << "error: The file is not .wav format.\n" << endl; }
		return -1;
	}


	// wavread() must be called after GetAudioLength().
	int fs, nbit;
	double *x = new double[x_length];
	WavRead(argv[1], &fs, &nbit, x);

	DisplayInformation(fs, nbit, x_length);

	//---------------------------------------------------------------------------
	// Analysis part
	//---------------------------------------------------------------------------
	WorldParameters world_parameters;
	// You must set fs and frame_period before analysis/synthesis.
	world_parameters.fs = fs;
	// 5.0 ms is the default value.
	world_parameters.frame_period = 5.0;
	world_parameters.melceps_order = 59;

	// F0 estimation

	// Harvest
	auto f1 = get_time_now();
	F0EstimationHarvest(x, x_length, &world_parameters);
	auto f2 = get_time_now();
	cout << "\t F0 compute:\t" << get_elapsed_msec(f1, f2) << " [msec]" << endl;

	double f0_sum = 0;
	for (int i = 0; i < world_parameters.f0_length; i++) {
		f0_sum += world_parameters.f0[i];
	}

	// Spectral envelope estimation
	auto s1 = get_time_now();
	SpectralEnvelopeEstimation(x, x_length, &world_parameters);
	auto s2 = get_time_now();
	cout << "\t Spec compute:\t" << get_elapsed_msec(s1, s2) << " [msec]" << endl;

	double spec_sum = 0;
	for (int i = 0; i < world_parameters.f0_length; ++i)
		for (int j = 0; j < world_parameters.fft_size / 2 + 1; ++j)
		{ spec_sum +=  world_parameters.spectrogram[i][j]; }

	// Aperiodicity estimation by D4C
	auto a1 = get_time_now();
	AperiodicityEstimation(x, x_length, &world_parameters);
	auto a2 = get_time_now();
	cout << "\t Aperiodicty compute:\t" << get_elapsed_msec(a1, a2) << " [msec]" << endl;

	double ap_sum = 0;
	for (int i = 0; i < world_parameters.f0_length; ++i)
		for (int j = 0; j < world_parameters.fft_size / 2 + 1; ++j)
		{ ap_sum +=  world_parameters.aperiodicity[i][j]; }

	// Summarize the Parameters
	cout << "fs: " << world_parameters.fs << endl;
	cout << "frame_period: " << world_parameters.frame_period << endl;
	cout << "f_length: " << world_parameters.f0_length << endl;
	cout << "fft_size: " << world_parameters.fft_size << endl;
	cout << "spec_sum: " << spec_sum << endl;
	cout << "f0_sum: " << spec_sum << endl;


	/* Note that F0 must not be changed until all parameters are estimated.
	ParameterModification(argc, argv, fs, world_parameters.f0_length,
						  world_parameters.fft_size, world_parameters.f0,
						  world_parameters.spectrogram);
	*/

	// code aperiodicity
	auto c_a1 = get_time_now();
	InitAndCodeAperiodicity(&world_parameters);
	auto c_a2 = get_time_now();
	cout << "\t Coded aperiodicty compute:\t" << get_elapsed_msec(c_a1, c_a2) << " [msec]" << endl;

	double coded_ap_sum = 0;
	for (int i = 0; i < world_parameters.f0_length; ++i)
		for (int j = 0; j < world_parameters.num_coded_aperidicites; ++j)
		{ coded_ap_sum +=  world_parameters.coded_aperiodicity[i][j]; }
	cout << "coded_ap_sum: " << coded_ap_sum << endl;

	// Spectrogram to Mel cepstrum
	auto mc1 = get_time_now();
	CalculateBarkAlpha(&world_parameters);
	Spectrogram2MelCepstrum(&world_parameters);
	auto mc2 = get_time_now();
	cout << "\t Mel cepstrum compute:\t" << get_elapsed_msec(mc1, mc2) << " [msec]" << endl;


	//---------------------------------------------------------------------------
	// Synthesis part
	// There are three samples in speech synthesis
	// 1: Conventional synthesis
	// 2: Example of real-time synthesis
	// 3: Example of real-time synthesis (Ring buffer is efficiently used)
	//---------------------------------------------------------------------------
	char filename[100];
	// The length of the output waveform
	int y_length = static_cast<int>((world_parameters.f0_length - 1) *
									world_parameters.frame_period / 1000.0 * fs) + 1;
	double *y = new double[y_length]();


	// Synthesis 1 (conventional synthesis)
	auto syn1 = get_time_now();
	WaveformSynthesis1(&world_parameters, fs, y_length, y);
	if (argc == 2){
		sprintf(filename, "output_1.wav");
	}
	else{
		sprintf(filename, "%s_1.wav", argv[2]);
	}
	
	WavWrite(y, y_length, fs, 16, filename);
	auto syn2 = get_time_now();
	cout << "\t Synthesis compute:\t" << get_elapsed_msec(syn1, syn2) << " [msec]" << endl;

	double synt_sum_1 = 0;
	for (int ii = 0; ii < y_length; ii++) {
		synt_sum_1 += y[ii];
	}

	delete[] y;
	delete[] x;
	DestroyMemory(&world_parameters);

	cout << "complete." << endl;
	return 0;
}
