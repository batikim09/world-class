#include <cmath>
#include <iostream>


// WORLD core functions.
#include "world_matlabfunctions.hpp"
#include "harvest.hpp"
#include "cheaptrick.hpp"
#include "d4c.hpp"
#include "synthesis.hpp"
#include "codec.hpp"
#include "conversion.hpp"

using namespace std;
using namespace world_class;


//-----------------------------------------------------------------------------
// struct for WORLD
// This struct is an option.
// Users are NOT forced to use this struct.
//-----------------------------------------------------------------------------
typedef struct {
	double frame_period;
	int fs;

	double *f0;
	double *time_axis;
	int f0_length;

	double **spectrogram;
	double **aperiodicity;
	double **coded_aperiodicity; //bap
	double **mel_cepstrum;

	int fft_size;
	int num_coded_aperidicites;
	int spec_order;
	int melceps_order;

	double alpha;

} WorldParameters;


void DisplayInformation(int fs, int nbit, int x_length)
{
	cout << "File information" << endl;
	cout << "Sampling : " << fs << " [Hz] " << nbit << " [Bit]" << endl;
	cout << "Length " << x_length << " [sample]" << endl;
	cout << "Length " << (double)x_length / fs << " [sec]" << endl;
}

void CalculateBarkAlpha(WorldParameters *world_parameters)
{
		world_parameters->alpha =  0.8517 * sqrt( atan(0.06583 * world_parameters->fs /1000.0) ) - 0.1916;		
}

void F0EstimationHarvest(
	double *x, int x_length,  WorldParameters *world_parameters)
{
	cout << "\nF0 estimation (Harvest)" << endl;

	// You can change the frame period.
	// But the estimation is carried out with 1-ms frame shift.
	HarvestOption option;
	option.frame_period = world_parameters->frame_period;

	// You can set the f0_floor below world::kFloorF0.
	option.f0_floor = 40.0;

	// You can use a Cosine table for fast computation,
	// but the computation is not exactly same as the original algorithm.
	// This option is not included int the original codes.

	// option.use_cos_table = true;


	Harvest harvest = Harvest(world_parameters->fs, option);


	// Parameters setting and memory allocation.
	world_parameters->f0_length = harvest.getSamples(world_parameters->fs, x_length);
	world_parameters->f0 = new double[world_parameters->f0_length];
	world_parameters->time_axis = new double[world_parameters->f0_length];


	harvest.compute(x, x_length, world_parameters->time_axis, world_parameters->f0);


}


void SpectralEnvelopeEstimation(
	double *x, int x_length, WorldParameters *world_parameters)
{
	cout << "\nSpectral envelope estimation (CheapTrick)" << endl;

	// Important notice (2017/01/02)
	// You can set the fft_size.
	// Default is GetFFTSizeForCheapTrick(world_parameters->fs, &option);
	// When fft_size changes from default value,
	// a replaced f0_floor will be used in CheapTrick().
	// The lowest F0 that WORLD can work as expected is determined
	// by the following : 3.0 * fs / fft_size.
	CheapTrickOption option;
	option.f0_floor = 71.0;
	// We can directly set fft_size.
	// option.fft_size = 1024;

	// Default value was modified to -0.15.
	// option.q1 = -0.15;


	CheapTrick cheaptrick = CheapTrick(world_parameters->fs, option);


	// Parameters setting and memory allocation.
	int fft_size = cheaptrick.getFFTSizeForCheapTrick(world_parameters->fs, option.f0_floor);
	world_parameters->fft_size = fft_size;
	world_parameters->spec_order = fft_size / 2 + 1;
	world_parameters->spectrogram = new double *[world_parameters->f0_length];

	for (int i = 0; i < world_parameters->f0_length; ++i)
		world_parameters->spectrogram[i] =
			new double[world_parameters->spec_order];


	cheaptrick.compute(x, x_length, world_parameters->time_axis,
					   world_parameters->f0, world_parameters->f0_length,
					   world_parameters->spectrogram);


}


void AperiodicityEstimation(
	double *x, int x_length, WorldParameters *world_parameters)
{
	cout << "\nAperiodicity estimation (D4C)" << endl;

	// Parameters setting and memory allocation.
	world_parameters->aperiodicity = new double *[world_parameters->f0_length];
	for (int i = 0; i < world_parameters->f0_length; ++i) {
		world_parameters->aperiodicity[i] = new double[world_parameters->fft_size / 2 + 1];
	}

	// Comment was modified because it was confusing (2017/12/10).
	// It is used to determine the aperiodicity in whole frequency band.
	// D4C identifies whether the frame is voiced segment even if it had an F0.
	// If the estimated value falls below the threshold,
	// the aperiodicity in whole frequency band will set to 1.0.
	// If you want to use the conventional D4C, please set the threshold to 0.0.
	D4COption option;
	option.threshold = 0.85;


	D4C d4c = D4C(world_parameters->fs, option);


	d4c.compute(x, x_length, world_parameters->time_axis,
				world_parameters->f0, world_parameters->f0_length,
				world_parameters->fft_size, world_parameters->aperiodicity);


}

//-----------------------------------------------------------------------------
// InitAndCodeAperiodicity codes the aperiodicity. The number of dimensions is
// determined by fs.
//
// Input: world_parameters that contains:
//   aperiodicity       : Aperiodicity before coding
//   f0_length          : Length of F0 contour
//   fs                 : Sampling frequency
//   fft_size           : FFT size
//
// Output: in world_parameters
//   coded_aperiodicity : Coded aperiodicity
//-----------------------------------------------------------------------------
void InitAndCodeAperiodicity(WorldParameters *world_parameters){

	int num_coded_aperidicites = GetNumberOfAperiodicities(world_parameters->fs);
	cout << "num_coded_aperidicites: " << num_coded_aperidicites << endl;
	// Parameters setting and memory allocation.
	world_parameters->coded_aperiodicity = new double *[world_parameters->f0_length];
	world_parameters->num_coded_aperidicites = num_coded_aperidicites;

	for (int i = 0; i < world_parameters->f0_length; ++i) {
		world_parameters->coded_aperiodicity[i] = new double[num_coded_aperidicites];
	}


	CodeAperiodicity(world_parameters->aperiodicity, world_parameters->f0_length,
						world_parameters->fs, world_parameters->fft_size, world_parameters->coded_aperiodicity);


}

void Spectrogram2MelCepstrum(WorldParameters *world_parameters){

	world_parameters->mel_cepstrum = new double *[world_parameters->f0_length];

	for (int i = 0; i < world_parameters->f0_length; ++i) {
		world_parameters->mel_cepstrum[i] = new double[world_parameters->melceps_order + 1];
	}

	sp2mc(world_parameters->spectrogram, world_parameters->mel_cepstrum, world_parameters->f0_length, world_parameters->spec_order, world_parameters->melceps_order, world_parameters->alpha);

}

void ParameterModification(
	int argc, char *argv[], int fs, int f0_length,
	int fft_size, double *f0, double **spectrogram)
{
	// F0 scaling
	if (argc >= 4) {
		double shift = atof(argv[3]);
		for (int i = 0; i < f0_length; ++i) f0[i] *= shift;
	}

	if (argc < 5) return;

	// Spectral stretching
	double ratio = atof(argv[4]);
	double *freq_axis1 = new double[fft_size];
	double *freq_axis2 = new double[fft_size];
	double *spectrum1 = new double[fft_size];
	double *spectrum2 = new double[fft_size];

	for (int i = 0; i <= fft_size / 2; ++i) {
		freq_axis1[i] = ratio * i / fft_size * fs;
		freq_axis2[i] = static_cast<double>(i) / fft_size * fs;
	}

	for (int i = 0; i < f0_length; ++i) {
		for (int j = 0; j <= fft_size / 2; ++j)
			spectrum1[j] = log(spectrogram[i][j]);
		interp1(freq_axis1, spectrum1, fft_size / 2 + 1, freq_axis2,
				fft_size / 2 + 1, spectrum2);
		for (int j = 0; j <= fft_size / 2; ++j)
			spectrogram[i][j] = exp(spectrum2[j]);
		if (ratio >= 1.0) continue;
		for (int j = static_cast<int>(fft_size / 2.0 * ratio);
			 j <= fft_size / 2; ++j)
			spectrogram[i][j] =
				spectrogram[i][static_cast<int>(fft_size / 2.0 * ratio) - 1];
	}

	delete[] spectrum1;
	delete[] spectrum2;
	delete[] freq_axis1;
	delete[] freq_axis2;
}

void WaveformSynthesis1(
	WorldParameters *world_parameters, int fs, int y_length, double *y)
{
	cout << "\nSynthesis 1 (conventional algorithm)" << endl;


	Synthesis synthesis = Synthesis(fs, world_parameters->fft_size, world_parameters->frame_period);


	synthesis.compute(world_parameters->f0, world_parameters->f0_length,
					  world_parameters->spectrogram, world_parameters->aperiodicity,
					  y_length, y);


}

void DestroyMemory(WorldParameters *world_parameters) {
	delete[] world_parameters->time_axis;
	delete[] world_parameters->f0;
	for (int i = 0; i < world_parameters->f0_length; ++i) {
		delete[] world_parameters->spectrogram[i];
		delete[] world_parameters->aperiodicity[i];
		delete[] world_parameters->coded_aperiodicity[i];
		delete[] world_parameters->mel_cepstrum[i];
	}
	delete[] world_parameters->spectrogram;
	delete[] world_parameters->aperiodicity;
	delete[] world_parameters->coded_aperiodicity;
	delete[] world_parameters->mel_cepstrum;
}
