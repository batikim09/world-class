#include <cmath>
#include <iostream>


#include "SPTK.h"

using namespace std;

double ** sp2lsp(double **spectrogram, double **logperiodogram, int length, int dim){

  for (int i = 0; i < length; ++i)
    for (int j = 0; j < dim; ++j)
      logperiodogram[i][j] = log(spectrogram[i][j]);

  return logperiodogram;
}
/*Convert spectrum envelope to mel-cepstrum
This is a simplified implementation of ``mcep`` for input type
is 4.

Parameters
----------
powerspec : array
    Power spectrum

order : int
    Order of mel-cepstrum

alpha : float
    All-pass constant.

Returns
-------
mc : array, shape(``order+1``)
    mel-cepstrum
*/
double** sp2mc(double **powerspec, double **melspec, int length, int sp_dim, int mc_dim, double alpha)
{
    int ret = 0;

    double **logperiodogram = new double* [length];

    for (int i = 0; i < length; ++i) {
      logperiodogram[i] = new double[sp_dim];
    }

    double **cepstrum = new double* [length];
    int cep_dim = (sp_dim - 1) * 2;

    for (int i = 0; i < length; ++i) {
      cepstrum[i] = new double[cep_dim];
    }

    //|X(ω)|² -> log(|X(ω)²|)
    sp2lsp(powerspec, logperiodogram, length, sp_dim);
    cerr  << "LOG SPEC" << endl;

    cerr << "cep dim: " << cep_dim << endl;
    //transform log-periodogram to real cepstrum, log-periodogram becomes pseudo time dimension
    //log(|X(ω)|²) -> c(m), #IFFT points is decided by the number of dimensions
    for(int i = 0; i < length; i++){
       ret = ifftr(logperiodogram[i], cepstrum[i], cep_dim);
       if(ret != 0)
       {
         cerr << "ifftr error: " << ret << endl;
        }
    }
    cerr  << "ifftr" << endl;

    for(int d = 0; d < cep_dim; d++)
    {
      cerr  << "cep: " << d << endl;
      cepstrum[0][d] /= 2.0;
    }

    int src_dim = cep_dim - 1;
    cerr  << "ceps" << endl;

    //c(m) -> cₐ(m)    return freqt(c, order, alpha)
    for(int i = 0; i < length; i++){
      cerr  << "cep: " << i << endl;

     freqt(cepstrum[i], src_dim, melspec[i], mc_dim, alpha);
    }

    cerr  << "freqt" << endl;

    delete [] logperiodogram;
    delete [] cepstrum;
}
