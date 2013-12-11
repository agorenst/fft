#ifndef _FFT_H_
#define _FFT_H_

#include <complex>
#include <vector>
#include <string>

namespace fft {

    // remove?
    typedef std::vector<std::complex<double>> poly;

    // keep internal?
    size_t next_power_of_two(size_t x);

    // keep internal?
    std::complex<double> omega(double n, bool inverse = false);

    // move resizing, omega etc. to inside class definiton?
    // the "main" fft implementation.
    template<typename V=std::complex<double>>
        std::vector<V> rec_fft(std::vector<V> a, const V wn);

    // inverse fft---just a wrapper over some scaling logic etc.
    template<typename V=std::complex<double>>
        std::vector<V> ifft(std::vector<V> a, const V wn);

    // use the FFT algorithm to actually implement a polynomial
    // multipilcation algorithm
    template<typename V=std::complex<double>>
        std::vector<V> polymult(std::vector<V> a, std::vector<V> b); 

    // outputs a matlab-interperable string of a solution.
    // Idea is to reveal the computation as it happens...
    std::vector<std::string> str_fft(std::vector<std::string> a, bool inv = false);

}


#endif
