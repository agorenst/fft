#include "fft.h"

using namespace std;

size_t fft::next_power_of_two(size_t x) {
    size_t res = 1;
    while (res < x) {
        res *= 2;
    }
    return res;
}


complex<double> fft::omega(double n, bool inverse) {
    const double e = 2.718281;
    const double pi = 3.141592;
    const complex<double> i{0,1};
    if (!inverse) {
        return pow(e,(2.0*pi*i)/n);
    }
    else {
        return pow(e,(-2.0*pi*i)/n);
    }
}


template<typename V=complex<double>>
vector<V> fft::rec_fft(vector<V> a, const V wn) {
    if (a.size()==1) {
        return a;
    }

    vector<V> E(a.size()/2);
    vector<V> O(a.size()/2);
    for (unsigned i = 0; i < a.size(); ++i) {
        if (i % 2 == 0) {
            E[i/2] = a[i];
        }
        else {
            O[(i-1)/2] = a[i];
        }
    }

    auto ERES = rec_fft(E, wn*wn);
    auto ORES = rec_fft(O, wn*wn);

    vector<V> ARES(a.size());
    V w{1,0};
    for (unsigned i = 0; i < a.size()/2; ++i) {
        ARES[i]              = ERES[i] + (w*ORES[i]);
        ARES[i+(a.size()/2)] = ERES[i] - (w*ORES[i]);
        w *= wn;
    }
    return ARES;
}

template<typename V=complex<double>>
vector<V> fft::polymult(vector<V> a, vector<V> b) {
    auto size = max(a.size(), b.size());
    auto newsize = next_power_of_two(size)*2;
    a.resize(newsize,0);
    b.resize(newsize,0);
    auto w = omega(newsize);
    auto winv = omega(newsize,true);

    auto aeval = rec_fft(a,w);
    auto beval = rec_fft(b,w);
    auto r = vector<V>(newsize);
    for (unsigned i = 0; i < newsize; ++i) {
        r[i] = aeval[i]*beval[i];
    }
    r = ifft(r,winv);
    //for (unsigned i = 0; i < newsize; ++i) {
    //    if (imag(r[i]) < threshold) {
    //        r[i] = {real(r[i]), 0};
    //    }
    //}
    return r;
}

template<typename V=std::complex<double>>
std::vector<V> fft::ifft(std::vector<V> a, const V wn) {
    auto b = rec_fft(a, wn);
    for (unsigned i = 0; i < b.size(); ++i) {
        b[i] = b[i]/(std::complex<double>(b.size(),0));
    }
    return b;
}

template fft::poly fft::polymult(poly a, poly b);
template fft::poly fft::rec_fft(poly a, complex<double> wn);

string str_omega(int count, int base) {
    if (count == base || count == 0) {
        return "1";
    }
    if (count * 2 == base) {
        return "-1";
    }
    string cs, bs;
    stringstream csmaker(cs);
    stringstream bsmaker(bs);
    csmaker << count;
    bsmaker << base;
    return "exp(2*i*pi*"+csmaker.str()+"/"+bsmaker.str()+")";
}

vector<string> fft::str_fft(vector<string> a, bool inv) {
    if (a.size() == 1) { return a; }
    vector<string> e(a.size()/2);
    vector<string> o(a.size()/2);
    int omegabase = a.size();
    int omegacount = 0;
    for (unsigned int i = 0; i < a.size(); ++i) {
        if (i % 2 == 0) {
            e[i/2] = a[i];
        }
        else {
            o[(i-1)/2] = a[i];
        }
    }
    vector<string> E = str_fft(e, inv);
    vector<string> O = str_fft(o, inv);
    vector<string> A(a.size());
    for (unsigned int i = 0; i < a.size()/2; ++i) {
        A[i]       = "(" + E[i] + "+(" + str_omega(omegacount, omegabase) + "*" + O[i] + "))";
        A[i+(a.size()/2)] = "(" + E[i] + "-(" + str_omega(omegacount, omegabase) + "*" + O[i] + "))";
        if (inv) { omegacount--; }
        else { omegacount++; }
    }
    return A;
}

//ostream& operator<<(ostream& o, vector<complex<double>> a) {
//    for (auto it = a.begin(); it != a.end(); ++it) {
//        o << *it << " ";
//    }
//    return o;
//}

//void matrix_printer(const int n) {
//    for (int i = 0; i < n; ++i) {
//        for (int j = 0; j < n; ++j) {
//            // omega to the j times omega to the i
//            cout << (j*i) % n << " ";
//        }
//        cout << endl;
//    }
//}
