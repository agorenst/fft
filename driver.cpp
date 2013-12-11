#include <complex>
#include <iostream>
#include <vector>

#include "fft.h"

using namespace std;

using namespace fft;
std::ostream& operator<<(ostream& o, vector<complex<double>> a) {
    for (auto it = a.begin(); it != a.end(); ++it) {
        o << *it << " ";
    }
    return o;
}

vector<complex<double>> parse_poly(string s) {
    stringstream sstr(s);
    vector<complex<double>> coefs;
    while (!sstr.eof()) {
        double real, imag;
        sstr >> real;
        sstr >> imag;
        coefs.push_back(complex<double>(real,imag));
    }
    auto extend_size = next_power_of_two(coefs.size());
    coefs.resize(extend_size, 0);
    return coefs;
}

void evaluate_poly(bool inv = false) {
    string s;
    getline(cin, s);
    auto coefs = parse_poly(s);
    vector<complex<double>> res;
    if (!inv) {
        res = rec_fft(coefs, omega(coefs.size(), inv));
    }
    else {
        res = ifft(coefs, omega(coefs.size(), inv));
    }
    for (auto c : res) {
        cout << c << "\t";
    }
    cout << endl;
}

void mult_poly() {
    string s1, s2;
    cerr << "yo"  << endl;
    getline(cin, s1);
    cerr << s1 << endl;
    getline(cin, s2);
    cerr << s2 << endl;
    auto p1 = parse_poly(s1);
    auto p2 = parse_poly(s2);
    auto res = polymult(p1,p2);
    cout << res << endl;
}

void evaluate_string() {
    string s, w;
    getline(cin, s);
    vector<string> coefs;
    stringstream parser(s);
    while(!parser.eof()) {
        parser >> w;
        coefs.push_back(w);
    }
    auto res = str_fft(coefs);
    cout << "[";
    for (auto c : res) {
        cout << c << ",";
    }
    cout << "]" << endl;
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        evaluate_poly();
    }
    if (argc == 2) {
        if (string(argv[1]) == "inv") {
            evaluate_poly(true);
        }
        if (string(argv[1]) == "mult") {
            mult_poly();
        }
        if (string(argv[1]) == "str") {
            evaluate_string();
        }
    }
}
