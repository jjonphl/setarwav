#include <Rcpp.h>

using namespace Rcpp;

static NumericVector _upsample(NumericVector x, int by) {
    int len = x.size();
    int outlen =  (len-1)*(by+1)+1;   

    NumericVector out(outlen, 0.0);
    double *xp = x.begin();
    double *outp = out.begin();
    int i, j;
    for (i = 0, j = 0; i < outlen; i += (by+1), j++) {
        outp[i] = xp[j];
    }

    return out;
}

static NumericVector _convolve(NumericVector x, NumericVector y) {
    int xlen = x.size();
    int ylen = y.size();
    int outlen = xlen + ylen - 1;

    NumericVector out(outlen, 0.0);

    double *xp = x.begin();
    double *yp = y.begin();
    double *outp = out.begin();

    int i, j;
    for (i = 0; i < xlen; i++) {
        for (j = 0; j < ylen; j++) {
            outp[i+j] += xp[i] * yp[j];
        }
    }

    return out;
}

static NumericMatrix _betahat(IntegerVector idx, int Lj, int j) {
    int idxlen = idx.size();
    NumericMatrix out(Lj, Lj);
    double *outp = out.begin();

    int *idxp = idx.begin();
    double acc;     // accumulator
    int i, k, l;

    for (i = 0; i < Lj; i++) {
        for (k = i; k < Lj; k++) {
            for (acc = 0.0, l = Lj-1; l < idxlen; l++) {
                acc += idxp[l-i] * idxp[l-k];
            } 
            outp[k*Lj+i] = outp[i*Lj+k] = acc;   // column-major
        }
    }

    Function warning("warning");
    int mj = idxlen - Lj + 1;
    std::stringstream msg;

    for (i = 0; i < Lj; i++) {
        for (k = i; k < Lj; k++) {
            acc = outp[k*Lj+i];
            if (acc != 0.0) {
                acc = mj / acc;
            } else {
                msg.str(std::string());    // clear string
                msg << "beta^{-1} = 0 at Lj=" << Lj 
                    << ", J=" << j << "(" << (i+1) << "," << (k+1) << ")"
                    << "; consider lowering number of levels (j)";
                warning(msg.str());
                
                acc = 1.0;
            }
            outp[k*Lj+i] = outp[i*Lj+k] = acc;
        }
    }

    return out;
}

// up-sample
RcppExport SEXP upsample(SEXP _x, SEXP _by) {
    NumericVector x(_x);
    IntegerVector by(_by);
    return wrap(_upsample(x, by[0]));
}

// convolve
RcppExport SEXP convolve(SEXP _x, SEXP _y) {
    NumericVector x(_x);
    NumericVector y(_y);
    return wrap(_convolve(x, y));
}

RcppExport SEXP betahat(SEXP _idx, SEXP _Lj, SEXP _j) {
    IntegerVector idx(_idx);
    IntegerVector Lj(_Lj);
    IntegerVector j(_j);

    return _betahat(idx, Lj[0], j[0]);
}

// x = time series
// h = wavelet filter
// g = scaling filter
// levels = # of levels
// idx = int array of 0/1, which elements are missing/available
RcppExport SEXP modwt2missing(SEXP _x, SEXP _h, SEXP _g, SEXP _levels, SEXP _idx) {
    NumericVector x(_x);
    NumericVector h(_h);
    NumericVector g(_g);
    IntegerVector levels(_levels);
    IntegerVector idx(_idx);

    //int taps = g.size();   // filter length
    int len = x.size();
    NumericVector gj, hj, gjtmp;   // level j scaling & wavelet filters
    double *hjp;

    int *idxp = idx.begin();
    double *xp = x.begin();

    int j, J = levels[0], tauj;  // tauj = tau_j = 2^(j-1)
    int lj, mj, t, l1, l2;

    NumericMatrix b;
    double acc;

    std::stringstream colname;
    List out;

    for (j = 0, tauj=1; j < J; j++, tauj *= 2) {
        if (j == 0) {       // level j+1
            gj = g;
            hj = h;
        } else {
            gjtmp = gj;
            gj = _convolve(gjtmp, _upsample(g, tauj-1));
            hj = _convolve(gjtmp, _upsample(h, tauj-1));
        }

        hjp = hj.begin();
        lj = hj.size();
        mj = len - lj + 1;
        NumericMatrix b = _betahat(idx, lj, j);
        double *bp = b.begin();

        NumericVector wj(mj, 0.0);
        double *wjp = wj.begin();

        for (t = lj-1; t < len; t++) {
            acc = 0.0;
            for (l1 = 0; l1 < lj; l1++) {
                // term == 0 when l1 == l2
                for (l2 = l1+1; l2 < lj; l2++) {
                    acc += (bp[l1*lj + l2] + bp[l2*lj + l1]) *
                           idxp[t-l1] * idxp[t-l2] * 
                           hjp[l1] * hjp[l2] * 
                           (xp[t-l1] - xp[t-l2]) * (xp[t-l1] - xp[t-l2]);
                }
            }

            wjp[t-lj+1] = -0.5 * acc;
        }

        colname.str(std::string());  // clear
        colname << "d" << (j+1);
        out[colname.str()] = wj;
    }

    return out;
}
