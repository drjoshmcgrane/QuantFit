#include <algorithm>
#include <cmath>
#include <Rcpp.h>
#include <sstream>
#include <vector>

using Rcpp::clone;
using Rcpp::IntegerMatrix;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::runif;

using std::max;
using std::min;
using std::stringstream;
using std::vector;

// [[Rcpp::export]]
List CCIterateSingle(int nIter, NumericMatrix old, NumericMatrix old_ll,
                     int burn, IntegerMatrix N, NumericMatrix n)
{
    // Single cancellation only - works for any matrix size
    List chain;
    int nrow = old.nrow();
    int ncol = old.ncol();

    for (int I = 2; I <= nIter; ++I)
    {
        for (int i = 0; i < nrow; ++i)
        {
            for (int j = 0; j < ncol; ++j)
            {
                // Single cancellation constraints only
                double lh1 = (j == 0) ? 0.0 : old(i, j-1);
                double lh2 = (i == 0) ? 0.0 : old(i-1, j);
                double rh1 = (j == ncol-1) ? 1.0 : old(i, j+1);
                double rh2 = (i == nrow-1) ? 1.0 : old(i+1, j);

                double lh = max(lh1, lh2);
                double rh = min(rh1, rh2);
                if (rh < lh) rh = 1.0;

                double draw = runif(1, lh, rh)(0);
                double ar = 2;
                double new_ll = n(i,j)*log(draw) + (N(i,j) - n(i,j))*log(1.0 - draw);

                if (old(i,j) != 1.0 && old(i,j) != 0.0)
                    ar = exp(new_ll - old_ll(i,j));
                if (ar > runif(1)(0))
                {
                    old(i,j) = draw;
                    old_ll(i,j) = new_ll;
                }
            }
        }

        if (I > burn && I % 4 == 0)
        {
            stringstream hoop;
            hoop << I;
            chain[hoop.str()] = clone(old);
        }
    }
    return chain;
}

// [[Rcpp::export]]
List CCIterateDouble(int nIter, NumericMatrix old, NumericMatrix old_ll,
                     int burn, IntegerMatrix N, NumericMatrix n)
{
    // Single + Double cancellation for 3x3 matrices
    List chain;

    for (int I = 2; I <= nIter; ++I)
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                // Single cancellation constraints
                double lh1 = (j == 0) ? 0.0 : old(i, j-1);
                double lh2 = (i == 0) ? 0.0 : old(i-1, j);
                double rh1 = (j == 2) ? 1.0 : old(i, j+1);
                double rh2 = (i == 2) ? 1.0 : old(i+1, j);

                // Double cancellation constraints
                double lh3 = 0.0;
                double rh3 = 1.0;
                bool test1 = old(1,0) < old(0,1);
                bool test2 = old(2,1) < old(1,2);
                if (test1 && test2)
                {
                    if (i == 0 && j == 2) lh3 = old(2,0);
                    else if (i == 2 && j == 0) rh3 = old(0,2);
                }
                else if (!test1 && !test2)
                {
                    if (i == 2 && j == 0) lh3 = old(0,2);
                    else if (i == 0 && j == 2) rh3 = old(2,0);
                }

                double lh = max(max(lh1, lh2), lh3);
                double rh = (rh3 > lh) ? min(min(rh1, rh2), rh3) : min(rh1, rh2);
                if (rh < lh) rh = 1.0;

                double draw = runif(1, lh, rh)(0);
                double ar = 2;
                double new_ll = n(i,j)*log(draw) + (N(i,j) - n(i,j))*log(1.0 - draw);

                if (old(i,j) != 1.0 && old(i,j) != 0.0)
                    ar = exp(new_ll - old_ll(i,j));
                if (ar > runif(1)(0))
                {
                    old(i,j) = draw;
                    old_ll(i,j) = new_ll;
                }
            }
        }

        if (I > burn && I % 4 == 0)
        {
            stringstream hoop;
            hoop << I;
            chain[hoop.str()] = clone(old);
        }
    }
    return chain;
}

// [[Rcpp::export]]
List CCIterateTriple(int nIter, NumericMatrix old, NumericMatrix old_ll,
                     int burn, IntegerMatrix N, NumericMatrix n)
{
    // Single + Double + Triple cancellation for 4x4 matrices
    List chain;

    // Define all 16 embedded 3x3 submatrices for double cancellation
    // Row combinations: {0,1,2}, {0,1,3}, {0,2,3}, {1,2,3} (0-indexed)
    // Col combinations: same
    int row_combos[4][3] = {{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}};
    int col_combos[4][3] = {{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}};

    // Define all 14 coherent triple cancellation tests from Kyngdon & Richards (2006)
    // Each test: ant1(r1,c1,r2,c2), ant2(...), ant3(...), conseq(r1,c1,r2,c2)
    // Using 0-indexed positions
    int tc_tests[14][16] = {
        {1,0,0,1, 2,1,1,2, 3,2,2,3, 3,0,0,3},
        {1,0,0,1, 2,1,1,2, 3,1,2,2, 3,0,0,2},
        {1,0,0,1, 2,0,1,1, 3,2,2,3, 3,0,0,3},
        {1,0,0,1, 2,0,1,1, 3,1,2,2, 3,0,0,2},
        {1,0,0,2, 2,2,1,3, 3,1,2,2, 3,0,0,3},
        {1,0,0,2, 2,2,1,3, 3,0,2,1, 3,1,0,3},
        {1,0,0,2, 2,1,1,2, 3,2,2,3, 3,0,0,3},
        {1,0,0,2, 2,1,1,2, 3,0,2,1, 3,1,0,3},
        {1,0,0,2, 2,0,1,1, 3,2,2,3, 3,0,0,3},
        {1,1,0,2, 2,2,1,3, 3,0,2,1, 3,0,0,3},
        {1,1,0,2, 2,0,1,1, 3,2,2,3, 3,0,0,3},
        {1,1,0,2, 2,0,1,1, 3,1,2,2, 3,0,0,3},
        {1,0,0,1, 2,2,1,3, 3,1,2,2, 3,0,0,3},
        {1,0,0,1, 2,1,1,3, 3,2,2,3, 3,0,0,3}
    };

    for (int I = 2; I <= nIter; ++I)
    {
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                // Single cancellation constraints
                double lh1 = (j == 0) ? 0.0 : old(i, j-1);
                double lh2 = (i == 0) ? 0.0 : old(i-1, j);
                double rh1 = (j == 3) ? 1.0 : old(i, j+1);
                double rh2 = (i == 3) ? 1.0 : old(i+1, j);

                // Double cancellation constraints across all 16 embedded 3x3 submatrices
                double lh3 = 0.0;
                double rh3 = 1.0;
                for (int ri = 0; ri < 4; ++ri)
                {
                    for (int ci = 0; ci < 4; ++ci)
                    {
                        int r0 = row_combos[ri][0];
                        int r1 = row_combos[ri][1];
                        int r2 = row_combos[ri][2];
                        int c0 = col_combos[ci][0];
                        int c1 = col_combos[ci][1];
                        int c2 = col_combos[ci][2];

                        bool test1 = old(r1,c0) < old(r0,c1);
                        bool test2 = old(r2,c1) < old(r1,c2);
                        if (test1 && test2)
                        {
                            if (i == r0 && j == c2) lh3 = max(lh3, old(r2,c0));
                            if (i == r2 && j == c0) rh3 = min(rh3, old(r0,c2));
                        }
                        else if (!test1 && !test2)
                        {
                            if (i == r2 && j == c0) lh3 = max(lh3, old(r0,c2));
                            if (i == r0 && j == c2) rh3 = min(rh3, old(r2,c0));
                        }
                    }
                }

                // Triple cancellation constraints
                double lh4 = 0.0;
                double rh4 = 1.0;
                for (int t = 0; t < 14; ++t)
                {
                    int a1r1 = tc_tests[t][0], a1c1 = tc_tests[t][1];
                    int a1r2 = tc_tests[t][2], a1c2 = tc_tests[t][3];
                    int a2r1 = tc_tests[t][4], a2c1 = tc_tests[t][5];
                    int a2r2 = tc_tests[t][6], a2c2 = tc_tests[t][7];
                    int a3r1 = tc_tests[t][8], a3c1 = tc_tests[t][9];
                    int a3r2 = tc_tests[t][10], a3c2 = tc_tests[t][11];
                    int cr1 = tc_tests[t][12], cc1 = tc_tests[t][13];
                    int cr2 = tc_tests[t][14], cc2 = tc_tests[t][15];

                    bool t1 = old(a1r1,a1c1) < old(a1r2,a1c2);
                    bool t2 = old(a2r1,a2c1) < old(a2r2,a2c2);
                    bool t3 = old(a3r1,a3c1) < old(a3r2,a3c2);

                    if (t1 && t2 && t3)
                    {
                        if (i == cr2 && j == cc2) lh4 = max(lh4, old(cr1,cc1));
                        if (i == cr1 && j == cc1) rh4 = min(rh4, old(cr2,cc2));
                    }
                    else if (!t1 && !t2 && !t3)
                    {
                        if (i == cr1 && j == cc1) lh4 = max(lh4, old(cr2,cc2));
                        if (i == cr2 && j == cc2) rh4 = min(rh4, old(cr1,cc1));
                    }
                }

                // Combine all constraints
                double lh = max(max(max(lh1, lh2), lh3), lh4);
                double rh_dc = min(rh3, rh4);
                double rh = (rh_dc > lh) ? min(min(rh1, rh2), rh_dc) : min(rh1, rh2);
                if (rh < lh) rh = 1.0;

                double draw = runif(1, lh, rh)(0);
                double ar = 2;
                double new_ll = n(i,j)*log(draw) + (N(i,j) - n(i,j))*log(1.0 - draw);

                if (old(i,j) != 1.0 && old(i,j) != 0.0)
                    ar = exp(new_ll - old_ll(i,j));
                if (ar > runif(1)(0))
                {
                    old(i,j) = draw;
                    old_ll(i,j) = new_ll;
                }
            }
        }

        if (I > burn && I % 4 == 0)
        {
            stringstream hoop;
            hoop << I;
            chain[hoop.str()] = clone(old);
        }
    }
    return chain;
}
