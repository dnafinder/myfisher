[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/myfisher)

📌 Overview
MyFisher is a unified MATLAB toolkit for Fisher’s Exact Test on RxC contingency tables. It includes exact algorithms for 2×2, 2×3, 2×4 and 3×3 tables, plus a Monte Carlo conditional method for larger matrices, providing accurate two-tailed and Mid-p values.

📂 Contents
myfisher.m General RxC Fisher’s exact / conditional test with automatic method selection and Monte Carlo sampling for larger tables.
myfisher22.m Exact Fisher’s test for 2×2 tables, including left- and right-tail p-values, two-tailed p-value, Mid-p correction, power analysis and recommended sample size.
myfisher23.m Exact Fisher’s test for 2×3 tables with full enumeration of all tables consistent with the observed margins.
myfisher24.m Exact Fisher’s test for 2×4 tables using vectorized enumeration for computational efficiency.
myfisher33.m Exact Fisher’s test for 3×3 tables with efficient generation of all feasible tables under fixed margins.

🧮 Methods
For small tables (2×2, 2×3, 2×4, 3×3), MyFisher enumerates all possible contingency tables with the same row and column sums as the observed one, computes their exact probabilities in log-space using gammaln, and sums the probabilities of tables as or more extreme than the observed table. For larger RxC tables, MyFisher uses a Monte Carlo conditional test that preserves margins by shuffling individual-level column labels.

🔧 Usage
Place all .m files on your MATLAB path, then call:

P = myfisher(X);

for a general RxC table X, or:

P = myfisher22(X); % 2×2
P = myfisher23(X); % 2×3
P = myfisher24(X); % 2×4
P = myfisher33(X); % 3×3

In the general RxC function, optional arguments control Monte Carlo accuracy:

P = myfisher(X, delta); % set half-width of the Monte Carlo error
P = myfisher(X, delta, alpha); % set error and confidence level

📊 Output
The exact routines return the two-tailed Fisher p-value (and, where implemented, Mid-p correction and tail-specific p-values). The 2×2 routine further reports power and suggested sample sizes based on an asymptotic normal approximation with continuity correction. The general RxC Monte Carlo routine prints the number of simulated tables, the estimated p-value and the guaranteed accuracy interval.

📎 Citation
If you use MyFisher in scientific work, please cite it as:

Cardillo G. (2007–2010). MyFisher: a compact and complete family of routines for Fisher’s Exact Test on RxC contingency tables. GitHub: https://github.com/dnafinder/myfisher

✉️ Contact
Author: Giuseppe Cardillo
E-mail: giuseppe.cardillo.75@gmail.com
