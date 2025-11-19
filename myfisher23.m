function Pvalue = myfisher23(x)
%MYFISHER23 Fisher's Exact Probability Test on a 2x3 contingency table.
%
%   P = MYFISHER23(X) computes the two-tailed Fisher's exact test p-value
%   for a 2x3 contingency table X, where X is a 2-by-3 matrix of
%   nonnegative integers.
%
%   The function:
%     - Computes all possible 2x3 tables with the same row and column sums
%       (margins) as the observed table.
%     - For each possible table, computes its exact probability under the
%       null hypothesis of independence, using factorials in log-space.
%     - Sums the probabilities of all tables that are as or more extreme
%       (in probability) than the observed table to obtain a two-tailed
%       p-value.
%
%   Example:
%
%                A   B   C
%           -------------------
%      X         0   3   2
%           -------------------
%      Y         6   5   1
%           -------------------
%
%   Calling:
%       myfisher23([0 3 2; 6 5 1])
%
%   Output:
%
%   ------------------------------------------------------------------------
%   2x3 matrix Fisher's exact test
%   ------------------------------------------------------------------------
%        Tables    two_tails_p_value    Mid_p_correction
%        ______    _________________    ________________
%
%        18        0.088235             0.074661
%
%   Created by Giuseppe Cardillo
%   giuseppe.cardillo.75@gmail.com
%
%   To cite this file, a suitable format is:
%   Cardillo G. (2007) MyFisher23: a very compact routine for Fisher's exact
%   test on 2x3 matrix.
%   GitHub: https://github.com/dnafinder/myfisher
%
%   ----------------------------------------------------------------------
%   Algorithm overview:
%   ----------------------------------------------------------------------
%   Let X be a 2x3 table with row sums R1, R2 and column sums C1, C2, C3.
%   A 2x3 table with fixed margins has 2 degrees of freedom.
%
%   We parameterize the first row as:
%
%       X(1,1) = a
%       X(1,2) = b
%       X(1,3) = R1 - a - b
%
%   The second row is then completely determined by the column sums:
%
%       X(2,j) = Cj - X(1,j),  j = 1,2,3
%
%   We enumerate all integer pairs (a, b) that satisfy:
%     - Nonnegativity of all cells
%     - Row and column margins exactly
%
%   Constraints:
%     1) 0 <= a <= min(R1, C1)
%     2) For a fixed a, b is bounded by:
%          - Row 1:           a + b <= R1         -> b <= R1 - a
%          - Column 2:        b <= C2
%          - Nonnegativity of X(2,3):
%                X(2,3) = C3 - (R1 - a - b) >= 0
%                => C3 >= R1 - a - b
%                => b >= R1 - a - C3
%
%        Overall:
%          b_min(a) = max(0, R1 - a - C3)
%          b_max(a) = min(C2, R1 - a)
%
%        If b_min(a) > b_max(a), there are no feasible tables for that a.
%
%   Each table's probability is:
%
%       P(Y) = [prod_i R_i! * prod_j C_j!] / [ N! * prod_{i,j} Y(i,j)! ]
%
%   We work in log-space using gammaln to avoid overflow:
%       log(n!) = gammaln(n+1).
%
%   The two-tailed p-value is the sum of probabilities of all tables with
%   probability less than or equal to that of the observed table.
%   A Mid-p correction is also reported.
%   ----------------------------------------------------------------------

%% Step 1: Input validation
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x, {'numeric'}, ...
    {'real','finite','integer','nonnegative','nonnan','size',[2 3]}));
parse(p,x);
clear p

%% Step 2: Basic margins
Rs = sum(x, 2);   % Row sums: [R1; R2]
Cs = sum(x, 1);   % Column sums: [C1 C2 C3]
N  = sum(Rs);     % Total count

%% Step 3: Optional reordering of rows and columns
% For numerical stability and a more compact enumeration, we sort:
%   - columns by increasing column sums
%   - rows by increasing row sums
% This does NOT change the test, only the enumeration order.
if ~issorted(Cs)
    [Cs, ind] = sort(Cs);
    x         = x(:, ind);
    clear ind
end
if ~issorted(Rs)
    [Rs, ind] = sort(Rs);
    x         = x(ind, :);
    clear ind
end

% After sorting:
R1 = Rs(1);
C1 = Cs(1);
C2 = Cs(2);
C3 = Cs(3);

%% Step 4: Determine the range for X(1,1) = a

% The first cell in the top row (column 1) can range from 0 up to the
% minimum of its row margin and its column margin:
%   0 <= a <= min(R1, C1)
A = 0:1:min(R1, C1);  % all possible values of X(1,1)

%% Step 5: For each a, determine the valid range of X(1,2) = b

% For each a, the maximum value of b is limited by:
%   - row 1 margin: a + b <= R1        -> b <= R1 - a
%   - column 2 margin: b <= C2
%
% So the upper bound is:
%   b_max(a) = min(C2, R1 - a)

b_max = min(C2, R1 - A);

% The lower bound of b is controlled by nonnegativity of X(2,3):
%
%   X(1,3) = R1 - a - b
%   X(2,3) = C3 - X(1,3) = C3 - (R1 - a - b)
%          = C3 - R1 + a + b
%
% We require X(2,3) >= 0:
%
%   C3 - R1 + a + b >= 0
%   b >= R1 - a - C3
%
% and of course b >= 0.
%
% Therefore:
%   b_min(a) = max(0, R1 - a - C3)

b_min = max(R1 - A - C3, 0);

% Number of feasible b values for each a:
%   count(a) = max(0, b_max(a) - b_min(a) + 1)
count_per_a = b_max - b_min + 1;
count_per_a(count_per_a < 0) = 0;

% Total number of feasible tables
et = sum(count_per_a);
if et <= 0
    error('No valid tables can be generated with the given margins.');
end

%% Step 6: Allocate the table matrix and compute index ranges

% Tables will hold all feasible 2x3 tables, flattened row-wise:
%   Tables(:,1:3) = first row  [X(1,1) X(1,2) X(1,3)]
%   Tables(:,4:6) = second row [X(2,1) X(2,2) X(2,3)]
Tables = zeros(et, 6);

% Compute start/stop row indices in Tables for each a
stop_idx  = cumsum(count_per_a);
start_idx = [1, stop_idx(1:end-1) + 1];

%% Step 7: Fill X(1,1) = a and X(1,2) = b for all tables

% For the first value of a, all entries in column 1 are already zero,
% so we only need to fill column 2.
current_a = A(1);
Tables(start_idx(1):stop_idx(1), 2) = b_min(1):1:b_max(1);

% For subsequent values of a, fill column 1 and column 2
for I = 2:length(A)
    if count_per_a(I) == 0
        % No feasible b for this a, skip
        continue
    end

    current_a = A(I);
    i_start   = start_idx(I);
    i_stop    = stop_idx(I);

    % X(1,1) = a
    Tables(i_start:i_stop, 1) = current_a;

    % X(1,2) = all feasible b between b_min(I) and b_max(I)
    Tables(i_start:i_stop, 2) = b_min(I):1:b_max(I);
end

%% Step 8: Complete X(1,3) and the second row

% Complete X(1,3) using top row margin:
%   X(1,3) = R1 - X(1,1) - X(1,2)
Tables(:, 3) = R1 - sum(Tables(:, 1:2), 2);

% Complete second row from column margins:
%   X(2,j) = Cj - X(1,j)
Tables(:, 4:6) = repmat(Cs, et, 1) - Tables(:, 1:3);

% For the 2x3 case, by construction all entries are >= 0
% and satisfy the given margins exactly.

%% Step 9: Compute exact probabilities for all tables

% Use log-factorials via gammaln:
%   log(n!) = gammaln(n+1)
zf = gammaln(Tables + 1);

% Constant factor:
%   K = log( prod(Rs!)*prod(Cs!) / N! )
K  = sum(gammaln([Rs' Cs] + 1)) - gammaln(N + 1);

% Probability of each table:
np = exp(K - sum(zf, 2));

%% Step 10: Locate the observed table among generated tables

% Match the first row of the observed table (after any sorting)
% with Tables(:,1:3).
[~, obt] = ismember(x(1,:), Tables(:, 1:3), 'rows');

if obt == 0
    error('Observed table not found among generated tables. Check margins.');
end

%% Step 11: Compute two-tailed p-value and Mid-p correction

% Two-tailed p-value:
%   sum of probabilities of all tables with probability <= that of observed
P = sum(np(np <= np(obt)));

% Mid-p correction:
%   MidP = 0.5 * P(observed) + sum_{P(Y) < P(observed)} P(Y)
MidP = 0.5*np(obt) + sum(np(np < np(obt)));

%% Step 12: Display results

tr = repmat('-', 1, 80);
disp(tr)
disp('2x3 matrix Fisher''s exact test')
disp(tr)
disp(array2table([et, P, MidP], ...
    'VariableNames', {'Tables','two_tails_p_value','Mid_p_correction'}));

%% Step 13: Return output (if requested)

if nargout
    Pvalue = P;
end

end
