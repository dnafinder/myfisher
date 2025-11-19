function Pvalue = myfisher33(x)
%MYFISHER33 Fisher's Exact Probability Test on a 3x3 contingency table.
%
%   P = MYFISHER33(X) computes the two-tailed Fisher's exact test p-value
%   for a 3x3 contingency table X, where X is a 3-by-3 matrix of
%   nonnegative integers.
%
%   The function:
%     - Computes all possible 3x3 tables with the same row and column sums
%       (margins) as the observed table.
%     - For each possible table, computes its exact probability under the
%       null hypothesis of independence, using factorials in log-space.
%     - Sums the probabilities of all tables that are as or more extreme
%       (in probability) than the observed table to obtain a two-tailed
%       p-value.
%     - Also reports the Mid-p correction.
%
%   Syntax:
%       p = myfisher33(x)
%
%   Inputs:
%       X - 3x3 data matrix (nonnegative integers)
%
%   Outputs:
%       P - two-tailed Fisher exact p-value
%
%   Example:
%
%                A   B   C
%           -------------------
%      X         4   1   1
%           -------------------
%      Y         1   3   1
%           -------------------
%      Z         0   1   5
%           -------------------
%
%   Calling:
%       myfisher33([4 1 1; 1 3 1; 0 1 5])
%
%   Example output:
%
%   ------------------------------------------------------------------------
%   3x3 matrix Fisher's exact test
%   ------------------------------------------------------------------------
%        Tables    two_tails_p_value    Mid_p_correction
%        ______    _________________    ________________
%
%        301       0.033687             0.025239
%
%   Created by Giuseppe Cardillo
%   giuseppe.cardillo.75@gmail.com
%
%   To cite this file, a suitable format is:
%   Cardillo G. (2007) MyFisher33: a very compact routine for Fisher's exact
%   test on 3x3 matrix.
%   GitHub: https://github.com/dnafinder/myfisher
%
%   ----------------------------------------------------------------------
%   Algorithm overview:
%   ----------------------------------------------------------------------
%   Let X be a 3x3 table with row sums R1,R2,R3 and column sums C1,C2,C3.
%   A 3x3 table with fixed margins has 4 degrees of freedom.
%
%   We construct all feasible tables in two stages:
%
%   1) PRIMARY TABLES (first row):
%      We choose:
%        I = X(1,1)
%        J = X(1,2)
%        X(1,3) = R1 - I - J
%
%      subject to:
%        0 <= I <= min(R1, C1)
%        0 <= J <= min(R1 - I, C2)
%
%      Each choice of (I,J) defines a possible first row:
%        [X(1,1) X(1,2) X(1,3)] = [I J R1 - I - J]
%
%      This set of first rows is stored in m1 (et x 3).
%
%   2) SECONDARY TABLES (second row):
%      For each primary first row m1(k,:), we compute the residual column
%      sums after subtracting that row:
%
%        Cs2 = Cs - m1(k,:)
%
%      Then we choose the second row:
%        X(2,1) = a
%        X(2,2) = b
%        X(2,3) = R2 - a - b
%
%      subject to:
%        0 <= a <= min(R2, Cs2(1))          (column 1 and row 2)
%
%        For each a, b is constrained by:
%          - Nonnegativity of X(2,2): b >= 0
%          - Column 2 margin (residual):  b <= Cs2(2)
%          - Nonnegativity of X(2,3):     R2 - a - b >= 0  => b <= R2 - a
%          - Column 3 residual margin:    R2 - a - b <= Cs2(3)
%                => b >= R2 - a - Cs2(3)
%
%        So:
%          b_min = max(0, R2 - a - Cs2(3))
%          b_max = min(Cs2(2), R2 - a)
%
%        If b_min > b_max, no feasible b for that (a, first row).
%
%      The third row is then determined uniquely by the remaining column
%      margins:
%
%        X(3,j) = Cj - X(1,j) - X(2,j),  j = 1,2,3
%
%      For all feasible (I,J,a,b) combinations, we get a full 3x3 table
%      that satisfies all margins and nonnegativity constraints.
%
%   For each table Y, the probability under the null is:
%
%      P(Y) = [prod_i Ri! * prod_j Cj!] / [ N! * prod_{i,j} Y(i,j)! ]
%
%   We compute log-factorials via GAMMALN to avoid overflow:
%      log(n!) = gammaln(n+1).
%
%   The two-tailed p-value is the sum of probabilities of all tables with
%   probability <= P(observed).
%
%   The Mid-p correction is:
%      MidP = 0.5*P(observed) + sum_{P(Y) < P(observed)} P(Y)
%   ----------------------------------------------------------------------


%% Step 1: Input validation
p = inputParser;
addRequired(p, 'x', @(x) validateattributes(x, {'numeric'}, ...
    {'real','finite','integer','nonnegative','nonnan','size',[3 3]}));
parse(p, x);
clear p

%% Step 2: Basic margins
Rs = sum(x, 2);   % row sums: [R1; R2; R3]
Cs = sum(x, 1);   % column sums: [C1 C2 C3]
N  = sum(Rs);     % total number of observations

%% Step 3: Reorder rows and columns (if needed)
% For more compact enumeration and (sometimes) better numerical behavior,
% we sort:
%   - columns by increasing column sums
%   - rows by increasing row sums
% This does not affect the p-value since Fisher's test conditions on
% the margins only.
if ~issorted(Cs)
    [Cs, ind] = sort(Cs);
    x         = x(:, ind);
end
if ~issorted(Rs)
    [Rs, ind] = sort(Rs);
    x         = x(ind, :);
end

%% Step 4: Probability of the observed table

% Constant factor:
%   Kf = log( prod(Ri!)*prod(Cj!) / N! )
Kf = sum(gammaln([Rs' Cs] + 1)) - gammaln(N + 1);

% Log-factorials of the observed table entries:
zf_obs = gammaln(x + 1);

% Probability of the observed matrix:
op = exp(Kf - sum(zf_obs(:)));

%% Step 5: Enumerate all possible primary tables (first row)

% First row (row 1) sum is Rs(1).
% Let:
%   I = X(1,1)
%   J = X(1,2)
%
% Constraints:
%   0 <= I <= min(Rs(1), Cs(1))
%   0 <= J <= min(Rs(1) - I, Cs(2))
%
% Then X(1,3) = Rs(1) - I - J.

I = 0:1:min(Rs(1), Cs(1));        % all possible values of X(1,1)
J = min(Rs(1) - I, Cs(2));        % for each I, max X(1,2)

% Total number of primary tables:
et = sum(J + 1);

% m1 will store all possible first rows: [X(1,1) X(1,2) X(1,3)]
m1 = zeros(et, 3);

% Index vectors to fill m1 efficiently
idxAstop  = cumsum(J + 1);
idxAstart = [1, idxAstop(1:end-1) + 1];

% Build m1: enumerate all (I,J) pairs allowed
for K = 1:length(I)
    % For fixed I(K), J ranges from 0 to J(K).
    rows = idxAstart(K):idxAstop(K);
    m1(rows, 1) = I(K);
    m1(rows, 2) = 0:1:J(K);
end

% Complete X(1,3) so each row of m1 is a valid first row
m1(:,3) = Rs(1) - sum(m1(:,1:2), 2);


%% Step 6: Compute residual column sums after first row (Cs2)

% m2 is a replicated version of Cs to match the number of primary tables.
m2  = repmat(Cs, et, 1);

% Cs2(k,:) = remaining column sums after subtracting row 1 = m1(k,:)
% i.e., sums of rows 2+3 in each column.
Cs2 = m2 - m1;


%% Step 7: Enumerate all possible secondary tables (second row)

% For each primary table, we now enumerate all possible second rows.
% Let:
%   a = X(2,1)
%   b = X(2,2)
%   X(2,3) = Rs(2) - a - b
%
% For each primary row k, the range of a is:
%
%   0 <= a <= L(k)
%   where L(k) = min(Rs(2), Cs2(k,1))
%
L  = min(Rs(2), Cs2(:,1));
KK = max(L) + 1;    % we add 1 since a ranges from 0 to L(k)

% Preallocate matrices for upper/lower bounds of b
UpB   = zeros(et, KK);    % upper bounds for b
LoB   = zeros(et, KK);    % lower bounds for b
InT   = NaN(et, KK);      % length (count) of interval for b for each (k,a)
sectab = zeros(1, et);    % number of secondary tables per primary table

% Compute upper bounds for b for each possible a index (K-1)
% Here K-1 represents 'a', common for all primary tables in that column.
for K = 1:KK
    a_val = K - 1;  % a ranges from 0 to max(L)
    % For each primary row k, we have:
    %   b <= Cs2(k,2) (residual column 2)
    %   b <= Rs(2) - a_val (nonnegativity of X(2,3))
    UpB(:,K) = min(Cs2(:,2), Rs(2) - a_val);
end

% Compute lower bounds, interval lengths, and number of secondary tables
for K = 1:et
    % For primary row K:
    %   a ranges from 0 to L(K).
    % For each a, we define b_min and b_max:
    %
    %   b_min = max(0, Rs(2) - a - Cs2(K,3))
    %   b_max = UpB(K, a+1)
    %
    % So we fill LoB(K,1:(L(K)+1)) accordingly.
    z  = L(K) + 1;
    a  = 0:1:L(K);
    LoB(K,1:z) = max(0, Rs(2) - a - Cs2(K,3));
    % Length of each b interval:
    %   InT(K,j) = max(0, b_max - b_min + 1)
    InT(K,1:z) = UpB(K,1:z) - LoB(K,1:z) + 1;
    % Total number of secondary tables for this primary row:
    sectab(K) = sum(InT(K,1:z), 2);
end

% Total number of complete 3x3 tables:
et2 = sum(sectab);

% Preallocate matrix M to store all full tables (flattened 3x3 → 1x9)
M = zeros(et2, 9);


%% Step 8: Build full set of secondary tables and complete 3x3 tables

% Prepare linear index vectors for intervals of b for each (primary row, a)
InT  = InT';
InT(isnan(InT)) = [];
idxAstop  = cumsum(InT(:))';
idxAstart = [1, idxAstop(1:end-1) + 1];

% Index ranges for secondary tables per primary row
idxBstop  = cumsum(sectab);
idxBstart = [1, idxBstop(1:end-1) + 1];

% 'z' here is an index into the (a,b) interval slices
z_idx = 1;

for K = 1:et
    % Expand primary table: repeat the first row sectab(K) times
    rowsK = idxBstart(K):idxBstop(K);
    M(rowsK, 1:3) = repmat(m1(K,:), sectab(K), 1);

    % For this primary row, a ranges from 0 to L(K):
    a_vals = 0:1:L(K);

    % For each possible a, set X(2,1)=a and X(2,2)=b in the correct block
    for J = 1:length(a_vals)
        if idxAstart(z_idx) > idxAstop(z_idx)
            % No feasible b for this (primary row, a), skip
            z_idx = z_idx + 1;
            continue
        end
        rows_a = idxAstart(z_idx):idxAstop(z_idx);

        % X(2,1) = a_vals(J)
        M(rows_a, 4) = a_vals(J);

        % X(2,2) runs from LoB(K,J) to UpB(K,J)
        M(rows_a, 5) = LoB(K,J):1:UpB(K,J);

        z_idx = z_idx + 1;
    end
end

% Complete the second row in column 3:
%   X(2,3) = Rs(2) - X(2,1) - X(2,2)
M(:,6) = Rs(2) - sum(M(:,4:5), 2);

% Complete the third row from column margins:
%   X(3,j) = Cj - X(1,j) - X(2,j), j = 1,2,3
M(:,7:9) = repmat(Cs, et2, 1) - M(:,1:3) - M(:,4:6);


%% Step 9: Compute probabilities for all possible tables

% Compute log-factorials for all entries of all tables
zf_all = gammaln(M + 1);

% Probability of each table:
np = exp(Kf - sum(zf_all, 2));

% Two-tailed p-value: sum of probabilities <= that of the observed table
P = sum(np(np <= op));

% Mid-p correction:
MidP = 0.5*op + sum(np(np < op));


%% Step 10: Display results

tr = repmat('-', 1, 80);
disp(tr)
disp('3x3 matrix Fisher''s exact test')
disp(tr)
disp(array2table([et2, P, MidP], ...
    'VariableNames', {'Tables','two_tails_p_value','Mid_p_correction'}));


%% Step 11: Return value (if requested)

if nargout
    Pvalue = P;
end

end
