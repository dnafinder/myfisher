function Pout = myfisher(x, varargin)
%MYFISHER Fisher's Exact / Conditional Test for a general RxC contingency table.
%
%   P = MYFISHER(X) performs Fisher's exact test on a contingency table X.
%   For small tables (2x2, 2x3, 2x4, 3x3), it calls optimized, dedicated
%   functions. For larger tables, it uses a conventional Monte Carlo
%   conditional test to estimate the two-tailed p-value.
%
%   P = MYFISHER(X, DELTA) uses DELTA to control the precision of the
%   Monte Carlo estimate (for large tables).
%
%   P = MYFISHER(X, DELTA, ALPHA) also specifies the confidence level
%   (1 - ALPHA) for the Monte Carlo accuracy guarantee.
%
%   ------------------------------------------------------------------
%   Strategy:
%   ------------------------------------------------------------------
%   1) If X is:
%        - 2x2 → use MYFISHER22
%        - 2x3 → use MYFISHER23
%        - 2x4 → use MYFISHER24
%        - 3x3 → use MYFISHER33
%      these specialized routines (written by the same author) compute
%      the exact Fisher p-value by enumerating all feasible tables.
%
%   2) Otherwise, MYFISHER uses a Monte Carlo conditional test:
%        - Expand the RxC table into an N×2 “individual-level” matrix,
%          where each row is (row_index, column_index).
%        - Randomly shuffle the column indices to generate new tables
%          with the same row and column margins (Fisher–Yates shuffle).
%        - For each simulated table, compute its probability and count
%          how many are as or more extreme than the observed one.
%        - The p-value is estimated as:
%               p_hat = (# of P_sim <= P_obs) / #simulations
%
%      The number of simulations is chosen so that:
%           P( |p_hat - p_true| <= DELTA ) ≥ 1 - ALPHA
%      (Psychometrika 1979; Vol. 44: 75–83).
%
%   Syntax:
%       P = myfisher(X)
%       P = myfisher(X, DELTA)
%       P = myfisher(X, DELTA, ALPHA)
%
%   Inputs:
%       X     - RxC data matrix (nonnegative integers)
%       DELTA - (optional) half-width of the Monte Carlo error interval,
%               i.e. desired |p_hat − p_true| ≤ DELTA.
%               Default: 0.01
%       ALPHA - (optional) significance level for the accuracy guarantee:
%               P( |p_hat − p_true| ≤ DELTA ) ≥ 1 − ALPHA.
%               Default: 0.05
%
%   Output:
%       P     - two-tailed p-value (exact for small tables; Monte Carlo
%               estimate for larger tables).
%
%   ------------------------------------------------------------------
%   Author and citation:
%   ------------------------------------------------------------------
%   Created by:  Giuseppe Cardillo
%   E-mail:      giuseppe.cardillo.75@gmail.com
%
%   To cite this file:
%   Cardillo G. (2010) MyFisher: the definitive function for the Fisher's exact
%   and conditional test for any RxC matrix.
%   GitHub: https://github.com/dnafinder/myfisher
%
%   ------------------------------------------------------------------


%% Step 1: Input error handling and parsing
p = inputParser;
addRequired(p, 'x', @(x) validateattributes(x, {'numeric'}, ...
    {'real','finite','integer','nonnegative','nonnan','2d'}));
addOptional(p, 'delta', 0.01, @(x) validateattributes(x, {'numeric'}, ...
    {'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p, 'alpha', 0.05, @(x) validateattributes(x, {'numeric'}, ...
    {'scalar','real','finite','nonnan','>',0,'<',1}));
parse(p, x, varargin{:});
delta = p.Results.delta;
alpha = p.Results.alpha;
clear p


%% Step 2: Check if we can use a specialized exact function
[rows, columns] = size(x);

if rows == 2 && columns == 2
    % 2x2: use myfisher22 (exact Fisher + power)
    assert(exist('myfisher22.m', 'file') ~= 0, ...
        ['You must download myfisher22.m from ', ...
         'https://github.com/dnafinder/myfisher22']);
    P_vec = myfisher22(x, alpha);  % returns [left right two_tails]
    P     = P_vec(3);              % two-tailed p-value

elseif rows == 2 && columns == 3
    % 2x3: use myfisher23 (exact Fisher)
    assert(exist('myfisher23.m', 'file') ~= 0, ...
        ['You must download myfisher23.m from ', ...
         'https://github.com/dnafinder/myfisher23']);
    P = myfisher23(x);             % two-tailed p-value

elseif rows == 2 && columns == 4
    % 2x4: use myfisher24 (exact Fisher)
    assert(exist('myfisher24.m', 'file') ~= 0, ...
        ['You must download myfisher24.m from ', ...
         'https://github.com/dnafinder/myfisher24']);
    P = myfisher24(x);             % two-tailed p-value

elseif rows == 3 && columns == 3
    % 3x3: use myfisher33 (exact Fisher)
    assert(exist('myfisher33.m', 'file') ~= 0, ...
        ['You must download myfisher33.m from ', ...
         'https://github.com/dnafinder/myfisher33']);
    P = myfisher33(x);             % two-tailed p-value

else
    %% Step 3: General RxC case – Monte Carlo conditional test

    % Compute margins and the constant factor K in log-space:
    C = sum(x, 1);   % column sums
    R = sum(x, 2);   % row sums
    N = sum(x(:));   % total count

    % Kf = log( prod(R!)*prod(C!) / N! )
    Kf = sum(gammaln([R' C] + 1)) - gammaln(N + 1);

    % Log-factorials of observed table:
    zf = gammaln(x + 1);

    % Probability of the observed table:
    op = exp(Kf - sum(zf(:)));

    % ------------------------------------------------------------------
    % Step 3.1: Transform RxC into an N×2 “individual-level” matrix
    % ------------------------------------------------------------------
    %
    % We construct an N×2 matrix "table" whose rows are (row_index, col_index).
    % Each cell (i,j) with count x(i,j) contributes x(i,j) rows equal to [i,j].
    %
    % This is done using built-in NDGRID + REPELEM (no external functions).

    [row_ids, col_ids] = ndgrid(1:size(x,1), 1:size(x,2));
    row_ids = row_ids(:);
    col_ids = col_ids(:);
    counts  = x(:);

    % REPELEM will skip entries where counts == 0
    table = [repelem(row_ids, counts), repelem(col_ids, counts)];

    % (By construction, size(table,1) == N.)

    clear R C row_ids col_ids counts zf


    % ------------------------------------------------------------------
    % Step 3.2: Monte Carlo sampling
    % ------------------------------------------------------------------
    %
    % We repeatedly shuffle the second column of "table" to generate new
    % RxC tables with the same margins as the original.
    %
    % For each shuffled table:
    %   - reconstruct the RxC table g
    %   - compute its probability P(g)
    %   - count how many times P(g) <= P(observed)
    %
    % The simulation size "tbs" is chosen to ensure that the Monte Carlo
    % estimate of the p-value is within DELTA units of the true p-value
    % with (1 - ALPHA)*100% confidence (Psychometrika 1979; Vol. 44: 75–83).

    tbs = round(((-realsqrt(2) * erfcinv(2 - alpha)) / (2 * delta))^2);
    MCC = 0;  % Monte Carlo counter

    for I = 1:tbs
        % Fisher–Yates shuffle (Sattolo-like) on the second column.
        % This is O(N) and usually faster than RANDPERM (O(N log N)).

        for J = N:-1:2
            s = ceil((J - 1) .* rand);
            tmp          = table(s, 2);
            table(s, 2)  = table(J, 2);
            table(J, 2)  = tmp;
        end

        % Reconstruct the RxC table g from the shuffled "table" matrix
        g = zeros(size(x));
        % Explicit loop: for each individual, increment the corresponding cell
        for J = 1:N
            g(table(J,1), table(J,2)) = g(table(J,1), table(J,2)) + 1;
        end

        % Compute log-factorials and probability of the new table
        zf  = gammaln(g + 1);
        gpv = exp(Kf - sum(zf(:)));

        % Update counter if the new table is as or more extreme
        if gpv <= op
            MCC = MCC + 1;
        end
    end

    % Monte Carlo estimate of the p-value
    P = MCC / tbs;

    % Display Monte Carlo summary
    tr = repmat('-', 1, 80);
    disp(tr)
    disp('Fisher''s test - Conventional Monte Carlo Method')
    disp(tr)
    disp(array2table([tbs, P], ...
        'VariableNames', {'Tables','p_value'}))
    fprintf('p-value is within %0.4f units of the true one with %0.4f%% confidence\n', ...
            delta, (1 - alpha) * 100);
    disp(tr)
end

%% Step 4: Return output (if requested)
if nargout
    Pout = P;
end

end
