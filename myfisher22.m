function Pvalue = myfisher22(x, varargin)
%MYFISHER22 Fisher's Exact Test on a 2x2 contingency table + power analysis.
%
%   P = MYFISHER22(X, ALPHA, PLTS) performs Fisher's exact test on a 2x2
%   contingency table X, computes one-tailed and two-tailed p-values, and
%   optionally computes power and sample size using an asymptotic normal
%   approximation with continuity correction.
%
%   The procedure:
%     - Enumerates all 2x2 tables with the same margins as X.
%     - Computes the exact probability of each table using factorials
%       in log-space (via GAMMALN).
%     - Computes left-tail, right-tail, and two-tailed p-values.
%     - Computes Mid-p correction.
%     - Computes power and recommended sample sizes using the method
%       of Sahai & Khurshid (Statistics in Medicine, 1996).
%     - Optionally plots the distribution of the Wald statistic.
%
%   Syntax:
%       myfisher22(X)
%       myfisher22(X, ALPHA)
%       myfisher22(X, ALPHA, PLTS)
%
%   Inputs:
%       X     - 2x2 data matrix of nonnegative integers:
%               [ A  B
%                 C  D ]
%
%       ALPHA - Significance level (default = 0.05).
%
%       PLTS  - Flag to control plotting of the Wald statistic distribution:
%               0 = no plot (default)
%               1 = show plot
%
%   Outputs:
%       P     - 1x3 vector (if requested) with:
%                 P(1) = left-tail p-value
%                 P(2) = right-tail p-value
%                 P(3) = two-tailed p-value
%               The function also prints:
%                 - number of possible tables,
%                 - Mid-p correction,
%                 - power (one- and two-tailed),
%                 - recommended sample sizes for power ≈ 0.80.
%
%   Example:
%
%                                    Vaccine
%                               Yes           No
%                            ---------------------
%                    Yes         7            12
% Infectious status
%                     No         8             3
%                            ---------------------
%
%   Calling:
%       myfisher22([7 12; 8 3])
%
%   Example output (abridged):
%
%   ------------------------------------------------------------------------
%   2x2 matrix Fisher's exact test
%   ------------------------------------------------------------------------
%        Tables    left_tail    right_tail    two_tails    Mid_p_correction
%        ______    _________    __________    _________    ________________
%
%        11        0.032884     0.99635       0.050175     0.035557
%
%   Power Computation (Asymptotic normal method)
%   ------------------------------------------------------------------------
%   alpha = 0.0500  Sum of row 1 (n1) = 19  Sum of row 2 (n2) = 11
%   ------------------------------------------------------------------------
%                    one_tail    two_tails
%                    ________    _________
%
%          Z1_b      0.65726     0.29123
%          Power     0.74449     0.61456
%
%   To achieve a recommended Power=0.80
%              one_tail    two_tails
%              ________    _________
%
%       n1     15          29
%       n2     18          34
%
%   ------------------------------------------------------------------------
%
%   Created by Giuseppe Cardillo
%   giuseppe.cardillo.75@gmail.com
%
%   To cite this file, a suitable format is:
%   Cardillo G. (2007) MyFisher22: a very compact routine for Fisher's exact
%   test on 2x2 matrix.
%   GitHub: https://github.com/dnafinder/myfisher
%
%   ----------------------------------------------------------------------
%   Algorithm overview for the 2x2 case:
%   ----------------------------------------------------------------------
%   Let the 2x2 table be:
%
%           |  A   B  |  R1
%       X = |  C   D  |  R2
%           ----------
%            C1  C2     N
%
%   with row sums R1, R2; column sums C1, C2; total N.
%
%   Fisher's exact test conditions on the margins (R1, R2, C1, C2). The only
%   free cell is A, and all other cells are determined:
%
%       A      = a
%       C      = R2 - C1 + a
%       B      = R1 - a
%       D      = C1 - a
%
%   The feasible values of a are:
%
%       0 <= a <= min(R1, C1)
%
%   For each feasible a, we obtain a valid 2x2 table with the same margins.
%
%   The probability of a given table is:
%
%       P = [R1! R2! C1! C2!] / [ N! A! B! C! D! ]
%
%   We compute:
%       - The first probability explicitly in log-space using GAMMALN.
%       - All others via a recursive ratio:
%
%         P(a+1) = P(a) * [ B(a)*C(a) ] / [ A(a+1)*D(a+1) ]
%
%     This is vectorized using logs and CUMSUM to avoid loops.
%
%   Tail p-values:
%       - left-tail:  sum of P over tables with A <= A_obs
%       - right-tail: sum of P over tables with A >= A_obs
%       - two-tails:  sum of P over tables with P <= P(A_obs)
%
%   Mid-p correction:
%       MidP = 0.5*P(A_obs) + sum_{P < P(A_obs)} P
%
%   Power & sample size:
%       Implemented using the asymptotic normal approximation with
%       continuity correction as in Sahai & Khurshid (1996).
%   ----------------------------------------------------------------------


%% Step 1: Input error handling & parsing
p = inputParser;
addRequired(p, 'x', @(x) validateattributes(x, {'numeric'}, ...
    {'real','finite','integer','nonnegative','nonnan','size',[2 2]}));
addOptional(p, 'alpha', 0.05, @(x) validateattributes(x, {'numeric'}, ...
    {'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p, 'plts', 0, @(x) isnumeric(x) && isreal(x) && isfinite(x) ...
    && isscalar(x) && (x==0 || x==1));
parse(p, x, varargin{:});
alpha = p.Results.alpha;
plts  = p.Results.plts;
clear p

%% Step 2: Basic margins (row sums, column sums, total)
Rs = sum(x, 2);  % row sums: [R1; R2]
Cs = sum(x, 1);  % column sums: [C1 C2]
N  = sum(Rs);    % total count

%% Step 3: Rearrange matrix if necessary
% We sort rows and columns by their margins to have a canonical ordering.
% This simplifies the enumeration of feasible A values.
% We keep track of how many flips we do with the 'flip' flag:
%   - If exactly one flip is done (flip == 1), we switch the direction
%     of left/right tails later to maintain correct interpretation.

flip = 0; % flip flag

% Sort rows by increasing row sums (if needed)
if ~issorted(Rs)
    x  = flipud(x);        % flip rows
    Rs = sort(Rs);
    flip = flip + 1;
end

% Sort columns by increasing column sums (if needed)
if ~issorted(Cs)
    x  = fliplr(x);        % flip columns
    Cs = sort(Cs);
    flip = flip + 1;
end

% At this point:
%   x is the (possibly reordered) table.
%   Rs, Cs are sorted row/column sums.


%% Step 4: Enumerate all feasible 2x2 tables via A

% Layout reminder:
%          ___________
%         |  A  |  B  | Rs(1)
%         |_____|_____|
%         |  C  |  D  | Rs(2)
%         |_____|_____|____
%           Cs(1) Cs(2)  N
%
% A 2x2 table has only 1 degree of freedom: A.
%
% All feasible values for A:
A = 0:1:min(Rs(1), Cs(1));  % A ranges from 0 to min(R1, C1)

% For each A, we build the corresponding (A, C, B, D) in vector form.
%
% From margin constraints:
%   A      = A
%   C      = Rs(2) - Cs(1) + A
%   B      = Rs(1) - A
%   D      = Cs(1) - A
%
% Collect them in rows:
%   z(1,:) = A
%   z(2,:) = C
%   z(3,:) = B
%   z(4,:) = D
z = [A;
     Rs(2) - Cs(1) + A;
     Rs(1) - A;
     Cs(1) - A];

% Number of feasible tables is simply length(A)
numTables = length(A);


%% Step 5: Compute probabilities of all tables

% We compute:
%   P(A(1)) exactly using GAMMALN (log-factorials),
%   then use a recursive ratio to get all other P(A(i)).

np = zeros(1, numTables);  % p-values (probabilities) for each table

% Precompute log(z) for recursive ratios.
% Note: log(0) = -Inf is fine here because of subtraction rules; we do not
% explicitly exponentiate these logs except in stable combinations.
lz = log(z);

% First probability P(A(1)):
%
% For A(1), we know the corresponding C, B, D from z(:,1).
% The general formula:
%   P = [Rs(1)! Rs(2)! Cs(1)! Cs(2)!] / [ N! A! B! C! D! ]
%
% For A(1) = 0, the formula simplifies, but here we keep a direct form
% in terms of log-factorials to avoid overflow.
np(1) = sum(gammaln([Rs(2)+1, Cs(2)+1])) ...
         - gammaln(N+1) - gammaln(z(2,1)+1);
np(1) = exp(np(1));

% Now use the recursive relation:
%
%   P(i+1) = P(i) * [ B(i)*C(i) ] / [ A(i+1)*D(i+1) ]
%
% In log-space:
%   log P(i+1) = log P(i) + [log B(i) + log C(i)] - [log A(i+1) + log D(i+1)]
%
% Define:
%   f(i) = log B(i) + log C(i) - log A(i+1) - log D(i+1)
%
% Then:
%   log P(j) = log P(1) + sum_{i=1}^{j-1} f(i)
%
% which we vectorize via CUMSUM.
f  = sum(lz(3:4, 1:end-1), 1) - sum(lz(1:2, 2:end), 1); % 1-by-(numTables-1)
lf = log(np(1)) + cumsum(f);                            % log P(2:end)
np(2:end) = exp(lf);


%% Step 6: Compute one-tailed and two-tailed p-values

% Observed A in the (possibly flipped) table is x(1,1).
% Its index in A is W = A_obs + 1.
W = x(1) + 1;

% left/right tails depend on whether we flipped an odd number of times
% (flip == 1) vs 0 or 2 times.
if flip ~= 1
    % No single flip: "left" is A <= A_obs, "right" is A >= A_obs.
    left_tail  = sum(np(1:W));
    right_tail = sum(np(W:end));
else
    % One flip: direction is reversed.
    left_tail  = sum(np(W:end));
    right_tail = sum(np(1:W));
end

% Two-tailed p-value: all tables with probability <= that of the observed.
two_tails = sum(np(np <= np(W)));

% Mid-p correction:
MidP = 0.5*np(W) + sum(np(np < np(W)));

P  = [left_tail, right_tail, two_tails];


%% Step 7: Display main Fisher test results

tr = repmat('-', 1, 80);
disp(tr)
disp('2x2 matrix Fisher''s exact test')
disp(tr)
disp(array2table([numTables, P, MidP], ...
    'VariableNames', {'Tables','left_tail','right_tail','two_tails','Mid_p_correction'}));


%% Step 8: Power computation (Asymptotic normal method)

% This part implements the asymptotic normal method with continuity
% correction proposed by Sahai & Khurshid (1996).
%
% Za: critical values for one-tailed and two-tailed tests
%   Za(1) → one-tailed (alpha)
%   Za(2) → two-tailed (alpha/2)
Za = -realsqrt(2) .* erfcinv(2 .* [1-alpha, 1-alpha/2]);

% Zb corresponds to Beta = 0.2 (Power = 0.8), since:
%   1 - Beta = 0.8  →  Beta = 0.2
Zb = -realsqrt(2) * erfcinv(1.6);  % 1 - Beta = 0.8

% Proportions in column 1:
%   p(1) = A / R1
%   p(2) = C / R2
p = x(:,1) ./ Rs;
d = abs(diff(p));         % absolute difference in proportions
k = Rs(2) / Rs(1);        % ratio of row sizes (n2/n1)
q = 1 - p;

% Pooled proportion under H0:
%   pm = (p1 + k*p2) / (1 + k)
pm = (p(1) + k*p(2)) / (k + 1);
qm = 1 - pm;

% Z1_b is a vector with one-tailed and two-tailed components:
% The formula is derived from the asymptotic distribution of the
% Wald statistic under H1 and includes continuity correction.
Z1_b = abs( ( realsqrt(Rs(1)*d^2) - Za .* realsqrt((1 + 1/k)*pm*qm) ) ...
            ./ realsqrt(p(1)*q(1) + p(2)*q(2)/k) );

% Power for one-tail and two-tails:
pwr = 0.5 .* erfc(-Z1_b ./ realsqrt(2));

fprintf('Power Computation (Asymptotic normal method)\n')
disp(tr)
fprintf('alpha = %0.4f  Sum of row 1 (n1) = %d  Sum of row 2 (n2) = %d\n', ...
        alpha, Rs);
disp(tr)
disp(array2table([Z1_b; pwr], ...
    'VariableNames', {'one_tail','two_tails'}, ...
    'RowNames',      {'Z1_b','Power'}))

%% Step 9: Recommended sample sizes for power ≈ 0.80 (if needed)

if any(pwr < 0.8)
    % Sample size using modified asymptotic normal method with
    % continuity correction (Sahai & Khurshid).
    nstar = ( Za .* realsqrt(pm*qm*(1 + 1/k)) ...
              + Zb .* realsqrt(p(1)*q(1) + p(2)*q(2)/k) ).^2 ./ d^2;

    % Recommended sample sizes for row 1 and row 2:
    n1 = round(nstar./4 .* (1 + realsqrt(1 + 2*(k+1)./(k*d.*nstar))).^2);
    n2 = round(k .* n1);

    disp(' ')
    disp('To achieve a recommended Power=0.80')
    disp(array2table([n1' n2'], ...
        'VariableNames', {'one_tail','two_tails'}, ...
        'RowNames',      {'n1','n2'}))
end
disp(tr)


%% Step 10: Optional plot of Wald statistic distribution

if plts
    % D is the overall proportion in column 1 across rows:
    %   D(i) = (A(i) + C(i)) / N
    % built from z(1,:) and z(3,:)
    D = sum(z([1 3],:), 1) ./ N;

    % Wald statistic T(X) for each table:
    %   TX = (A/C1 - C/C2) / sqrt( D*(1-D) * (1/C1 + 1/C2) )
    TX = ( z(1,:)./Cs(1) - z(3,:)./Cs(2) ) ...
         ./ realsqrt( D.*(1-D) .* sum(1./Cs) );

    % We build a custom bar-like plot using FILL, because the built-in
    % BAR function has limitations for this specific overlay/legend logic.
    hold on
    Wh = 2 * abs(TX(2) - TX(1)) / 5;  % bar width based on tick spacing

    % Tables with probability <= observed (one color)
    leTX = TX(np <= np(W));
    lenp = np(np <= np(W));
    MX1  = repmat(leTX, 4, 1);
    MY1  = repmat(lenp, 4, 1);
    MX1([1 2],:) = MX1([1 2],:) - Wh;
    MX1([3 4],:) = MX1([3 4],:) + Wh;
    MY1([1 4],:) = 0;
    H1 = fill(MX1, MY1, 'b');
    H1Group = hggroup;              % group these bars
    set(H1, 'Parent', H1Group)
    set(get(get(H1Group, 'Annotation'), 'LegendInformation'), ...
        'IconDisplayStyle', 'on');  % include in legend

    % Observed table (highlighted)
    eTX  = TX(W);
    enp  = np(W);
    MX2  = repmat(eTX, 4, 1);
    MY2  = repmat(enp, 4, 1);
    MX2([1 2],:) = MX2([1 2],:) - Wh;
    MX2([3 4],:) = MX2([3 4],:) + Wh;
    MY2([1 4],:) = 0;
    fill(MX2, MY2, 'g');

    % Tables with probability > observed (another color, if any)
    if any(np > np(W))
        gTX  = TX(np > np(W));
        gnp  = np(np > np(W));
        MX3  = repmat(gTX, 4, 1);
        MY3  = repmat(gnp, 4, 1);
        MX3([1 2],:) = MX3([1 2],:) - Wh;
        MX3([3 4],:) = MX3([3 4],:) + Wh;
        MY3([1 4],:) = 0;
        H3 = fill(MX3, MY3, 'r');
        H3Group = hggroup;          % group these bars
        set(H3, 'Parent', H3Group)
        set(get(get(H3Group,'Annotation'), 'LegendInformation'), ...
            'IconDisplayStyle','on');  % include in legend

        legend('P<=P_o_b_s_e_r_v_e_d _t_a_b_l_e', ...
               'P_o_b_s_e_r_v_e_d _t_a_b_l_e', ...
               'P>P_o_b_s_e_r_v_e_d _t_a_b_l_e');
    else
        legend('P<=P_o_b_s_e_r_v_e_d _t_a_b_l_e', ...
               'P_o_b_s_e_r_v_e_d _t_a_b_l_e');
    end
    hold off

    axis square
    title('Distribution of Wald Statistic for Fisher''s Exact test', ...
          'FontName','Arial','FontSize',12,'FontWeight','Bold');
    xlabel('T(X)', 'FontName','Arial','FontSize',12,'FontWeight','Bold');

    txt = ['P[T(X)|X\in\Gamma(' num2str(Rs(1)) ')]'];
    ylabel(txt, 'FontSize',12,'FontWeight','Bold');
end


%% Step 11: Return value (if requested)
if nargout
    Pvalue = P;
end

end
