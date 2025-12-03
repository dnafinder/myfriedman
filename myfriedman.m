function [stats, mc] = myfriedman(x, varargin)
%MYFRIEDMAN Friedman test for nonparametric two-way ANOVA with post-hoc.
%
%   [STATS, MC] = MYFRIEDMAN(X)
%   [STATS, MC] = MYFRIEDMAN(X, 'Alpha', ALPHA, 'Reps', REPS, ...
%                                'PostHoc', POSTHOC, 'Display', DISPLAY, ...
%                                'Exact', EXACT)
%
%   The Friedman test is a non-parametric statistical test for detecting
%   differences between treatments across multiple blocks in a two-way
%   balanced, complete block design. It is the rank-based analogue of
%   repeated-measures ANOVA.
%
%   This implementation extends MATLAB's FRIEDMAN by:
%     - using exact critical values from tables for small designs,
%     - using both chi-square and F approximations for larger designs,
%     - optionally performing post-hoc multiple comparisons when the global
%       test is significant.
%
%   Inputs:
%     X        - data matrix (Blocks x Treatments).
%
%   Name–Value pair arguments:
%     'Alpha'  - significance level (scalar in (0,1), default = 0.05).
%
%     'Reps'   - number of replicate observations per block–treatment cell
%                (positive integer, default = 1). If REPS > 1, rows in X
%                are interpreted as stacked replicates for each block.
%
%     'PostHoc'- logical-like flag to enable post-hoc multiple comparisons:
%                true  -> perform multiple comparisons when global H0 is rejected
%                false -> only perform the global Friedman test
%                (default = true)
%
%     'Display'- logical-like flag controlling command-window output:
%                true  -> print tables and messages (default)
%                false -> run silently, only return outputs
%
%     'Exact'  - logical-like flag to use exact critical values for small
%                designs when available:
%                true  -> use tables from myfriedmantables.mat (default)
%                false -> always use chi-square/F approximations
%
%   Outputs:
%     STATS - structure with test results (see README for details).
%     MC    - structure with post-hoc multiple comparison results (or empty).
%
%   Example:
%
%     x = [115 142  36  91  28;
%           28  31   7  21   6;
%          220 311 108  51 117;
%           82  56  24  46  33;
%          256 298 124  46  84;
%          294 322 176  54  86;
%           98  87  55  84  25];
%
%     [stats, mc] = myfriedman(x);
%
%   Created by Giuseppe Cardillo
%   giuseppe.cardillo.75@gmail.com
%
%   To cite this file, this would be an appropriate format:
%   Cardillo G. (2009). MYFRIEDMAN: Friedman test for non parametric
%   two way ANalysis Of VAriance. Available on GitHub:
%   https://github.com/dnafinder/myfriedman

% -------------------------------------------------------------------------
% Input parsing
% -------------------------------------------------------------------------
p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'x', @(x) validateattributes(x, {'numeric'}, ...
    {'2d','real','finite','nonnan','nonempty'}, mfilename, 'X', 1));
addParameter(p, 'Alpha',  0.05, @(a) validateattributes(a, {'numeric'}, ...
    {'scalar','real','finite','nonnan','>',0,'<',1}, mfilename, 'Alpha'));
addParameter(p, 'Reps',   1,    @(r) validateattributes(r, {'numeric'}, ...
    {'scalar','real','finite','nonnan','positive','integer'}, mfilename, 'Reps'));
addParameter(p, 'PostHoc', true,  @validateLogicalLike);
addParameter(p, 'Display', true,  @validateLogicalLike);
addParameter(p, 'Exact',   true,  @validateLogicalLike);
parse(p, x, varargin{:});

x           = p.Results.x;
alpha       = p.Results.Alpha;
reps        = p.Results.Reps;
postHocFlag = logical(normalizeLogicalLike(p.Results.PostHoc));
displayFlag = logical(normalizeLogicalLike(p.Results.Display));
exactFlag   = logical(normalizeLogicalLike(p.Results.Exact));

assert(size(x,1) > 1, 'X must be a matrix with more than one row.');

[b, k] = size(x);        % b = blocks (rows), k = treatments (columns)

% Prepare output containers
stats = struct();
mc    = struct([]);  % verrà riempito come mc(1).field nel post-hoc

% -------------------------------------------------------------------------
% Rank transformation within blocks (with possible replicates)
% -------------------------------------------------------------------------
R        = zeros(b, k);
tiesFlag = false;  % any ties across blocks?

z = 1;
for I = 1:reps:b
    % Keep REPS rows and transform them into a vector
    S = reshape(x(I:I+reps-1, :), 1, k*reps);
    % Rank the values with tie-handling
    [Sr, ts] = tiedrank(S);
    % Reshape back into REPS x k and assign into R
    R(I:I+reps-1, :) = reshape(Sr, reps, k);
    if ts
        tiesFlag = true;
    end
    z = z + 1;
end

T  = sum(R);           % observed sum of ranks for each treatment
A  = sum(sum(R.^2));   % sum of squared ranks
Te = b * (k*reps + 1) / 2;      % expected rank sum per treatment under H0
Tx = sum((T - Te).^2);          % Friedman statistic (unscaled form)

% Common header
trLine = repmat('-', 1, 90);

if displayFlag
    disp('FRIEDMAN TEST FOR IDENTICAL TREATMENT EFFECTS: TWO-WAY BALANCED, COMPLETE BLOCK DESIGNS');
    disp(trLine);
end

% -------------------------------------------------------------------------
% Decide exact vs approximation
% -------------------------------------------------------------------------
stats.blocks     = b;
stats.treatments = k;
stats.reps       = reps;
stats.sumRanks   = T(:);
stats.sumSqRanks = A;
stats.Fr         = Tx;
stats.ties       = tiesFlag;
stats.tieCorrection = [];
stats.alpha      = alpha;
stats.method     = '';
stats.chi2       = struct('value', NaN, 'df', NaN, 'p', NaN);
stats.F          = struct('value', NaN, 'df_num', NaN, 'df_denom', NaN, 'p', NaN);
stats.criticalValue = NaN;
stats.rejectNull    = false;

flag = false;  % whether global null is rejected (for post-hoc)

% Small-sample exact region conditions (as original)
useExact = exactFlag && ( ...
    (k == 3 && b < 16) || ...
    (k == 4 && b < 16) || ...
    (k == 5 && b < 11) || ...
    (k == 6 && b < 11) );

% Special case: very few blocks for k=3
if k == 3 && b < 3
    if displayFlag
        disp('You must increase the number of blocks');
        disp(trLine);
    end
    stats.method = 'Insufficient blocks for k=3';
    if nargout == 0
        clear stats mc
    end
    return;
end

if useExact
    % Exact Friedman distribution via critical values
    if displayFlag
        disp('Exact Friedman distribution for small size samples');
    end
    
    % Load critical values
    S   = load('myfriedmantables.mat', 'friedman');
    tab = S.friedman;
    
    if k == 3
        critvalstab = tab.A;
    elseif k == 4
        critvalstab = tab.B;
    elseif k == 5
        critvalstab = tab.C;
    elseif k == 6
        critvalstab = tab.D;
    else
        critvalstab = [];
    end
    
    [cv, alphaUsed] = friedmanExactCritical(critvalstab, b, alpha);
    stats.alpha         = alphaUsed;
    stats.criticalValue = cv;
    stats.method        = 'Exact Friedman';
    
    if displayFlag
        Texact = array2table([b, k, reps, A, Tx, alphaUsed, cv], ...
            'VariableNames', {'Blocks','Treatments','Replicates', ...
                              'Sum_of_Squared_Ranks','Fr','alpha','cv'});
        disp(Texact);
    end
    
    if Tx > cv
        stats.rejectNull = true;
        flag = true;
        if displayFlag
            fprintf('The %i treatments have not identical effects\n', k);
        end
    else
        if displayFlag
            fprintf('The %i treatments have identical effects\n', k);
        end
    end
    
    % Degrees of freedom for post-hoc (used later)
    dfd = (k - 1) * (b - 1);
    
else
    % Large-sample approximations: chi-square and F
    N = b * k / reps;
    
    if tiesFlag
        % Chi-square approximation with ties
        C = N * reps^2 * (k*reps + 1)^2 / 4;
        T1 = (k - 1) * Tx / (A - C);
        stats.tieCorrection = C;
    else
        C  = [];
        T1 = 12 * Tx / (N * reps^2 * (k*reps + 1));
        stats.tieCorrection = [];
    end
    
    df  = k - 1;           % chi-square degrees of freedom
    P1  = 1 - chi2cdf(T1, df);
    
    db  = b - 1;
    T2  = db * T1 / (b*df - T1);   % transform chi-square into F
    dfd = df * db;                 % denominator degrees of freedom
    P2  = 1 - fcdf(T2, df, dfd);
    
    stats.method      = 'Chi-square/F approximation';
    stats.chi2.value  = T1;
    stats.chi2.df     = df;
    stats.chi2.p      = P1;
    stats.F.value     = T2;
    stats.F.df_num    = df;
    stats.F.df_denom  = dfd;
    stats.F.p         = P2;
    
    if displayFlag
        disp(cell2table({b, k, reps, A, tiesFlag, C}, ...
            'VariableNames', {'Blocks','Treatments','Replicates', ...
                              'Sum_of_Squared_Ranks','Ties','Correction_factor'}));
        disp(trLine);
        fprintf('Chi-square approximation (the most conservative)\n');
        disp(trLine);
        disp(table(T1, df, P1, ...
            'VariableNames', {'Chi_square','df','two_tailed_p_value'}));
        fprintf('F-statistic approximation (the less conservative)\n');
        disp(trLine);
        disp(table(T2, df, dfd, P2, ...
            'VariableNames', {'F','df_num','df_denom','two_tailed_p_value'}));
    end
    
    if P2 > alpha
        if displayFlag
            fprintf('The %i treatments have identical effects\n', k);
        end
        stats.rejectNull = false;
    else
        if displayFlag
            fprintf('The %i treatments have not identical effects\n', k);
        end
        stats.rejectNull = true;
        flag = true;
    end
end

% -------------------------------------------------------------------------
% Post-hoc multiple comparisons (if requested and global H0 rejected)
% -------------------------------------------------------------------------
mc = struct([]); % verrà riempito solo se flag && postHocFlag

if flag && postHocFlag
    if displayFlag
        disp(' ');
        disp('POST-HOC MULTIPLE COMPARISONS');
        disp(trLine);
    end
    
    % Generate matrix of absolute differences among rank sums
    tmp   = repmat(T, k, 1);
    Rdiff = tril(abs(tmp - tmp'), -1);
    
    % Decide which method to use for post-hoc
    if all(Rdiff(:) == fix(Rdiff(:)))
        % Exact pairwise comparisons (Bioinformatics 2017 18:68)
        mask   = tril(true(size(Rdiff)), -1);
        d      = unique(Rdiff(mask)');
        pvalue = zeros(size(Rdiff));
        clear mask
        
        % The following follows the original algorithm:
        % a, h, B are intermediate combinatorial weights
        a = repmat(k, 1, k+1);
        h = 0:1:k;
        B = exp(gammaln(a + 1) - gammaln(h + 1) - gammaln(a - h + 1)) ...
            .* ((1 / (1 - b)^k) ./ (b.^h));
        clear a h
        
        for I = 1:length(d)
            if d(I) == 0
                p = 1;
            else
                Aacc = 0;
                for h = 0:1:k
                    ss = ceil((d(I) + h) / b);
                    if h >= ss
                        E = 2*h + 1;
                        G = d(I) + h;
                        s = ss:1:h;
                        D = b.*s - G;
                        Fv = h + s;
                        Aacc = Aacc + B(h+1) .* sum( ...
                            (-1).^s .* exp(-gammaln(Fv + 1) - gammaln(E - Fv) + ...
                            gammaln(D + E) - gammaln(D + 1)));
                    end
                end
                p = 2 * Aacc;
            end
            pvalue(tril(Rdiff == d(I), -1)) = p;
            clear Aacc D E Fv G ss s p
        end
        
        if displayFlag
            disp('Absolute difference among mean ranks');
            disp(Rdiff);
            disp('p-values');
            disp(pvalue);
            disp('p-values < alpha');
            disp(tril(pvalue < alpha, -1));
        end
        
        mc(1).method      = 'Exact pairwise (Bioinformatics 2017)';
        mc(1).Rdiff       = Rdiff;
        mc(1).pvalue      = pvalue;
        mc(1).cv          = [];
        mc(1).significant = tril(pvalue < alpha, -1);
    else
        % Conover-type nonparametric LSD
        cv = tinv(1 - alpha/2, dfd) * sqrt(2 * (b*A - sum(T.^2)) / dfd);
        mcMatrix = Rdiff > cv;
        
        if displayFlag
            fprintf('Critical value: %0.4f\n', cv);
            disp('Absolute difference among mean ranks');
            disp(Rdiff);
            disp('Absolute difference > Critical Value');
            disp(tril(mcMatrix));
        end
        
        mc(1).method      = 'Conover-type LSD';
        mc(1).Rdiff       = Rdiff;
        mc(1).pvalue      = [];
        mc(1).cv          = cv;
        mc(1).significant = tril(mcMatrix, -1);
    end
end

% -------------------------------------------------------------------------
% Output behaviour
% -------------------------------------------------------------------------
if nargout == 0
    clear stats mc
end

end

% -------------------------------------------------------------------------
% Local helper: validate logical-like values
% -------------------------------------------------------------------------
function tf = validateLogicalLike(x)
    try
        normalizeLogicalLike(x);
        tf = true;
    catch
        tf = false;
    end
end

function y = normalizeLogicalLike(x)
    if islogical(x)
        y = x;
    elseif isnumeric(x) && isscalar(x)
        y = (x ~= 0);
    elseif ischar(x) || (isstring(x) && isscalar(x))
        s = lower(char(x));
        if any(strcmp(s, {'true','on','yes'}))
            y = true;
        elseif any(strcmp(s, {'false','off','no'}))
            y = false;
        else
            error('Invalid logical-like value: %s', s);
        end
    else
        error('Invalid type for logical-like option.');
    end
end

% -------------------------------------------------------------------------
% Local helper: exact Friedman critical value
% -------------------------------------------------------------------------
function [cv, alphaUsed] = friedmanExactCritical(critvalstab, b, alpha)
%FRIEDMANEXACTCRITICAL Selects the exact critical value for given b and alpha.
    alphacol = [0.1 0.05 0.025 0.01 0.005 0.001];
    
    % Find row matching number of blocks
    idxRow = find(critvalstab(:,1) == b, 1, 'first');
    if isempty(idxRow)
        error('No exact critical values available for b = %d.', b);
    end
    critvalsrow = critvalstab(idxRow, 2:end);
    
    % Require alpha to be one of the tabulated values
    [tf, idxCol] = ismember(alpha, alphacol);
    if ~tf
        error('For exact Friedman critical values, alpha must be one of: %s', ...
            num2str(alphacol));
    end
    
    cv = critvalsrow(idxCol);
    if isnan(cv)
        error('No exact critical value available for b = %d and alpha = %g.', b, alpha);
    end
    
    alphaUsed = alpha;
end
