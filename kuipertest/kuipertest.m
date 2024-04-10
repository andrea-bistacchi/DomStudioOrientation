function [H, pValue, KuiperStat, criticalValue] = kuipertest(x, varargin)
%KUIPERTEST Single Sample Kuiper Goodness-Of-Fit Hypothesis Test.  
%   H = KUIPERTEST(X) performs a Kuiper test to determine if the random
%   sample X could have the hypothesized continuous cumulative distribution
%   function F. 
%   Null Hypothesis Ho: "The Sample is taken form a population with 
%   cumulative distribution F(X) for all X". This is a 2-sided test.
%   H indicates the result of the hypothesis test:
%      H = 0 => Do not reject the null hypothesis at the 5% significance
%      level. 
%      H = 1 => Reject the null hypothesis at the 5% significance
%      level.
%
%   The Kuiper test is similar to Kolmogorov-Smirnov (K-S) test, but K-S
%   test tend to be most sensitive around median value of the distribution
%   and less sensitive at the distribution tails. 
%
%   X is a vector representing a random sample from some underlying
%   distribution, with cumulative distribution function F. Missing 
%   observations in X, indicated by NaNs (Not-a-Number), are ignored.
%
%   [H,P] = KUIPERTEST(...) also returns the asymptotic P-value P.
%
%   [H,P,KSTAT] = KUIPERTEST(...) also returns the Kuiper test statistic
%   KSTAT defined above.
%
%   [H,P,KSSTAT,CV] = KUIPERTEST(...) returns the critical value of the
%   test CV.
%
%   [...] = KUIPERTEST(X,'PARAM1',val1,'PARAM2',val2,...) specifies one or
%   more of the following parameter name/value pairs to specify whether to
%   perform a simple or composite hypothesis test, to specify the
%   distribution being tested for, to control for the significance level,
%   and to specify whether to perform the test using Monte-Carlo
%   simulations:
%
%   Parameter       Value
%
%   'alpha'         A value ALPHA between 0 and 1 specifying the
%                   significance level. Default is 0.05 for 5% significance.
%
%   'CDF'           CDF is the c.d.f. under the null hypothesis.  It can
%                   be specified either as a ProbabilityDistribution object
%                   or as a two-column matrix. CDF must be completely
%                   specified. If missing, the default is the standard
%                   normal with sample mean and variance.
%  
%   Test statistics is KSTAT = max(S(x) - F(x)) + max(F(x) - S(x))
%   where S(x) is the empirical cumulative distribution function
%
%   Calculation are based on the asymtotic approximation of pValue.
%
%   In the matrix version of CDF, column 1 contains the x-axis data and
%   column 2 the corresponding y-axis c.d.f data. Since the Kuiper test
%   statistic will occur at one of the observations in X, the calculation
%   is most efficient when CDF is only specified at the observations in X.
%   When column 1 of CDF represents x-axis points independent of X, CDF is
%   're-sampled' at the observations found in the vector X via
%   interpolation. In this case, the interval along the x-axis (the column
%   1 spread of CDF) must span the observations in X for successful
%   interpolation.
%
%   See also KSTEST, CMTEST, UTEST, ADTEST.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Copyright (c) 20 March 2015 by Ahmed Ben Saïda           %
%                 LaREMFiQ Laboratory, IHEC Sousse - Tunisia             %
%                       Email: ahmedbensaida@yahoo.com                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%   Reference:
%   Press, W.H., et. al., "Numerical Recipes in C++", 3rd Ed.
%      Cambridge University Press, 2007.
%   Stephens (1965), "The goodness-of-fit statistic Vn. distribution and
%       significance points", Biometrika, Vol. 52, No. 3-4, pp.309-321.
%

% Parse arguments and check if parameter/value pairs are valid
paramNames = {'Alpha', 'CDF'};
dflts  =     {  0.05,   []  };

[valAlpha, CDF] = parseArgs(paramNames, dflts, varargin{:});

% Ensure the sample data is a real vector.
if ~isvector(x) || ~isreal(x)
    error('Sample data X must be a real vector.');
end

% Ensure the significance level, ALPHA, is a scalar between 0 and 1.
if ~isscalar(valAlpha) || ~(valAlpha > 0 && valAlpha < 1)
    error('Significance level ALPHA must be a scalar between 0 and 1.');
else
    alpha = valAlpha;
end

try
    x = x(~isnan(x));
    % n is the sample size. 
    % Calculate it now because ecdf removes duplicates.
    n = length(x);
    [sampleCDF,x] = ecdf(x);
catch ME
    if isequal(ME.identifier,'stats:ecdf:VectorRequired')
        error('Input X must be a vector.');
    elseif isequal(ME.identifier,'stats:ecdf:NotEnoughData')
        error('Input sample has no valid data (all missing values).');
    else
       error('Error calculating CDF of input ''X''.');
    end
end

if n == 1
    error('Cannot conduct the Kuiper test on a single observation.')
end

x = x(2:end);

% Check & scrub the hypothesized CDF specified under the null hypothesis.
% If CDF has been specified, remove any rows with NaN's in them and sort
% x-axis data found in column 1 of CDF. If CDF has not been specified, then
% allow the convenience of x ~ N(0,1) under the null hypothesis.
if (isa(CDF,'ProbDist') || isa(CDF,'prob.ProbabilityDistribution'))
   xCDF = x;
   yCDF = cdf(CDF,x);
elseif ~isempty(CDF)
    if size(CDF,2) ~= 2
      error('Hypothesized CDF matrix must have 2 columns.');
    end

    CDF  =  CDF(~isnan(sum(CDF,2)),:);

    if size(CDF,1) == 0
      error('Hypothesized CDF matrix must have at least 1 valid row.');
    end

    [xCDF,i] =  sort(CDF(:,1));    % Sort the theoretical CDF.
    yCDF = CDF(i,2);

    ydiff = diff(yCDF);
    if any(ydiff<0)
      error('CDF must define an increasing function of X.');
    end

    % Remove duplicates, but it's an error if they are not consistent
    dups = find(diff(xCDF) == 0);
    if ~isempty(dups)
      if ~all(ydiff(dups) == 0)
         error('CDF must not have duplicate X values.');
      end
      xCDF(dups) = [];
      yCDF(dups) = [];
    end
    
else
    xCDF = x;
    yCDF = normcdf(x,mean(x),std(x));
end

% If CDF's x-axis values have been specified by the data sample X, then just
% assign column 2 to the null CDF; if not, then we interpolate subject to the
% check that the x-axis interval of CDF must bound the observations in X. Note
% that both X and CDF have been sorted and have had NaN's removed.
if isequal(x,xCDF)
   nullCDF  =  yCDF;   % CDF has been specified at the observations in X.
else
   if (x(1) < xCDF(1)) || (x(end) > xCDF(end))
     error('Hypothesized CDF matrix must span the observations interval in X.');
   end
   nullCDF  =  interp1(xCDF, yCDF, x, 'linear');
end

% Compute the test statistic of interest.

% max (F(x) - S(x))
delta1    =  nullCDF - sampleCDF(1:end-1);   % Vertical difference at jumps approaching from the LEFT.
delta2    =  nullCDF - sampleCDF(2:end);     % Vertical difference at jumps approaching from the RIGHT.
deltaCDF  =  [delta1 ; delta2];
Kminus    =  max(deltaCDF);

% max(S(x) - F(x))
delta1    =  sampleCDF(1:end-1) - nullCDF;   % Vertical difference at jumps approaching from the LEFT.
delta2    =  sampleCDF(2:end)   - nullCDF;   % Vertical difference at jumps approaching from the RIGHT.
deltaCDF  =  [delta1 ; delta2];
Kplus     =  max(deltaCDF);  

KuiperStat = Kplus + Kminus;

% Now we compute the Asymtotic p-Value 
lambda  =  max((sqrt(n) + 0.155 + 0.24/sqrt(n)) * KuiperStat, 0); % This max is useless if CDF in input is correct.

if lambda < 0.4 % Useless to compute pValue (pValue equals 1 at the 7th decimal
                % For small value of lambda some problems may also arise with sum convergence
    pValue = 1; % KuiperStat very small for this sample length: never reject H0
    H = 0;
    return
end

% Use the approximation in Press et al. (2007, p. 739)
fun = @(j) 2 * (4 * j.^2 * lambda^2 - 1) .* exp(-2 * j.^2 * lambda^2);
j   = 1:100; % j tends to infinity.
pValue   = sum(fun(j));

pValue(pValue < 0) = 0;
pValue(pValue > 1) = 1;

% "H = 0" implies that we "Do not reject the null hypothesis at the
% significance level of alpha," and "H = 1" implies that we "Reject null
% hypothesis at significance level of alpha."
H = (pValue < alpha);

% Compute the critical value only when asked.
if nargout > 3
   % Calculate the critical value which 'Kuiperstatistic' must exceed for the null
   % hypothesis to be rejected. If the sample size 'n' is greater than 100, use
   % a reverse search; otherwise interpolate into his 'exact' table.

   if n <= 12        % Small sample exact values.
      % For n <= 12, Stephens tabularized the exact critical values by solving
      % the nth order polynomial. These exact values for n <= 12 are hard-coded
      % into the matrix 'exact' shown below. Rows 2:12 correspond to sample
      % sizes n = 2:12.
      a     =  [0.001	0.005	0.010	0.025	0.050	0.100	0.150	0.850	0.900	0.95	0.975	0.990	0.995	0.999]';  % 2-sided significance level

      exact =  [0.9995	0.9975	0.995	0.9875	0.975	0.950	0.925	0.575	0.550	0.525	0.513	0.505	0.503	0.501
                0.982	0.959	0.942	0.909	0.871	0.817	0.776	0.491	0.462	0.425	0.398	0.374	0.362	0.346
                0.937	0.892	0.864	0.816	0.768	0.714	0.683	0.434	0.411	0.378	0.351	0.325	0.309	0.285
                0.881	0.822	0.789	0.740	0.700	0.652	0.619	0.388	0.370	0.343	0.320	0.296	0.280	0.254
                0.824	0.762	0.732	0.687	0.646	0.601	0.571	0.356	0.337	0.314	0.295	0.273	0.260	0.234
                0.775	0.716	0.686	0.641	0.604	0.561	0.532	0.333	0.315	0.290	0.273	0.255	0.243	0.219
                0.734	0.676	0.647	0.605	0.569	0.528	0.501	0.313	0.296	0.274	0.256	0.239	0.228	0.207
                0.699	0.642	0.614	0.574	0.539	0.500	0.475	0.296	0.281	0.259	0.243	0.225	0.215	0.196
                0.668	0.613	0.586	0.547	0.514	0.477	0.452	0.282	0.267	0.247	0.231	0.214	0.204	0.187
                0.641	0.587	0.562	0.524	0.492	0.456	0.432	0.270	0.256	0.237	0.221	0.205	0.195	0.178
                0.617	0.565	0.540	0.503	0.471	0.437	0.415	0.259	0.245	0.227	0.213	0.197	0.188	0.170];

      criticalValue  =  spline(a , exact(n-1,:)' , alpha);
   
   elseif (n > 12) && (n <= 100) 
      a     =  [0.005	0.010	0.025	0.050	0.100	0.150	0.850	0.900	0.95	0.975	0.990	0.995];  % 2-sided significance level
      
      exact =  [0.527	0.503	0.469	0.439	0.408	0.386	0.241	0.228	0.211	0.198	0.184	0.175   %n=14
                0.496	0.473	0.441	0.414	0.384	0.363	0.227	0.215	0.198	0.186	0.173	0.165   %n=16
                0.470	0.448	0.417	0.392	0.363	0.343	0.215	0.203	0.188	0.176	0.163	0.156   %n=18
                0.447	0.427	0.397	0.372	0.346	0.326	0.204	0.193	0.179	0.168	0.156	0.148   %n=20
                0.369	0.352	0.328	0.307	0.285	0.269	0.169	0.160	0.147	0.138	0.128	0.122   %n=30
                0.322	0.307	0.286	0.268	0.248	0.235	0.147	0.139	0.129	0.121	0.112	0.107   %n=40
                0.289	0.276	0.256	0.241	0.223	0.211	0.141	0.125	0.116	0.108	0.101	0.095   %n=50
                0.264	0.252	0.235	0.220	0.204	0.193	0.121	0.115	0.106	0.099	0.093	0.088   %n=60
                0.245	0.234	0.218	0.204	0.189	0.179	0.112	0.107	0.098	0.092	0.086	0.082   %n=70
                0.230	0.219	0.204	0.191	0.178	0.168	0.105	0.100	0.092	0.086	0.081	0.077   %n=80
                0.206	0.197	0.183	0.172	0.159	0.151	0.095	0.090	0.083	0.078	0.073	0.069]; %n=100
      
      sampleSizes = [14 16  18  20  30  40  50  60  70  80  100];
      
      [OneOverSampleSizes, LogAlpha] = meshgrid(1./sampleSizes, log(a));
      criticalValue = interp2(OneOverSampleSizes, LogAlpha, exact', 1/n, log(alpha));

   else                % Large sample approximate values.
      criticalValue = fzero(@(x) kuiperpvalue(x, n) - alpha, KuiperStat);
   end
end
    
%----------------Subfunctions--------------------------------------------%
function varargout = parseArgs(paramNames, dflts, varargin)

varargout = dflts;
for i = 1:length(dflts)
    varargout{i} = dflts{i};
end

n = length(varargin);

% Check the input argument list:
if rem(n,2) ~= 0
    error('Wrong input variables.')
else
    parameters = varargin(1:2:n-1);
    values = varargin(2:2:n);
end

for j = 1:length(parameters)
    % Ensure that all parameter names are strings:
    if ~ischar(parameters{j})
       error('Parameter name at input argument number %s must be a string.', num2str(2*j - 1))
    end
    % Ensure that all parameter names represent valid field names:
    matches = find(strncmpi(parameters{j},paramNames,length(parameters{j})));
    if isempty(matches)
       % No matches found:
       error('Unrecognized parameter name ''%s''.', parameters{j})
    elseif length(matches) > 1
       % More than one match:
       exacts = find(strcmpi(parameters{j},paramNames));
       if length(exacts) == 1
        % An exact match:
          matches = exacts;
       else
          % The parameter name is ambiguous; additional information needed:
          error('Ambiguous parameter name ''%s''. Additional information needed.', parameters{j})
       end
    end
    varargout{matches} = values{j};
end

%------------------------------------------------------------------------%
function pval = kuiperpvalue(KuiperStat, n)

lambda  =  (sqrt(n) + 0.155 + 0.24/sqrt(n)) * KuiperStat;

% Use the approximation in Press et al. (2007, p. 739)
fun = @(j) 2 * (4 * j.^2 * lambda^2 - 1) .* exp(-2 * j.^2 * lambda^2);
j   = 1:100; % j tends to infinity.
pval   = sum(fun(j));