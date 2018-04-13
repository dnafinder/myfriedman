function myfriedman(x,varargin)
% The Friedman test is a non-parametric statistical test developed by the U.S.
% economist Milton Friedman. Similar to the parametric repeated measures ANOVA,
% it is used to detect differences in treatments across multiple test attempts.
% The procedure involves ranking each row (or block) together, then considering
% the values of ranks by columns. Applicable to complete block designs, it is
% thus a special case of the Durbin test. The Friedman test is used for two-way
% repeated measures analysis of variance by ranks. In its use of ranks it is
% similar to the Kruskal-Wallis one-way analysis of variance by ranks. 
% When the number of blocks or treatments is large the probability
% distribution can be approximated by chi-square or F distribution. If n or
% k is small, the  approximation to chi-square becomes  poor and the
% p-value should be obtained from tables of Q specially prepared 
% for the Friedman test. The MatLab function FRIEDMAN only uses the chi-square
% approximation. On the contrary, MYFRIEDMAN uses the exact distribution for
% small size samples and chi-square and F distribution for large sample
% size. If the p-value is significant, a post-hoc multiple comparisons
% tests is performed. 
%
% Syntax: 	myfriedman(X,ALPHA,REPS)
%      
%     Inputs:
%           X - data matrix
%           ALPHA - Significance levele (default=0.05)
%           REPS - If there is more than one observation per row-column pair,
%           then the argument REPS indicates the number of observations per
%           "cell". A cell contains REPS number of rows (default=1).
%     Outputs:
%           - Used Statistic
%           - Multiple comparisons (eventually)
%
%      Example: 
%
% x=[115 142 36 91 28; 28 31 7 21 6; 220 311 108 51 117; 82 56 24 46 33; 256 298 124 46 84; 294 322 176 54 86; 98 87 55 84 25];
%
%           Calling on Matlab the function: myfriedman(x)
%
%           Answer is:
%
% FRIEDMAN TEST FOR IDENTICAL TREATMENT EFFECTS:
% TWO-WAY BALANCED, COMPLETE BLOCK DESIGNS
% ------------------------------------------------------------------------------------------
% Exact Friedman distribution for small size samples
%     Blocks    Treatments    Replicates    Sum_of_Squared_Ranks    Fr     alpha     cv  
%     ______    __________    __________    ____________________    ___    _____    _____
% 
%     7         5             1             385                     378    0.05     9.143
% 
% The 5 treatments have not identical effects
%  
% POST-HOC MULTIPLE COMPARISONS
% ------------------------------------------------------------------------------------------
% Critical value: 6.3053
% Absolute difference among mean ranks
%      0     0     0     0     0
%      3     0     0     0     0
%     15    18     0     0     0
%     15    18     0     0     0
%     18    21     3     3     0
% 
% Absolute difference > Critical Value
%    0   0   0   0   0
%    0   0   0   0   0
%    1   1   0   0   0
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% See also durbin, quadetest
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2009). MYFRIEDMAN: Friedman test for non parametric two way ANalysis Of VAriance
% http://www.mathworks.com/matlabcentral/fileexchange/25882

%Input Error handling
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty'}));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p,'reps',1, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','positive','integer'}));
parse(p,x,varargin{:});
assert(size(x,1)>1,'Warning: X must be a matrix, not a vector')
alpha=p.Results.alpha; reps=p.Results.reps;
clear p

[b,k]=size(x);
R=zeros(b,k); ties=zeros(b/reps,1); z=1; %array preallocation

for I=1:reps:b
    % Keep REPS rows and transform them into an array
    S=reshape(x(I:I+reps-1,:),1,k*reps); 
    % Rank the values
    [Sr,ts]=tiedrank(S);
    % Reshape the S array into REPSxk matrix and assign it to the proper R
    % slice.
    R(I:I+reps-1,:)=reshape(Sr,reps,k);
    if ts % check for ties
        ties(z)=1;
    end
    z=z+1;
end

T=sum(R); %The observed sum of ranks for each treatment
A=sum(sum(R.^2)); %sum of squared ranks
Te=b*(k*reps+1)/2; %The expected value of ranks sum under Ho
Tx=sum((T-Te).^2); % The Friedman statistic
clear S Sr ts z Te I

%display results
tr=repmat('-',1,90); %set the divisor
disp('FRIEDMAN TEST FOR IDENTICAL TREATMENT EFFECTS: TWO-WAY BALANCED, COMPLETE BLOCK DESIGNS')
disp(tr)
flag=0;
if k==3 && b<16
    if b<3
        disp('You must increase the number of blocks')
        disp(tr)
        return
    else
        load('myfriedmantables.mat','friedman')
        critvalstab=friedman.A; clear friedman
        exactdist
    end
elseif k==4 && b<16
        load('myfriedmantables.mat','friedman')
        critvalstab=friedman.B; clear friedman
        exactdist
elseif k==5 && b<11
        load('myfriedmantables.mat','friedman')
        critvalstab=friedman.C; clear friedman
        exactdist
elseif k==6 && b<11
        load('myfriedmantables.mat','friedman')
        critvalstab=friedman.D; clear friedman
        exactdist
else
    N=b*k/reps;
    %T1 is the chi square approximation...
    if any(ties) %...with ties
        C=N*reps^2*(k*reps+1)^2/4;
        T1=(k-1)*Tx/(A-C);
    else %...without ties
        C=false;
        T1=12*Tx/(N*reps^2*(k*reps+1));
    end
    disp(cell2table({b,k,reps,A,any(ties),C},'VariableNames',{'Blocks','Treatments','Replicates','Sum_of_Squared_Ranks','Ties','Correction_factor'}))
    clear C
    df=k-1; %chi-square degrees of freedom
    P1=1-chi2cdf(T1,df);  %probability associated to the Chi-squared-statistic.
    db=b-1;
    T2=db*T1/(b*df-T1); %Transform chi-square into F
    dfd=df*db; %denominator degrees of freedom
    P2=1-fcdf(T2,df,dfd);  %probability associated to the F-statistic.
    fprintf('Chi-square approximation (the most conservative)\n')
    disp(tr)
    disp(table(T1,df,P1,'VariableNames',{'Chi_square','df','two_tailed_p_value'}))
    fprintf('F-statistic approximation (the less conservative)\n')
    disp(tr)
    disp(table(T2,df,dfd,P2,'VariableNames',{'F','df_num','df_denom','two_tailed_p_value'}))
    if P2>alpha
        fprintf('The %i treatments have identical effects\n',k)
    else
        fprintf('The %i treatments have not identical effects\n',k)
        flag=1;
    end
    clear T1 T2 P1 P2 N df db
end

if flag
    clear ties reps R Tx
%when the test is significant myfriedman computes multiple comparisons
%between the individual samples. 
    disp(' ')
    disp('POST-HOC MULTIPLE COMPARISONS')
    disp(tr)
    tmp=repmat(T,k,1); Rdiff=tril(abs(tmp-tmp'),-1); %Generate a matrix with the absolute differences among ranks
    clear tmp tr
    if all(Rdiff==fix(Rdiff))
        %Bioinformatics 2017 18:68
        mask=tril(true(size(Rdiff)),-1);
        d=unique(Rdiff(mask)'); clear mask
        pvalue=zeros(size(Rdiff));
        clear T tmp ties R reps 

        a=repmat(k,1,k+1);
        h=0:1:k;
        B=exp(gammaln(a+1)-gammaln(h+1)-gammaln(a-h+1)).*((1/(1-b)^k)./(b.^h));
        clear a h 

        for I=1:length(d)
            if d(I)==0
                p=1;
            else
                A=0; 
                for h=0:1:k
                    ss=ceil((d(I)+h)/b);
                    if h>=ss
                        E=2*h+1; G=d(I)+h;
                        s=ss:1:h; D=b.*s-G; F=h+s;
                        A=A+B(h+1).*sum((-1).^s.*exp(-gammaln(F+1)-gammaln(E-F)+gammaln(D+E)-gammaln(D+1)));
                    end
                end
                p=2*A;
            end
            pvalue(tril(Rdiff==d(I),-1))=p;
            clear A D E F G ss s p
        end

        clear I h p d B b k
        disp('Absolute difference among mean ranks')
        disp(Rdiff)
        disp('p-values')
        disp(pvalue)
        disp('p-values < alpha')
        disp(tril(pvalue<alpha,-1))
    else
        %These comparisons are performed for all possible contrasts. A
        %contrast is considered significant if the following inequality is
        %satisfied:
        %|Rj-Ri|>tinv(1-alpha/2,dfd)*realsqrt(2*(b*A-sum(T.^2))/((b-1)*(k-1)))
        %where t is a quantile from the Student t distribution on (b-1)(k-1)
        %degrees of freedom.
        %This method is a nonparametric equivalent to Fisher's least significant
        %difference method as described in: Pratical Nonparametric Statistics by W.J. Conover
        cv=tinv(1-alpha/2,dfd)*realsqrt(2*(b*A-sum(T.^2))/dfd); %critical value
        clear dfd b k A T alpha
        mc=Rdiff>cv; %Find differences greater than critical value
        %display results
        fprintf('Critical value: %0.4f\n',cv)
        disp('Absolute difference among mean ranks')
        disp(Rdiff)
        disp('Absolute difference > Critical Value')
        disp(tril(mc))
    end
end

function exactdist
    %critical values from: The Canadian Journal of Statistics - 1993 - 21(1):39-43
    idx=find(critvalstab(:,1)==b,1,'first');
    critvalsrow=critvalstab(idx,2:end); clear critvalstab idx
    alphacol=[0.1 0.05 0.025 0.01 0.005 0.001];
    alphacolcell={'0.100','0.050','0.025','0.010','0.005','0.001'};
    [~,idx]=ismember(alpha,alphacol);
    if idx==0
        idx=listdlg('PromptString','Please, select an alpha value:','ListSize',[300 150],...
            'Name','Disposable values', 'SelectionMode','single',...
            'ListString',alphacolcell);
        alpha=alphacol(idx);
    end
    cv=critvalsrow(idx);
    if isnan(cv)
        idx=listdlg('PromptString','Please, select again an alpha value:','ListSize',[300 150],...
            'Name','Disposable values', 'SelectionMode','single',...
            'ListString',alphacolcell(1:length(critvalsrow(~isnan(critvalsrow)))));
        cv=critvalsrow(idx);
        alpha=alphacol(idx);
    end
    clear alphacol* idx critvalsrow
    disp('Exact Friedman distribution for small size samples')
    disp(array2table([b,k,reps,A,Tx,alpha,cv],'VariableNames',{'Blocks','Treatments','Replicates','Sum_of_Squared_Ranks','Fr','alpha','cv'}))
    if Tx>cv
        flag=1;
        fprintf('The %i treatments have not identical effects\n',k)
        dfd=(k-1)*(b-1);
    else
        fprintf('The %i treatments have identical effects\n',k)
    end
    clear cv
end
end

