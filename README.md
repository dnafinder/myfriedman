# myfriedman
Friedman test for non parametric two way ANalysis Of VAriance<br/>
The Friedman test is a non-parametric statistical test developed by the U.S.
economist Milton Friedman. Similar to the parametric repeated measures ANOVA,
it is used to detect differences in treatments across multiple test attempts.
The procedure involves ranking each row (or block) together, then considering
the values of ranks by columns. Applicable to complete block designs, it is
thus a special case of the Durbin test. The Friedman test is used for two-way
repeated measures analysis of variance by ranks. In its use of ranks it is
similar to the Kruskal-Wallis one-way analysis of variance by ranks. 
When the number of blocks or treatments is large the probability
distribution can be approximated by chi-square or F distribution. If n or
k is small, the  approximation to chi-square becomes  poor and the
p-value should be obtained from tables of Q specially prepared 
for the Friedman test. The MatLab function FRIEDMAN only uses the chi-square
approximation. On the contrary, MYFRIEDMAN uses the exact distribution for
small size samples and chi-square and F distribution for large sample
size. If the p-value is significant, a post-hoc multiple comparisons.

           Created by Giuseppe Cardillo
           giuseppe.cardillo-edta@poste.it

 To cite this file, this would be an appropriate format:
 Cardillo G. (2009). MYFRIEDMAN: Friedman test for non parametric two way ANalysis Of VAriance
 http://www.mathworks.com/matlabcentral/fileexchange/25882
