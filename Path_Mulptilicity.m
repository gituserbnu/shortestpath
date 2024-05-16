function [H_matrix,PHD,PHI]=Path_Mulptilicity(A)
% Path_Mulptilicity is a fast algorithm for the number of shortest paths
%    Inputs:     A,          the adjacency matrix of graph G(V, E)
%               
%
%    Outputs:    H_matrix,   the Path Hesitation Matrix
%                PHD,        the Path Hesitation Degree of all nodes
%                PHI,        the Path Hesitation Index
%    Example:
%       [H_matrix,PHD,PHI]=Path_Mulptilicity(A)
%
%    Source: https://wujunpla.net/uploads/Path-Mulptilicity.zip.
%
% Version 1.0   (2023 May)
% Version 1.0.1 (2023 September)
% Copyright (C) 2023 Ye Deng (Beijing Normal University)
% Distributed under GPL 2.0
% http://www.gnu.org/copyleft/gpl.html
% Notes:
% 
[n,m]=size(A);
H_matrix=zeros(n,n);
C_distance=zeros(n,n);
i=1;
while sum(sum(~(H_matrix+eye(n,n))))~=0
A_temp=A^i;
  A_temp=A_temp-diag(diag(A_temp));
  A_temp=(~H_matrix).*A_temp;
if sum(sum(A_temp))==0
break
else
  H_matrix=H_matrix+A_temp;
A_temp(A_temp~=0)=1;
C_distance=C_distance+A_temp*i;
i=i+1;
end
end
count=n*(n-1)/2;
%%%%calculate the value of PHD%%%%
PHD=sum(H_matrix,2)/(n-1);
%%%%calculate the value of PHI%%%%
PHI=sum(sum(triu(H_matrix)))/count;


