
def Path_Mulptilicity(A):
% Path_Mulptilicity is a fast algorithm for the number of shortest paths
% 
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
    n, m = A.shape
    H_matrix = np.zeros((n, n))
    i = 1
    while np.sum(np.where((C_num+np.eye(n)) != 0, 0, 1))  != 0:
        A_temp=np.linalg.matrix_power(A, i)
        A_temp=A_temp-np.diag(np.diag(A_temp))
        A_temp=(np.where(C_num != 0, 0, 1))*A_temp
        H_matrix = H_matrix + A_temp
        A_temp[A_temp != 0] = 1
        i=i+1
    count = n * (n - 1) / 2
    %%%%calculate the value of PHD%%%%
    PHD = np.sum(H_matrix, axis=1) / (n-1)
    %%%%calculate the value of PHI%%%%
    PHI = np.sum(np.triu(H_matrix )) / count
   

    return H_matrix, PHD, PHI




