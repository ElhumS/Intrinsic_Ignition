% WRITTEN BY: GUSTAVO DECO
% CODE SHOULD NOT BE DISTRIBUTED AND ANALYSIS SHOULD NOT BE CONDUCTED WITHOUT PRIOR CONSENT FROM  THE AUTHOR
% If you would like to use this software for publication please contact

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

function x=demean(x,dim)

% DEMEAN(X) 
% Removes the Average or mean value.
%
% DEMEAN(X,DIM)
% Removes the mean along the dimension DIM of X. 

if(nargin==1),
   dim = 1;
   if(size(x,1) > 1)
      dim = 1;
   elseif(size(x,2) > 1)
      dim = 2;
   end;
end;

dims = size(x);
dimsize = size(x,dim);
dimrep = ones(1,length(dims));
dimrep(dim) = dimsize;

x = x - repmat(mean(x,dim),dimrep);
