% WRITTEN BY: GUSTAVO DECO
% CODE SHOULD NOT BE DISTRIBUTED AND ANALYSIS SHOULD NOT BE CONDUCTED WITHOUT PRIOR CONSENT FROM  THE AUTHOR
% If you would like to use this software for publication please contact

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

function c=adif(a,b)
 if abs(a-b)>pi
  c=2*pi-abs(a-b);
 else
  c=abs(a-b);
 end
