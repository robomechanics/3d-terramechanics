function A = ellipse_int(r, beta, h)
%calculates the area of half of a partial ellipse
%        _____
%     /####|    \           |
%   /###A##| h    \         |
%  /#######|       \        | r
% |                 |       |
% |                 |       |
% |                 |
% |                 | 
%  \               /
%   \             /
%     \  _____  /
%       
% ---2*r*sin(beta)---

syms z
A = r*sin(beta)*int((1-((z-h+r)/r)^2)^.5, z, [0 h]);
A = eval(A);
end
