% Author: Shunsuke Ono (ono@isl.titech.ac.jp)
% Last version: May. 1, 2017

function[u] = ProjFroball(u, f, epsilon)

temp = u-f;
radius = norm(temp(:),'fro');
if radius > epsilon
    u = f + (epsilon/radius)*temp;
end