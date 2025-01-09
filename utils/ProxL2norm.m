% Author: Shunsuke Ono (ono@isl.titech.ac.jp)
% Last version: May. 1, 2017

function[out] = ProxL2norm(x, gamma)

out = x/(2*gamma+1);
