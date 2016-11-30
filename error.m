1;

function [y, dy] = exactResult(f, T, x)
  y = (f / (2*T)) .* x .* (1 - x);
  dy = (f / (2*T)) * (1 - 2 .*x);
endfunction