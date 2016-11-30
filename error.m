1;

function y = exactResult(f, T, x)
  y = (f / (2*T)) .* x .* (1 - x);
endfunction