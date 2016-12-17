function y = dirac(x)

%  x
  eps = 1e-1;
  if (abs(x) <= eps)
    y = 1 / (2*eps)
  else
    y = 0;
  endif
  
endfunction