function adapt
  
  clc
  close all
  
  for i = 1:100
    clc
    refine(2*i) = 0;
    refine(2*i + 1) = 0; 
    [L2, H1, ElemL2, ElemH1] = main(2, 2, 0.0001, refine);
%    H1
%    ElemH1
    [dummy, refine(2*i)] = max(ElemH1);
    ElemH1(refine(2*i)) = 0;
    [dummy, refine(2*i+1)] = max(ElemH1);
  endfor
  
endfunction