function adapt
  
  clc
  close all
  
  numRefine = 1;
  refine(1) = 0;
  
  for i = 1:10
%    clc
    [L2, H1, ElemL2, ElemH1] = main(2, 2, 0.01, refine);
%    H1
%    ElemH1
    avg = mean(ElemH1);
    for j = 1:numel(ElemH1)
      if (ElemH1(j) > avg)
        refine(numRefine) = j;
        numRefine = numRefine + 1;
      endif
    endfor
    
    info(1) = numel(ElemH1);
    info
    
  endfor
  
endfunction