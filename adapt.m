function adapt
  
  clc
  close all
  
  
  for d = 1:1
  
    numRefine = 1;
    refine(1, d + 1) = 0;
    
    eps = 0.0001;
    nElemStart = 2;
    pdeg = d + 1;
  
    for i = 1:15
      [L2, H1, ElemL2, ElemH1] = main(nElemStart, pdeg, eps, refine(:, d + 1));
      avg = mean(ElemH1);
      max_elem = max(ElemH1);
      std_elem = std(ElemH1) / avg;
      multiplier = 1;
      if (std_elem > 0.4)
        multiplier = 1.0;
      endif
      for j = 1:numel(ElemH1)
        if (ElemH1(j) > avg * multiplier)
          refine(numRefine, d + 1) = j;
          numRefine = numRefine + 1;
        endif
      endfor
      
      info(i,d*4+1) = numel(ElemH1) + numel(ElemH1)*d;
      info(i,d*4+2) = H1;
      info(i,d*4+3) = L2;
      info(i,d*4+4) = numel(ElemH1);
      
    endfor
  endfor
  
  info
  
  figure(3)
  hold on
%  loglog((info(:,1)), (info(:,2)), 'k--*')
%  loglog((info(:,1)), (info(:,3)), 'b-*')
  
  semilogy((info(:,5)), (info(:,6)), 'r--o')
  semilogy((info(:,5)), (info(:,7)), 'm-o')
  
  
  xlabel('# degrés de liberté')
  ylabel('||e||')
    
  title(['epsilon = ', num2str(eps)])
%  legend('Norme H1, linéaire', 'Norme L2, linéaire', 'Norme H1, quadratique', 'Norme L2, quadratique')
  legend('Norme H1', 'Norme L2')
  
endfunction