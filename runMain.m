function runMain

  clc
  close all

  epsilons = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001];
  %epsilons = [1, 0.1];
  logElem = [2:10];
  pDegrees = [1:2];

  entryId = 1;
  
  figureId = 2;
  
  for k = 1:numel(epsilons)
    for j = 1:numel(pDegrees)
      for i = 1:numel(logElem)
        [L2(i, j), H1(i, j)] = main(2^logElem(i), pDegrees(j), epsilons(k));
 
        elemSize2 = 1 / (2^logElem(i));

 
        results(entryId, 3) = logElem(i);
        results(entryId, 2) = pDegrees(j);
        results(entryId, 1) = epsilons(k);
        
        results(entryId, 4) = elemSize2;
        
        results(entryId, 5) = log10(L2(i, j));
        results(entryId, 7) = log10(H1(i, j));

 
        
        if (i > 1)
          elemSize1 = 1 / (2^logElem(i - 1));
          
          alphaL2 = log(L2(i, j) / L2(i-1, j)) / log(elemSize2 / elemSize1);
          alphaH1 = log(H1(i, j) / H1(i-1, j)) / log(elemSize2 / elemSize1);
          
          results(entryId, 6) = alphaL2;
          results(entryId, 8) = alphaH1;
        endif
        
        entryId = entryId + 1;

      end
    end
    
    start1 = entryId - 2*numel(logElem);
    stop1 = start1 + numel(logElem) - 1;
    
    start2 = entryId - numel(logElem);
    stop2 = start2 + numel(logElem) - 1;
    
    
    figure(figureId)
    hold on
    
    plot(results(start1:stop1,3), results(start1:stop1,5), 'r--o')
    plot(results(start1:stop1,3), results(start1:stop1,7), 'b-o')
    plot(results(start2:stop2,3), results(start2:stop2,5), 'k--*')
    plot(results(start2:stop2,3), results(start2:stop2,7), 'c-*')
      
    xlabel('log(1/h)')
    ylabel('log(error norm)')
    
    title(['epsilon = ', num2str(epsilons(k))])
    legend('L2, linéaire', 'H1, linéaire', 'L2, quadratique', 'H1, quadratique')
   
     figureId = figureId + 1;
 
  end
  
  results
%  L2

  
endfunction