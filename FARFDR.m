function [FAR,FDR] = FARFDR(res)
    times = length(res);
    FAR = sum(res(:,1:times/2) >1)/(times/2) ;
    FDR = 1-sum(res(:,times/2+1:end)<1)/(times/2) ;
end