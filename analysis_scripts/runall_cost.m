% script for looping over all hier data for calculations
cd('/Volumes/HelmsDeep/IBAMR/peri-gPC/peri-pinch')
load('allpara.txt');

n = length(allpara);
%n=5;
cost = zeros(n,2);

for i =1:n
    
    if allpara(i,3)<=0.7 
        endtime=350000;
    elseif allpara(i,3)<1.0
        endtime=300000;
    else
        endtime=250000;
    end
    
    [cost(i,1),cost(i,2)] = compute_Work_and_Cost_Of_Transport(i,endtime);
    
    cd('/Volumes/HelmsDeep/IBAMR/peri-gPC/peri-pinch')
    
end

csvwrite(['cost_of_transport_',num2str(n),'_',date,'.csv'],cost)