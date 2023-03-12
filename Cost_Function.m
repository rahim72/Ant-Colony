
function Cost=Cost_Function(tour,sys)
    Cost=0;
    numbertour=numel(tour);
    tour=[tour tour(1)];
    
    for i=1:numbertour
        Cost=Cost+sys.D(tour(i),tour(i+1));
    end

end