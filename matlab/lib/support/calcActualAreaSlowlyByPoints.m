%{
function calcActualAreaSlowlyByPoints: calculate the "actual" by checking
    point by point whether it is covered or not.
    The ratio [qty of covered points] / [total of points] * field's area
    is equal to the actual covered area
        
input: outputFile: file path to save the outcomes.

output: the "actual" covered area used to compare the results 
    obtained by calcArea function.

@author Israel C. S. Rocha (israel.c.rocha@gmail.com)

%}

function actualArea = calcActualAreaSlowlyByPoints(outputFile, areaByAlgorithm)

    disp('Calculando área por pontos...')
    global robots
    field = getField(robots);

    delta = 0.05;
    
    total = 0;
    actualArea = 0;
    [m, ~] = size(robots);
    
    tic
    for(i = field(1) : delta : field(2))        
        for(j = field(3) : delta : field(4))
            total = total + 1;
            
            distances = pdist2([i, j], [robots(:,1), robots(:,2)]) - robots(:,3)';
            distances = min(distances);
            
            if(lesserOrEqualThan(distances, 0))
                actualArea = actualArea + 1;
            end
                        
        end        
    end
    toc
    
    if(~exist('./saves'))
        mkdir ./saves;
    end
        
    actualArea = actualArea / total;
    actualArea = actualArea * (field(4) - field(3)) * (field(2) - field(1))
    
    error = (areaByAlgorithm - actualArea)
    
    disp(['Salvando em: saves/ ' outputFile])
    save(['saves/ ' outputFile], 'actualArea', 'areaByAlgorithm', 'error', 'field', 'robots')
end