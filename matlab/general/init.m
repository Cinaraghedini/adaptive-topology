function init()
    format long
    clear -global inner outer outerPolygons innerPolygons robots field sectors centroids tools actualIndexes
        
    addpath('./', './lib/core', './lib/support', './test');
    
    global tools
    tools = Tools();    
   
end