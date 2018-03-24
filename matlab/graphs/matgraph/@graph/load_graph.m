function load_graph(g,filename)
% load(g,filename) --- read a saved graph on disk
% reads in a graph saved by a previous call to save(g,filename)


fileID= fopen(filename,'r');

if fileID == -1
    error(['File "', filename, '" cannot be opened for reading']);
end

delimiter = '\t';
formatSpec = '%s%s%[^\n\r]';

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

fclose(fileID);

elista = dataArray{:, 1};

nl = length(elista);

for k=1:nl
    line = elista{k};
    eval(line);
end

resize(g,0);

if (sp == 1)
    sparse(g)
else
    full(g)
end

resize(g,nverts);
add(g,elist);

if (length(xy) > 0)
    embed(g,xy);
end

n = nv(g);

if (length(labs)>0)
    for k=1:n
        label(g,k,labs{k});
    end
end