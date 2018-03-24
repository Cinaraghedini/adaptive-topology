
function [data] = create_Ninf_g(matrix)

data=[];
for i=1:size(matrix,2)
    ddsN=[];
    listNi=kneighbors(matrix,i,1);
    for j=listNi
        ddsNj=[];
        listNj=kneighbors(matrix,j,1);
        for k=listNj
            if k~=i
                ddsNj=[ddsNj; k];
            end
        end
        ddsN=[ddsN; {i} {j} {ddsNj}];
    end
    data=[data; {ddsN}];
end