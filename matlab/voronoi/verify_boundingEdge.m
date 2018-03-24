function boundingSet = verify_boundingEdge(v1,center,r)

boundingSet=1;

if pdist([v1; center],'euclidean') > r
   boundingSet(1,1)=-1;
end
   