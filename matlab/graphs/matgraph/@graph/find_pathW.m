function p = find_pathW(g,u,v)
% find_pathM(g,u,v) --- find a shortest path from u to v
% wv matriz de pesos

n = nv(g);



if (u<1) | (u>n) | (v<1) | (v>n)
    p = [];
    return
end

if u==v
    p = u;
    return
end

q_init(n+1); 
track = zeros(1,n);
track(v) = v;

q_push(v); 
c=0
while(q_size > 0)
    t = q_pop_front;
    if t==u
        c=c+1  
        p (c,:)= [];
        last = u;
        while (last ~= v)
            p = [p,last];
            last = track(last);
        end
        p = [p,v];
    end
    push_list = neighbors(g,t);
    for s = push_list
        if (track(s) == 0)
            track(s) = t;
            q_push(s);
        end
    end
end

if track(u) == 0
    p = [];
    return
end

p = [];
last = u;
while (last ~= v)
    p = [p,last];
    last = track(last);
end
p = [p,v];


q_init(1);
