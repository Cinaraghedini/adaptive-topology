% Subject: Programming data structures in Matlab 14 (HOWTO)
% From: sturla.molden (https://www.mathworks.com/matlabcentral/newsreader/author/45479)
% Date: 24 Oct, 2004 02:11:59
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/79044


function f = Stack()
    
    f.push = @push;
    f.pop = @pop;
    f.top = @top;
    f.bottom = @bottom;
    f.isEmpty = @isEmpty;
    f.getSize = @getSize;
    f.flush = @flush;
    s = [];
    bot = [];
    size = 0;
    
    function s = getSize()
        s = size;
    end

    function flush()
        s = [];
        bot = [];
        size = 0;
    end    

    function push(data)      
      wasEmpty = isEmpty();
        
      tmp = Node;
      tmp.setData(data);
      tmp.setNext(s);
      s = tmp;
      size = size + 1;
      
      disp(['Empilhando ' num2str(data)]);
      
      if(wasEmpty)
          bot = top();          
      end
      
    end

    function b = bottom()
        b = bot;
    end

    function y = pop()
      if isempty(s)
         y = [];
      else
         y = top();
         s = s.getNext();
         size = size - 1;
      end
      disp(['Desempilhando ' num2str(y)]);
    end

    function y = top()
         if isempty(s)
             y = [];
         else
            y = s.getData();
         end
    end

    function empty = isEmpty()
        empty = isempty(s);
    end

end

function f = Node()

   d = []; next = [];
   f.getData = @getData;
   f.setData = @setData;
   f.getNext = @getNext;
   f.setNext = @setNext;

   function y = getData
      y = d;
   end

   function setData(data)
      d = data;
   end

   function y = getNext
      y = next;
   end

   function setNext(x)
      next = x;
   end

end