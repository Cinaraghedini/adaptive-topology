classdef dlnode < handle
% DLNODE  A class to represent a doubly-linked list node.
% Multiple dlnode objects may be linked together to create linked listes.
% Each node contains a piece of data and provides access to the next
% and previous nodes.
   properties
      Data
   end
   properties(SetAccess = private)
      Next
      Prev
   end
    
   methods
      function node = dlnode(Data)
      % DLNODE  Constructs a dlnode object.
         if nargin > 0
            node.Data = Data;
         end
      end
      
      function insertAfter(newNode, nodeBefore)
      % insertAfter  Inserts newNode after nodeBefore.
         disconnect(newNode);
         newNode.Next = nodeBefore.Next;
         newNode.Prev = nodeBefore;
         if ~isempty(nodeBefore.Next)
            nodeBefore.Next.Prev = newNode;
         end
         nodeBefore.Next = newNode;
      end
      
      function insertBefore(newNode, nodeAfter)
      % insertBefore  Inserts newNode before nodeAfter.
         disconnect(newNode);
         newNode.Next = nodeAfter;
         newNode.Prev = nodeAfter.Prev;
         if ~isempty(nodeAfter.Prev)
             nodeAfter.Prev.Next = newNode;
         end
         nodeAfter.Prev = newNode;
      end 

      function disconnect(node)
      % DISCONNECT  Removes a node from a linked list.  
      % The node can be reconnected or moved to a different list.
         Prev = node.Prev; %#ok<*PROP>
         Next = node.Next;
         if ~isempty(Prev)
             Prev.Next = Next;
         end
         if ~isempty(Next)
             Next.Prev = Prev;
         end
         node.Next = [];
         node.Prev = [];
      end
      
      function delete(node)
      % DELETE  Deletes a dlnode from a linked list.
         disconnect(node);
      end        
      function disp(node)
      % DISP  Displays a link node.
         disp('Doubly-linked list node with data:')
         disp(node.Data);
      end
   end % methods
end % classdef

    
    