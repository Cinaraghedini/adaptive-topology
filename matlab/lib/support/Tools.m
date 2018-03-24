function t = Tools()

    t.toRelativeIndex = @toRelativeIndex;
    t.toActualIndex = @toActualIndex;

    %{
    function toRelativeIndex: convert real robots id to relative 1 to n index.

    input: 
        actualIndex: true not sequential index.
        actualIndexes: GLOBAL vector mapping all actual indexes.

    output: 
        relative index 1 to n index.
    %}
    function k = toRelativeIndex(actualIndex)
        global actualIndexes;
        k = find(actualIndexes == actualIndex);
    end

    %{
    function toRealIndex: convert relative 1 to n index to actual index

    input: 
        relativeIndex: relative 1 to n sequential index.
        actualIndexes: GLOBAL vector mapping all actual indexes.

    output: 
        actual index.
    %}

    function k = toActualIndex(relativeIndex)
        global actualIndexes;
        
        if(isempty(actualIndexes))
            k = relativeIndex;
            return;
        end
        k = actualIndexes(relativeIndex);
    end




end