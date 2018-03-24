function plotGrid(robots, field)

    if(nargin == 1 || ~exist('field', 'var') || isempty(field))
        field = getField(robots);
    end
        
    axis(field);
    axis('square');
    color('b');

    grid on;

    [m, ~] = size(robots);

    for i = 1 : m

        viscircles([robots(i, 1), robots(i, 2)], robots(i, 3), 'EdgeColor', 'g');
        viscircles([robots(i, 1), robots(i, 2)], min(robots(i, 3)/100, 10), 'EdgeColor', 'g');

    end

end