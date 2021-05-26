function [data] = sym_to_data_ang(symbol,N)
%SYM_TO_DATA returns an N-sample upchirp of data symbol

    data = [];
    accumulator = 0;
    
    for j = symbol
        phase = -pi + ((j-1)*(2*pi/(N)));
        temp = [];

        for i = 1:N
            accumulator = accumulator + phase;
            polar_radius = 1;

            [x, y] = pol2cart(accumulator, polar_radius);

            temp(i) = complex(x, y);

            phase = phase + (2*pi/(N));
        end
        data = [data temp];
    end
end

