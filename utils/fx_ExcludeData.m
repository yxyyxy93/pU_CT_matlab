function inph = fx_ExcludeData(inph, varargin)
% exclude the extra data before and after
% input (inph, 'f', front_I, 'r', rear_I)

[lx, ly, ~] = size(inph);

if (varargin{1} == 'f')
    front_I = varargin{2};
    for i=1:lx
        for j=1:ly    
            inph(i, j, 1:round(front_I(i, j))) = 0;
        end
    end
end
    
if (varargin{3} == 'r')
    rear_I = varargin{4};
    for i=1:lx
        for j=1:ly    
            inph(i, j, round(rear_I(i, j) + 1):end) = 0;
        end
    end
end

end

