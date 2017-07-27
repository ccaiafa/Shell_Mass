function [ scan ] = clear_flags( scan )
P = size(scan,2);
for n = 1:P
    scan{n}.flag = 0;
end
end

