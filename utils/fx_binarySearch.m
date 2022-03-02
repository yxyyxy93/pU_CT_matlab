function [index] = fx_binarySearch(A, n, num)

%--------------------------------------------------------------------------
% Syntax:       [index] = binarySearch(A, n, num);
%               
% Inputs:       A: Array (sorted) that you want to search
%               n: Length of array A
%               num: Number you want to search in array A
%               
% Outputs:      index: Return position in A that A(index) == num
%                      or the index of the first number that larger than
%                      num
%               
% Description:  This function find number in array (sorted) using binary
%               search
%               
% Complexity:   O(1)    best-case performance
%               O(log_2 (n))    worst-case performance
%               O(1)      auxiliary space
%               
% Author:       Trong Hoang Vo
%               hoangtrong2305@gmail.com
%               modified by XiaoyuYang 
%               
% Date:         March 31, 2016
%--------------------------------------------------------------------------

left  = 1;
right = n;

while left <= right
    mid = ceil((left + right) / 2);
    if A(mid) == num
        index = mid;
        return;
    elseif A(mid) > num
        right = mid - 1;
    else
        left = mid + 1;
    end
end

index = left; 

end

