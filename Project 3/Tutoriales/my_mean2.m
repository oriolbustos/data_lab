function n = my_mean2(x)
    % MY_MEAN calculates the mean of a vector x
    sum_x = sum(x);       % Sum the elements of the vector x
    count_x = length(x);  % Count the number of elements in x
    n = sum_x / count_x;  % Calculate the mean by dividing the sum by the count
end
