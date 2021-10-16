include("./bezier-fourier.jl")


function chebfit(x, y)
    n = length(x_data)
    xmax = max(x_data)
    xmin = min(x_data)
    x_data = (2*x_data - xmax - xmin)/(xmax - xmin)
    T = zeros(n, n)
    T[:,1] = ones(n,1)
    T[:,2] = x_data
    for j = 3:n
    T[:,j] = 2*x_data.*T[:,j-1] - T[:,j-2]
    end
    b = T \ f_data
    println(b)
    x = (2*x - xmax - xmin)/(xmax - xmin)
    y = zeros(size(x))
    for j = 1:n
        y = y + b[j]*cos( (j-1)*acos(x) )
    end
end
