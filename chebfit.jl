include("./bezier-fourier.jl")


function chebfit(x, y)
    n = length(x)
    xmax = maximum(x)
    xmin = minimum(x)
    x = (2*x .- (xmax + xmin))./(xmax - xmin)
    T = zeros(n, n)
    T[:,1] = ones(n,1)
    T[:,2] = x
    for j = 3:n
        T[:,j] = 2*x.*T[:,j-1] - T[:,j-2]
    end
    b = T \ y
    println(b)
end

function chebeval(x, b)
    n = length(x)
    xmax = maximum(x)
    xmin = minimum(x)
    x = (2*x - xmax - xmin)/(xmax - xmin)
    y = zeros(size(x))
    for j = 1:n
        y = y + b[j]*cos( (j-1)*acos(x) )
    end
end


ctrl_pts = [(-1.0, -0.0),
            (-0.5, -0.5),
            (0.5, 0.5),
            (1.0, 0.0)]::Vector{Tuple{Float64, Float64}}
spline = bezier_curve(ctrl_pts, 0.01)
#println(spline)
x = collect(I.map(x -> x[1], spline))
y = collect(I.map(x -> x[2], spline))

chebfit(x, y)
