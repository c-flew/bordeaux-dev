include("./bezier-fourier.jl")

using Plots


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
    return b
end

function chebeval(x, b)
    n = length(x)
    xmax = maximum(x)
    xmin = minimum(x)
    x = (2*x .- (xmax + xmin))./(xmax - xmin)
    y = zeros(size(x))
    println(length(x))
    println(size(x))
    for j = 1:n
        y = y + b[j]*cos.( (j-1)*acos.(x) )
    end
    return y
end


ctrl_pts = [(-1.0, -0.0),
            (-0.5, -0.5),
            (0.5, 0.5),
            (1.0, 0.0)]::Vector{Tuple{Float64, Float64}}
spline = bezier_curve(ctrl_pts, 0.1)
#println(spline)
x = collect(I.map(x -> x[1], spline))
y = collect(I.map(x -> x[2], spline))

b = chebfit(x, y)

ŷ = chebeval(x, b)
# ŷ = chebeval(x, [-2.51865E-9, 0.0999683,
# -4.88679E-9,
# -0.0900597,
# -6.12917E-10,
# -0.00868024,
# 1.83171E-9,
# -0.00109985,
# 3.13157E-9,
# -0.000128458,
# 3.05507E-9])
print(ŷ)

plot(x, ŷ)
