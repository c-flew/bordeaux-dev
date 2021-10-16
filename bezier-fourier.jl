using FFTW
using ApproxFun


const I = Iterators
#const Point = Tuple{Real, Real}
const Point = Tuple{Float64, Float64}
const Spline = Vector{Vector{Float64}}


omega(degree::Integer, k::Integer)::ComplexF64 = exp((-2*pi*im) / (degree+1))^k

transform_ctrl_pt(ctrl_pt::Point) = (ctrl_pt[1] + ctrl_pt[2], ctrl_pt[1] - ctrl_pt[2])

function bezier_curve(ctrl_pts::Vector{Point}, percision::Float64 = 0.01, position=[])::Spline
    degree = length(ctrl_pts) - 1

    if percision >= 1 || percision <= 0
        error("percision must be in range 0..1")
    end

    omegas = omega.(degree, 0:degree)

    ctrl_pts = transform_ctrl_pt.(ctrl_pts)

    Q = ifft(reshape(collect(I.flatten(ctrl_pts)), (2, length(ctrl_pts))))

    if isempty(position)
        spline = Spline(undef, convert(Int, 1/percision)+1)

        for s in 0:convert(Int, 1/percision)
            spline[s+1] = sum((x -> real(Q[1:2, x] * (1 + (s*percision) * (omegas[x] - 1))^degree)).(1:degree+1))
        end

        spline
    else
        (s -> sum((x -> real(Q[1:2, x] * (1 + (s*percision) * (omegas[x] - 1))^degree)).(1:degree+1))).(position)
    end

end

function curvature(ctrl_pts::Vector{Point}, percision::Float64 = 0.01)::Vector{Float64}

    degree = length(ctrl_pts) - 1

    if percision >= 1 || percision <= 0
        error("percision must be in range 0..1")
    end

    d1_ctrl_pts = d_ctrl_pts(ctrl_pts, degree)
    d1 = bezier_curve(d1_ctrl_pts, percision)::Spline

    d2_ctrl_pts = d_ctrl_pts(d1_ctrl_pts, degree-1)
    d2 = bezier_curve(d2_ctrl_pts, percision)::Spline

    curvature = Vector{Float64}(undef, length(d1))
    for i in 1:length(d1)
        curvature[i] = (d1[i][1] * d2[i][2] - d1[i][2] * d2[i][1]) / (sqrt(d1[i][1]^2 + d1[i][2]^2)^3)
    end

    curvature
end

function choose_waypts(pt1, pt2, t1, t2)
    waypts = fill((0.0, 0.0), 4)
    waypts[1] = Tuple(pt1)
    waypts[2] = Tuple(pt1 + (1/5)*t1)
    waypts[3] = Tuple(pt2 - (1/5)*t2)
    waypts[4] = Tuple(pt2)

    return waypts
end

# function cheby_bezier(ctrl_pts::Vector{Point}, percision::Float64 = 0.01)
#     fx = Fun(x->bezier_curve(ctrl_pts, percision, x)[1][1], 0..1)
#     fy = Fun(x->bezier_curve(ctrl_pts, percision, x)[1][2], 0..1)

#     (fx, fy)
# end

# eval_cheby_bezier(f, t)::Point = (f[1](t), f[2](t))
