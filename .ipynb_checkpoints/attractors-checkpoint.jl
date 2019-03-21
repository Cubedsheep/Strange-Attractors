using Random
using Distributed


#---------------------------------------------------------------------------
# functions to deal with color


"""
    two_schemes(n, cutoff, n_max, scheme1, scheme2)

calculates the color for a given pixel with value n with n_max the max value for a pixel and cuttof the value to switch from scheme1 to scheme2
"""
function two_schemes(n::Int64, cutoff::Int64, n_max::Int64, scheme1, scheme2; gamma1=1, gamma2=1)
    if n < cutoff
        ratio = 1-(1-n/cutoff)^gamma1
        return(get(scheme1, ratio))
    end
    ratio = ((n-cutoff)/(n_max-cutoff))^gamma2
    return(get(scheme2, ratio))
end


"""
function for converting n to a greyscale value, default colormap for
attractors
"""
function map_gray(n, alpha)
    return UInt8(floor(UInt8, (clamp(n^alpha * typemax(UInt8)^(1-alpha), typemin(UInt8), typemax(UInt8)))))
end


function linear_color(n, n_max, alpha, rgb)
    color = rgb * n^alpha/n_max^alpha
    return RGB(color[1], color[2], color[3])
end


function linear_color_white(n, n_max, alpha, rgb)
    if n == 0
        color = [1., 1., 1.]
    else
        color = rgb * n^alpha/n_max^alpha
    end
    return RGB(color[1], color[2], color[3])
end
    
    
function color_map_background(n, n_max, gamma, rgb, base)
    color = [red(rgb), green(rgb), blue(rgb)]
    if n == 0
        color = base*color
    else
        color = (n/n_max)^gamma*color
    end
    return RGB(color[1], color[2], color[3])
end

        
#---------------------------------------------------------------------------
# functions to calculate the attractors

"""
    next_DeJong(X::Array{Float64,1}, params::Array{Float64,1})

calculate the next point of the Peter De Jong attractor starting from X with 
parameter set params = [a, b, c, d]
"""
function next_DeJong(X::Array{Float64,1}, params::Array{Float64,1})
    return [sin(params[1]*X[2])-cos(params[2]*X[1]), sin(params[3]*X[1])-cos(params[4]*X[2])]
end

        
"""
    calc_hist(grid, params, MAX, X, x, y)

Calculate the histogram of a Peter De Jong attractor for a grid with dimensions (x, y),
using X as initial value until a pixel reaches MAX hits.

"""
function calc_hist(grid::Array{Int64,2}, params::Array{Float64,1}, MAX::Int64, X::Array{Float64,1}, x::Int64, y::Int64)
    maxi = 0
    # calculate the points
    while maxi <= MAX
        # calculate the next point
        X = next_DeJong(X, params)
        # add 1 to the cell in the grid containing this point
        i = floor(Int64, div(y * (X[1]+2), 4)) + 1
        j = floor(Int64, div(x * (X[2]+2), 4)) + 1
        grid[i,j] += 1
        if grid[i,j] > maxi
            maxi = grid[i,j]
        end
    end
    return grid
end


"""
    calc_n_points(grid, params, N, X, x, y)

calculates N points of the Peter De Jong attractor with parameters params starting
from the point X on a grid of size x, y
"""
function calc_n_points(grid::Array{Int64,2}, params::Array{Float64,1}, N::Int64,
    X::Array{Float64,1}, x::Int64, y::Int64)
    for i=1:N
        X = next_DeJong(X, params)
        i = floor(Int64, div(y * (X[1]+2), 4)) + 1
        j = floor(Int64, div(x * (X[2]+2), 4)) + 1
        grid[i,j] += 1
    end
    return grid
end
                
                
"""
    calc_hist(grid, params, MAX, X, x, y, n_max)

Calculate the histogram of a Peter De Jong attractor for a grid with dimensions (x, y),
using X as initial value until a pixel reaches MAX hits or n_max iterations are reached.

"""
function calc_hist(grid::Array{Int64,2}, params::Array{Float64,1}, MAX::Int64, X::Array{Float64,1}, x::Int64, y::Int64, n_max)
    maxi = 0
    n = 0
    # calculate the points
    while (maxi <= MAX) & (n <= n_max)
        # calculate the next point
        X = next_DeJong(X, params)
        # add 1 to the cell in the grid containing this point
        i = floor(Int64, div(y * (X[1]+2), 4)) + 1
        j = floor(Int64, div(x * (X[2]+2), 4)) + 1
        grid[i,j] += 1
        if grid[i,j] > maxi
            maxi = grid[i,j]
        end
        n += 1
    end
    return grid
end


"""
    peter_de_jong(params, MAX, X, Y, dX, dY, spacing="centered", N_max=-1, threads=1)

calculate the histogram for the peter de jong attractor

# Arguments
- `params::Array{Float64,1}`: array containing the for parameters that define the attractor
- `MAX::Int64`: the max amount of points in 1 pixel, when this is reached the function stops
- `X::Int64`: the amount of pixels in the horizontal direction (including borders)
- `Y::Int64`: the amount of pixels in the vertical direction(including borders)
- `dX::Int64`: the amount of pixels to use as horizontal border
- `dY::Int64`: the amount of pixels to use as vertical border
- `spacing::String=centered`: one of "stretched", "center", "top-right" or "bottom-left". stretched will stretch the attractor to fit exactly in the given space, all the other options respect aspect ratios.
- `N_max::Int64=-1`: the maximum amount of iterations to use (per thread). If set to -1 this value is ignored
- `threads::Int64=1`: the amount of threads to use to do the calculation
"""
function peter_de_jong(params::Array{Float64,1}, MAX::Int64, X::Int64, Y::Int64, dX::Int64, dY::Int64; spacing::String="centered", N_max::Int64=-1, threads::Int64=1)
    # determine the canvas to draw on (size - borders)
    width = X - 2*dX
    height = Y - 2*dY
    side = 0
    grids = []
    # initialize the grid to calculate the attractor on
    if spacing == "stretched"
        grids = fill(zeros(Int64, height, width), threads)
    elseif spacing in ["centered", "top-right", "bottom-left"]
        side = min(width, height)
        grids = fill(zeros(Int64, side, side), threads)
    else
        error("""unknown option for spacing, use on of "stretched", "centered", "top-right" or "bottom-left" """)
    end
            
    # initialize arrays for the parralel map
    PARAMS = fill(params, threads)
    MAX = fill(ceil(Int64, MAX/threads), threads)
    WIDTH = fill(width, threads)
    HEIGHT = fill(height, threads)
    if spacing in ["centered", "top-right", "bottom-left"]
        WIDTH = fill(side, threads)
        HEIGHT = fill(side, threads)
    end
    N_MAX = fill(N_max, threads)
    rng = MersenneTwister(12345)
    X0 = [rand(rng, Float64, 2)*4 .- 2 for i=1:threads]
    
    # calculate the attractor
    if N_max == -1
        grids = pmap(calc_hist, grids, PARAMS, MAX, X0, WIDTH, HEIGHT)
    else
        grids = pmap(calc_hist, grids, PARAMS, MAX, X0, WIDTH, HEIGHT, N_MAX)
    end
    
    # add all grids together
    grid = grids[1]
    if threads > 1
        for i=1:threads
            grid += grids[i]
        end
    end
    
    # stitch the borders to the grid
    if spacing == "centered"
        left = zeros(Int64, side, floor(Int64, (X-width)/2))
        right = zeros(Int64, side, ceil(Int64, (X-width)/2))
        top = zeros(Int64, floor(Int64, (Y-side)/2), X)
        bottom = zeros(Int64, ceil(Int64, (Y-side)/2), X)
    elseif spacing == "top-right"
        left = zeros(Int64, side, dX+(width-side))
        right = zeros(Int64, side, dX)
        top = zeros(Int64, dY, X)
        bottom = zeros(Int64, dY+(height-side), X)
    elseif spacing == "bottom-left"
        left = zeros(Int64, side, dX)
        right = zeros(Int64, side, dX+(width-side))
        top = zeros(Int64, dY+(height-side), X)
        bottom = zeros(Int64, dY, X)
    else
        left = zeros(Int64, height, dX)
        right = zeros(Int64, height, dX)                        
        top = zeros(Int64, dY, X)                        
        bottom = zeros(Int64, dY, X)
    end
                            
    grid = hcat(left, grid, right)
    grid = vcat(top, grid, bottom)
    return grid
end
                        
                        
function peter_de_jong(params::Array{Float64}, MAX::Int64, SIDE::Int64; threads::Int64=1)
    return peter_de_jong(params, MAX, SIDE, SIDE, 0, 0, threads=threads)
end
