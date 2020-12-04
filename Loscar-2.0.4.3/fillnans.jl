function lininterp(datain::Vector{Float64},ind::Int,maxgap::Int)
	xl,xr,yl,yr = prepwindow(datain,ind);
	if xr - xl <= maxgap+1 # make sure the window is not too long
		return (yr-yl)/(xr-xl)*(ind-xl) + yl; # (slope)*distance + offset
	else
		return NaN;
	end
end

function prepwindow(datain::Vector{Float64},ind::Int)
	li = datain[1:ind]; # values left of NaN
	ri = datain[ind:end]; # values right of NaN
	rl = .!isnan.(li); # find all valid values on left
	rr = .!isnan.(ri); # ------------------------ right
	xl = 1:ind |> x -> x[rl][end]; # get only left x coodinate closest to NaN
	xr = ind:length(datain) |> x -> x[rr][1];#right--------------------
	yl = li[rl][end]; # y coordinate on left
	yr = ri[rr][1];	  # ----------------right
	return xl,xr,yl,yr
end

function fillnans(datavec::Vector{Float64},maxgap::Int)
    nanlines = .!isnan.(datavec);
    dataout = copy(datavec);
    for i in 2:length(datavec)-1
        if nanlines[i]==false # only for NaNs
            # Get maximum possible data range (the maxgap length will be checked inside lininterp fce)
            starti = i <= maxgap ? 1 : i - maxgap
            stopi  = i >= length(datavec)-maxgap ? length(datavec) : i + maxgap
            # Interpolate only if valid data exist at both sides of the NaN
            if any(nanlines[starti:i]) && any(nanlines[i:stopi])
                dataout[i] = lininterp(datavec[starti:stopi],i-starti+1,maxgap);
            end
        end
    end
    return dataout
end