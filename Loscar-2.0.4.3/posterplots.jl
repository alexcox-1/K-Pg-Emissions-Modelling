emisstime = (1:5000:1495001) .- 500000;
p1 = plot(emisstime,co2vals,lw=3,framestyle=:box,xflip=true,label="",size=(800,300),xtickfontsize=10,ytickfontsize=10)
plot!(emisstime,svals,lw=3,label="",xflip=true)
ylims!(-0.005,0.09)
figpath = "/Users/alexcox/Documents/DartmouthF2020/deccan/figs";
savefig(p1,figpath*"/sampleemiss.pdf")

p2 = plot(timev .- 500000,temp,lw=3,framestyle=:box,xflip=true,label="",size=(800,300),xtickfontsize=10,ytickfontsize=10)
plot!(timev .- 500000,temp .+ temperror,lw=1,label="",xflip=true)
plot!(timev .- 500000,temp .- temperror,lw=1,label="",xflip=true)
plot!(timev .- 500000,mu,lw=3,linestyle=:dot,xflip=true,label="")
savefig(p2,figpath*"/sampleoutput.pdf")

p3 = plot(emisstime,co2vals,lw=1,framestyle=:box,xflip=true,label="",size=(800,300),xtickfontsize=10,ytickfontsize=10,linestyle=:dash)
plot!(emisstime,svals,lw=1,label="",xflip=true,linestyle=:dash)
newco2vals = zeros(300);
newco2vals[1:35] .= 0.055;
newco2vals[36:75] .= 0.048;
newco2vals[76:125] .= 0.003;
newco2vals[106:120] .= 0.063;
newco2vals[121:150] .= 0.054;
newco2vals[151:200] .= 0.0;
newco2vals[201:250] .= 0.070;
newco2vals[251:270] .= 0.032;
newco2vals[271:300] .= 0.042;
newsvals = zeros(300) .+ 0.01;
newsvals[75:125] .= 0.025;
newsvals[150:180] .= 0.005
newsvals[251:300] .= 0.015;
plot!(emisstime,newco2vals,lw=3,label="",xflip=true)
plot!(emisstime,newsvals,lw=3,label="",xflip=true)
ylims!(-0.005,0.09)
savefig(p3,figpath*"/newemiss.pdf")

