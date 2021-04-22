p1 = plot(timev,temp,lw=3,framestyle=:box,size=(900,300),label="",xflip=true
,xtickfont=14,ytickfont=14,ribbon=(temperror))
plot!(timev,temps2e4[:,51200],lw=3,label="")
savefig(p1,"/Users/alexcox/Documents/manuscripts/loscarletter/tempcurve.pdf")

p2 = plot(timev,exp.(co22e4[:,51200]),lw=3,framestyle=:box,size=(925,300),
label="", ribbon=(0.005),ytickfont=14,xtickfont=14,xflip=true)
plot!(timev,exp.(so22e4[:,51200]),lw=3,ribbon=(0.002))
savefig(p2,"/Users/alexcox/Documents/manuscripts/loscarletter/emissionscurve.pdf")