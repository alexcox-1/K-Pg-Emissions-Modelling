# processing loscar solution
plot(d13csolncores[:,200,:],label="")

plot!(d13cvals,ribbon=d13cerror,lw=2,label="")


plot(tempsolncores[:,200,:],label="")

plot!(temp,ribbon=temperror,lw=2,label="actual value")

plot(meanslosc,ribbon=stdslosc)

plot!(meanclosc,ribbon=stdclosc)