function genEmaxAll(Domain_set::OrderedDict, T)
    println("\n Backward induction \n")
    println("\n Solving Exact Model \n")
    println("== Iteration t=$T ==\n")
    fEmax = EmaxT(Domain_set[T])
    sendto(procs(),fEmax=fEmax)
    Emaxall = OrderedDict(T => fEmax)
    for t = reverse(2:T-1)
        println("== Iteration t=$t ==\n")
        fEmax = Emaxt(Domain_set[t])
        sendto(procs(),fEmax=fEmax)
        tempDict = OrderedDict(t => fEmax)
        Emaxall = merge(Emaxall,tempDict)
    end
    return Emaxall
end

function EmaxT(SST::Array)
    Emax = pmap_EMAXTi(SST)
    return fEmax = OrderedDict(zip(SST,Emax))
end

# Calculates Emax at t=2,...,T-1
# SSt = State Space at t
function Emaxt(SSt::Array)
    Emax = pmap_EMAXti(SSt)
    return fEmaxt = OrderedDict(zip(SSt,Emax))
end
