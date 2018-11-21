function genEmaxAll(Domain_set::OrderedDict, MC_ϵ::Array, T)
    global fEmax
    println("\n Backward induction \n")
    println("\n Solving Exact Model \n")
    println("== Iteration t=$T ==\n")
    fEmax=EmaxT(Domain_set[T])
    Emaxall = OrderedDict(T => fEmax)
    for t = reverse(2:T-1)
        global fEmax
        println("== Iteration t=$t ==\n")
        fEmax= Emaxt(Domain_set[t])
        tempDict = OrderedDict(t => fEmax)
        Emaxall = merge(Emaxall,tempDict)
    end
    return Emaxall
end

function EmaxT(SST::Array)
    Emax = [EMAXTi(x) for x in SST]
    return fEmax = OrderedDict(zip(SST,Emax))
end

function Emaxt(SSt::Array)
    Emax = [EMAXti(x) for x in SSt]
    return fEmaxt = OrderedDict(zip(SSt,Emax))
end
