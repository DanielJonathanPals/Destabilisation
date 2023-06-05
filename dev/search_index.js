var documenterSearchIndex = {"docs":
[{"location":"matrices/","page":"Matrices","title":"Matrices","text":"CurrentModule = Destabilisation","category":"page"},{"location":"matrices/#Matrices","page":"Matrices","title":"Matrices","text":"","category":"section"},{"location":"matrices/","page":"Matrices","title":"Matrices","text":"Here some usefull functions are presented which return matrices which are for instance needed to fit a model to the recored data.","category":"page"},{"location":"matrices/","page":"Matrices","title":"Matrices","text":"slice_traj\nY\nZ_t\nZ\nId\nJ","category":"page"},{"location":"matrices/#Destabilisation.slice_traj","page":"Matrices","title":"Destabilisation.slice_traj","text":"slice_traj(data::Matrix; p::Int64=1, T::Int64=100)\n\nSlices out the last T+p columns from the data argument. Here p denotes the number of presamples and T encodes the number of actual samples. \n\n\n\n\n\n","category":"function"},{"location":"matrices/#Destabilisation.Y","page":"Matrices","title":"Destabilisation.Y","text":"Y(traj::Matrix; p::Int64=1, T::Union{Int64,Nothing}=nothing)\n\nThis function turns a trajectory traj, which is a matrix with K rows into a matrix, Y as defined     in Lütkepohls book in equation (3.2.1) by slicing the last T colomns from the traj argument.      If T===nothing then T is set to the number of     columns of traj minus p i.e. the traj argument is interpreted in such a way that it     consists of all the presamples plus all the actual samples. If on the other hand T is given     and T+p is not equal to the length of the trajectory, i.e. traj has less than T+p columns,     an error is throuwn.     As always p denotes the number of presamples and T encodes the number of actual samples. \n\n\n\n\n\n","category":"function"},{"location":"matrices/#Destabilisation.Z_t","page":"Matrices","title":"Destabilisation.Z_t","text":"Z_t(traj::Matrix, t; p::Int64=1)\n\nReturns the vector Zₜ as defined e.g. in equation (3.2.1) in Lütkepohls book. Here it is assumed     that the trajectory traj includes the presample datapoints i.e. yₜ from the book corresponds     to the column of traj with index t+p.     As always p denotes the number of presamples.\n\n\n\n\n\n","category":"function"},{"location":"matrices/#Destabilisation.Z","page":"Matrices","title":"Destabilisation.Z","text":"Z(traj; p::Int64=1, T::Union{Int64,Nothing}=nothing)\n\nReturns Z as defined in Lütkepohls book in equation (3.2.1). If the number of columns of     the trajectory traj does not equal T+p then an error is throuwn.      As always p denotes the number of presamples and T encodes the number of actual samples. \n\n\n\n\n\n","category":"function"},{"location":"matrices/#Destabilisation.Id","page":"Matrices","title":"Destabilisation.Id","text":"Id(d)\n\nReturns a d by d sparse identity matrix\n\n\n\n\n\n","category":"function"},{"location":"matrices/#Destabilisation.J","page":"Matrices","title":"Destabilisation.J","text":"J(K::Int64, p::Int64)\n\nReturns the matrix J as defined in equation (2.1.11)\n\n\n\n\n\n","category":"function"},{"location":"home/","page":"Home","title":"Home","text":"CurrentModule = Destabilisation","category":"page"},{"location":"home/#FiveBoxModel","page":"Home","title":"FiveBoxModel","text":"","category":"section"},{"location":"home/","page":"Home","title":"Home","text":"Documentation for Destabilisation.","category":"page"},{"location":"home/","page":"Home","title":"Home","text":"This model can be use to destabilize stochastic dynamical systems.","category":"page"},{"location":"formatTests/","page":"Format Test","title":"Format Test","text":"CurrentModule = Destabilisation","category":"page"},{"location":"formatTests/#Format-Test","page":"Format Test","title":"Format Test","text":"","category":"section"},{"location":"formatTests/","page":"Format Test","title":"Format Test","text":"These functions can be used to verify if a given object has indeed the desired format needed for the respecive application.","category":"page"},{"location":"formatTests/","page":"Format Test","title":"Format Test","text":"check_traj\ncheck_DynamicalSystem","category":"page"},{"location":"formatTests/#Destabilisation.FormatTests.check_traj","page":"Format Test","title":"Destabilisation.FormatTests.check_traj","text":"check_traj(traj)\n\nChecks if traj has the format of a trajectory i.e. if traj is a two dimentional matrix with all entries     given by numbers. Further it is always assumed that each column of trajectory represents a datapoint     of the time series, i.e. each row denotes the time evolution of a single variable.\n\n\n\n\n\n","category":"function"},{"location":"formatTests/#Destabilisation.FormatTests.check_DynamicalSystem","page":"Format Test","title":"Destabilisation.FormatTests.check_DynamicalSystem","text":"check_DynamicalSystem(progressor,observable,x_init,p_init,p_bounds,x_length,p_length,obs_length)\n\nCompatibility test for initialisation of an object of DynamicalSystem.\n\n\n\n\n\ncheck_DynamicalSystem(progressor,observable,x_init,p_init,p_bounds)\n\nCompatibility test for initialisation of an object of DynamicalSystem.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Index","title":"Index","text":"CurrentModule = Destabilisation","category":"page"},{"location":"#Destabilisation","page":"Index","title":"Destabilisation","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"Documentation for Destabilisation.","category":"page"},{"location":"","page":"Index","title":"Index","text":"","category":"page"},{"location":"coupling/","page":"Coupling to dynamical system","title":"Coupling to dynamical system","text":"CurrentModule = Destabilisation","category":"page"},{"location":"coupling/#Coupling-the-destabilataions-package-to-the-dynamical-system","page":"Coupling to dynamical system","title":"Coupling the destabilataions package to the dynamical system","text":"","category":"section"},{"location":"coupling/","page":"Coupling to dynamical system","title":"Coupling to dynamical system","text":"The struct DynamicalSystem allows acces to the parameter dependent dynamics of the underlying system and also gives information about how to comupte the observables of interest.","category":"page"},{"location":"coupling/","page":"Coupling to dynamical system","title":"Coupling to dynamical system","text":"DynamicalSystem\nDynamicalSystem(progressor,observable,x_init,p_init;p_bounds=nothing)","category":"page"},{"location":"coupling/#Destabilisation.DynamicalSystem","page":"Coupling to dynamical system","title":"Destabilisation.DynamicalSystem","text":"DynamicalSystem\n\nAn instance of the struct DynamicalSystem contains all the relevant information of the stochastic dynamical      system which we want to destabilize.\n\nFields\n\nprogressor::Function: The Argument progressor should contains a function that takes two arguments, each of type   Vector{Flaot64}, representing the current state variables x and the current parameter values p and returns   the state variables as a Vector{Flaot64} for the next timestep. I.e. this function describes the progression   of the dynamical system for one time step\nobservable::Function: This function should take the current state variables x and the current parameter values p,   both as Vector{Float64}, as input and return an output of type Vector{Float64} containing all the observables   of the system which are relevant for further analysis.\nx_init::Vector{Float64}: Initial values for the state space variables.\np_init::Vector{Float64}: Initial parameter values.\np_bounds::Union{Vector{Tuple},Nothing}: In case the parameter values are bounded to certain intervalls, the   bounds of each of these intervalls can be specified in this argument. I.e. the i-th tuple in p_bounds    represent the interval bounds for the i-th parameter. If p_bounds = nothing then it is assumes that there   are no restictions to the parameter values\nx_length::Int64: Number of state space variables\np_length::Int64: Number of parameters\nobs_length::Int64: Number of observables\n\n\n\n\n\n","category":"type"},{"location":"coupling/#Destabilisation.DynamicalSystem-NTuple{4, Any}","page":"Coupling to dynamical system","title":"Destabilisation.DynamicalSystem","text":"DynamicalSystem(progressor,observable,x_init,p_init;p_bounds=nothing)\n\nAlternative initialisation of an element of type DynamicalSystem where x_length, p_length and      obs_length are automatically computed.\n\n\n\n\n\n","category":"method"}]
}
