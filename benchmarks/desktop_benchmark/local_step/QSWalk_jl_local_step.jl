using QSWalk
using SparseArrays
using LinearAlgebra
using Glob
using HDF5
using MatrixMarket

include("./../rusage.jl")

function natural(x, y)

        k(x) = [occursin(r"\d+", s) ? parse(Int, s) : s

        for s in split(replace(x, r"\d+" => s->" $s "))]

		A = k(x); B= k(y)

        	for (a, b) in zip(A, B)

        	        if !isequal(a, b)
				return typeof(a) <: typeof(b) ? isless(a, b) :
        	                isa(a,Int) ? true : false
			end
        	end

        	return length(A) < length(B)
end

function run_test(omega, files, files_sym, results_file, times_file, max_sim_time)

	total_sim_time = 0

	for (file, file_sym) in zip(files, files_sym)

		G = sparse(MatrixMarket.mmread(file))
	    	G_sym = sparse(MatrixMarket.mmread(file_sym))

		n = size(G)

		total_step_time = time()

		lind_local = local_lind(G)
                rho =  Matrix{Float64}(I, n[1], n[1])/n[1]

		so_time = time()
		generator = evolve_generator(G_sym,lind_local,omega)
		so_time = time() - so_time

		step_time = time()
		rhot=evolve(generator, rho, 100.0)
		step_time = time() - step_time

		peak_memory = get_vmsize()*(1024^(-2))
		total_step_time = time() - total_step_time
		total_sim_time += total_step_time

		outputName = splitext(basename(file))[1]

		if isfile("$results_file.h5")

			save_time = time()
			h5open("$results_file.h5", "r+") do file
				write(file, "$outputName/Re", real(rhot))
				write(file, "$outputName/Im", imag(rhot))
			end
			save_time = time() - save_time

		else

			save_time = time()
			h5open("$results_file.h5", "w") do file
				write(file, "$outputName/Re", real(rhot))
				write(file, "$outputName/Im", imag(rhot))
			end
			save_time = time() - save_time
		end

		nnz = count(!iszero, G)
		so_nnz = count(!iszero, generator)
		dim = n[1]

		open("$times_file.csv", "a") do file
			write(file, "$outputName, $dim, $nnz, $so_nnz, $so_time, $step_time, $total_step_time, $peak_memory, $save_time\n")
		end

		if total_sim_time > max_sim_time
			break
		end

	end
end

times_file = "Results/QSWalk_jl_local_step"
max_sim_time = parse(Float64, ARGS[4])

if ! isfile("$times_file.csv")
	open("$times_file.csv", "w") do file
		write(file, "name,dim,nnz,SO_nnz,SO_time,step_time,total_time,peak_memory,save_time\n")
	end
end

omega = 0.1

files=sort(glob(string(ARGS[1],"*.mtx")), lt=natural)
files_sym=sort(glob(string(ARGS[2],"*.mtx")), lt=natural)
run_test(omega, files, files_sym, string(ARGS[3],"_QSWalk_jl"), times_file, max_sim_time)
