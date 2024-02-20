include("../src/sampling.jl")
include("../src/lambda_guess.jl")
include("../src/lambda_cyclic.jl")

Z4 = FiniteCyclicGroup([4])
H4 = generateSunletArr(Z4)

lambda_guess = lambdaGuess(Z4)
M_guess = getMatrix(lambda_guess, Z4)
rank(M_guess)
keys(lambda_guess)
mueta_guess = lambda_to_mueta(lambda_guess, Z4)
I_guess = getArrangementInequality(mueta_guess, H4)

S4 = sample(Z4, H4, 2^20)
R4 = collect(keys(S4))
mueta_max = S4[R4[1]]
I_max = R4[1]
for I in collect(keys(S4))
    if (getRank(S4[I], Z4) == 16) & (sum(I .⊻ I_guess) == 1)
        I_max = I
        mueta_max = S4[I]
    end
end
for r=10:16
    println(r, " ", sum([getRank(S4[I], Z4) == r for I in R4]))
end
areNeighbors(I_guess, I_max)
I_max .⊻ I_guess
lambda_max = to_lambda(mueta_max[1:4], mueta_max[5:end], Z4)
M_max = getMatrix(lambda_max, Z4)

low_peak = collect(keys(S4))[1]
for s in keys(S4)
    upper_found = false
    r = getRank(S4[s], Z4)

    if r == 16
        upper_found = true
    end

    for t in keys(S4) 
        if sum(s .⊻ t) == 1
            if upper_found
                break
            end

            if getRank(S4[t], Z4) > r
                upper_found = true
            end
        end
    end

    if upper_found == false
        print(s)
        low_peak = s
        break
    end 
end

getRank(S4[low_peak], Z4)

for s in keys(S4)
    if sum(s .⊻ low_peak) == 1
        print(getRank(S4[s], Z4))
    end
end

H4

function matrix_to_latex(M, headers)
    to_return = "\\[\\begin{array}{*{16}c}\n\t"
    for h=eachindex(headers)
        to_return = Base.string(to_return, "\\begin{turn}{90} $(headers[h]) \\end{turn}")
        if h < length(headers)
            to_return = Base.string(to_return, " & ")
        else
            to_return = Base.string(to_return, " \\\\\n")
        end
    end

    for i=1:20
        to_return = Base.string(to_return, "\t")
        for j=1:16
            to_return = Base.string(to_return, M[i,j])
            if j < 16
                to_return = Base.string(to_return, " & ")
            end
        end
        to_return = Base.string(to_return, " \\\\\n")
    end
    to_return = Base.string(to_return, "\\end{array}\\]")
    return to_return
end

write("data/Z4_matrices/max.txt", matrix_to_latex(M_max, getValidGroupTriples(Z4)))
write("data/Z4_matrices/guess.txt", matrix_to_latex(M_guess, getValidGroupTriples(Z4)))

mu_eta = S4[R4[3]]
lambda_ex = to_lambda(mu_eta[1 : 4], mu_eta[(4 + 1):end], Z4)
M_ex = getMatrix(lambda_ex, Z4)
rank(M_ex)

Z3 = FiniteCyclicGroup([3])
H3 = generateSunletArr(Z3)
S3 = sample(Z3, H3, 2^10)
R3 = collect(keys(S3))
ranks = [getRank(S3[I], Z3) for I in R3]
R3[1]
for r=7:9
    println(sum([rank == r for rank in ranks]))
end
S3[R3[1]]