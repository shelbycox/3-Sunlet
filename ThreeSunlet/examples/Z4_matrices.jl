include("../src/sampling.jl")
include("../src/lambda_guess.jl")
include("../src/lambda_cyclic.jl")

Z4 = FiniteCyclicGroup([4])
H4 = generateSunletArr(Z4)
LinearAlgebra.dot(H4[1], mueta_guess)

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

function matrix_to_latex(M, headers)
    print("\\begin{bmatrix}\n\t")
    for h=eachindex(headers)
        print("$(headers[h])")
        if h < length(headers)
            print(" & ")
        else
            print(" \\\\\n")
        end
    end

    for i=eachindex(M)
        print("\t")
        for j=eachindex(M[1])
            print(M[i][j])
            if j < length(M[1])
                print(" & ")
            end
        end
        print("\\\\\n")
    end
    print("\\end{bmatrix}")
end

matrix_to_latex(M_max, getValidGroupTriples(Z4))

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