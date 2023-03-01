function sample(group::FiniteCyclicGroup, H, N::Int)
    sample = Dict()
    for i=1:N
        new_point = getSamplePoint(group)
        sample[ineq(new_point, H)] = new_point
    end
    return sample
end

function getSamplePoint(group::FiniteCyclicGroup)
    groupSize = getGroupSize(group)
    return rand(Float64, 2*groupSize - 1)
end