"""
    GeneralFieldModule

Module for general fields.
"""
module GeneralFieldModule

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import ..FieldModule.AbstractField
import ..FieldModule.@add_Field_fields
import ..FieldModule.FieldDOFData

"""
    GeneralField{T} <: AbstractField{T}

General field, meaning the entities can be anything.
"""
mutable struct GeneralField{T} <: AbstractField{T}
    @add_Field_fields(T)
end


"""
    GeneralField(data::FMat{T}=[]) where {T<:Number}

Constructor of general field.  The values of the field are given by the array
on input, `data`. This array needs to have as many rows as there are entities,
and as many columns as there are degrees of freedom per entities.
"""
function GeneralField(datamatrix::FMat{T}=[]) where {T<:Number}
    data = [[FieldDOFData(datamatrix[idx, j], zero(T), 0, false) for j in 1:size(datamatrix, 2)] for idx in 1:size(datamatrix, 1)]
    nfreedofs = 0
    return GeneralField(data, nfreedofs)
end

"""
    GeneralField(data::FVec{T}) where {T<:Number}

Constructor of general field.  The values of the field are given by the vector
on input, `data`. This vector needs to have as many rows as there are entities.
"""
function GeneralField(data::FVec{T}) where {T<:Number}
    return GeneralField(reshape(data, length(data), 1))
end

end
