"""
    NodalFieldModule

Module for nodal fields.
"""
module NodalFieldModule

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import ..FieldModule: AbstractField, nents
import ..FieldModule.@add_Field_fields
import ..FieldModule.FieldDOFData

"""
    NodalField{T} <: AbstractField{T}

Nodal field, meaning the entities are the finite element nodes.
"""
mutable struct NodalField{T} <: AbstractField{T}
    @add_Field_fields(T)
end

"""
    NodalField(datamatrix::FMat{T}=[]) where {T<:Number}

Constructor of nodal field. The values of the field are given by the array
on input, `datamatrix`. This array needs to have as many rows as there are nodes,
and as many columns as there are degrees of freedom per node.
"""
function NodalField(datamatrix::FMat{T}=[]) where {T<:Number}
    data = [[FieldDOFData(datamatrix[idx, j], zero(T), 0, false) for j in 1:size(datamatrix, 2)] for idx in 1:size(datamatrix, 1)]
    nfreedofs = 0
    return NodalField(data, nfreedofs)
end

"""
    NodalField(data::FVec{T}) where {T<:Number}

Constructor of nodal field. The values of the field are given by the vector
on input, `data`. This vector needs to have as many entries as there are nodes;
there is just one degree of freedom per nodes.
"""
function NodalField(data::FVec{T}) where {T<:Number}
    return NodalField(reshape(data, length(data), 1))
end

"""
    nnodes(self::NodalField)::FInt = nents(self)

Provide the number of nodes  in the nodal field.
"""
nnodes(self::NodalField)::FInt = nents(self)

end
