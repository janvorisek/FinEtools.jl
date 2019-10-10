"""
    ElementalFieldModule

Module for elemental fields.
"""
module ElementalFieldModule

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import ..FieldModule: AbstractField, nents
import ..FieldModule.@add_Field_fields
import ..FieldModule.FieldDOFData

"""
    ElementalField{T} <: AbstractField{T}

Elemental field, meaning the entities are finite elements.

The values in the field are indexed by the element number.  This means  that
there needs to be one field per finite element set.
"""
mutable struct ElementalField{T} <: AbstractField{T}
	@add_Field_fields(T)
end

"""
   ElementalField(datamatrix::FMat{T}=[]) where {T<:Number}

Constructor of elemental field. The values of the field are given by the array
on input, `datamatrix`. This array needs to have as many rows as there are elements,
and as many columns as there are degrees of freedom per element.
"""
function ElementalField(datamatrix::FMat{T}=[]) where {T<:Number}
	data = [[FieldDOFData(datamatrix[idx, j], zero(T), 0, false) for j in 1:size(datamatrix, 2)] for idx in 1:size(datamatrix, 1)]
    nfreedofs = 0
	return ElementalField(data, nfreedofs)
end

"""
    ElementalField(data::FVec{T}) where {T<:Number}

Constructor of elemental field. The values of the field are given by the vector
on input, `data`. This vector needs to have as many entries as there are elements;
there is just one degree of freedom per element.
"""
function ElementalField(data::FVec{T}) where {T<:Number}
    return ElementalField(reshape(data, length(data), 1))
end

"""
    nelems(self::ElementalField)::FInt = nents(self)

Provide the number of elements  in the elemental field.
"""
nelems(self::ElementalField)::FInt = nents(self)

end
