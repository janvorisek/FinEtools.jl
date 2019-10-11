"""
    FieldModule

Module for abstract fields.
"""
module FieldModule

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import Base.copyto!

"""
    fixed_dofnum

    This is how we  recognize fixed degrees of freedom: the degree of freedom
    number is zero
"""
const fixed_dofnum = 0

"""
    FieldDOFData{T}

Type of data stored  per degree of freedom for each entity in a field.
"""
mutable struct FieldDOFData{T}
    dofvalue::T;
    fixeddofvalue::T;
    dofnum::FInt;
    isfixed::Bool;
end

"""
	AbstractField{T}

Abstract field.

Expected  attributes:
  + `data::Vector{Vector{FieldDOFData{T}}}`: Vector of vectors of degree of freedom data per entity
  + `nfreedofs::FInt`: Total number of free degrees of freedom

See also: [`@add_Field_fields()`](@ref) .
"""
abstract type AbstractField{T<:Number} end

"""
    add_Field_fields(T)

Generate the attributes (i. e. fields) of a `Field`. The methods defined for
the abstract type depend on these attributes to be present.
"""
macro add_Field_fields(T)
    return esc(:(
    data::Vector{Vector{FieldDOFData{$T}}};
    nfreedofs::FInt;
    )
    )
end

"""
    ndofs(self::AbstractField{T}) where {T}

How many degrees of freedom per entity?

When the field stores no entities, zero (0) is returned.
"""
ndofs(self::AbstractField{T}) where {T}  = length(self.data) > 0 ? length(self.data[1]) : 0

"""
    nents(self::AbstractField{T}) where {T}

Number of nodes associated with the field.
"""
nents(self::AbstractField{T}) where {T} = length(self.data)

"""
    copyto!(DEST::F,  SRC::F) where {T, F<:AbstractField{T}}

Copy data from one field to another.
"""
function copyto!(DEST::F,  SRC::F) where {T, F<:AbstractField{T}}
    copyto!(DEST.data, SRC.data)
    DEST.nfreedofs = SRC.nfreedofs
    return  DEST
end

"""
    wipe!(self::AbstractField{T}) where {T}

Wipe all the data from the field.

This includes values, prescribed values, degree of freedom numbers, and "is
fixed" flags. The number of free degrees of freedom is set to zero.
"""
function wipe!(self::AbstractField{T}) where {T}
    z = zero(T)
    nent, ndof = nents(self), ndofs(self)
    for i in 1:nent
        for j in 1:ndof
            self.data[i][j].dofvalue = z
            self.data[i][j].fixeddofvalue = z
            self.data[i][j].dofnum = 0
            self.data[i][j].isfixed = false
        end
    end
    self.nfreedofs = 0
    return self
end

"""
    gathersysvec(self::F) where {T, F<:AbstractField{T}}

Gather values from the field for the whole system vector.
"""
function gathersysvec(self::F) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    vec = zeros(T, self.nfreedofs)
    for i in 1:nent
        for j in 1:ndof
            en = self.data[i][j].dofnum
            if (en > 0) && (en <= self.nfreedofs)
                vec[en] = self.data[i][j].dofvalue
            end
        end
    end
    return vec
end

"""
    gathersysvec!(self::F, vec::FVec{T}) where {T, F<:AbstractField{T}}

Gather values from the field for the whole system vector.

The writing to `vec` occurs in-place.
"""
function gathersysvec!(self::F, vec::FVec{T}) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    @assert length(vec) == self.nfreedofs
    for i in 1:nent
        for j in 1:ndof
            en = self.data[i][j].dofnum
            if (en > 0) && (en <= self.nfreedofs)
                vec[en] = self.data[i][j].dofvalue
            end
        end
    end
    return vec
end

"""
    gathervalues_asvec!(self::F, dest::AbstractArray{T, 1}, conn::CC) where {T, F<:AbstractField{T}, CC}

Gather values from the field into a vector.

The order is: for each node  in the connectivity, copy into the buffer all the
degrees of freedom,  then the next node and so on.

`dest` = destination buffer: overwritten  inside,  must be preallocated
in the correct size
"""
function gathervalues_asvec!(self::F, dest::AbstractArray{T, 1}, conn::CC) where {T, F<:AbstractField{T}, CC}
    nent, ndof = nents(self), ndofs(self)
    en::FInt = 1;
    for i in 1:length(conn)
        for j in 1:ndof
            dest[en] = self.data[conn[i]][j].dofvalue;
            en = en + 1;
        end
    end
    return dest
end

"""
    gathervalues_asmat!(self::F, dest::AbstractArray{T, 2}, conn::CC) where {T, F<:AbstractField{T}, CC}

Gather values from the field into a two-dimensional array.

The order is: for each node  in the connectivity, copy into the corresponding
row of the buffer all the degrees of freedom,  then the next node into the next
row and so on.

`dest` = destination buffer: overwritten  inside,  must be preallocated
in the correct size
"""
function gathervalues_asmat!(self::F, dest::AbstractArray{T, 2}, conn::CC) where {T, F<:AbstractField{T}, CC}
    nent, ndof = nents(self), ndofs(self)
    for i in 1:length(conn)
        for j in 1:ndof
            dest[i, j] = self.data[conn[i]][j].dofvalue;
        end
    end
    return dest
end

"""
    gatherfixedvalues_asvec!(self::F, dest::AbstractArray{T, 1}, conn::CC) where {T, F<:AbstractField{T}, CC}

Gather FIXED values from the field into a vector.

The order is: for each node  in the connectivity, copy into the buffer all the
fixed degrees of freedom,  then the next node and so on. If a degree of freedom
is NOT fixed, the corresponding entry is  set to zero.

`dest` = destination buffer: overwritten  inside,  must be preallocated
in the correct size
"""
function gatherfixedvalues_asvec!(self::F, dest::AbstractArray{T, 1}, conn::CC) where {T, F<:AbstractField{T}, CC}
    nent, ndof = nents(self), ndofs(self)
    en::FInt = 1;
    for i in 1:length(conn)
        for j in 1:ndof
            if self.data[conn[i]][j].isfixed # fixed degree of freedom
                dest[en] = self.data[conn[i]][j].fixeddofvalue
            else
                dest[en] = zero(T)
            end
            en = en + 1;
        end
    end
    return dest
end

"""
    gatherfixedvalues_asmat!(self::F, dest::AbstractArray{T, 2}, conn::CC) where {T, F<:AbstractField{T}, CC}

Gather FIXED values from the field into a two-dimensional array.

The order is: for each node  in the connectivity, copy into the corresponding
row of the buffer all the degrees of freedom,  then the next node into the next
row and so on.  If a degree of freedom
is NOT fixed, the corresponding entry is  set to zero.

`dest` = destination buffer: overwritten  inside,  must be preallocated
in the correct size
"""
function gatherfixedvalues_asmat!(self::F, dest::AbstractArray{T, 2}, conn::CC) where {T, F<:AbstractField{T}, CC}
    nent, ndof = nents(self), ndofs(self)
    for i in 1:length(conn)
        for j in 1:ndof
            if self.data[conn[i]][j].isfixed # fixed degree of freedom
                dest[i, j] = self.data[conn[i]][j].fixeddofvalue
            else
                dest[i, j] = zero(T)
            end
        end
    end
    return dest
end

"""
    anyfixedvaluenz(self::F, conn::CC) where {T, F<:AbstractField{T}, CC}

Is any degree of freedom value fixed (prescribed) and nonzero?

Degrees of freedom for all entities listed in the `conn` argument are checked.
Returns a Boolean.
"""
function anyfixedvaluenz(self::F, conn::CC) where {T, F<:AbstractField{T}, CC}
    nent, ndof = nents(self), ndofs(self)
    for i in 1:length(conn)
        for j in 1:ndof
            if self.data[conn[i]][j].isfixed # fixed degree of freedom
                if  abs(self.data[conn[i]][j].fixeddofvalue) > zero(T)
                    return true
                end
            end
        end
    end
    return false
end

"""
    gatherdofnums!(self::F, dest::A, conn::CC) where {T, F<:AbstractField{T}, A, CC}

Gather degree of freedom numbers from the field.

Degrees of freedom for all entities listed in the `conn` argument are collected.
"""
function gatherdofnums!(self::F, dest::A, conn::CC) where {T, F<:AbstractField{T}, A, CC}
    nent, ndof = nents(self), ndofs(self)
    en::FInt = 1;
    for i in 1:length(conn)
        for j in 1:ndof
            dest[en] = self.data[conn[i]][j].dofnum;
            en = en+1;
        end
    end
    return dest
end

"""
    numberdofs!(self::F) where {T, F<:AbstractField{T}}

Number the degrees of freedom.

The free components in the field are numbered consecutively. No effort is
made to optimize the numbering in any way. If you'd like to optimize the
numbering of the degrees of freedom, use the form that sets the
permutation of the degrees of freedom, or the permutation of the nodes.
"""
function numberdofs!(self::F) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    self.nfreedofs::FInt =0
    for i in 1:nent
        for j in 1:ndof
            if !self.data[i][j].isfixed # free degree of freedom
                self.nfreedofs = self.nfreedofs + 1
                self.data[i][j].dofnum = self.nfreedofs
            else # fixed degree of freedom: no equation
                self.data[i][j].dofnum = fixed_dofnum
            end
        end
    end
    return  self
end

"""
    resetdofnums!(self::F) where {T, F<:AbstractField{T}}

Set all degree of freedom numbers to an invalid number.
"""
function resetdofnums!(self::F) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    if self.nfreedofs != 0 # We assume that if this is zero, the zeroing out doesn't need to happen.
        for i in 1:nent
            for j in 1:ndof
                self.data[i][j].dofnum = fixed_dofnum
            end
        end
    end
    self.nfreedofs = 0
end

"""
    resetfixeddofvalues!(self::F) where {T, F<:AbstractField{T}}

Set the fixed values for all degree of freedom numbers to zero.
"""
function resetfixeddofvalues!(self::F) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    for i in 1:nent
        for j in 1:ndof
            self.data[i][j].fixeddofvalue = zero(T)
        end
    end
end

"""
    setebc!(self::F, entids::FIntVec, isfixed::Bool, comp::FInt, val::FVec{T}) where {T, F<:AbstractField{T}}

Set the EBCs (essential boundary conditions).

`entids`         - array of N entity identifiers
`isfixed` = scaler Boolean: are the degrees of freedom being fixed (true)
             or released (false),
`comp` = integer, which  degree of freedom (component),
`val` = array of N values of type T

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::F, entids::FIntVec, isfixed::Bool, comp::FInt, val::FVec{T}) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    @assert 1 <= comp <= ndof "Requested  nonexistent  degree of freedom"
    @assert maximum(entids) <= nent "Requested nonexistent entity"
    @assert size(entids) == size(val) "Arrays of mismatched sizes"
    for  j in 1:length(entids)
        self.data[entids[j]][comp].isfixed = isfixed;
        if self.data[entids[j]][comp].isfixed
            self.data[entids[j]][comp].fixeddofvalue = val[j];
        else
            self.data[entids[j]][comp].fixeddofvalue = zero(T)
        end
    end
    resetdofnums!(self)
    return  self
end

"""
    setebc!(self::F, entids::FIntVec, isfixed::Bool, comp::FInt, val::T) where {T, F<:AbstractField{T}}

Set the EBCs (essential boundary conditions).

`entids`         - array of N node identifiers
`isfixed` = scaler Boolean: are the degrees of freedom being fixed (true)
             or released (false),
`comp` = integer, which  degree of freedom (component),
`val` = scalar of type T

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::F, entids::FIntVec, isfixed::Bool, comp::FInt, val::T) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    @assert 1 <= comp <= ndof "Requested  nonexistent  degree of freedom"
    @assert maximum(entids) <= nent "Requested nonexistent entity"
    for  j in 1:length(entids)
        self.data[entids[j]][comp].isfixed = isfixed;
        if self.data[entids[j]][comp].isfixed
            self.data[entids[j]][comp].fixeddofvalue = val;
        else
            self.data[entids[j]][comp].fixeddofvalue = zero(T)
        end
    end
    resetdofnums!(self)
    return  self
end

"""
    setebc!(self::F, entids::FIntVec, comp::FInt, val::FVec{T}) where {T, F<:AbstractField{T}}

Set the EBCs (essential boundary conditions).

`entids`         - array of N node identifiers
`comp` = integer, which  degree of freedom (component),
`val` = array of N values of type T

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::F, entids::FIntVec, comp::FInt, val::FVec{T}) where {T, F<:AbstractField{T}}
    return setebc!(self, entids, true, comp, val)
end


"""
    setebc!(self::F, entids::FIntVec, comp::FInt; val::T=0.0) where {T, F<:AbstractField{T}}

Set the EBCs (essential boundary conditions).

`entids`         - array of N node identifiers
`comp` = integer, which  degree of freedom (component),
`val` = scalar of type T

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::F, entids::FIntVec, comp::FInt; val::T=zero(T)) where {T, F<:AbstractField{T}}
    return setebc!(self, entids, true, comp, val)
end

"""
    setebc!(self::F, entids::FIntVec, isfixed::Bool, comp::FInt; val::T=zero(T)) where {T, F<:AbstractField{T}}
    j = comp

Set the EBCs (essential boundary conditions).

`entids`         - array of N node identifiers
`comp` = integer, which  degree of freedom (component),
`val` = scalar of type T

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::F, entids::FIntVec, isfixed::Bool, comp::FInt; val::T=zero(T)) where {T, F<:AbstractField{T}}
    j = comp
    @assert (j >= 1) && (j <= ndofs(self))
    setebc!(self, entids, isfixed, j, val)
    return self
end

"""
    setebc!(self::F, entids::FIntVec) where {T, F<:AbstractField{T}}

Set the EBCs (essential boundary conditions).

Suppress all degrees of freedom at the given nodes.

`entids`         - array of N node identifiers

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::F, entids::FIntVec) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    Zer = zero(T)
    for comp in 1:ndof
        setebc!(self, entids, true, comp, Zer)
    end
    return self
end

"""
    setebc!(self::F, fenid::FInt) where {T, F<:AbstractField{T}}

Set the EBCs (essential boundary conditions).

Suppress all degrees of freedom at the given node.

`fenid`         - One integer as a node identifier

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::F, fenid::FInt) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    Zer = zero(T)
    for comp in 1:ndof
        setebc!(self, [fenid], true, comp, Zer)
    end
    return self
end

"""
    setebc!(self::F) where {T, F<:AbstractField{T}}

Set the EBCs (essential boundary conditions).

All essential boundary conditions are CLEARED.

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`nfreedofs = 0`; and sets all degrees of freedom to an invalid number.
"""
function setebc!(self::F) where {T, F<:AbstractField{T}}
    resetdofnums!(self)
    return  self
end

"""
    applyebc!(self::F) where {T, F<:AbstractField{T}}

Apply EBCs (essential boundary conditions).

The prescribed values of the degrees of freedom are copied to the degree of
freedom values for fixed (prescribed) degrees of freedom.
"""
function applyebc!(self::F) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    for i in 1:nent
        for j in 1:ndof
            if self.data[i][j].isfixed
                self.data[i][j].dofvalue = self.data[i][j].fixeddofvalue
            end
        end
    end
    return  self
end

"""
    scattersysvec!(self::F, vec::FVec{T}) where {T, F<:AbstractField{T}}

Scatter values to the field from a system vector.
"""
function scattersysvec!(self::F, vec::FVec{T}) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    for i in 1:nent
        for j in 1:ndof
            dn = self.data[i][j].dofnum
            if (dn > 0) && (dn <= self.nfreedofs)
                self.data[i][j].dofvalue = vec[dn];
            end
        end
    end
    return  self
end

"""
    incrscattersysvec!(self::F, vec::FVec{T}) where {T, F<:AbstractField{T}}

Increment values of the field by scattering a system vector.
"""
function incrscattersysvec!(self::F, vec::FVec{T}) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    for i in 1:nent
        for j in 1:ndof
            dn = self.data[i][j].dofnum
            if (dn > 0) && (dn <= self.nfreedofs)
                self.data[i][j].dofvalue += vec[dn];
            end
        end
    end
    return  self
end

"""
    prescribeddofs(uebc::F, u::F) where {T, F<:AbstractField{T}}

Find which degrees of freedom are prescribed and to which values.

`uebc` = field which defines the constraints (is the dof fixed and to which value?),
`u` = field which does not have the constraints applied, and serves as the source of equation numbers,
`uebc` and `u` may be one and the same field.
"""
function prescribeddofs(uebc::F, u::F) where {T, F<:AbstractField{T}}
    @assert nents(uebc) == nents(u)
    @assert ndofs(uebc) == ndofs(u)
	dofnums = FInt[]
	prescribedvalues = T[]
	nent, ndof = nents(uebc), ndofs(uebc)
	for i in 1:nent
		for j in 1:ndof
			if uebc.data[i][j].isfixed
			    push!(prescribedvalues, uebc.data[i][j].fixeddofvalue);
			    dn = u.data[i][j].dofnum
			    @assert  (dn > 0) && (dn <= u.nfreedofs)
			    push!(dofnums, dn)
			end
		end
	end
	return dofnums, prescribedvalues
end

"""
    dofvaluesasarray(self::F) where {T, F<:AbstractField{T}}

Return the degree of freedom values as an array.

Return the degree of freedom values as a matrix, where the number of rows
matches the number of entities in the field.
"""
function dofvaluesasarray(self::F) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    values = fill(zero(T), nent, ndof)
    for i in 1:nent
        for j in 1:ndof
            values[i, j] = self.data[i][j].dofvalue
        end
    end
    return values
end

"""
    fixeddofvaluesasarray(self::F) where {T, F<:AbstractField{T}}

Return the fixed (prescribed) degree of freedom values as an array.

Return the fixed (prescribed) as a matrix, where the number of rows
matches the number of entities in the field.
"""
function fixeddofvaluesasarray(self::F) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    values = fill(zero(T), nent, ndof)
    for i in 1:nent
        for j in 1:ndof
            values[i, j] = self.data[i][j].fixeddofvalue
        end
    end
    return values
end

"""
    isfixedasarray(self::F) where {T, F<:AbstractField{T}}

Return the is-fixed (prescribed) Boolean flags as an array.

Return the Boolean flags whether a degree of freedom is fixed or free as a
matrix, where the number of rows matches the number of entities in the field.
"""
function isfixedasarray(self::F) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    values = fill(false, nent, ndof)
    for i in 1:nent
        for j in 1:ndof
            values[i, j] = self.data[i][j].isfixed
        end
    end
    return values
end

"""
    setdofvalue!(self::F, ent, dof, val::T) where {T, F<:AbstractField{T}}

Set the value of the degree of freedom at a given entity.
"""
function setdofvalue!(self::F, ent, dof, val::T) where {T, F<:AbstractField{T}}
    self.data[ent][dof].dofvalue = val
    return self
end

"""
    setdofvalues!(self::F, val::T) where {T, F<:AbstractField{T}}

Set the value of all the degree of freedom at all entities.
"""
function setdofvalues!(self::F, val::T) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    for i in 1:nent
        for j in 1:ndof
            setdofvalue!(self, i, j, val)
        end
    end
    return self
end

"""
    setdofvalues!(self::F, ent, vals) where {T, F<:AbstractField{T}}

Set the value of all the degree of freedom of a given entity.
"""
function setdofvalues!(self::F, ent, vals) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    for j in 1:ndof
        setdofvalue!(self, ent, j, vals[j])
    end
    return self
end

"""
    getdofvalue(self::F, ent, dof) where {T, F<:AbstractField{T}}

Get the value of the degree of freedom at a given entity.
"""
function getdofvalue(self::F, ent, dof) where {T, F<:AbstractField{T}}
    return self.data[ent][dof].dofvalue
end

"""
    getdofnum(self::F, ent, dof) where {T, F<:AbstractField{T}}

Get the degree of freedom number at a given entity and for a given degree of freedom.
"""
function getdofnum(self::F, ent, dof) where {T, F<:AbstractField{T}}
    return self.data[ent][dof].dofnum
end

"""
    getdofvalues!(self::F, vals, ent) where {T, F<:AbstractField{T}}

Get the values of all the degrees of freedom at a given entity.
"""
function getdofvalues!(self::F, vals, ent) where {T, F<:AbstractField{T}}
    nent, ndof = nents(self), ndofs(self)
    for j in 1:ndof
        vals[j] = getdofvalue(self, ent, j)
    end
    return vals
end

"""
    getdofsvector(self::F, ent) where {T, F<:AbstractField{T}}

Provide access to the vector for all the degrees of freedom for a given entity.
"""
function getdofsvector(self::F, ent) where {T, F<:AbstractField{T}}
    return self.data[ent]
end

"""
    getdofvalue(v::Vector{FieldDOFData{T}}, dof) where {T}

Get the value of a degree of freedom from the vector of dofs.
"""
function getdofvalue(v::Vector{FieldDOFData{T}}, dof) where {T}
    return v[dof].dofvalue
end

"""
    setdofnum!(self::F, ent, dof, val::Int) where {T, F<:AbstractField{T}}

Set the degree of freedom number at a given entity.
"""
function setdofnum!(self::F, ent, dof, val::Int) where {T, F<:AbstractField{T}}
    self.data[ent][dof].dofnum = val
    return self
end

"""
    getdofnum(v::Vector{FieldDOFData{T}}, dof) where {T}

Get the degree of freedom number from the vector of dofs.
"""
function getdofnum(v::Vector{FieldDOFData{T}}, dof) where {T}
    return v[dof].dofnum
end

end
