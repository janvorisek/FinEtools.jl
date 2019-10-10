module momap2para13
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cross
function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

    Meshing = L2blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linearspace(0.0, A, nA + 1))
    fensc,fesc = Meshing(xs)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i = 1:count(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x = centroid
        valuesasarray(fc)[i, :] = sin.(2*x/A)
    end
    File = "momap2para12-coarse.vtk"
    MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", valuesasarray(fc))])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linearspace(0.0, A, nA + 1))
    fensf,fesf = Meshing(xs)
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i = 1:count(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x = centroid
        valuesasarray(referenceff)[i, :] = sin.(2*x/A)
    end
    File = "momap2para12-reference.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", valuesasarray(referenceff))])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para12-fine.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", valuesasarray(ff))])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    diffff = ElementalField(valuesasarray(referenceff) - valuesasarray(ff))
    femm  = FEMMBase(IntegDomain(fesf, GaussRule(1, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error/ref -  0.23860709149331033) < 1.0e-4
end
end
using .momap2para13
momap2para13.test()
