module mTEST13H_vib
using FinEtools
using FinEtools.AlgoDeforLinearModule: ssit
using DataFrames
using CSV
using Base.Test
function test()
    # Harmonic forced vibration problem is solved for a homogeneous square plate,
    # simply-supported on the circumference.
    # This is the TEST 13H from the Abaqus v 6.12 Benchmarks manual.
    # The test is recommended by the National Agency for Finite Element Methods and Standards (U.K.):
    # Test 13 from NAFEMS “Selected Benchmarks for Forced Vibration,” R0016, March 1993.
    #
    #
    # The plate is discretized with hexahedral solid elements. The simple support
    # condition is approximated by distributed rollers on the boundary.
    # Because only the out of plane displacements are prevented, the structure
    # has three rigid body modes in the plane of the plate.
    #
    #
    # The nonzero benchmark frequencies are (in hertz): 2.377, 5.961, 5.961,
    # 9.483, 12.133, 12.133, 15.468, 15.468 [Hz].

    println("""
    Homogeneous square plate, simply-supported on the circumference from
    the test 13 from NAFEMS “Selected Benchmarks for Forced Vibration,” R0016, March 1993.
    The nonzero benchmark frequencies are (in hertz): 2.377, 5.961, 5.961,
    9.483, 12.133, 12.133, 15.468, 15.468 [Hz].
    """)

    # t0 = time()

    E = 200*phun("GPa");# Young's modulus
    nu = 0.3;# Poisson ratio
    rho = 8000*phun("KG*M^-3");# mass density
    qmagn = 100.0*phun("Pa")
    L = 10.0*phun("M"); # side of the square plate
    t = 0.05*phun("M"); # thickness of the square plate
    nL = 16; nt = 4;
    tolerance = t/nt/100;
    # neigvs = 11;
    # OmegaShift = (2*pi*0.5) ^ 2; # to resolve rigid body modes
    frequencies = vcat(linspace(0,2.377,20), linspace(2.377,15,70))

    # Compute the parameters of Rayleigh damping. For the two selected
    # frequencies we have the relationship between the damping ratio and
    # the Rayleigh parameters
    # $\xi_m=a_0/\omega_m+a_1\omega_m$
    # where $m=1,2$.  Solving for the Rayleigh parameters $a_0,a_1$ yields:
    zeta1= 0.02; zeta2  =0.02;
    o1 =2*pi*2.377;  o2 =2*pi*15.468;
    Rayleigh_mass = 2*(o1*o2)/(o2^2-o1^2)*(o2*zeta1-o1*zeta2);# a0
    Rayleigh_stiffness = 2*(o1*o2)/(o2^2-o1^2)*(-1/o2*zeta1+1/o1*zeta2);# a1

    Rayleigh_mass = Rayleigh_mass;
    Rayleigh_stiffness = Rayleigh_stiffness;

    MR = DeforModelRed3D
    fens,fes  = H8block(L, L, t, nL, nL, nt)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(FCplxFlt, size(fens.xyz,1), 3)) # displacement field
    nl = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, 3)
    nl = selectnode(fens, box=[L L -Inf Inf -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, 3)
    nl = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, 3)
    nl = selectnode(fens, box=[-Inf Inf L L -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, 3)
    applyebc!(u)
    numberdofs!(u)
    println("nfreedofs = $(u.nfreedofs)")

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearMSH8(MR, GeoD(fes, GaussRule(3,2)), material)
    femm = associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(3,3)), material)
    M = mass(femm, geom, u)
    C = Rayleigh_mass*M + Rayleigh_stiffness*K

    # if true
    #     t0 = time()
    #     d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    #     d = d - OmegaShift;
    #     fs = real(sqrt.(complex(d)))/(2*pi)
    #     println("Reference Eigenvalues: $fs [Hz]")
    #     println("eigs solution ($(time() - t0) sec)")
    # end

    bdryfes = meshboundary(fes)
    topbfl = selectelem(fens, bdryfes, facing=true, direction=[0.0 0.0 1.0])
    el1femm =  FEMMBase(GeoD(subset(bdryfes,topbfl), GaussRule(2,2)))
    function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
        forceout .=  [0.0, 0.0, -qmagn]
        return forceout
    end
    fi = ForceIntensity(FFlt, 3, pfun);
    F = distribloads(el1femm, geom, u, fi, 2);

    U1 = zeros(FCplxFlt, u.nfreedofs, length(frequencies))
    for k = 1:length(frequencies)
        frequency = frequencies[k];
        omega = 2*pi*frequency;
        U1[:, k] = (-omega^2*M + 1im*omega*C + K)\F;
    end

    midpoint = selectnode(fens, box=[L/2 L/2 L/2 L/2 0 0], inflate=tolerance);
    midpointdof = u.dofnums[midpoint, 3]
    umidAmpl = abs(U1[midpointdof, :])/phun("MM")
    df = DataFrame(Frequency=vec(frequencies), umidAmpl=vec(umidAmpl))
    File = "mTEST13H_vib.CSV"
    CSV.write(File, df)
    @async run(`"paraview.exe" $File`)
    true
end
end
using mTEST13H_vib
mTEST13H_vib.test()
