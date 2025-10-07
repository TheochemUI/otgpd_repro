project(
    'gprd',
    'cpp',
    version: '0.1',
    default_options: ['warning_level=3', 'cpp_std=c++17'],
)

gprd_deps = []
gprd_args = []
gprd_incdirs = [include_directories('.')]

# Handle Eigen dependencies gracefully
eigen = dependency(
    'Eigen3',
    method: 'cmake',
    modules: ['Eigen3::Eigen'],
    required: false,
)
if eigen.found()
    eigen_dep = eigen
else
    eigen_proj = subproject(eigen)
    eigen_dep = eigen_proj.get_variable('eigen_dep')
endif

gprd_deps += eigen_dep

if get_option('use_highs')
    gprd_args += ['-DUSE_HIGHS=TRUE']
    highs_proj = subproject('highs')
    gprd_deps += highs_proj.get_variable('highs_dep')
    gprd_incdirs += highs_proj.get_variable('_incdirs')
endif

ceres_dep = dependency('Ceres', required: false, not_found_message: 'falling back to subproject')

if not ceres_dep.found()
    ceres_opts = cmake.subproject_options()
    ceres_opts.add_cmake_defines(
        {
            'BUILD_BENCHMARKS': false,
            'BUILD_EXAMPLES': false,
            'BUILD_TESTING': false,
            'MINIGLOG': true,
            'USE_CUDA': false,
            'GFLAGS': false,
        },
    )
    ceres_opts.set_install(false)
    ceres = cmake.subproject('ceres-solver', options: ceres_opts)
    ceres_dep = ceres.dependency('ceres')
endif

gprd_deps += ceres_dep

if get_option('use_hdf5')
    gprd_args += ['-DWITH_HDF5']
    gprd_deps += dependency(
        'HighFive',
        method: 'cmake',
        modules: ['HighFive_HighFive'],
        required: true,
    )
endif


if get_option('use_mkl')
    message('Be warned, MKL has not been tested extensively!')
    mkldep = dependency('mkl-dynamic-ilp64-iomp', required: true)
    add_project_arguments('-DEIGEN_USE_MKL_ALL', language: 'cpp')
    add_project_link_arguments(
        '-lmkl_intel_lp64',
        '-lmkl_core',
        '-lmkl_sequential',
        language: 'cpp',
    )
endif

libgprd = library(
    'gprdimer',
    'backend/DistributionFunctions.cpp',
    'gpr/auxiliary/Distance.cpp',
    'gpr/auxiliary/ProblemSetUp.cpp',
    'gpr/auxiliary/Gradient.cpp',
    'gpr/covariance_functions/ConstantCF.cpp',
    'gpr/covariance_functions/SexpatCF.cpp',
    'gpr/dimer/Dimer.cpp',
    'gpr/optimization/LBFGS.cpp',
    'gpr/ml/GaussianProcessRegression.cpp',
    'gpr/ml/SCG.cpp',
    'gpr/ml/ADAM.cpp',
    'gpr/ml/Ceres.cpp',
    'gpr/observation_models/LikGaussian.cpp',
    'gpr/prior/PriorGaussian.cpp',
    'gpr/prior/PriorLogUnif.cpp',
    'gpr/prior/PriorSqrtt.cpp',
    'gpr/prior/PriorT.cpp',
    'gpr/prior/PriorUnif.cpp',
    'gpr/potentials/EAMPotential.cpp',
    'gpr/AtomicDimer.cpp',
    'gpr/AtomicDimerInit.cpp',
    'managers/io/FileManager.cpp',
    dependencies: gprd_deps,
    cpp_args: gprd_args,
    install: get_option('do_install'),
)

libgprd_dep = declare_dependency(
    link_with: libgprd,
    include_directories: gprd_incdirs,
)


if get_option('use_capnp')
    gprd_args += ['-DUSE_CAPNP=TRUE']
    subdir('schema')
    capnp_fileloader = library(
        'capnp_fileloader',
        'managers/io/CapnProtoManager.cpp',
        dependencies: [gprd_deps, capnp_rpc_dep, libgprd_dep],
        cpp_args: gprd_args,
        link_with: param_capnp,
        install: get_option('do_install'),
    )
endif

if (not meson.is_subproject()) or get_option('with_tests')
    subdir('tests')

    add_languages('fortran', required: true)
    fc = meson.get_compiler('fortran')

    if get_option('use_capnp')
        executable(
            'gprd',
            'main.cpp',
            dependencies: [test_deps, gprd_deps, libgprd_dep, capnp_rpc_dep],
            link_with: capnp_fileloader,
            cpp_args: gprd_args,
            install: true,
        )
    else
        executable(
            'gprd',
            'main.cpp',
            dependencies: [test_deps, gprd_deps, libgprd_dep],
            cpp_args: gprd_args,
            install: true,
        )
    endif

    # For the tests
    executable(
        'mEAMCUH2',
        'gpr/eam_potential_fortran/main.f90',
        'gpr/eam_potential_fortran/eamroutines.f90',
    )
    # Needed to run the executable
    configure_file(input: 'pot.par', output: 'pot.par', copy: true)
    configure_file(input: 'tmp.con', output: 'tmp.con', copy: true)
    run_command(
        'bash',
        '-c',
        'cp -r $MESON_SOURCE_ROOT/input $MESON_BUILD_ROOT/input',
    )

endif  # if subproject
