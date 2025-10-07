test_args = gprd_args
test_deps = gprd_deps
# Handle GTest dependencies gracefully
gtest = dependency('gtest_main', required: false)
if gtest.found()
    test_deps += gtest
else
    gtest_proj = subproject('gtest')
    test_deps += gtest_proj.get_variable('gtest_dep')
    test_deps += gtest_proj.get_variable('gmock_dep')
    test_deps += gtest_proj.get_variable('gtest_main')
endif

test_array = [  #
    ['test_field', files('data_types/FieldTest.cpp')],
    ['test_distance', files('gpr/auxiliary/DistanceTest.cpp')],
    ['test_gradient', files('gpr/auxiliary/GradientTest.cpp')],
    ['test_problemsetup', 'gpr/auxiliary/ProblemSetUpTest.cpp'],
    ['test_sexpatcf', files('gpr/covariance_functions/SexpatCFTest.cpp')],
    ['test_dimer', files('gpr/dimer/DimerTest.cpp')],
    ['test_dimerrotrem', files('gpr/dimer/DimerRotRemTest.cpp')],
    ['test_gpr', files('gpr/ml/GaussianProcessRegressionTest.cpp')],
    ['test_scg', files('gpr/ml/SCGTest.cpp')],
    ['test_adam', files('gpr/ml/ADAMTest.cpp')],
    ['test_likgauss', files('gpr/observation_models/LikGaussianTest.cpp')],
    ['test_eampot', files('gpr/potentials/EAMPotentialTest.cpp')],
    ['test_prior_gauss', files('gpr/prior/PriorGaussianTest.cpp')],
    ['test_prior_logunif', files('gpr/prior/PriorLogUnifTest.cpp')],
    ['test_prior_sqrtt', files('gpr/prior/PriorSqrttTest.cpp')],
    ['test_prior_tt', files('gpr/prior/PriorTTest.cpp')],
    ['test_atomic_dimer', files('gpr/AtomicDimerTest.cpp')],
    ['test_managers_io', files('managers/io/FileManagerTest.cpp')],
]

foreach test : test_array
    test(
        test.get(0),
        executable(
            test.get(0),
            sources: [test.get(1)],
            dependencies: test_deps,
            include_directories: gprd_incdirs,
            link_with: libgprd,
            cpp_args: test_args,
        ),
        workdir: meson.project_source_root(),
    )
endforeach


if get_option('use_capnp')
    test(
        'test_capnp_params',
        executable(
            'test_capnp_params',
            sources: [files('managers/io/CapnpParamsTest.cpp')],
            dependencies: [test_deps, capnp_rpc_dep],
            include_directories: gprd_incdirs,
            link_with: [libgprd, capnp_fileloader],
            cpp_args: test_args,
        ),
        workdir: meson.project_source_root(),
    )
endif
