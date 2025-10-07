capnp_rpc_dep = dependency('capnp-rpc')
_mbproot = meson.project_build_root()
capnpc_ = find_program('capnpc')
gen_capnp = custom_target(
    'gen_capnp',
    input: ['gprd_params.capnp'],
    output: ['gprd_params.capnp.c++', 'gprd_params.capnp.h'],
    # TODO(rg): Can this be cleaner?
    command: [
        capnpc_,
        f'-oc++:@_mbproot@/schema',
        '@INPUT@',
        '--verbose',
        '--src-prefix=../schema',
    ],
)
param_capnp = library(
    'param_capnp',
    gen_capnp,
    dependencies: gprd_deps + capnp_rpc_dep,
    cpp_args: gprd_args,
    include_directories: gprd_incdirs,
    install: true,
)
