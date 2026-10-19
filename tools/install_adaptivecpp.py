"""Install the pinned AdaptiveCpp OpenMP SDK used by extension tests.

Requires git, CMake, Ninja, a C++ compiler, Boost and NUMA development headers.
This profile tests generated SYCL on CPU; it does not provide a GPU compiler.
"""
from argparse import ArgumentParser
from pathlib import Path
import subprocess
import tempfile

REVISION = '31db2e83903c1db9f50b03b620eb0a09224cd0eb'


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    prefix = args.prefix.resolve()
    with tempfile.TemporaryDirectory(prefix='symox-adaptivecpp-') as directory:
        source, build = Path(directory)/'source', Path(directory)/'build'
        commands = [
            ['git', 'init', str(source)],
            ['git', '-C', str(source), 'fetch', '--depth', '1',
             'https://github.com/AdaptiveCpp/AdaptiveCpp.git', REVISION],
            ['git', '-C', str(source), 'checkout', '--detach', 'FETCH_HEAD'],
            ['cmake', '-S', str(source), '-B', str(build), '-GNinja',
             f'-DCMAKE_INSTALL_PREFIX={prefix}', '-DCMAKE_BUILD_TYPE=Release',
             '-DACPP_COMPILER_FEATURE_PROFILE=none', '-DACPP_DEFAULT_TARGETS=omp',
             '-DWITH_CUDA_BACKEND=OFF', '-DWITH_ROCM_BACKEND=OFF',
             '-DWITH_OPENCL_BACKEND=OFF', '-DWITH_LEVEL_ZERO_BACKEND=OFF',
             '-DWITH_VULKAN_BACKEND=OFF', '-DCMAKE_POLICY_VERSION_MINIMUM=3.5'],
            ['cmake', '--build', str(build), '--target', 'install', '--parallel', str(args.jobs)],
        ]
        for command in commands:
            subprocess.run(command, check=True)
    print(f'AdaptiveCpp {REVISION} installed in {prefix}; use --acpp-targets=omp.')


if __name__ == '__main__':
    main()
