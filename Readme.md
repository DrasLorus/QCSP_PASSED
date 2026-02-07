# QCSP PASSED: Platform Agnostic Standalone Software Efficient Detector

QCSP detector surpassed all expectancies. Then it passed all the tests.

## Usage

Not much currently. However, classes can already be used ! The `CDetector` can be used with
single and double precision floating point to compute scores and detect frames for arbitrary
PN sequences, frequency span et frequency resolution (as a symbol rotation of π / step\_denominator).

## Compilation

It has been tested (and validated) on Archlinux, Debian 12 (Bullseye), Ubuntu 24.04
(and Fedora 40, VoidLinux, Gentoo, and so on, but unsupported).

You will need basics compilation tools (resp. `build-essential` on Debian and Ubuntu, `base-devel`
on Arch, and groups `"Development Tools" "Development Libraries"` on Fedora) and CMake. You will
also need Boost program-options (resp. `libboost-program-options-dev`, `boost`, `boost-devel`),
MatIO (resp. `libmatio-dev`, `libmatio`, `matio-devel`) and Catch2 (`catch2`), the latter can also
be automatically fetched by CMake if `ENABLE_SYSTEM_CATCH2` is set to `OFF`.

### Example (On archlinux)

```bash
sudo pacman -Syu base-devel cmake boost boost-libs libmatio catch2 python git --needed

git clone <url.git> QCSP_PASSED
cd QCSP_PASSED

cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j "$(nproc)"

```

### Example with automatic fetch (On Ubuntu 22.04)

```bash
sudo apt update && sudo apt upgrade
sudo apt install build-essential cmake libboost-program-options-dev libmatio-dev python3 git

git clone <url.git> QCSP_PASSED
cd QCSP_PASSED

mkdir build && pushd build
cmake .. -DCMAKE_BUILD_TYPE:STRING=Release -DENABLE_SYSTEM_CATCH2:BOOL=OFF
cmake --build build -j "$(nproc)" 

```

### Windows special case

Windows 10 (24H1) has also been tested and validated, using

- MinGW64 GCC and [MSYS2](https://www.msys2.org/),
- MSVC thanks to [vcpkg](https://vcpkg.io/en/index.html).

It requires to clone the repository with its submodule (using `--recursive`) or run `git submodule init && git submodule update` for an already cloned repo.
Then, `vcpkg` must be bootstrapped. A target `vcpkg` has been conveniently added to the CMakeLists.txt inside the [`windows` directory](./windows).
Go under the directory, run CMake on the local CMakeLists.txt, and build the target `vcpkg`.
Then everything should work the same, provided you give CMake the toolchain of vcpkg (`CMAKE_TOOLCHAIN_FILE=<repo root>/windows/vcpkg/scripts/buildsystems/vcpkg.cmake`)

In summary, in any case, in PowerShell (hit `win + x`, then select PowerShell), do:

```powershell
git clone --recursive <url.git> QCSP_PASSED
cd QCSP_PASSED
```

---

Then, to use MSVC with `vcpkg`, from the QCSP\_PASSED directory, do:

```powershell
cd windows

cmake -B build -S .
cmake --build build -t vcpkg

cd ..

.\windows\vcpkg\vcpkg.exe install matio[mat73]:x64-windows boost-program-options:x64-windows

cmake -B build -S . -DCMAKE_BUILD_TYPE:STRING=Release -DENABLE_SYSTEM_CATCH2:BOOL=OFF -DCMAKE_TOOLCHAIN_FILE=$(pwd)/windows/vcpkg/scripts/buildsystems/vcpkg.cmake

cmake --build build -j 4
```

If it fails, first ensure that you have git and CMake installed ([winget](https://github.com/microsoft/winget-cli)
is recommended for the matter) and that you have MSVC (if not, run `winget install Microsoft.VisualStudio.2022.BuildTools`, then install C/C++ development utilities).
Then, if vcpkg is successfully bootstrapped, but dependencies aren't added automatically,
[check your configuration](https://vcpkg.io/en/docs/users/buildsystems/cmake-integration.html), or [install dependencies manually](https://vcpkg.io/en/docs/examples/installing-and-using-packages.html).

---

To use MSYS2, from the QCSP\_PASSED directory, do:

```powershell
cmake -B build -S . -DCMAKE_BUILD_TYPE:STRING=Release -DENABLE_SYSTEM_CATCH2:BOOL=OFF -G Ninja -DCMAKE_CXX_COMPILER=<path/to/msys-root>/ucrt64/bin/g++ -DCMAKE_C_COMPILER=<path/to/msys-root>/ucrt64/bin/gcc

cmake --build build -j 4
```

If it fails, yet again, first ensure that you have git and CMake installed ([winget](https://github.com/microsoft/winget-cli)
is recommended for the matter) and that you have MSYS2 (if not, run `winget install msys2.msys2`, then follow the [installation procedure](https://www.msys2.org/), but be sure to use UCRT64).
Then, inside the MSYS2 UCRT64 shell, run

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-toolchain mingw-w64-ucrt-x86_64-ninja mingw-w64-ucrt-x86_64-boost mingw-w64-ucrt-x86_64-matio
```

It should work afterward.

Windows 11 should work the same.

## Testing

Testing can be done when two conditions are met:

1. Catch v2.13 or higher is available,
2. Test vectors are available.

The first should be already satisfied (either by installation or by auto-fetching). The second is
partially satisfied since one test vector is hosted in this repository (`test_data_w1_nofreq.mat`
in the [data directory](./data)), with corresponding parameters (`parameters_20210903.mat`). The
others can be generated using the MATLAB/Octave files (suffix `.m`) in the [data directory](./data).

The functions they define can be used within MATLAB (tested with R2021a) or
Octave (tested with v7.1.0). It results in slightly different test vectors for each, due to a
different randomness, to some implementation details (especially [here](./data/save_test_vectors.m))
and to the fact that Octave (at the time of writing) does not support MATLAB binary format v7.3
(while [MatIO](https://github.com/tbeu/matio) does).

However, they are valid and can be used to validate the build (using `ctest`, see lower), to
implement new tests or to explore insights of the algorithm without tinkering with C++.

### Test Vectors Generation

Two test vectors must be generated to perform the tests.

- One for a p\_omega of 34, a maximum rotation error of π / 8 (step\_denominator = 8) and a rotation
  span of 2π.
- One for a p\_omega of 4, a maximum rotation error of π / 4 (step\_denominator = 4) and a rotation
  span of π / 2.

A variety of others are needed depending on the LIST\_TARGET\_Q CMake list content.

Those can be automatically generated if the option `ENABLE\_DATA\_AUTOGENERATION` (available when testing is enabled) is activated, and a suitable mat\_interpreter has been found. If it is not the case, they can be generated by calling the auto-generated script `generator\_script.m` in the [data directory](./data) from within [MATLAB](https://fr.mathworks.com/) or [Octave](https://octave.org/).

#### Additional Test Vectors

Functions can be used to generate any test vector. They can also be looked at to generate
other values for more runs, different SNRs, and so on and so forth. Functions are self documented
in a way that should be sufficient to understand their usage. This means that running
`help <function>` should provide enough information to allow understanding and use of the function.

For example

```matlab
help save_test_vectors
```

results in

```text
  SAVE_TEST_VECTORS Generate a set of test vectors and saves it in a matfile.
    You must provide a PN sequence, a overmodulation sequence `best_N`, the rotation span (symbol
    rotation is in [-rotation_span/2, rotation_span/2[), the number of frequency hypothesis, the
    step_denominator (maximum of error is pi/step_denominator) and a seed for the randomness. Be
    aware that files grow with the length of PN and best_N, and with p_omega.

    Usage:
        fn = save_test_vectors(PN, best_N, rotation_span, p_omega, step_denominator, seed)

    Example:
        fn = save_test_vectors(PN, best_N, 2 * pi, 4, 4, 0) produce a file
        "test_data_w4_step4_span0.5.mat" that contains:
               cabs_max_sqr_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]
               cabs_max_sqr_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
              cabs_max_sqr_raw_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]
              cabs_max_sqr_raw_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
              cabs_max_sqrt_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]
              cabs_max_sqrt_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
             cabs_max_sqrt_raw_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]
             cabs_max_sqrt_raw_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
                    data_input_infdB_w4_q64_N60_1pi_2_n10: [192000x1 single]
                    data_input_m10dB_w4_q64_N60_1pi_2_n30: [576000x1 single]
                        deltas_infdB_w4_q64_N60_1pi_2_n10: [10x1     single]
                        deltas_m10dB_w4_q64_N60_1pi_2_n30: [30x1     single]
              iter_fcts_sqr_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]
              iter_fcts_sqr_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
             iter_fcts_sqr_raw_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]
             iter_fcts_sqr_raw_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
             iter_fcts_sqrt_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]
             iter_fcts_sqrt_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
            iter_fcts_sqrt_raw_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]
            iter_fcts_sqrt_raw_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
                  norms_sqr_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x1 single]
                  norms_sqr_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x1 single]
                 norms_sqrt_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x1 single]
                 norms_sqrt_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x1 single]
                     rotations_infdB_w4_q64_N60_1pi_2_n10: [10x1     single]
                     rotations_m10dB_w4_q64_N60_1pi_2_n30: [30x1     single]
                  score_sqr_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]
                  score_sqr_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
                 score_sqr_raw_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]
                 score_sqr_raw_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
                 score_sqrt_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]
                 score_sqrt_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
                score_sqrt_raw_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]
                score_sqrt_raw_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
        and that can be loaded either using load or matfile (if supported).

    See also LOAD, MATFILE
```

### Validating the build

The validation process use `ctest` (that comes bundled with `cmake`). Just run
`ctest -j <number of jobs>` in the build directory after the generation of the test vectors, and
the 32 tests should run and pass (the majority takes less than 1s to run, only those related to
p\_omega = 34 take up to 7s on a quad-core i5 8th gen).

