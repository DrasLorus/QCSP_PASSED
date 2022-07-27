# QCSP PASSED: Platform Agnostic Standalone Software Efficient Detector

QCSP detector surpassed all expectancies. Then it passed all the tests.

## Usage

Not much currently. However, classes can already be used ! The CDetector can be used with single precision floating point to compute scores and detect frames for arbitrary PN sequences, frequency span et frequency resolution (as a symbol rotation of π / step_denominator).

## Compilation

It has been tested (and validated) on ArchLinux, Debian 11 (Bullseye), Ubuntu 20.04, Ubuntu 22.04 and Fedora 36 (and VoidLinux but unsupported).

You will need basics compilation tools (resp. `build-essential` on Debian and Ubuntu, `base-devel` on Arch, and groups `"Development Tools" "Development Libraries"` on Fedora) and CMake. You will also need Boost program-options (resp. `libboost-program-options-dev`, `boost`, `boost-devel`), MatIO (resp. `libmatio-dev`, `libmatio`, `matio-devel`) and Catch2 (`catch2`), the latter can also be automatically fetched by CMake if `ENABLE_SYSTEM_CATCH2` is set to `OFF`.

### Exemple (On Archlinux)

```bash
sudo pacman -Syu base-devel cmake boost boost-libs libmatio catch2 git --needed

git clone <url.git> QCSP_PASSED
cd QCSP_PASSED

cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4

```

### Exemple with automatic fetch (On Ubuntu 20.04)

```bash
sudo apt update
sudo apt upgrade -Syu build-essential cmake libboost-program-options-dev libmatio-dev git -y

git clone <url.git> QCSP_PASSED
cd QCSP_PASSED

# CMake version on Ubuntu 20.04 does not support `-B` flag
mkdir build && pushd build
cmake . -DCMAKE_BUILD_TYPE:STRING=Release -DENABLE_SYSTEM_CATCH2:BOOL=OFF
cmake --build build -j 4

```

## Testing

Testing can be done when two conditions are met:

1. Catch v2.13 or higher is available,
2. Test vectors are available.

The first should be already satisfied (either by installation or by auto-fetching). The second is partially satisfied since one test vector is hosted in this repository (`test_data_w1_nofreq.mat` in the [data directory](./data)), with corresponding parameters (`parameters_20210903.mat`). The others can be generated using the MATLAB/Octave files (suffix `.m`) in the [data directory](./data).

The functions they define can be used within MATLAB (tested with R2021a) or Octave (tested with v7.1.0), and results in slightly different test vectors for both, due to a different randomness, to some implementation details (especially [here](./data/save_test_vectors.m)) and to the fact that Octave (at the time of writing) does not support MATLAB binary format v7.3 (while [MatIO](https://github.com/tbeu/matio) does).

However, they are valid and can be used to validate the build (using `ctest`, see lower), to implement new tests or to explore insights of the algorithm without tinkering with C++.

### Test Vectors Generation

Two test vectors must be generated to perform the tests.

- One for a p_omega of 34, a maximum rotation error of π / 8 (step_denominator = 8) and a rotation span of 2π.
- One for a p_omega of 4, a maximum rotation error of π/4 (step_denominator = 4) and a rotation span of π / 2.

#### With MATLAB

Assuming that the `matlab` executable is on your path (or that an alias has been set), you must run

```bash
pushd data
matlab -nosplash -nodesktop -batch 'load parameters_20210903.mat;
                                    display("Parameters loaded.");
                                    fn_34_8 = save_test_vectors(PN64, best_N, pi * 2, 34, 8, 0);
                                    fprintf("File %s written\n", fn_34_8);
                                    fn_4_4 = save_test_vectors(PN64, best_N, pi * 2, 4, 4, 0);
                                    fprintf("File %s written\n", fn_4_4);'
popd
```

#### With Octave

Assuming that the `octave` executable is on your path (or that an alias has been set), you must run

```bash
pushd data
octave --no-gui --eval 'load parameters_20210903.mat;
                        display("Parameters loaded.");
                        fn_34_8 = save_test_vectors(PN64, best_N, pi * 2, 34, 8, 0);
                        fprintf("File %s written\n", fn_34_8);
                        fn_4_4 = save_test_vectors(PN64, best_N, pi * 2, 4, 4, 0);
                        fprintf("File %s written\n", fn_4_4);'
popd
```

#### Additional Test Vectors

The last function can be used to generated any test vector. It can also be looked at to generate other values for more runs, different SNRs, and so on and so forth. Functions are self documented in a way that should be sufficient to understand their usage. This means that running `help <function>` should provide enough information to allow understanding and use of the function.

For example

```Octave
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

The validation process use `ctest` (that comes bundled with `cmake`). Just run `ctest -j <number of jobs>` in the build directory after the generation of the test vectors, and the 32 tests should run and pass (the majority takes less than 1s to run, only those related to p_omega = 34 take up to 7s on a quad-core i5 8th gen).
