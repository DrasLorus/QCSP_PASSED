# QCSP PASSED: Platform Agnostic Standalone Software Efficient Detector

QCSP detector surpassed all expectancies. Then it passed all the tests.

## Usage

Not much currently. However, classes can already be used ! The CDetector can be used with single precision floating point to compute scores and detect frames for arbitrary PN sequences, frequency span et frequency resolution (as a symbol rotation of pi / step_denominator).

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
pushd build && ctest -j 4 && popd

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
ctest -j 4 && popd

```
