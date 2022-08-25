#include "./CGenerator.hpp"

using namespace QCSP::StandaloneDetector::Simulators;

namespace nw = QCSP::StandaloneDetector::NbldpcWrapper;

void CGenerator::nbldpc_encode(const vector<int32_t> & message, vector<int32_t> & codeword) {
    nw::Encoding(code, table, message, codeword);
}
