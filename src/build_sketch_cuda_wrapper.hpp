#pragma once
#include <vector>
#include <string>

// sketch_type: 0 = hll, 1 = smh
void launch_sketch_files_cuda(const std::vector<std::string> &files, uint k, int sketch_type, std::vector<void*> &sketches);
