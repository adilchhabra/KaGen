#include "kagen/kagen.h"

#include <mpi.h>

namespace kagen {
void GenerateStreamed(PGeneratorConfig config, MPI_Comm comm);

void GenerateVertexStreamed(PGeneratorConfig config, MPI_Comm comm);
}
