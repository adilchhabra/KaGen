#pragma once

#include "kagen/generators/generator.h"
#include "kagen/io/graph_format.h"

namespace kagen {
class FileGraphFactory : public GeneratorFactory {
public:
    std::unique_ptr<Generator> Create(const PGeneratorConfig& config, PEID rank, PEID size) const final;
};

class FileGraphGenerator : public Generator {
public:
    FileGraphGenerator(const PGeneratorConfig& config, const PEID rank, const PEID size);

protected:
    void GenerateEdgeList() final;

    void GenerateCSR() final;

    void FinalizeEdgeList(MPI_Comm comm) final;

    void FinalizeCSR(MPI_Comm comm) final;

private:
    void GenerateImpl(GraphRepresentation representation);

    bool CheckDeficit(ReaderDeficits deficit) const;

    bool RequiresPostprocessing() const;

    bool Output(MPI_Comm comm) const;

    const PGeneratorConfig& config_;

    PEID rank_;
    PEID size_;

    int deficits_;
};
} // namespace kagen
