#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "kagen/context.h"
#include "kagen/io/seq_graph_writer.h"
#include "kagen/tools/statistics.h"

namespace kagen {
SequentialGraphWriter::SequentialGraphWriter(const OutputGraphConfig& config, Graph& graph, MPI_Comm comm)
    : GraphWriter(config, graph, comm) {}

void SequentialGraphWriter::Write(const bool report_progress) {
    PEID rank, size;
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &size);
    const bool output = report_progress && rank == ROOT;

    const std::string filename = config_.distributed ? config_.filename + "." + std::to_string(rank) : config_.filename;

    const bool requires_sorted_edges      = Requirements() & Requirement::SORTED_EDGES;
    const bool requires_coordinates       = Requirements() & Requirement::COORDINATES;
    const bool requires_coordinates2d     = Requirements() & Requirement::COORDINATES_2D;
    const bool requires_coordinates3d     = Requirements() & Requirement::COORDINATES_3D;
    const bool supports_no_vertex_weights = Requirement() & Requirement::NO_VERTEX_WEIGHTS;
    const bool supports_no_edge_weights   = Requirement() & Requirement::NO_EDGE_WEIGHTS;
    const bool has_coordinates2d          = coordinates_.first.size() == vertex_range_.second - vertex_range_.first;
    const bool has_coordinates3d          = coordinates_.second.size() == vertex_range_.second - vertex_range_.first;

    // Check if edges have to be sorted
    auto cmp_from = [](const auto& lhs, const auto& rhs) {
        return std::get<0>(lhs) < std::get<0>(rhs);
    };

    if (requires_sorted_edges && !std::is_sorted(edges_.begin(), edges_.end(), cmp_from)) {
        if (!edge_weights_.empty()) {
            std::vector<SSInt> indices(edges_.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](const auto& lhs, const auto& rhs) {
                return cmp_from(edges_[lhs], edges_[rhs]);
            });
            for (std::size_t e = 0; e < edges_.size(); ++e) {
                indices[e] = edge_weights_[indices[e]];
            }
            std::swap(edge_weights_, indices);
        }

        std::sort(edges_.begin(), edges_.end(), cmp_from);
    }

    // Check if coordinates are required
    const bool local_coordinates_ok = (!requires_coordinates2d || has_coordinates2d)
                                      && (!requires_coordinates3d || has_coordinates3d)
                                      && (!requires_coordinates || has_coordinates2d || has_coordinates3d);
    bool global_coordinates_ok;
    MPI_Allreduce(&local_coordinates_ok, &global_coordinates_ok, 1, MPI_CXX_BOOL, MPI_LAND, comm_);
    if (!global_coordinates_ok) {
        if (rank == ROOT) {
            std::cerr << "Output format requires coordinates, but the graph was generated without coordinates\n";
            std::cerr << "This may happen because coordinates were disabled or the graph format does not support "
                         "coordinates\n";
        }
        std::exit(1);
    }

    // Warn if we have vertex or edge weights and the output format does not support them
    if (supports_no_vertex_weights && HasVertexWeights() && output) {
        std::cout
            << "Warning: graph was generated with vertex weights, but the output format does not support vertex weights"
            << std::endl;
    }
    if (supports_no_edge_weights && HasEdgeWeights() && output) {
        std::cout
            << "Warning: graph was generated with edge weights, but the output format does not support edge weights"
            << std::endl;
    }

    // Everything OK, write graph
    CreateFile(filename);

    if (!config_.distributed) {
        const SInt n = FindNumberOfGlobalNodes(vertex_range_, comm_);
        const SInt m = FindNumberOfGlobalEdges(edges_, comm_);

        if (output) {
            std::cout << "Writing graph to " << filename << " ..." << std::endl;
        }

        if (rank == ROOT && config_.header != OutputHeader::NEVER) {
            AppendHeaderTo(filename, n, m);
        }

        for (PEID pe = 0; pe < size; ++pe) {
            if (output) {
                std::cout << "  Writing subgraph of PE " << pe << " ... " << std::flush;
            }
            if (rank == pe) {
                AppendTo(filename);
            }
            MPI_Barrier(comm_);
            if (output) {
                std::cout << "OK" << std::endl;
            }
        }

        if (rank == ROOT && config_.header != OutputHeader::NEVER) {
            AppendFooterTo(filename);
        }
    } else {
        const SInt n = vertex_range_.second - vertex_range_.first;
        const SInt m = edges_.size();
        const bool write_header_footer =
            (rank == ROOT && config_.header == OutputHeader::ROOT) || config_.header == OutputHeader::ALWAYS;

        if (output) {
            std::cout << "Writing graph to [" << config_.filename << ".0";
            if (size > 2) {
                std::cout << ", ...";
            }
            if (size > 1) {
                std::cout << ", " << config_.filename << "." << size - 1;
            }
            std::cout << "] ..." << std::endl;
        }

        if (write_header_footer) {
            AppendHeaderTo(filename, n, m);
        }
        AppendTo(filename);
        if (write_header_footer) {
            AppendFooterTo(filename);
        }
    }
}

void SequentialGraphWriter::AppendFooterTo(const std::string&) {}

int SequentialGraphWriter::Requirements() const {
    return Requirement::NONE;
}

void SequentialGraphWriter::CreateFile(const std::string& filename) {
    std::ofstream ofs(filename);
}
} // namespace kagen