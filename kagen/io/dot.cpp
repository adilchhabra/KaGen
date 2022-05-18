#include "kagen/io/dot.h"

#include "kagen/io/buffered_writer.h"

namespace kagen {
DotWriter::DotWriter(EdgeList& edges, const VertexRange vertex_range, Coordinates& coordinates, MPI_Comm comm)
    : SequentialGraphWriter(edges, vertex_range, coordinates, comm) {}

std::string DotWriter::DefaultExtension() const {
    return "dot";
}

void DotWriter::AppendHeaderTo(const std::string& filename, SInt, SInt) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteString("graph G{\n").Flush();
}

void DotWriter::AppendTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);

    if (!coordinates_.first.empty()) {
        auto& coordinates = coordinates_.first; // 2D
        for (SInt node = vertex_range_.first; node < vertex_range_.second; ++node) {
            const auto& [x, y] = coordinates[node];
            out.WriteInt(node + 1)
                .WriteString("[pos=\"")
                .WriteFloat(x * 10)
                .WriteChar(',')
                .WriteFloat(y * 10)
                .WriteString("!\"]\n")
                .Flush();
        }
    }

    for (const auto& [from, to]: edges_) {
        if (from < to) { // need edges only once
            out.WriteInt(from + 1).WriteString("--").WriteInt(to + 1).WriteChar('\n').Flush();
        }
    }
}

void DotWriter::AppendFooterTo(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);
    out.WriteString("}\n").Flush();
}
} // namespace kagen