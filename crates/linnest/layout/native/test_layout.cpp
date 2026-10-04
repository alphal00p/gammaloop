#include "seed.hpp"
#include "spqr.hpp"
#include <fstream>
#include <iostream>

int main(int argc, char **argv) {
  try {
    ec::require(argc == 3, "Expected constraint and seed fixtures");
    std::ifstream constraints_file(argv[1]);
    const auto constraints = ec::Json::parse(constraints_file);
    for (const auto &fixture : constraints) {
      const auto &request = fixture.at("request");
      ec::OrderedMap<ec::Ends> edges;
      for (const auto &edge : request.at("edges"))
        edges[edge.at("id")] = {edge.at("source"), edge.at("target")};
      ec::Graph graph(request.at("nodes").get<ec::Set>(), edges);
      graph.rotation = request.at("rotation").get<std::map<ec::Id, ec::Ids>>();
      ec::Expansion expansion(graph, request.at("constraints"));
      auto embedded = ec::try_embed(expansion.graph);
      const auto &expected = fixture.at("expected");
      ec::require(bool(embedded) == expected.at("planar").get<bool>(),
                  "Native feasibility differs from exhaustive rotation oracle");
      if (!embedded)
        continue;
      auto collapsed = expansion.collapse(*embedded);
      ec::require(collapsed.rotation ==
                      expected.at("rotation").get<std::map<ec::Id, ec::Ids>>(),
                  "Native constrained rotation differs from reference");
      for (const auto &[edge, endpoints] : collapsed.edges)
        ec::require(endpoints == expected.at("edges").at(edge).get<ec::Ends>(),
                    "Constraint expansion changed physical incidence");
    }
    std::ifstream seed_file(argv[2]);
    const auto seeds = ec::Json::parse(seed_file);
    for (const auto &fixture : seeds) {
      const auto result = ec::initialize(fixture.at("diagram"));
      const auto &expected = fixture.at("expected");
      for (const auto &key :
           {"routes", "node_ids", "external_ids", "edge_endpoints"})
        ec::require(result.at(key) == expected.at(key),
                    fixture.at("name").get<std::string>() + ": " + key);
      ec::require(result.at("positions").size() ==
                      expected.at("positions").size(),
                  "Seed changed point identities");
      for (const auto &[point, coordinate] : expected.at("positions").items())
        for (size_t axis = 0; axis < 2; ++axis)
          ec::require(
              std::abs(result.at("positions").at(point).at(axis).get<double>() -
                       coordinate.at(axis).get<double>()) <= 1e-12,
              fixture.at("name").get<std::string>() + ": seed coordinate");
      ec::require(result.at("crossings").size() ==
                      expected.at("crossings").size(),
                  "Seed changed crossing count");
      for (size_t i = 0; i < expected.at("crossings").size(); ++i)
        for (const auto &key : {"id", "edges", "route_points"})
          ec::require(result.at("crossings")[i].at(key) ==
                          expected.at("crossings")[i].at(key),
                      "Seed changed crossing provenance");
      for (const auto &key : {"contacts", "overlaps", "degenerate"})
        ec::require(result.at("report").at("geometry").at(key).empty(),
                    "Seed has geometric degeneracy");
    }
    std::cout << constraints.size() << " constrained embeddings and "
              << seeds.size() << " seed fixtures passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
