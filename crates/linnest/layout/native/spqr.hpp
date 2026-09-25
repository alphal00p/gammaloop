#pragma once

#include <nlohmann/json.hpp>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/decomposition/StaticPlanarSPQRTree.h>
#include <ogdf/planarity/BoyerMyrvold.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace ec {
using Json = nlohmann::ordered_json;
using json = Json;

class Decomposition {
  ogdf::Graph graph;
  ogdf::NodeArray<std::string> node_ids{graph};
  ogdf::EdgeArray<std::string> edge_ids{graph};
  std::unordered_map<std::string, ogdf::node> nodes;
  std::unordered_map<std::string, ogdf::edge> edges;

  static std::string skeleton_id(ogdf::node component, ogdf::edge edge) {
    return std::to_string(component->index()) + ":" +
           std::to_string(edge->index());
  }

public:
  explicit Decomposition(const json &input) {
    for (const auto &value : input.at("nodes")) {
      auto id = value.get<std::string>();
      if (nodes.count(id)) {
        throw std::invalid_argument("duplicate node ID: " + id);
      }
      const auto node = graph.newNode();
      node_ids[node] = id;
      nodes.emplace(id, node);
    }
    for (const auto &value : input.at("edges")) {
      auto id = value.at("id").get<std::string>();
      const auto source = value.at("source").get<std::string>();
      const auto target = value.at("target").get<std::string>();
      if (edges.count(id)) {
        throw std::invalid_argument("duplicate edge ID: " + id);
      }
      if (!nodes.count(source) || !nodes.count(target)) {
        throw std::invalid_argument("unknown endpoint for edge: " + id);
      }
      if (source == target) {
        throw std::invalid_argument(
            "self-loop must be separated before SPQR: " + id);
      }
      const auto edge = graph.newEdge(nodes.at(source), nodes.at(target));
      edge_ids[edge] = id;
      edges.emplace(id, edge);
    }
  }

  json run(const json &input) {
    ogdf::BoyerMyrvold planarity;
    const bool planar = planarity.planarEmbed(graph);
    const bool biconnected = ogdf::isBiconnected(graph);
    json output = {{"planar", planar},
                   {"biconnected", biconnected},
                   {"root", nullptr},
                   {"components", json::array()},
                   {"rotation_cw", json::object()}};
    if (!planar || !biconnected) {
      return output;
    }
    if (graph.numberOfNodes() < 2 || graph.numberOfEdges() < 3) {
      throw std::invalid_argument(
          "SPQR requires a block with at least three edges");
    }
    ogdf::edge root_edge = graph.firstEdge();
    if (input.contains("root_edge")) {
      const auto id = input.at("root_edge").get<std::string>();
      if (!edges.count(id)) {
        throw std::invalid_argument("unknown root edge: " + id);
      }
      root_edge = edges.at(id);
    }
    ogdf::StaticPlanarSPQRTree tree(graph, root_edge, true);
    output["root"] = std::to_string(tree.rootNode()->index());
    for (const auto node : graph.nodes) {
      std::vector<std::string> rotation;
      for (const auto adj : node->adjEntries) {
        rotation.push_back(edge_ids[adj->theEdge()]);
      }
      std::reverse(rotation.begin(), rotation.end());
      output["rotation_cw"][node_ids[node]] = rotation;
    }
    for (const auto component : tree.tree().nodes) {
      const auto &skeleton = tree.skeleton(component);
      const auto type = tree.typeOf(component);
      json item = {
          {"id", std::to_string(component->index())},
          {"type", type == ogdf::SPQRTree::NodeType::SNode   ? "S"
                   : type == ogdf::SPQRTree::NodeType::PNode ? "P"
                                                             : "R"},
          {"reference_edge", skeleton_id(component, skeleton.referenceEdge())},
          {"nodes", json::array()},
          {"edges", json::array()},
          {"rotation_cw", json::object()}};
      for (const auto node : skeleton.getGraph().nodes) {
        const auto &id = node_ids[skeleton.original(node)];
        item["nodes"].push_back(id);
        std::vector<std::string> rotation;
        for (const auto adj : node->adjEntries) {
          rotation.push_back(skeleton_id(component, adj->theEdge()));
        }
        std::reverse(rotation.begin(), rotation.end());
        item["rotation_cw"][id] = rotation;
      }
      for (const auto edge : skeleton.getGraph().edges) {
        json record = {{"id", skeleton_id(component, edge)},
                       {"source", node_ids[skeleton.original(edge->source())]},
                       {"target", node_ids[skeleton.original(edge->target())]},
                       {"real_edge", nullptr},
                       {"twin", nullptr}};
        if (skeleton.isVirtual(edge)) {
          const auto other = skeleton.twinTreeNode(edge);
          record["twin"] = {
              {"component", std::to_string(other->index())},
              {"edge", skeleton_id(other, skeleton.twinEdge(edge))}};
        } else {
          record["real_edge"] = edge_ids[skeleton.realEdge(edge)];
        }
        item["edges"].push_back(record);
      }
      output["components"].push_back(item);
    }
    return output;
  }
};

inline Json decompose(const Json &input) {
  Decomposition decomposition(input);
  return decomposition.run(input);
}
} // namespace ec
