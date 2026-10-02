#pragma once
// Gutwenger–Klein–Mutzel EC expansion and two-sided SPQR insertion.
// Edge-identity rotations are clockwise. Container order matches the reference
// implementation: ties in shortest paths use the original incidence order.
#include <algorithm>
#include <array>
#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <nlohmann/json.hpp>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ec {
using Json = nlohmann::ordered_json;
using Id = std::string;
using Ids = std::vector<Id>;
using Set = std::set<Id>;
using Ends = std::array<Id, 2>;
using Dart = std::pair<Id, Id>;
struct NoPath : std::runtime_error {
  using std::runtime_error::runtime_error;
};
inline void require(bool valid, const std::string &message) {
  if (!valid)
    throw std::runtime_error(message);
}
// Python dictionaries preserve insertion order, including for repeated edge
// insertions. A sorted map here would silently change the certified seed.
template <class T> struct OrderedMap {
  std::vector<std::pair<Id, T>> values;
  auto begin() { return values.begin(); }
  auto end() { return values.end(); }
  auto begin() const { return values.begin(); }
  auto end() const { return values.end(); }
  auto find(const Id &id) {
    return std::find_if(begin(), end(),
                        [&](const auto &p) { return p.first == id; });
  }
  auto find(const Id &id) const {
    return std::find_if(begin(), end(),
                        [&](const auto &p) { return p.first == id; });
  }
  bool contains(const Id &id) const { return find(id) != end(); }
  size_t size() const { return values.size(); }
  bool empty() const { return values.empty(); }
  T &operator[](const Id &id) {
    auto it = find(id);
    if (it == end()) {
      values.emplace_back(id, T{});
      return values.back().second;
    }
    return it->second;
  }
  const T &at(const Id &id) const {
    auto it = find(id);
    require(it != end(), "Unknown identity: " + id);
    return it->second;
  }
  T &at(const Id &id) {
    auto it = find(id);
    require(it != end(), "Unknown identity: " + id);
    return it->second;
  }
  void erase(const Id &id) {
    auto it = find(id);
    if (it != end())
      values.erase(it);
  }
  Set keys() const {
    Set s;
    for (const auto &[k, v] : values)
      s.insert(k);
    return s;
  }
  void update(const OrderedMap &other) {
    for (const auto &[k, v] : other)
      (*this)[k] = v;
  }
};
inline Set intersection(const Set &a, const Set &b) {
  Set out;
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                        std::inserter(out, out.end()));
  return out;
}
inline Set difference(const Set &a, const Set &b) {
  Set out;
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));
  return out;
}
inline Set as_set(const Ids &v) { return {v.begin(), v.end()}; }
inline size_t index_of(const Ids &v, const Id &k) {
  auto it = std::find(v.begin(), v.end(), k);
  require(it != v.end(), "Missing incidence " + k);
  return it - v.begin();
}
inline bool cyclic_equal(const Ids &a, const Ids &b) {
  if (a.size() != b.size())
    return false;
  if (a.empty())
    return true;
  auto it = std::find(b.begin(), b.end(), a[0]);
  if (it == b.end())
    return false;
  size_t j = it - b.begin();
  for (size_t i = 0; i < a.size(); ++i)
    if (a[i] != b[(i + j) % b.size()])
      return false;
  return true;
}
inline Ids after(const Ids &v, const Id &k) {
  size_t i = index_of(v, k);
  Ids r(v.begin() + i + 1, v.end());
  r.insert(r.end(), v.begin(), v.begin() + i);
  return r;
}
inline void replace(Ids &v, const Id &old, const Ids &replacement) {
  size_t i = index_of(v, old);
  v.erase(v.begin() + i);
  v.insert(v.begin() + i, replacement.begin(), replacement.end());
}

struct Graph {
  Set nodes, protected_edges, wheel_hubs;
  OrderedMap<Ends> edges;
  std::map<Id, Ids> rotation, hubs;
  Graph() = default;
  Graph(Set ns, OrderedMap<Ends> es)
      : nodes(std::move(ns)), edges(std::move(es)) {
    for (const auto &n : nodes)
      rotation[n] = {};
    for (const auto &[e, p] : edges) {
      require(nodes.count(p[0]) && nodes.count(p[1]) && p[0] != p[1],
              "Invalid graph endpoints");
      rotation[p[0]].push_back(e);
      rotation[p[1]].push_back(e);
    }
  }
  Id other(const Id &e, const Id &n) const {
    const auto &p = edges.at(e);
    if (p[0] == n)
      return p[1];
    require(p[1] == n, "Vertex not incident to edge");
    return p[0];
  }
  Graph subgraph(const Set &chosen, bool keep_nodes = false) const {
    Graph out;
    if (keep_nodes)
      out.nodes = nodes;
    for (const auto &[e, p] : edges)
      if (chosen.count(e)) {
        out.edges[e] = p;
        out.nodes.insert(p.begin(), p.end());
      }
    for (const auto &n : out.nodes) {
      Ids r;
      for (const auto &e : rotation.at(n))
        if (chosen.count(e))
          r.push_back(e);
      out.rotation[n] = r;
    }
    out.protected_edges = intersection(protected_edges, chosen);
    for (const auto &[n, r] : hubs)
      if (out.nodes.count(n))
        out.hubs[n] = r;
    out.wheel_hubs = intersection(wheel_hubs, out.nodes);
    return out;
  }
  struct Faces {
    std::vector<std::vector<Dart>> boundaries;
    std::map<Dart, size_t> face_of;
  };
  Faces faces() const {
    std::map<Dart, Id> previous;
    std::map<Id, Set> incidence;
    for (const auto &[e, p] : edges)
      for (const auto &n : p)
        incidence[n].insert(e);
    for (const auto &n : nodes) {
      const auto &r = rotation.at(n);
      require(as_set(r) == incidence[n] && r.size() == incidence[n].size(),
              "Invalid rotation at " + n);
      for (size_t i = 0; i < r.size(); ++i)
        previous[{r[i], n}] = r[(i + r.size() - 1) % r.size()];
    }
    Faces out;
    for (const auto &[e, p] : edges)
      for (const auto &n : p) {
        Dart start{e, n}, d = start;
        if (out.face_of.count(start))
          continue;
        std::vector<Dart> boundary;
        while (!out.face_of.count(d)) {
          out.face_of[d] = out.boundaries.size();
          boundary.push_back(d);
          Id target = other(d.first, d.second);
          d = {previous.at({d.first, target}), target};
        }
        require(d == start, "Invalid face permutation");
        out.boundaries.push_back(boundary);
      }
    return out;
  }
  bool valid_embedding() const {
    auto fs = faces();
    Set seen;
    std::map<Id, Ids> adj;
    for (const auto &[e, p] : edges) {
      adj[p[0]].push_back(p[1]);
      adj[p[1]].push_back(p[0]);
    }
    for (const auto &root : nodes) {
      if (seen.count(root))
        continue;
      Set component;
      Ids pending{root};
      while (!pending.empty()) {
        Id n = pending.back();
        pending.pop_back();
        if (!component.insert(n).second)
          continue;
        auto &v = adj[n];
        pending.insert(pending.end(), v.begin(), v.end());
      }
      seen.insert(component.begin(), component.end());
      long ne = 0, nf = 0;
      for (const auto &[e, p] : edges)
        ne += component.count(p[0]);
      for (const auto &f : fs.boundaries)
        nf += !f.empty() && component.count(f[0].second);
      if (ne && long(component.size()) - ne + nf != 2)
        return false;
    }
    for (const auto &[hub, expected] : hubs)
      if (!cyclic_equal(rotation.at(hub), expected))
        return false;
    for (const auto &f : fs.boundaries) {
      bool wheel = false, protected_face = true;
      for (const auto &[e, n] : f) {
        wheel |= wheel_hubs.count(n);
        protected_face &= protected_edges.count(e);
      }
      if (wheel && (f.size() != 3 || !protected_face))
        return false;
    }
    return true;
  }
  void validate() const {
    require(valid_embedding(), "Invalid constrained embedding");
  }
};
// Implemented by the pinned OGDF primitive, shared by native and WASM callers.
Json decompose(const Json &request);
inline Json spqr(const Graph &g) {
  Json req = {{"nodes", g.nodes}, {"edges", Json::array()}};
  for (const auto &[e, p] : g.edges)
    req["edges"].push_back({{"id", e}, {"source", p[0]}, {"target", p[1]}});
  auto out = decompose(req);
  if (out.contains("error"))
    throw std::runtime_error(out["error"].get<Id>());
  if (!out.at("planar").get<bool>())
    return out;
  require(out.value("biconnected", true), "SPQR requires a biconnected block");
  return out;
}
inline Id skeleton_id(const Id &component, const Id &edge,
                      const Set &reserved) {
  Id id = "@spqr:" + component + ":" + edge;
  while (reserved.count(id))
    id = "@" + id;
  return id;
}
inline std::optional<OrderedMap<Graph>> feasible(const Graph &g,
                                                 const Json &tree) {
  OrderedMap<Graph> out;
  for (const auto &c : tree.at("components")) {
    Id cid = c.at("id");
    std::map<Id, Id> names, real;
    for (const auto &e : c.at("edges")) {
      Id eid = e.at("id");
      names[eid] = e.at("real_edge").is_null()
                       ? skeleton_id(cid, eid, g.edges.keys())
                       : e.at("real_edge").get<Id>();
      if (!e.at("real_edge").is_null())
        real[e.at("real_edge")] = names[eid];
    }
    Graph s;
    s.nodes = c.at("nodes").get<Set>();
    for (const auto &e : c.at("edges")) {
      Ends p;
      if (e.at("real_edge").is_null()) {
        p = {e.at("source"), e.at("target")};
        std::sort(p.begin(), p.end());
      } else
        p = g.edges.at(e.at("real_edge"));
      s.edges[names.at(e.at("id"))] = p;
    }
    for (const auto &[n, r] : c.at("rotation_cw").items())
      for (const auto &e : r)
        s.rotation[n].push_back(names.at(e));
    for (const auto &e : g.protected_edges)
      if (real.count(e))
        s.protected_edges.insert(real.at(e));
    s.wheel_hubs = intersection(g.wheel_hubs, s.nodes);
    std::set<bool> orientation;
    for (const auto &[hub, expected] : g.hubs)
      if (s.nodes.count(hub)) {
        Ids mapped;
        for (const auto &e : expected)
          mapped.push_back(real.at(e));
        s.hubs[hub] = mapped;
        require(c.at("type") == "R", "Wheel hub is outside R skeleton");
        auto actual = s.rotation.at(hub);
        if (cyclic_equal(actual, mapped))
          orientation.insert(false);
        else {
          std::reverse(actual.begin(), actual.end());
          require(cyclic_equal(actual, mapped), "SPQR changed wheel order");
          orientation.insert(true);
        }
      }
    if (orientation.size() > 1)
      return std::nullopt;
    if (orientation == std::set<bool>{true})
      for (auto &[n, r] : s.rotation)
        std::reverse(r.begin(), r.end());
    out[cid] = s;
  }
  return out;
}
inline std::vector<Set> blocks(const Graph &g) {
  std::map<Id, std::vector<std::pair<Id, Id>>> adj;
  for (const auto &[e, p] : g.edges) {
    adj[p[0]].push_back({e, p[1]});
    adj[p[1]].push_back({e, p[0]});
  }
  std::map<Id, size_t> discovery, low;
  Ids stack;
  std::vector<Set> result;
  std::function<void(const Id &, const Id &)> visit = [&](const Id &u,
                                                          const Id &parent) {
    low[u] = discovery[u] = discovery.size();
    for (const auto &[e, v] : adj[u]) {
      if (e == parent)
        continue;
      if (!discovery.count(v)) {
        stack.push_back(e);
        visit(v, e);
        low[u] = std::min(low[u], low[v]);
        if (low[v] >= discovery[u]) {
          Set component;
          while (!stack.empty()) {
            Id current = stack.back();
            stack.pop_back();
            component.insert(current);
            if (current == e)
              break;
          }
          result.push_back(component);
        }
      } else if (discovery[v] < discovery[u]) {
        stack.push_back(e);
        low[u] = std::min(low[u], discovery[v]);
      }
    }
  };
  for (const auto &n : g.nodes)
    if (!discovery.count(n))
      visit(n, "");
  return result;
}
inline Graph splice(const Graph &parent, const Id &pe, const Graph &child,
                    const Id &ce) {
  Set poles(parent.edges.at(pe).begin(), parent.edges.at(pe).end());
  require(poles == Set(child.edges.at(ce).begin(), child.edges.at(ce).end()),
          "Virtual edges have different poles");
  require(difference(intersection(parent.edges.keys(), child.edges.keys()),
                     {pe, ce})
              .empty(),
          "SPQR edge identities overlap");
  require(intersection(parent.nodes, child.nodes) == poles,
          "SPQR expansions share non-pole vertices");
  Graph out = parent;
  out.nodes.insert(child.nodes.begin(), child.nodes.end());
  out.edges.erase(pe);
  for (const auto &[e, p] : child.edges)
    if (e != ce)
      out.edges[e] = p;
  for (const auto &n : child.nodes)
    if (poles.count(n))
      replace(out.rotation[n], pe, after(child.rotation.at(n), ce));
    else
      out.rotation[n] = child.rotation.at(n);
  out.protected_edges.erase(pe);
  for (const auto &e : child.protected_edges)
    if (e != ce)
      out.protected_edges.insert(e);
  out.hubs.insert(child.hubs.begin(), child.hubs.end());
  out.wheel_hubs.insert(child.wheel_hubs.begin(), child.wheel_hubs.end());
  return out;
}
struct Tree {
  OrderedMap<Json> raw;
  OrderedMap<Graph> skeletons;
  std::map<Id, OrderedMap<Ends>> adjacency;
  Id root;
  Tree(const Json &s, OrderedMap<Graph> graphs)
      : skeletons(std::move(graphs)), root(s.at("root")) {
    Set reserved;
    for (const auto &c : s.at("components")) {
      raw[c.at("id")] = c;
      for (const auto &e : c.at("edges"))
        if (!e.at("real_edge").is_null())
          reserved.insert(e.at("real_edge"));
    }
    for (const auto &[cid, c] : raw)
      for (const auto &e : c.at("edges"))
        if (!e.at("twin").is_null()) {
          auto t = e.at("twin");
          Id other = t.at("component");
          adjacency[cid][other] = {skeleton_id(cid, e.at("id"), reserved),
                                   skeleton_id(other, t.at("edge"), reserved)};
        }
  }
  Graph expand(const Id &start, const std::set<Set> &blocked = {},
               const Graph *replacement = nullptr) const {
    std::function<Graph(const Id &, const Id &)> visit = [&](const Id &n,
                                                             const Id &parent) {
      Graph out = (n == start && replacement) ? *replacement : skeletons.at(n);
      auto it = adjacency.find(n);
      if (it != adjacency.end())
        for (const auto &[next, p] : it->second)
          if (next != parent && !blocked.count({n, next}))
            out = splice(out, p[0], visit(next, n), p[1]);
      return out;
    };
    return visit(start, "");
  }
  Ids allocation_path(const Id &start, const Id &end) const {
    Ids starts;
    Set ends;
    for (const auto &[k, g] : skeletons) {
      if (g.nodes.count(start))
        starts.push_back(k);
      if (g.nodes.count(end))
        ends.insert(k);
    }
    std::sort(starts.begin(), starts.end());
    std::map<Id, std::optional<Id>> pred;
    std::deque<Id> queue;
    for (const auto &k : starts) {
      pred[k] = std::nullopt;
      queue.push_back(k);
    }
    while (!queue.empty()) {
      Id n = queue.front();
      queue.pop_front();
      if (ends.count(n)) {
        Ids path;
        std::optional<Id> cur = n;
        while (cur) {
          path.push_back(*cur);
          cur = pred.at(*cur);
        }
        std::reverse(path.begin(), path.end());
        return path;
      }
      auto it = adjacency.find(n);
      if (it != adjacency.end())
        for (const auto &[next, p] : it->second)
          if (!pred.count(next)) {
            pred[next] = n;
            queue.push_back(next);
          }
    }
    throw std::runtime_error("Missing SPQR endpoint");
  }
};
struct Dual {
  const Graph &graph;
  Graph::Faces faces;
  std::vector<std::vector<std::pair<size_t, Id>>> arcs;
  explicit Dual(const Graph &g)
      : graph(g), faces(g.faces()), arcs(faces.boundaries.size()) {
    for (const auto &[e, p] : g.edges)
      if (!g.protected_edges.count(e)) {
        auto a = faces.face_of.at({e, p[0]}), b = faces.face_of.at({e, p[1]});
        if (a != b) {
          arcs[a].push_back({b, e});
          arcs[b].push_back({a, e});
        }
      }
  }
  std::set<size_t> at(const Id &n) const {
    std::set<size_t> out;
    for (const auto &e : graph.rotation.at(n))
      out.insert(faces.face_of.at({e, n}));
    return out;
  }
  std::set<size_t> side(const Id &e, size_t s) const {
    return {faces.face_of.at({e, graph.edges.at(e)[s]})};
  }
  std::optional<Ids> shortest(const std::set<size_t> &sources,
                              const std::set<size_t> &targets) const {
    std::map<size_t, std::optional<std::pair<size_t, Id>>> pred;
    std::deque<size_t> queue;
    for (auto n : sources) {
      pred[n] = std::nullopt;
      queue.push_back(n);
    }
    while (!queue.empty()) {
      size_t n = queue.front();
      queue.pop_front();
      if (targets.count(n)) {
        Ids path;
        while (pred.at(n)) {
          auto p = *pred.at(n);
          path.push_back(p.second);
          n = p.first;
        }
        std::reverse(path.begin(), path.end());
        return path;
      }
      for (const auto &[next, e] : arcs[n])
        if (!pred.count(next)) {
          pred[next] = std::make_pair(n, e);
          queue.push_back(next);
        }
    }
    return std::nullopt;
  }
};
struct Witness {
  std::vector<Dart> darts;
  Id start_corner, end_corner;
};
inline Witness witness(const Graph &g, const Id &start, const Id &end,
                       const Ids &crossings) {
  Dual dual(g);
  auto targets = dual.at(end);
  for (size_t source : dual.at(start)) {
    size_t face = source;
    std::vector<Dart> darts;
    bool valid = true;
    for (const auto &e : crossings) {
      const auto &p = g.edges.at(e);
      auto a = dual.faces.face_of.at({e, p[0]}),
           b = dual.faces.face_of.at({e, p[1]});
      if (face == a) {
        darts.push_back({e, p[0]});
        face = b;
      } else if (face == b) {
        darts.push_back({e, p[1]});
        face = a;
      } else {
        valid = false;
        break;
      }
    }
    if (valid && targets.count(face)) {
      Id a, b;
      for (const auto &[e, n] : dual.faces.boundaries[source])
        if (n == start) {
          a = e;
          break;
        }
      for (const auto &[e, n] : dual.faces.boundaries[face])
        if (n == end) {
          b = e;
          break;
        }
      return {darts, a, b};
    }
  }
  throw std::runtime_error("Embedding does not realize insertion path");
}
inline Id free_corner(const Graph &g, const Id &n) {
  auto fs = g.faces();
  for (const auto &e : g.rotation.at(n)) {
    bool wheel = false;
    for (const auto &[edge, v] : fs.boundaries[fs.face_of.at({e, n})])
      wheel |= g.wheel_hubs.count(v);
    if (!wheel)
      return e;
  }
  throw std::runtime_error("Cannot attach block inside wheel");
}
inline Graph join_vertex(const Graph &parent, const Graph &child, const Id &n,
                         const Id &pc, const Id &cc) {
  require(intersection(parent.nodes, child.nodes) == Set{n},
          "Blocks must meet at one cut vertex");
  Graph out = parent;
  out.nodes.insert(child.nodes.begin(), child.nodes.end());
  out.edges.update(child.edges);
  for (const auto &[v, r] : child.rotation)
    if (v != n)
      out.rotation[v] = r;
  Ids order = after(child.rotation.at(n), cc);
  order.push_back(cc);
  auto &r = out.rotation[n];
  r.insert(r.begin() + index_of(r, pc) + 1, order.begin(), order.end());
  out.protected_edges.insert(child.protected_edges.begin(),
                             child.protected_edges.end());
  out.hubs.insert(child.hubs.begin(), child.hubs.end());
  out.wheel_hubs.insert(child.wheel_hubs.begin(), child.wheel_hubs.end());
  return out;
}
inline Graph unite(const Graph &a, const Graph &b) {
  require(intersection(a.nodes, b.nodes).empty(),
          "Disjoint union shares vertices");
  Graph out = a;
  out.nodes.insert(b.nodes.begin(), b.nodes.end());
  out.edges.update(b.edges);
  out.rotation.insert(b.rotation.begin(), b.rotation.end());
  out.hubs.insert(b.hubs.begin(), b.hubs.end());
  out.protected_edges.insert(b.protected_edges.begin(),
                             b.protected_edges.end());
  out.wheel_hubs.insert(b.wheel_hubs.begin(), b.wheel_hubs.end());
  return out;
}
inline Graph glue(const Graph &original, std::vector<Graph> pending) {
  Graph out;
  while (!pending.empty()) {
    bool attached = false;
    for (size_t i = 0; i < pending.size(); ++i) {
      auto shared = intersection(out.nodes, pending[i].nodes);
      if (shared.empty())
        continue;
      require(shared.size() == 1, "Blocks share multiple cut vertices");
      Id n = *shared.begin();
      out = join_vertex(out, pending[i], n, free_corner(out, n),
                        free_corner(pending[i], n));
      pending.erase(pending.begin() + i);
      attached = true;
      break;
    }
    if (!attached) {
      out = unite(out, pending.front());
      pending.erase(pending.begin());
    }
  }
  out.nodes.insert(original.nodes.begin(), original.nodes.end());
  for (const auto &n : original.nodes)
    out.rotation.try_emplace(n, Ids{});
  return out;
}
inline std::optional<Graph> try_embed(const Graph &g) {
  std::vector<Graph> embedded;
  for (const auto &es : blocks(g)) {
    Graph block = g.subgraph(es);
    if (block.edges.size() > 2) {
      Json s = spqr(block);
      if (!s.at("planar").get<bool>())
        return std::nullopt;
      auto skeletons = feasible(block, s);
      if (!skeletons)
        return std::nullopt;
      Tree tree(s, *skeletons);
      block = tree.expand(tree.root);
    }
    embedded.push_back(block);
  }
  Graph out = glue(g, embedded);
  if (!out.valid_embedding())
    return std::nullopt;
  return out;
}
inline Graph embed(const Graph &g) {
  auto out = try_embed(g);
  require(bool(out), "Expected an EC-planar graph");
  return *out;
}
struct Expansion {
  Graph original, graph;
  Json constraints;
  std::map<Id, Set> members;
  std::map<Id, std::map<Id, Id>> leaf_ports;
  Set reserved;
  size_t counter = 0;
  static Ids leaves(const Json &t) {
    if (t.is_string())
      return {t.get<Id>()};
    require(t.is_object() && t.size() == 2 && t.contains("kind") &&
                t.contains("children"),
            "Invalid constraint tree");
    Id kind = t.at("kind");
    require(kind == "group" || kind == "mirror" || kind == "oriented",
            "Unknown constraint kind");
    const auto &children = t.at("children");
    require(children.is_array() && children.size() >= 2,
            "Constraint node needs two children");
    Ids out;
    for (const auto &child : children) {
      auto ls = leaves(child);
      out.insert(out.end(), ls.begin(), ls.end());
    }
    return out;
  }
  Id fresh(const Id &prefix) {
    for (;;) {
      Id id = "__ec_" + prefix + "_" + std::to_string(++counter);
      if (reserved.insert(id).second)
        return id;
    }
  }
  Id node(const Id &owner, const Id &prefix) {
    Id n = fresh(prefix);
    graph.nodes.insert(n);
    graph.rotation[n] = {};
    members[owner].insert(n);
    return n;
  }
  Id edge(const Id &e, const Id &a, const Id &b, bool protect = false) {
    graph.edges[e] = {a, b};
    graph.rotation[a].push_back(e);
    graph.rotation[b].push_back(e);
    if (protect)
      graph.protected_edges.insert(e);
    return e;
  }
  void expand(const Id &owner, const Json &tree,
              std::optional<Id> parent = std::nullopt) {
    if (tree.is_string()) {
      leaf_ports[owner][tree.get<Id>()] = node(owner, "leaf");
      return;
    }
    const auto &children = tree.at("children");
    size_t degree = children.size() + bool(parent);
    Ids ports;
    if (tree.at("kind") == "group") {
      Id center = node(owner, "group");
      ports.assign(degree, center);
    } else {
      Id hub = node(owner, "hub");
      graph.wheel_hubs.insert(hub);
      Ids rim, spokes;
      for (size_t i = 0; i < 2 * degree; ++i)
        rim.push_back(node(owner, "rim"));
      for (size_t i = 0; i < rim.size(); ++i) {
        spokes.push_back(edge(fresh("spoke"), hub, rim[i], true));
        edge(fresh("rimedge"), rim[i], rim[(i + 1) % rim.size()], true);
      }
      if (tree.at("kind") == "oriented")
        graph.hubs[hub] = spokes;
      for (size_t i = 0; i < rim.size(); i += 2)
        ports.push_back(rim[i]);
    }
    if (parent)
      edge(fresh("treeedge"), *parent, ports.back(), true);
    for (size_t i = 0; i < children.size(); ++i)
      if (children[i].is_string())
        leaf_ports[owner][children[i].get<Id>()] = ports[i];
      else
        expand(owner, children[i], ports[i]);
  }
  Expansion(const Graph &g, Json c)
      : original(g), constraints(std::move(c)), reserved(g.nodes) {
    auto es = g.edges.keys();
    reserved.insert(es.begin(), es.end());
    for (const auto &n : g.nodes) {
      members[n] = {n};
      leaf_ports[n] = {};
      if (!constraints.contains(n)) {
        graph.nodes.insert(n);
        graph.rotation[n] = {};
      }
    }
    graph.hubs = g.hubs;
    graph.wheel_hubs = g.wheel_hubs;
    for (const auto &[n, t] : constraints.items()) {
      require(g.nodes.count(n) && !g.wheel_hubs.count(n),
              "Invalid constrained vertex");
      Set incidence;
      for (const auto &[e, p] : g.edges)
        if (p[0] == n || p[1] == n)
          incidence.insert(e);
      auto ls = leaves(t);
      require(ls.size() == as_set(ls).size() && as_set(ls) == incidence,
              "Constraint must contain all incidences once");
      members[n].clear();
      expand(n, t);
    }
    for (const auto &[e, p] : g.edges) {
      auto a = leaf_ports[p[0]].find(e), b = leaf_ports[p[1]].find(e);
      edge(e, a == leaf_ports[p[0]].end() ? p[0] : a->second,
           b == leaf_ports[p[1]].end() ? p[1] : b->second);
    }
    graph.protected_edges.insert(g.protected_edges.begin(),
                                 g.protected_edges.end());
  }
  static bool admissible(const Json &tree, const Ids &rotation) {
    auto expected = leaves(tree);
    if (rotation.size() != expected.size() ||
        as_set(rotation) != as_set(expected))
      return false;
    std::function<bool(const Json &, const Ids &)> linear =
        [&](const Json &t, const Ids &order) {
          if (t.is_string())
            return order == Ids{t.get<Id>()};
          const auto &children = t.at("children");
          std::map<Id, size_t> owner;
          for (size_t i = 0; i < children.size(); ++i)
            for (const auto &e : leaves(children[i]))
              owner[e] = i;
          std::vector<std::pair<size_t, Ids>> runs;
          for (const auto &e : order) {
            size_t i = owner.at(e);
            if (runs.empty() || runs.back().first != i)
              runs.push_back({i, {}});
            runs.back().second.push_back(e);
          }
          std::set<size_t> unique;
          for (const auto &r : runs)
            unique.insert(r.first);
          if (runs.size() != children.size() ||
              unique.size() != children.size())
            return false;
          bool forward = true, reverse = true;
          for (size_t i = 0; i < runs.size(); ++i) {
            forward &= runs[i].first == i;
            reverse &= runs[i].first == runs.size() - 1 - i;
          }
          if (t.at("kind") == "oriented" && !forward)
            return false;
          if (t.at("kind") == "mirror" && !forward && !reverse)
            return false;
          for (const auto &[i, r] : runs)
            if (!linear(children[i], r))
              return false;
          return true;
        };
    for (size_t i = 0; i < rotation.size(); ++i) {
      Ids r(rotation.begin() + i, rotation.end());
      r.insert(r.end(), rotation.begin(), rotation.begin() + i);
      if (linear(tree, r))
        return true;
    }
    return false;
  }
  Graph collapse(const Graph &embedded) const {
    std::map<Id, Id> owner;
    for (const auto &[n, ns] : members)
      for (const auto &m : ns)
        owner[m] = n;
    auto of = [&](const Id &n) {
      auto it = owner.find(n);
      return it == owner.end() ? n : it->second;
    };
    Set ns;
    for (const auto &n : embedded.nodes)
      ns.insert(of(n));
    OrderedMap<Ends> edges;
    for (const auto &[e, p] : embedded.edges)
      if (of(p[0]) != of(p[1]))
        edges[e] = {of(p[0]), of(p[1])};
    Graph out(ns, edges);
    out.protected_edges = intersection(embedded.protected_edges, edges.keys());
    out.hubs = original.hubs;
    out.wheel_hubs = original.wheel_hubs;
    for (const auto &n : embedded.nodes)
      if (!owner.count(n))
        out.rotation[n] = embedded.rotation.at(n);
    for (const auto &n : original.nodes) {
      const auto &ms = members.at(n);
      Set internal;
      for (const auto &[e, p] : embedded.edges)
        if (ms.count(p[0]) && ms.count(p[1]))
          internal.insert(e);
      Ids incident;
      for (const auto &[e, p] : edges)
        if (p[0] == n || p[1] == n)
          incident.push_back(e);
      if (incident.empty()) {
        out.rotation[n] = {};
        continue;
      }
      Id edge = incident.front();
      const auto &p = embedded.edges.at(edge);
      Id point = ms.count(p[0]) ? p[0] : p[1];
      std::set<Dart> visited;
      Ids boundary;
      while (visited.insert({edge, point}).second) {
        const auto &r = embedded.rotation.at(point);
        edge = r[(index_of(r, edge) + 1) % r.size()];
        if (internal.count(edge))
          point = embedded.other(edge, point);
        else
          boundary.push_back(edge);
      }
      require(boundary.size() == incident.size() &&
                  as_set(boundary) == as_set(incident),
              "Expansion boundary lost incidences");
      out.rotation[n] = boundary;
    }
    out.validate();
    for (const auto &[n, t] : constraints.items())
      if (as_set(out.rotation.at(n)) == as_set(leaves(t)))
        require(admissible(t, out.rotation.at(n)),
                "Collapsed constraint violated");
    return out;
  }
};
struct Insertion {
  Graph embedding;
  Ids crossings;
  Id start, end;
  Json states = Json::array();
};
inline Graph p_embedding(const Tree &tree, const Id &cid, const Id &previous,
                         const Id &following, size_t side) {
  Graph g = tree.skeletons.at(cid);
  Id in = tree.adjacency.at(cid).at(previous)[0],
     out = tree.adjacency.at(cid).at(following)[0];
  const auto poles = g.edges.at(out);
  Ids order{out}, middle;
  for (const auto &e : g.rotation.at(poles[0]))
    if (e != in && e != out)
      middle.push_back(e);
  if (side == 0)
    order.push_back(in);
  order.insert(order.end(), middle.begin(), middle.end());
  if (side == 1)
    order.push_back(in);
  g.rotation[poles[0]] = order;
  std::reverse(order.begin(), order.end());
  g.rotation[poles[1]] = order;
  return g;
}
inline Insertion block_insertion(const Graph &g, const Id &start,
                                 const Id &end) {
  if (g.edges.size() <= 2)
    return {g, {}, start, end};
  Json raw = spqr(g);
  require(raw.at("planar").get<bool>(), "Insertion requires planar graph");
  auto skeletons = feasible(g, raw);
  require(bool(skeletons), "Insertion constraints infeasible");
  Tree tree(raw, *skeletons);
  Ids path = tree.allocation_path(start, end);
  std::set<Set> boundaries;
  for (size_t i = 1; i < path.size(); ++i)
    boundaries.insert({path[i - 1], path[i]});
  const size_t inf = std::numeric_limits<size_t>::max() / 4;
  std::array<size_t, 2> costs{0, 0};
  struct Choice {
    size_t entering = 0;
    Ids crossed;
    Graph graph;
    bool mirrored = false;
  };
  std::vector<std::array<Choice, 2>> choices;
  Json history = Json::array();
  for (size_t i = 0; i < path.size(); ++i) {
    const Id &cid = path[i];
    std::optional<Id> previous =
                          i ? std::optional<Id>(path[i - 1]) : std::nullopt,
                      following = i + 1 < path.size()
                                      ? std::optional<Id>(path[i + 1])
                                      : std::nullopt;
    Id kind = tree.raw.at(cid).at("type");
    std::array<Choice, 2> local;
    if (kind == "P" && path.size() > 1) {
      require(previous && following, "Shortest path ends at P node");
      for (size_t side = 0; side < 2; ++side) {
        auto s = p_embedding(tree, cid, *previous, *following, side);
        local[side] = {side, {}, tree.expand(cid, boundaries, &s), false};
      }
    } else {
      std::array<size_t, 2> best{inf, inf};
      for (bool mirrored : {false, true}) {
        if (mirrored && (kind != "R" || !tree.skeletons.at(cid).hubs.empty()))
          continue;
        Graph skeleton = tree.skeletons.at(cid);
        if (mirrored)
          for (auto &[n, r] : skeleton.rotation)
            std::reverse(r.begin(), r.end());
        Graph part = tree.expand(cid, boundaries, &skeleton);
        std::optional<Id> in =
                              previous
                                  ? std::optional<Id>(
                                        tree.adjacency.at(cid).at(*previous)[0])
                                  : std::nullopt,
                          out =
                              following
                                  ? std::optional<Id>(tree.adjacency.at(cid).at(
                                        *following)[0])
                                  : std::nullopt;
        if (in)
          part.protected_edges.insert(*in);
        if (out)
          part.protected_edges.insert(*out);
        Dual dual(part);
        for (size_t entering = 0; entering < 2; ++entering) {
          auto sources = in ? dual.side(*in, 1 - entering) : dual.at(start);
          for (size_t leaving = 0; leaving < 2; ++leaving) {
            auto targets = out ? dual.side(*out, leaving) : dual.at(end);
            auto crossed = dual.shortest(sources, targets);
            if (crossed && costs[entering] != inf) {
              size_t cost = costs[entering] + crossed->size();
              if (cost < best[leaving]) {
                best[leaving] = cost;
                local[leaving] = {entering, *crossed, part, mirrored};
              }
            }
          }
        }
      }
      costs = best;
    }
    choices.push_back(local);
    Json cs = Json::array();
    for (auto c : costs)
      cs.push_back(c == inf ? Json(nullptr) : Json(c));
    history.push_back({{"component", cid}, {"type", kind}, {"costs", cs}});
  }
  if (std::min(costs[0], costs[1]) == inf)
    throw NoPath("No EC insertion avoids protected edges");
  require(costs[0] == costs[1], "Terminal insertion costs differ");
  std::vector<Choice> selected(path.size());
  size_t side = 0;
  for (size_t i = path.size(); i-- > 0;) {
    selected[i] = choices[i][side];
    history[i]["mirrored"] = selected[i].mirrored;
    history[i]["entering"] = selected[i].entering;
    history[i]["leaving"] = side;
    side = selected[i].entering;
  }
  Graph embedding = selected[0].graph;
  Ids crossings = selected[0].crossed;
  for (size_t i = 1; i < path.size(); ++i) {
    auto pair = tree.adjacency.at(path[i]).at(path[i - 1]);
    embedding = splice(selected[i].graph, pair[0], embedding, pair[1]);
    crossings.insert(crossings.end(), selected[i].crossed.begin(),
                     selected[i].crossed.end());
  }
  embedding.protected_edges =
      intersection(embedding.protected_edges, g.protected_edges);
  witness(embedding, start, end, crossings);
  return {embedding, crossings, start, end, history};
}
inline Insertion optimal_insertion(const Graph &g, const Id &start,
                                   const Id &end) {
  require(start != end && g.nodes.count(start) && g.nodes.count(end),
          "Insertion needs two distinct vertices");
  std::vector<Graph> subgraphs;
  std::map<Id, std::vector<size_t>> incidence;
  for (const auto &es : blocks(g)) {
    Graph b = g.subgraph(es);
    for (const auto &n : b.nodes)
      incidence[n].push_back(subgraphs.size());
    subgraphs.push_back(b);
  }
  using Key = std::pair<char, Id>;
  Key source{'v', start}, target{'v', end};
  std::deque<Key> queue{source};
  std::map<Key, std::optional<Key>> pred{{source, std::nullopt}};
  while (!queue.empty() && !pred.count(target)) {
    Key current = queue.front();
    queue.pop_front();
    std::vector<Key> adjacent;
    if (current.first == 'v') {
      for (auto i : incidence[current.second])
        adjacent.push_back({'b', std::to_string(i)});
    } else
      for (const auto &n : subgraphs[std::stoul(current.second)].nodes)
        adjacent.push_back({'v', n});
    for (const auto &n : adjacent)
      if (!pred.count(n)) {
        pred[n] = current;
        queue.push_back(n);
      }
  }
  require(pred.count(target), "Insertion endpoints are disconnected");
  std::vector<Key> sequence;
  std::optional<Key> cursor = target;
  while (cursor) {
    sequence.push_back(*cursor);
    cursor = pred.at(*cursor);
  }
  std::reverse(sequence.begin(), sequence.end());
  std::optional<Graph> embedding;
  Ids crossings;
  std::set<size_t> used;
  Json history = Json::array();
  for (size_t offset = 1; offset < sequence.size(); offset += 2) {
    size_t i = std::stoul(sequence[offset].second);
    Id x = sequence[offset - 1].second, y = sequence[offset + 1].second;
    auto step = block_insertion(subgraphs[i], x, y);
    if (!embedding)
      embedding = step.embedding;
    else {
      auto left = witness(*embedding, start, x, crossings),
           right = witness(step.embedding, x, y, step.crossings);
      embedding = join_vertex(*embedding, step.embedding, x, left.end_corner,
                              right.start_corner);
    }
    crossings.insert(crossings.end(), step.crossings.begin(),
                     step.crossings.end());
    used.insert(i);
    for (const auto &state : step.states)
      history.push_back(state);
    witness(*embedding, start, y, crossings);
  }
  std::set<size_t> pending;
  for (size_t i = 0; i < subgraphs.size(); ++i)
    if (!used.count(i))
      pending.insert(i);
  while (!pending.empty()) {
    bool progress = false;
    auto todo = pending;
    for (size_t i : todo) {
      const auto &current = subgraphs[i];
      auto shared = intersection(embedding->nodes, current.nodes);
      if (shared.empty())
        continue;
      require(shared.size() == 1, "Invalid block decomposition");
      Id n = *shared.begin();
      embedding =
          join_vertex(*embedding, current, n, free_corner(*embedding, n),
                      free_corner(current, n));
      pending.erase(i);
      progress = true;
    }
    if (!progress) {
      size_t i = *pending.begin();
      embedding = unite(*embedding, subgraphs[i]);
      pending.erase(i);
    }
  }
  embedding->nodes.insert(g.nodes.begin(), g.nodes.end());
  for (const auto &n : g.nodes)
    embedding->rotation.try_emplace(n, Ids{});
  witness(*embedding, start, end, crossings);
  return {*embedding, crossings, start, end, history};
}
struct Planarization {
  Graph embedding;
  OrderedMap<Ids> chains;
  OrderedMap<Ends> crossing_nodes;
  Json insertions = Json::array();
};
inline Id fresh_id(const Id &base, const Set &existing) {
  Id value = base;
  size_t serial = 0;
  while (existing.count(value))
    value = base + ":" + std::to_string(++serial);
  return value;
}
inline Planarization planarize_insertion(const Insertion &result,
                                         const Id &id) {
  const Graph &g = result.embedding;
  require(!g.edges.contains(id), "Inserted edge exists");
  require(result.crossings.size() == as_set(result.crossings).size(),
          "Path crosses edge twice");
  auto proof = witness(g, result.start, result.end, result.crossings);
  Planarization out;
  out.embedding = g;
  for (const auto &[e, p] : g.edges)
    out.chains[e] = {e};
  Ids nodes{result.start};
  std::vector<Ends> split;
  for (size_t i = 0; i < proof.darts.size(); ++i) {
    const auto &[e, origin] = proof.darts[i];
    require(!g.protected_edges.count(e), "Insertion crosses protected edge");
    Id n =
        fresh_id("@cross:" + id + ":" + std::to_string(i), out.embedding.nodes);
    out.embedding.nodes.insert(n);
    nodes.push_back(n);
    Ends p = out.embedding.edges.at(e);
    out.embedding.edges.erase(e);
    Id a = fresh_id(e + ":part:0", out.embedding.edges.keys());
    out.embedding.edges[a] = {p[0], n};
    Id b = fresh_id(e + ":part:1", out.embedding.edges.keys());
    out.embedding.edges[b] = {n, p[1]};
    replace(out.embedding.rotation[p[0]], e, {a});
    replace(out.embedding.rotation[p[1]], e, {b});
    out.chains[e] = {a, b};
    split.push_back(origin == p[0] ? Ends{a, b} : Ends{b, a});
    out.crossing_nodes[n] = {e, id};
  }
  nodes.push_back(result.end);
  Ids inserted;
  for (size_t i = 1; i < nodes.size(); ++i) {
    Id e = result.crossings.empty()
               ? id
               : fresh_id(id + ":part:" + std::to_string(i - 1),
                          out.embedding.edges.keys());
    out.embedding.edges[e] = {nodes[i - 1], nodes[i]};
    inserted.push_back(e);
  }
  for (size_t i = 1; i + 1 < nodes.size(); ++i)
    out.embedding.rotation[nodes[i]] = {split[i - 1][0], inserted[i],
                                        split[i - 1][1], inserted[i - 1]};
  std::array<std::tuple<Id, Id, Id>, 2> corners{
      {{result.start, proof.start_corner, inserted.front()},
       {result.end, proof.end_corner, inserted.back()}}};
  for (auto [v, c, e] : corners) {
    if (!out.embedding.edges.contains(c))
      for (const auto &part : out.chains.at(c)) {
        auto p = out.embedding.edges.at(part);
        if (p[0] == v || p[1] == v) {
          c = part;
          break;
        }
      }
    auto &r = out.embedding.rotation[v];
    r.insert(r.begin() + index_of(r, c) + 1, e);
  }
  out.chains[id] = inserted;
  return out;
}
inline Planarization planarize(const Graph &g) {
  if (auto initial = try_embed(g)) {
    Planarization result;
    result.embedding = *initial;
    for (const auto &[e, p] : g.edges)
      result.chains[e] = {e};
    return result;
  }
  Set chosen = g.protected_edges;
  Graph embedded = embed(g.subgraph(chosen, true));
  Ids removed;
  for (const auto &[e, p] : g.edges) {
    if (chosen.count(e))
      continue;
    Set candidate = chosen;
    candidate.insert(e);
    if (auto next = try_embed(g.subgraph(candidate, true))) {
      chosen = std::move(candidate);
      embedded = std::move(*next);
    } else {
      removed.push_back(e);
    }
  }
  Planarization out;
  out.embedding = embedded;
  for (const auto &e : chosen)
    out.chains[e] = {e};
  for (const auto &e : removed) {
    std::map<Id, Id> owners;
    for (const auto &[owner, chain] : out.chains)
      for (const auto &part : chain)
        owners[part] = owner;
    Json constraints = Json::object();
    for (const auto &[n, record] : out.crossing_nodes)
      constraints[n] = {{"kind", "mirror"},
                        {"children", out.embedding.rotation.at(n)}};
    Expansion preserved(out.embedding, constraints);
    Graph constrained = embed(preserved.graph);
    const auto &ends = g.edges.at(e);
    auto path = optimal_insertion(constrained, ends[0], ends[1]);
    auto step = planarize_insertion(path, e);
    out.embedding = preserved.collapse(step.embedding);
    for (const auto &[n, record] : step.crossing_nodes)
      out.crossing_nodes[n] = {owners.at(record[0]), e};
    for (auto &[owner, chain] : out.chains) {
      Ids next;
      for (const auto &part : chain) {
        const auto &replacement = step.chains.at(part);
        next.insert(next.end(), replacement.begin(), replacement.end());
      }
      chain = next;
    }
    out.chains[e] = step.chains.at(e);
    owners.clear();
    for (const auto &[owner, chain] : out.chains)
      for (const auto &part : chain)
        owners[part] = owner;
    for (const auto &[n, r] : out.crossing_nodes) {
      Ids incidence;
      for (const auto &part : out.embedding.rotation.at(n))
        incidence.push_back(owners.at(part));
      require(incidence.size() == 4 && incidence[0] == incidence[2] &&
                  incidence[1] == incidence[3] && incidence[0] != incidence[1],
              "Crossing alternation changed");
      require(as_set(incidence) == Set{r[0], r[1]},
              "Crossing provenance changed");
    }
    Ids crossed;
    for (const auto &[n, r] : out.crossing_nodes)
      if (step.crossing_nodes.contains(n))
        crossed.push_back(r[0]);
    out.insertions.push_back({{"edge", e},
                              {"crossed_edges", crossed},
                              {"crossings", path.crossings.size()},
                              {"states", path.states}});
  }
  require(out.chains.keys() == g.edges.keys(),
          "Planarization lost original edges");
  out.embedding.validate();
  return out;
}
} // namespace ec
