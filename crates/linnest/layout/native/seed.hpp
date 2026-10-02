#pragma once
#include "coordinates.hpp"
namespace ec {
inline Json geometry(const Positions &positions, const Routes &routes,
                     double tolerance = 1e-8) {
  struct Segment {
    Id edge;
    size_t index;
    Id a, b;
  };
  std::vector<Segment> segments;
  for (const auto &[e, path] : routes)
    for (size_t i = 1; i < path.size(); ++i)
      segments.push_back({e, i - 1, path[i - 1], path[i]});
  Json crossings = Json::array(), contacts = Json::array(),
       overlaps = Json::array(), degenerate = Json::array();
  auto orient = [](Point a, Point b, Point c) {
    return cross(subtract(b, a), subtract(c, a));
  };
  auto on = [&](Point a, Point b, Point p) {
    return std::abs(orient(a, b, p)) <= tolerance &&
           p[0] >= std::min(a[0], b[0]) - tolerance &&
           p[0] <= std::max(a[0], b[0]) + tolerance &&
           p[1] >= std::min(a[1], b[1]) - tolerance &&
           p[1] <= std::max(a[1], b[1]) + tolerance;
  };
  for (size_t i = 0; i < segments.size(); ++i) {
    const auto &first = segments[i];
    Point a = positions.at(first.a), b = positions.at(first.b);
    if (distance(a, b) <= tolerance)
      degenerate.push_back({first.edge, first.index});
    for (size_t j = i + 1; j < segments.size(); ++j) {
      const auto &second = segments[j];
      Point c = positions.at(second.a), d = positions.at(second.b);
      Json pair = Json::array({Json::array({first.edge, first.index}),
                               Json::array({second.edge, second.index})});
      std::array<double, 4> values{orient(a, b, c), orient(a, b, d),
                                   orient(c, d, a), orient(c, d, b)};
      bool collinear = true;
      for (double v : values)
        collinear &= std::abs(v) <= tolerance;
      if (collinear) {
        size_t axis = std::abs(a[0] - b[0]) >= std::abs(a[1] - b[1]) ? 0 : 1;
        double extent =
            std::min(std::max(a[axis], b[axis]), std::max(c[axis], d[axis])) -
            std::max(std::min(a[axis], b[axis]), std::min(c[axis], d[axis]));
        if (extent > tolerance) {
          overlaps.push_back(pair);
          continue;
        }
      }
      std::array<int, 4> signs;
      for (size_t k = 0; k < 4; ++k)
        signs[k] = std::abs(values[k]) <= tolerance ? 0
                   : values[k] > 0                  ? 1
                                                    : -1;
      if (signs[0] * signs[1] < 0 && signs[2] * signs[3] < 0)
        crossings.push_back(pair);
      else if (first.a != second.a && first.a != second.b &&
               first.b != second.a && first.b != second.b &&
               (on(a, b, c) || on(a, b, d) || on(c, d, a) || on(c, d, b)))
        contacts.push_back(pair);
    }
  }
  return {{"valid", crossings.empty() && contacts.empty() && overlaps.empty() &&
                        degenerate.empty()},
          {"crossings", crossings},
          {"contacts", contacts},
          {"overlaps", overlaps},
          {"degenerate", degenerate},
          {"segments", segments.size()}};
}
inline Json open_crossings(Positions &positions, Routes &routes,
                           const Set &nodes) {
  auto original_positions = positions;
  auto original_routes = routes;
  std::set<Dart> segments;
  for (const auto &[e, r] : routes)
    for (size_t i = 1; i < r.size(); ++i)
      segments.insert({r[i - 1], r[i]});
  std::map<std::pair<Id, size_t>, Ids> replacements;
  Json records = Json::array();
  for (const auto &n : nodes) {
    std::vector<std::pair<Id, size_t>> occurrences;
    for (const auto &[e, r] : original_routes)
      for (size_t i = 0; i < r.size(); ++i)
        if (r[i] == n)
          occurrences.push_back({e, i});
    require(occurrences.size() == 2, "Crossing vertex needs two route spans");
    Ids neighbors;
    for (const auto &[e, i] : occurrences) {
      const auto &r = original_routes.at(e);
      require(i > 0 && i + 1 < r.size(), "Crossing at route endpoint");
      neighbors.push_back(r[i - 1]);
      neighbors.push_back(r[i + 1]);
    }
    require(as_set(neighbors).size() == 4, "Crossing needs four neighbors");
    Point center = original_positions.at(n);
    double radius = std::numeric_limits<double>::infinity();
    for (const auto &[key, p] : original_positions)
      if (key != n)
        radius = std::min(radius, distance(center, p) / 8);
    for (const auto &[a, b] : segments)
      if (a != n && b != n)
        radius =
            std::min(radius, point_segment(center, original_positions.at(a),
                                           original_positions.at(b)) /
                                 4);
    require(std::isfinite(radius) && radius > 0, "Crossing lacks clearance");
    std::vector<Ids> knots;
    for (size_t occurrence = 0; occurrence < occurrences.size(); ++occurrence) {
      auto [e, i] = occurrences[occurrence];
      const auto &r = original_routes.at(e);
      Ids pair;
      for (size_t end = 0; end < 2; ++end) {
        Id neighbor = r[end == 0 ? i - 1 : i + 1];
        Point direction = subtract(original_positions.at(neighbor), center);
        double length = std::hypot(direction[0], direction[1]);
        Id key = "b:" + e + ":cross:" + n + ":" + std::to_string(occurrence) +
                 ":" + std::to_string(end);
        require(!positions.contains(key), "Crossing knot collision");
        positions[key] = {center[0] + radius * direction[0] / length,
                          center[1] + radius * direction[1] / length};
        pair.push_back(key);
      }
      replacements[{e, i}] = pair;
      knots.push_back(pair);
    }
    Point a = positions.at(knots[0][0]), b = positions.at(knots[0][1]),
          c = positions.at(knots[1][0]), d = positions.at(knots[1][1]);
    Point u = subtract(b, a), v = subtract(d, c), offset = subtract(c, a);
    double determinant = cross(u, v);
    require(determinant != 0, "Crossing chords are parallel");
    double t = cross(offset, v) / determinant,
           s = cross(offset, u) / determinant;
    require(t > 0 && t < 1 && s > 0 && s < 1,
            "Routes do not alternate at crossing");
    records.push_back(
        {{"id", n},
         {"edges", Ids{occurrences[0].first, occurrences[1].first}},
         {"point", Point{a[0] + t * u[0], a[1] + t * u[1]}},
         {"route_points", knots},
         {"clearance_radius", radius}});
  }
  for (auto &[e, r] : routes) {
    Ids result;
    for (size_t i = 0; i < r.size(); ++i) {
      auto it = replacements.find({e, i});
      if (it == replacements.end())
        result.push_back(r[i]);
      else
        result.insert(result.end(), it->second.begin(), it->second.end());
    }
    r = result;
  }
  for (const auto &n : nodes)
    positions.erase(n);
  return records;
}
struct Fan {
  Id flow, owner, spoke;
  Ids endpoints;
};
inline void spread_fans(Positions &positions, const Routes &routes,
                        const std::vector<Fan> &groups) {
  bool repeated = false;
  for (const auto &f : groups)
    repeated |= f.endpoints.size() > 1;
  if (!repeated)
    return;
  std::vector<Dart> segments;
  for (const auto &[e, r] : routes) {
    bool all = true;
    for (const auto &n : r)
      all &= positions.contains(n);
    if (all)
      for (size_t i = 1; i < r.size(); ++i)
        segments.push_back({r[i - 1], r[i]});
  }
  double radius = std::numeric_limits<double>::infinity();
  for (const auto &[a, b] : segments)
    radius = std::min(radius, distance(positions.at(a), positions.at(b)) / 4);
  for (const auto &f : groups) {
    if (f.endpoints.size() == 1)
      continue;
    Id endpoint = f.endpoints.front();
    Point a = positions.at(f.owner), b = positions.at(endpoint),
          u = subtract(b, a);
    for (const auto &[key, p] : positions)
      if (key != f.owner && key != endpoint)
        radius = std::min(radius, point_segment(p, a, b) / 4);
    for (const auto &[first, second] : segments) {
      if (Set{first, second} == Set{f.owner, endpoint})
        continue;
      Point c = positions.at(first), d = positions.at(second);
      if (f.owner == first || f.owner == second) {
        Point other = f.owner == first ? d : c, w = subtract(other, a);
        double det = cross(u, w), ul = std::hypot(u[0], u[1]),
               wl = std::hypot(w[0], w[1]);
        if (std::abs(det) <= 1e-12 * ul * wl) {
          require(u[0] * w[0] + u[1] * w[1] < 0, "Auxiliary spokes overlap");
          radius = std::min({radius, ul / 4, wl / 4});
        } else
          radius = std::min(
              radius, std::abs(det) / (4 * (std::abs(u[0]) + std::abs(w[0]))));
      } else {
        double gap = std::min({point_segment(a, c, d), point_segment(b, c, d),
                               point_segment(c, a, b), point_segment(d, a, b)});
        radius = std::min(radius, gap / 4);
      }
    }
  }
  require(std::isfinite(radius) && radius > 0, "External fan lacks clearance");
  for (const auto &f : groups)
    if (f.endpoints.size() > 1) {
      Point p = positions.at(f.endpoints.front());
      for (size_t i = 0; i < f.endpoints.size(); ++i)
        positions[f.endpoints[i]] = {
            p[0],
            p[1] +
                radius * (2 * double(i) / double(f.endpoints.size() - 1) - 1)};
    }
}
inline Json group(const Json &children) {
  return children.size() == 1 ? children[0]
                              : Json{{"kind", "group"}, {"children", children}};
}
inline std::map<Id, Id> side_hubs(Graph &graph, const Expansion &expansion,
                                  const Id &exterior,
                                  const std::map<Id, Ids> &sides) {
  const auto &ports = expansion.leaf_ports.at(exterior);
  std::map<Id, Id> hubs;
  for (const Id side : {"incoming", "outgoing"}) {
    const auto &spokes = sides.at(side);
    Id node = ports.at(spokes[0]);
    if (spokes.size() == 1) {
      Id spoke = spokes[0], owner = graph.other(spoke, node),
         hub = "__metric_" + side + "_hub__",
         ref = "__metric_" + side + "_reference__";
      require(!graph.nodes.count(hub) && !graph.edges.contains(ref),
              "Metric hub identity collision");
      graph.nodes.insert(hub);
      graph.edges[spoke] = {owner, hub};
      graph.edges[ref] = {hub, node};
      replace(graph.rotation[node], spoke, {ref});
      graph.rotation[hub] = {spoke, ref};
      graph.protected_edges.insert(ref);
      node = hub;
    } else
      for (const auto &s : spokes)
        require(ports.at(s) == node, "Side did not expand to one group");
    hubs[side] = node;
  }
  Set hub_ids;
  for (const auto &[side, id] : hubs)
    hub_ids.insert(id);
  auto centers = difference(expansion.members.at(exterior), hub_ids);
  require(centers.size() == 1, "Exterior has no unique center");
  Id center = *centers.begin();
  Ids refs = graph.rotation.at(center);
  require(refs.size() == 2 && graph.protected_edges.count(refs[0]) &&
              graph.protected_edges.count(refs[1]),
          "Exterior center is not protected degree two");
  Id bridge = "__metric_exterior_bridge__";
  for (const auto &e : refs) {
    Id other = graph.other(e, center);
    replace(graph.rotation[other], e, {bridge});
    graph.edges.erase(e);
    graph.protected_edges.erase(e);
  }
  graph.rotation.erase(center);
  graph.nodes.erase(center);
  graph.edges[bridge] = {hubs.at("incoming"), hubs.at("outgoing")};
  graph.protected_edges.insert(bridge);
  graph.validate();
  return hubs;
}
inline Json initialize(const Json &diagram, double scale = 2.4,
                       bool external_sides = true) {
  require(std::isfinite(scale) && scale > 0,
          "Spacing must be finite and positive");
  Set nodes;
  Json node_ids = Json::object(), external_ids = Json::object(),
       endpoints = Json::object();
  for (size_t i = 0; i < diagram.at("vertices").size(); ++i) {
    Id key = std::to_string(i), n = "v:" + key;
    node_ids[key] = n;
    nodes.insert(n);
  }
  OrderedMap<Ends> edges;
  Routes routes, segments;
  Ids incoming, outgoing;
  std::map<Id, Ids> spokes{{"incoming", {}}, {"outgoing", {}}};
  std::vector<Fan> fans;
  const Id exterior = "__ec_exterior__";
  Set unique;
  for (size_t index = 0; index < diagram.at("edges").size(); ++index) {
    const auto &record = diagram.at("edges")[index];
    Json e;
    if (record.is_array()) {
      e = record[0];
      e["id"] = index;
      auto ext = record[1].value("external", Json::object());
      if (ext.is_object())
        e["state"] = ext.value("state", Json(nullptr));
    } else
      e = record;
    Id id = e.contains("id")
                ? (e["id"].is_string() ? e["id"].get<Id>() : e["id"].dump())
                : std::to_string(index);
    require(unique.insert(id).second, "Duplicate physical edge ID");
    auto source = e.at("source"), target = e.at("target");
    require(!source.is_null() || !target.is_null(), "Edge has no endpoint");
    auto key = [](const Json &v) {
      return v.is_string() ? v.get<Id>() : v.dump();
    };
    for (const auto &v : {source, target})
      if (!v.is_null())
        require(node_ids.contains(key(v)), "Unknown edge endpoint");
    Id a = source.is_null() ? "x:" + id : node_ids.at(key(source)).get<Id>(),
       b = target.is_null() ? "x:" + id : node_ids.at(key(target)).get<Id>();
    endpoints[id] = {{"source", source}, {"target", target}};
    if (source.is_null() || target.is_null()) {
      Id endpoint = source.is_null() ? a : b, owner = source.is_null() ? b : a;
      Id flow = !e.contains("state") || e["state"].is_null()
                    ? (source.is_null() ? "incoming" : "outgoing")
                    : e["state"].get<Id>();
      require(flow == "incoming" || flow == "outgoing",
              "Unknown external flow");
      (flow == "incoming" ? incoming : outgoing).push_back(endpoint);
      external_ids[id] = endpoint;
      routes[id] = {a, b};
      auto it = std::find_if(fans.begin(), fans.end(), [&](const Fan &f) {
        return f.flow == flow && f.owner == owner;
      });
      if (it == fans.end()) {
        Id spoke = "spoke:" + flow + ":" + owner;
        fans.push_back({flow, owner, spoke, {endpoint}});
        spokes[flow].push_back(spoke);
        edges[spoke] = {owner, exterior};
      } else
        it->endpoints.push_back(endpoint);
    } else {
      Ids route{a};
      size_t bends = a == b ? 2 : 1;
      for (size_t i = 0; i < bends; ++i)
        route.push_back("b:" + id + ":" + std::to_string(i));
      route.push_back(b);
      routes[id] = route;
      nodes.insert(route.begin(), route.end());
      for (size_t i = 1; i < route.size(); ++i) {
        Id s = "route:" + id + ":" + std::to_string(i - 1);
        edges[s] = {route[i - 1], route[i]};
        segments[id].push_back(s);
      }
    }
  }
  Json constraints = Json::object();
  bool side_mode = external_sides && !incoming.empty() && !outgoing.empty();
  if (!fans.empty()) {
    nodes.insert(exterior);
    Json leaves = Json::array(), sides = Json::array();
    for (const Id side : {"incoming", "outgoing"}) {
      for (const auto &e : spokes[side])
        leaves.push_back(e);
      if (!spokes[side].empty())
        sides.push_back(group(spokes[side]));
    }
    constraints[exterior] = side_mode ? group(sides) : group(leaves);
  }
  Graph initial(nodes, edges);
  for (const auto &f : fans)
    initial.protected_edges.insert(f.spoke);
  Expansion expansion(initial, constraints);
  Planarization result = planarize(expansion.graph);
  Graph embedded = result.embedding;
  std::map<Id, Id> hubs;
  std::optional<Ids> outer;
  if (side_mode) {
    hubs = side_hubs(embedded, expansion, exterior, spokes);
    outer = Ids{hubs.at("incoming"), hubs.at("outgoing")};
  } else if (!fans.empty()) {
    Id hub = *expansion.members.at(exterior).begin();
    hubs["radial"] = hub;
    outer = Ids{hub};
  }
  Positions positions = draw(embedded, outer);
  for (const auto &[id, parts] : segments) {
    Id current = routes.at(id).front();
    Ids route{current};
    for (const auto &segment : parts)
      for (const auto &part : result.chains.at(segment)) {
        current = embedded.other(part, current);
        route.push_back(current);
      }
    require(current == routes.at(id).back(), "Planarization changed endpoint");
    routes[id] = route;
  }
  std::map<Id, Point> hub_positions;
  for (const auto &[side, hub] : hubs) {
    hub_positions[side] = positions.at(hub);
    positions.erase(hub);
  }
  std::map<Id, double> rails;
  if (side_mode) {
    Point left = hub_positions.at("incoming"),
          right = hub_positions.at("outgoing");
    double low = std::numeric_limits<double>::infinity(), high = -low;
    for (const auto &[n, p] : positions) {
      low = std::min(low, p[0]);
      high = std::max(high, p[0]);
    }
    require(left[0] < low && low <= high && high < right[0],
            "Exterior hub edge did not bound drawing");
    rails["incoming"] = (left[0] + low) / 2;
    rails["outgoing"] = (right[0] + high) / 2;
  }
  for (const auto &f : fans) {
    Point p = positions.at(f.owner),
          hub = hub_positions.at(side_mode ? f.flow : "radial");
    double fraction =
        side_mode ? (rails.at(f.flow) - hub[0]) / (p[0] - hub[0]) : .15;
    positions[f.endpoints.front()] = {hub[0] + fraction * (p[0] - hub[0]),
                                      hub[1] + fraction * (p[1] - hub[1])};
  }
  spread_fans(positions, routes, fans);
  Json crossings =
      open_crossings(positions, routes, result.crossing_nodes.keys());
  Point center{0, 0};
  for (size_t axis = 0; axis < 2; ++axis) {
    long double sum = 0;
    for (const auto &[n, p] : positions)
      sum += p[axis];
    center[axis] = positions.empty() ? 0 : double(sum / positions.size());
  }
  std::vector<double> lengths;
  for (const auto &[e, r] : routes) {
    double length = 0;
    for (size_t i = 1; i < r.size(); ++i)
      length += distance(positions.at(r[i - 1]), positions.at(r[i]));
    lengths.push_back(length);
  }
  std::sort(lengths.begin(), lengths.end());
  double factor = 1;
  if (!lengths.empty()) {
    size_t m = lengths.size() / 2;
    double median =
        lengths.size() % 2 ? lengths[m] : (lengths[m - 1] + lengths[m]) / 2;
    require(median > 0, "Zero median route length");
    factor = scale / median;
  }
  for (auto &[n, p] : positions)
    for (size_t axis = 0; axis < 2; ++axis)
      p[axis] = factor * (p[axis] - center[axis]);
  for (auto &c : crossings) {
    Point p = c.at("point");
    for (size_t axis = 0; axis < 2; ++axis)
      p[axis] = factor * (p[axis] - center[axis]);
    c["point"] = p;
    c["clearance_radius"] = c.at("clearance_radius").get<double>() * factor;
  }
  Json audit = geometry(positions, routes);
  require(audit["contacts"].empty() && audit["overlaps"].empty() &&
              audit["degenerate"].empty(),
          "EC coordinates have geometric degeneracy");
  std::map<Ends, size_t> actual, expected;
  for (const auto &p : audit["crossings"]) {
    Ends es{p[0][0].get<Id>(), p[1][0].get<Id>()};
    std::sort(es.begin(), es.end());
    ++actual[es];
  }
  for (const auto &c : crossings) {
    Ends es = c.at("edges");
    std::sort(es.begin(), es.end());
    ++expected[es];
  }
  require(actual == expected, "Geometric crossings disagree with provenance");
  Json ps = Json::object(), rs = Json::object();
  for (const auto &[n, p] : positions)
    ps[n] = p;
  for (const auto &[e, r] : routes)
    rs[e] = r;
  return {{"positions", ps},
          {"routes", rs},
          {"node_ids", node_ids},
          {"external_ids", external_ids},
          {"incoming_ids", incoming},
          {"outgoing_ids", outgoing},
          {"edge_endpoints", endpoints},
          {"crossings", crossings},
          {"report",
           {{"status", crossings.empty() ? "ec-planar" : "ec-planarized"},
            {"method", "Gutwenger–Klein–Mutzel EC expansion and optimal "
                       "individual edge insertion; Chrobak–Payne coordinates"},
            {"constraint_tree", constraints.value(exterior, Json(nullptr))},
            {"constraint_scope",
             side_mode ? "consecutive exterior incoming/outgoing groups; "
                         "metric rails constructed separately"
                       : "unconstrained radial/cofacial external order"},
            {"external_routes_straight", true},
            {"external_side_constraints",
             {{"requested", external_sides}, {"realized", side_mode}}},
            {"insertions", result.insertions},
            {"geometry", audit}}}};
}
} // namespace ec
