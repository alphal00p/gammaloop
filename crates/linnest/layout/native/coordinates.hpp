#pragma once
#include "embedding.hpp"
#include <cmath>
#include <numeric>

namespace ec {
using Point = std::array<double, 2>;
using Positions = OrderedMap<Point>;
using Routes = OrderedMap<Ids>;
inline double distance(Point a, Point b) {
  return std::hypot(a[0] - b[0], a[1] - b[1]);
}
inline double cross(Point a, Point b) { return a[0] * b[1] - a[1] * b[0]; }
inline Point subtract(Point a, Point b) { return {a[0] - b[0], a[1] - b[1]}; }
inline double point_segment(Point p, Point a, Point b) {
  Point d = subtract(b, a);
  double squared = d[0] * d[0] + d[1] * d[1];
  require(squared > 0, "Zero-length planar segment");
  double f = std::clamp(((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / squared,
                        0.0, 1.0);
  return distance(p, {a[0] + f * d[0], a[1] + f * d[1]});
}
// Integer-set iteration follows CPython's stable positive-integer hash table.
// The reference's canonical ordering selects with set.pop(); reproducing its
// tie order keeps the initial geometry unchanged across the native transition.
class IntSet {
  std::vector<int> table = std::vector<int>(8, -1);
  size_t used = 0, fill = 0, finger = 0;
  size_t slot(int value, bool insertion) const {
    size_t mask = table.size() - 1, i = size_t(value) & mask, perturb = value,
           free = table.size();
    for (;;) {
      size_t j = i, probes = i + 9 <= mask ? 9 : 0;
      do {
        int entry = table[j];
        if (entry == value)
          return j;
        if (entry == -1)
          return insertion && free < table.size() ? free : j;
        if (entry == -2)
          free = j;
        ++j;
      } while (probes--);
      perturb >>= 5;
      i = (i * 5 + 1 + perturb) & mask;
    }
  }
  void resize() {
    auto old = table;
    size_t size = 8;
    while (size <= used * 4)
      size <<= 1;
    table.assign(size, -1);
    fill = used;
    for (int value : old)
      if (value >= 0)
        table[slot(value, true)] = value;
  }

public:
  IntSet() = default;
  explicit IntSet(const std::vector<int> &values) {
    for (int v : values)
      add(v);
  }
  bool contains(int v) const { return table[slot(v, false)] == v; }
  void add(int v) {
    size_t i = slot(v, true);
    if (table[i] == v)
      return;
    bool empty = table[i] == -1;
    table[i] = v;
    ++used;
    if (empty && ++fill * 5 >= (table.size() - 1) * 3)
      resize();
  }
  void erase(int v) {
    size_t i = slot(v, false);
    if (table[i] == v) {
      table[i] = -2;
      --used;
    }
  }
  int pop() {
    require(used > 0, "No canonical-ordering candidate");
    size_t i = finger & (table.size() - 1);
    while (table[i] < 0)
      i = (i + 1) & (table.size() - 1);
    int v = table[i];
    table[i] = -2;
    --used;
    finger = i + 1;
    return v;
  }
  std::vector<int> iteration() const {
    std::vector<int> out;
    for (int v : table)
      if (v >= 0)
        out.push_back(v);
    return out;
  }
};
// Chrobak–Payne triangulation/canonical ordering translated from NetworkX 3.7.
// Copyright NetworkX Developers, BSD-3-Clause; see LICENSE.networkx.
struct Embedding {
  std::map<int, std::vector<int>> rotation;
  std::vector<int> nodes;
  void node(int n) {
    if (!rotation.count(n)) {
      rotation[n] = {};
      nodes.push_back(n);
    }
  }
  void add(int u, int v, std::optional<int> ref = std::nullopt,
           bool before = true) {
    node(u);
    node(v);
    auto &r = rotation[u];
    if (r.empty()) {
      require(!ref, "Reference on empty rotation");
      r.push_back(v);
      return;
    }
    require(bool(ref), "Missing cyclic insertion reference");
    auto it = std::find(r.begin(), r.end(), *ref);
    require(it != r.end(), "Missing cyclic neighbor");
    size_t i = it - r.begin();
    if (!before)
      ++i;
    r.insert(r.begin() + i, v);
  }
  int neighbor(int u, int v, int delta) const {
    const auto &r = rotation.at(u);
    auto it = std::find(r.begin(), r.end(), v);
    require(it != r.end(), "Missing halfedge");
    long i = it - r.begin();
    return r[(i + long(r.size()) + delta) % r.size()];
  }
  bool edge(int u, int v) const {
    const auto &r = rotation.at(u);
    return std::find(r.begin(), r.end(), v) != r.end();
  }
  std::vector<int> face(int u, int v) const {
    int start = u, next = v;
    std::vector<int> result;
    do {
      result.push_back(u);
      int w = neighbor(v, u, -1);
      u = v;
      v = w;
      require(result.size() <= rotation.size() * rotation.size() + 2,
              "Face walk did not close");
    } while (u != start || v != next);
    return result;
  }
  static Embedding from(const std::map<int, std::vector<int>> &data) {
    Embedding out;
    for (const auto &[u, r] : data) {
      std::optional<int> ref;
      for (auto it = r.rbegin(); it != r.rend(); ++it) {
        out.add(u, *it, ref, true);
        ref = *it;
      }
    }
    for (const auto &[u, r] : data)
      out.node(u);
    return out;
  }
  std::vector<int> biconnect(int start, int outgoing,
                             std::set<std::pair<int, int>> &counted) {
    if (!counted.insert({start, outgoing}).second)
      return {};
    int a = start, b = outgoing, c = neighbor(b, a, -1);
    std::vector<int> face{start};
    std::set<int> seen{start};
    while (b != start || c != outgoing) {
      require(a != b, "Invalid halfedge");
      if (seen.count(b)) {
        add(a, c, b, false);
        add(c, a, b, true);
        counted.insert({b, c});
        counted.insert({c, a});
        b = a;
      } else {
        seen.insert(b);
        face.push_back(b);
      }
      a = b;
      b = c;
      c = neighbor(b, a, -1);
      counted.insert({a, b});
    }
    return face;
  }
  void triangulate_face(int a, int b) {
    int c = neighbor(b, a, -1), d = neighbor(c, b, -1);
    if (a == b || a == c)
      return;
    while (a != d) {
      if (edge(a, c)) {
        a = b;
        b = c;
        c = d;
      } else {
        add(a, c, b, false);
        add(c, a, b, true);
        b = c;
        c = d;
      }
      d = neighbor(c, b, -1);
    }
  }
  std::vector<int> triangulate(bool full) {
    std::set<std::pair<int, int>> counted;
    std::vector<std::vector<int>> faces;
    size_t outer = 0;
    for (int u : nodes) {
      if (rotation[u].empty())
        continue;
      int first = rotation[u][0], v = first;
      do {
        auto face = biconnect(u, v, counted);
        if (!face.empty()) {
          faces.push_back(face);
          if (face.size() > faces[outer].size())
            outer = faces.size() - 1;
        }
        v = neighbor(u, v, 1);
      } while (v != first);
    }
    require(!faces.empty(), "No planar face");
    for (size_t i = 0; i < faces.size(); ++i)
      if (i != outer || full)
        triangulate_face(faces[i][0], faces[i][1]);
    auto result = faces[outer];
    if (full)
      result = {result[0], result[1], neighbor(result[1], result[0], -1)};
    return result;
  }
  std::vector<std::pair<int, std::vector<int>>>
  canonical(const std::vector<int> &outer) const {
    int v1 = outer[0], v2 = outer[1];
    std::map<int, int> chords, ccw, cw;
    std::set<int> marked;
    IntSet ready(outer);
    int prev = v2;
    for (size_t i = 2; i < outer.size(); ++i) {
      ccw[prev] = outer[i];
      prev = outer[i];
    }
    ccw[prev] = v1;
    prev = v1;
    for (size_t i = outer.size(); i-- > 1;) {
      cw[prev] = outer[i];
      prev = outer[i];
    }
    auto on = [&](int x) {
      return !marked.count(x) && (ccw.count(x) || x == v1);
    };
    auto adjacent = [&](int x, int y) {
      if (!ccw.count(x))
        return cw.at(x) == y;
      if (!cw.count(x))
        return ccw.at(x) == y;
      return ccw.at(x) == y || cw.at(x) == y;
    };
    for (int v : outer)
      for (int w : rotation.at(v))
        if (on(w) && !adjacent(v, w)) {
          ++chords[v];
          ready.erase(v);
        }
    std::vector<std::pair<int, std::vector<int>>> order(nodes.size());
    order[0] = {v1, {}};
    order[1] = {v2, {}};
    ready.erase(v1);
    ready.erase(v2);
    for (size_t k = nodes.size(); k-- > 2;) {
      int v = ready.pop();
      marked.insert(v);
      int wp = -1, wq = -1;
      for (int n : rotation.at(v)) {
        if (marked.count(n))
          continue;
        if (on(n)) {
          if (n == v1)
            wp = n;
          else if (n == v2)
            wq = n;
          else if (cw.at(n) == v)
            wp = n;
          else
            wq = n;
        }
        if (wp >= 0 && wq >= 0)
          break;
      }
      require(wp >= 0 && wq >= 0, "Canonical boundary not found");
      std::vector<int> contour{wp};
      int n = wp;
      while (n != wq) {
        int next = neighbor(v, n, -1);
        contour.push_back(next);
        cw[n] = next;
        ccw[next] = n;
        n = next;
        require(contour.size() <= nodes.size(),
                "Canonical contour did not close");
      }
      if (contour.size() == 2) {
        if (--chords[wp] == 0)
          ready.add(wp);
        if (--chords[wq] == 0)
          ready.add(wq);
      } else {
        std::vector<int> middle(contour.begin() + 1, contour.end() - 1);
        IntSet newly(middle);
        for (int w : newly.iteration()) {
          ready.add(w);
          for (int other : rotation.at(w))
            if (on(other) && !adjacent(w, other)) {
              ++chords[w];
              ready.erase(w);
              if (!newly.contains(other)) {
                ++chords[other];
                ready.erase(other);
              }
            }
        }
      }
      order[k] = {v, contour};
    }
    return order;
  }
  OrderedMap<Point>
  coordinates(std::optional<std::array<int, 2>> outer_edge = std::nullopt,
              std::optional<int> outer_vertex = std::nullopt) {
    OrderedMap<Point> out;
    if (nodes.size() < 4) {
      std::array<Point, 3> triangle{{{0, 0}, {2, 0}, {1, 1}}};
      auto sequence = nodes;
      if (nodes.size() == 3 && outer_edge) {
        sequence = {(*outer_edge)[0], (*outer_edge)[1]};
        for (int n : nodes)
          if (n != sequence[0] && n != sequence[1])
            sequence.push_back(n);
      }
      for (size_t i = 0; i < sequence.size(); ++i)
        out[std::to_string(sequence[i])] = triangle[i];
      return out;
    }
    auto outer = triangulate(bool(outer_edge) || bool(outer_vertex));
    if (outer_edge)
      outer = face((*outer_edge)[0], (*outer_edge)[1]);
    else if (outer_vertex)
      outer = face(*outer_vertex, rotation.at(*outer_vertex)[0]);
    auto order = canonical(outer);
    std::map<int, std::optional<int>> left, right;
    std::map<int, long long> dx, y;
    int a = order[0].first, b = order[1].first, c = order[2].first;
    dx[a] = 0;
    y[a] = 0;
    right[a] = c;
    left[a] = std::nullopt;
    dx[b] = 1;
    y[b] = 0;
    right[b] = std::nullopt;
    left[b] = std::nullopt;
    dx[c] = 1;
    y[c] = 1;
    right[c] = b;
    left[c] = std::nullopt;
    for (size_t k = 3; k < order.size(); ++k) {
      int v = order[k].first;
      const auto &r = order[k].second;
      int wp = r[0], wp1 = r[1], wq = r.back(), wq1 = r[r.size() - 2];
      ++dx[wp1];
      ++dx[wq];
      long long gap = 0;
      for (size_t i = 1; i < r.size(); ++i)
        gap += dx[r[i]];
      dx[v] = (-y[wp] + gap + y[wq]) / 2;
      y[v] = (y[wp] + gap + y[wq]) / 2;
      dx[wq] = gap - dx[v];
      if (r.size() > 2)
        dx[wp1] -= dx[v];
      right[wp] = v;
      right[v] = wq;
      if (r.size() > 2) {
        left[v] = wp1;
        right[wq1] = std::nullopt;
      } else
        left[v] = std::nullopt;
    }
    out[std::to_string(a)] = {0, double(y[a])};
    std::vector<int> pending{a};
    while (!pending.empty()) {
      int parent = pending.back();
      pending.pop_back();
      for (auto child : {left[parent], right[parent]})
        if (child) {
          out[std::to_string(*child)] = {out.at(std::to_string(parent))[0] +
                                             double(dx[*child]),
                                         double(y[*child])};
          pending.push_back(*child);
        }
    }
    return out;
  }
};
inline Positions draw(const Graph &graph,
                      std::optional<Ids> outer = std::nullopt) {
  std::vector<Id> names(graph.nodes.begin(), graph.nodes.end());
  std::map<Id, int> ids;
  for (size_t i = 0; i < names.size(); ++i)
    ids[names[i]] = i;
  std::map<int, std::vector<int>> rotations;
  for (const auto &n : names) {
    auto r = graph.rotation.at(n);
    if (!r.empty())
      std::rotate(r.begin(), std::min_element(r.begin(), r.end()), r.end());
    std::vector<int> neighbors;
    for (const auto &e : r)
      neighbors.push_back(ids.at(graph.other(e, n)));
    require(std::set<int>(neighbors.begin(), neighbors.end()).size() ==
                neighbors.size(),
            "Subdivide parallel drawing edges");
    rotations[ids[n]] = neighbors;
  }
  Embedding full = Embedding::from(rotations);
  std::set<int> seen;
  std::vector<std::set<int>> components;
  for (size_t root = 0; root < names.size(); ++root) {
    if (seen.count(root))
      continue;
    std::set<int> component;
    std::vector<int> pending{int(root)};
    while (!pending.empty()) {
      int n = pending.back();
      pending.pop_back();
      if (!component.insert(n).second)
        continue;
      const auto &r = rotations.at(n);
      pending.insert(pending.end(), r.begin(), r.end());
    }
    seen.insert(component.begin(), component.end());
    components.push_back(component);
  }
  std::sort(
      components.begin(), components.end(), [&](const auto &a, const auto &b) {
        bool ao = outer && !a.count(ids.at(outer->front())),
             bo = outer && !b.count(ids.at(outer->front()));
        return std::make_pair(ao, *a.begin()) < std::make_pair(bo, *b.begin());
      });
  Positions positions;
  double top = 0;
  for (const auto &component : components) {
    std::map<int, std::vector<int>> data;
    for (int n : component)
      data[n] = full.rotation.at(n);
    auto local = Embedding::from(data);
    std::optional<std::array<int, 2>> edge;
    std::optional<int> vertex;
    if (outer && component.count(ids.at(outer->front()))) {
      if (outer->size() == 2)
        edge = std::array<int, 2>{ids.at((*outer)[0]), ids.at((*outer)[1])};
      else
        vertex = ids.at((*outer)[0]);
    }
    auto coordinates = local.coordinates(edge, vertex);
    if (!positions.empty()) {
      double low = std::numeric_limits<double>::infinity(), high = -low,
             minimum = low, maximum = high;
      for (const auto &[n, p] : positions) {
        low = std::min(low, p[0]);
        high = std::max(high, p[0]);
      }
      for (const auto &[n, p] : coordinates) {
        minimum = std::min(minimum, p[0]);
        maximum = std::max(maximum, p[0]);
      }
      double factor = (high - low) / (2 * std::max(1.0, maximum - minimum)),
             center = (minimum + maximum) / 2;
      for (auto &[n, p] : coordinates)
        p = {(low + high) / 2 + (p[0] - center) * factor,
             top + 4 + p[1] * factor};
    }
    for (const auto &[n, p] : coordinates)
      positions[names[std::stoul(n)]] = p;
    for (const auto &[n, p] : positions)
      top = std::max(top, p[1]);
  }
  return positions;
}
} // namespace ec
