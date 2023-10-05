graphs = []
# To be imported after having the gammaloop UFO model named "scalars"

# Horizontal edge name convention: q1XY with X Y coordinate of the horizontal edge on a grid
# Vertical edge name convention: q2XY with X Y coordinate of the horizontal edge on a grid
# Vertex name convention: XY with X Y coordinate of the node on a grid

# Fishnet 2x2 massless
graphs.append(
    {
        "edges": {
            1001: {
                "name": "p1",
                "PDG": 1000,
                "type": "in",
                "momentum": "p1",
                "indices": (),
                "vertices": (101, 11)
            },
            1002: {
                "name": "p2",
                "PDG": 1000,
                "type": "in",
                "momentum": "p2",
                "indices": (),
                "vertices": (102, 31)
            },
            1003: {
                "name": "p3",
                "PDG": 1000,
                "type": "out",
                "momentum": "p3",
                "indices": (),
                "vertices": (13, 103)
            },
            1004: {
                "name": "p4",
                "PDG": 1000,
                "type": "out",
                "momentum": "p4",
                "indices": (),
                "vertices": (33, 104)
            },
            111: {
                "name": "q111",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (11, 12)
            },
            112: {
                "name": "q112",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (12, 13)
            },
            121: {
                "name": "q121",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (21, 22)
            },
            122: {
                "name": "q122",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (22, 23)
            },
            131: {
                "name": "q131",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (31, 32)
            },
            132: {
                "name": "q132",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (32, 33)
            },
            211: {
                "name": "q211",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (11, 21)
            },
            212: {
                "name": "q212",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (12, 22)
            },
            213: {
                "name": "q213",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (13, 23)
            },
            221: {
                "name": "q221",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (21, 31)
            },
            222: {
                "name": "q222",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (22, 32)
            },
            223: {
                "name": "q223",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (23, 33)
            },
        },
        "nodes": {
            101: {
                "PDGs": (1000,),
                "momenta": ("p1",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (1001,)
            },
            102: {
                "PDGs": (1000,),
                "momenta": ("p2",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (1002,)
            },
            103: {
                "PDGs": (1000,),
                "momenta": ("p3",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (1003,)
            },
            104: {
                "PDGs": (1000,),
                "momenta": ("p4",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (1004,)
            },
            11: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (1001, 211, 111)
            },
            12: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (111, 212, 112)
            },
            13: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (112, 213, 1003)
            },
            21: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (211, 121, 221)
            },
            22: {
                "PDGs": (1000, 1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 1,
                "edge_ids": (121, 212, 122, 222)
            },
            23: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (213, 122, 223)
            },
            31: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (221, 1002, 131)
            },
            32: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (131, 222, 132)
            },
            33: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (132, 1004, 223)
            },
        },
        "overall_factor": "1"
    }
)
