graphs = []
graphs.append(
    {
        "edges": {
            101: {"name": "p1", "PDG": 22, "type": "in", "momentum": "p1", "indices": (), "vertices": (101, 1)},
            102: {"name": "p2", "PDG": 22, "type": "in", "momentum": "p2", "indices": (), "vertices": (102, 2)},
            103: {"name": "p3", "PDG": 22, "type": "out", "momentum": "p3", "indices": (), "vertices": (3, 103)},
            104: {"name": "p4", "PDG": 22, "type": "out", "momentum": "p4", "indices": (), "vertices": (4, 104)},
            105: {"name": "p5", "PDG": 22, "type": "out", "momentum": "p5", "indices": (), "vertices": (5, 105)},
            106: {"name": "p6", "PDG": 22, "type": "out", "momentum": "p6", "indices": (), "vertices": (6, 106)},
            1: {"name": "q1", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (1, 2)},
            2: {"name": "q2", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (2, 21)},
            3: {"name": "q3", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (3, 11)},
            4: {"name": "q4", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (4, 5)},
            5: {"name": "q5", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (5, 22)},
            6: {"name": "q6", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (6, 12)},
            11: {"name": "q11", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (11, 4)},
            13: {"name": "q13", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (12, 1)},
            12: {"name": "q12", "PDG": 21, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (11, 12)},
            21: {"name": "q21", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (21, 3)},
            23: {"name": "q23", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (22, 6)},
            22: {"name": "q22", "PDG": 21, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (21, 22)},        
        },
        "nodes": {
            101: {"PDGs": (22,), "momenta": ("p1",), "indices": (), "vertex_id": -1, "edge_ids": (101,)},
            102: {"PDGs": (22,), "momenta": ("p2",), "indices": (), "vertex_id": -1, "edge_ids": (102,)},
            103: {"PDGs": (22,), "momenta": ("p3",), "indices": (), "vertex_id": -1, "edge_ids": (103,)},
            104: {"PDGs": (22,), "momenta": ("p4",), "indices": (), "vertex_id": -1, "edge_ids": (104,)},
            105: {"PDGs": (22,), "momenta": ("p5",), "indices": (), "vertex_id": -1, "edge_ids": (105,)},
            106: {"PDGs": (22,), "momenta": ("p6",), "indices": (), "vertex_id": -1, "edge_ids": (106,)},
            1: {"PDGs": (6, 22, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (13, 101, 1)},
            2: {"PDGs": (6, 22, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (1, 102, 2)},
            3: {"PDGs": (6, 22, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (2, 103, 3)},
            4: {"PDGs": (6, 22, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (11, 104, 4)},
            5: {"PDGs": (6, 22, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (4, 105, 5)},
            6: {"PDGs": (6, 22, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (5, 106, 6)},
            11: {"PDGs": (6, 21, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (3, 12, 11)},
            12: {"PDGs": (6, 21, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (6, 12, 13)},
            21: {"PDGs": (6, 21, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (2, 22, 21)},
            22: {"PDGs": (6, 21, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (5, 22, 23)},
        },
        "overall_factor": "1"
    }
)
