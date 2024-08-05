graphs = []
graphs.append(
    {
        "edges": {
            101: {"name": "p2", "PDG": 22, "type": "in", "momentum": "p1", "indices": (), "vertices": (101, 1)},
            102: {"name": "p1", "PDG": 22, "type": "in", "momentum": "p2", "indices": (), "vertices": (102, 2)},
            103: {"name": "p3", "PDG": 22, "type": "out", "momentum": "p3", "indices": (), "vertices": (9, 103)},
            104: {"name": "p4", "PDG": 22, "type": "out", "momentum": "p4", "indices": (), "vertices": (7, 104)},
            105: {"name": "p5", "PDG": 22, "type": "out", "momentum": "p5", "indices": (), "vertices": (6, 105)},
            106: {"name": "p6", "PDG": 22, "type": "out", "momentum": "p6", "indices": (), "vertices": (4, 106)},
            1: {"name": "q1", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (1, 2)},
            2: {"name": "q2", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (2, 4)},
            4: {"name": "q4", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (4, 6)},
            6: {"name": "q6", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (6, 7)},
            7: {"name": "q7", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (7, 9)},
            9: {"name": "q9", "PDG": 6, "type": "virtual", "momentum": "N/A", "indices": (), "vertices": (9, 1)},
        },
        "nodes": {
            101: {"PDGs": (22,), "momenta": ("p1",), "indices": (), "vertex_id": -1, "edge_ids": (101,)},
            102: {"PDGs": (22,), "momenta": ("p2",), "indices": (), "vertex_id": -1, "edge_ids": (102,)},
            103: {"PDGs": (22,), "momenta": ("p3",), "indices": (), "vertex_id": -1, "edge_ids": (103,)},
            104: {"PDGs": (22,), "momenta": ("p4",), "indices": (), "vertex_id": -1, "edge_ids": (104,)},
            105: {"PDGs": (22,), "momenta": ("p5",), "indices": (), "vertex_id": -1, "edge_ids": (105,)},
            106: {"PDGs": (22,), "momenta": ("p6",), "indices": (), "vertex_id": -1, "edge_ids": (106,)},
            1: {"PDGs": (6, 22, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (1, 101, 9)},
            2: {"PDGs": (6, 22, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (2, 102, 1)},
            4: {"PDGs": (6, 22, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (4, 106, 2)},
            6: {"PDGs": (6, 22, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (6, 105, 4)},
            7: {"PDGs": (6, 22, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (7, 104, 6)},
            9: {"PDGs": (6, 22, 6), "momenta": ("N/A", "N/A", "N/A"), "indices": (), "vertex_id": 0, "edge_ids": (9, 103, 7)},
        },
        "overall_factor": "1"
    }
)
