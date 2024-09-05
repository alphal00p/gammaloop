graphs = []

# Bubble
graphs.append(
    {
        "edges": {
            1001: {
                "name": "p1",
                "PDG": 22,
                "type": "in",
                "momentum": "p1",
                "indices": (),
                "vertices": (101, 11)
            }, 
            1002: {
                "name": "p2",
                "PDG": 22,
                "type": "out",
                "momentum": "p2",
                "indices": (),
                "vertices": (12, 102)
            },
            1: {
                "name": "q1",
                "PDG": 6,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (11, 12)
            },
            2: {
                "name": "q2",
                "PDG": 6,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (12, 11)
            }
        },
        "nodes": { 
            101: {
                "PDGs": (22,),
                "momenta": ("p1",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (1001,)
            },
            102: {
                "PDGs": (22,),
                "momenta": ("p2",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (1002,)
            },
            11: {
                "PDGs": (6, 22, 6),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (1,1001, 2)
            },
            12: {
                "PDGs": (6, 22, 6),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (2, 1002, 1)
            }
        },
        "overall_factor": "1",
    }
)
