graphs = []

# Bubble
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
                "type": "out",
                "momentum": "p2",
                "indices": (),
                "vertices": (12, 102)
            },
            1: {
                "name": "q1",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "k1",
                "indices": (),
                "vertices": (11, 12)
            },
            2: {
                "name": "q2",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "k1-p1",
                "indices": (),
                "vertices": (12, 11)
            }
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
            11: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("p1", "-k1", "k1-p1"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (1001, 1, 2)
            },
            12: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("-p2", "k1", "-k1+p1"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (1002, 2, 1)
            }
        },
        "overall_factor": "1",
    }
)
