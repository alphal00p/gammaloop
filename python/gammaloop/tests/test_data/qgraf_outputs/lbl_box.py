graphs = []
# To be imported after having the gammaloop UFO model named "scalars"

# Triangle massless
graphs.append(
    {
        "edges": {
            101: {
                "name": "p1",
                "PDG": 22,
                "type": "in",
                "momentum": "p1",
                "indices": (),
                "vertices": (101, 1)
            },
            102: {
                "name": "p2",
                "PDG": 22,
                "type": "out",
                "momentum": "p2",
                "indices": (),
                "vertices": (2, 102)
            },
            103: {
                "name": "p3",
                "PDG": 22,
                "type": "out",
                "momentum": "p3",
                "indices": (),
                "vertices": (3, 103)
            },
            104: {
                "name": "p4",
                "PDG": 22,
                "type": "out",
                "momentum": "p4",
                "indices": (),
                "vertices": (4, 104)
            },
            1: {
                "name": "q1",
                "PDG": 6,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (1, 2)
            },
            2: {
                "name": "q2",
                "PDG": 6,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (2, 3)
            },
            3: {
                "name": "q3",
                "PDG": 6,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (3, 4)
            },
            4: {
                "name": "q4",
                "PDG": 6,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (4, 1)
            }

        },
        "nodes": {
            101: {
                "PDGs": (22,),
                "momenta": ("p1",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (101,)
            },
            102: {
                "PDGs": (22,),
                "momenta": ("p2",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (102,)
            },
            103: {
                "PDGs": (22,),
                "momenta": ("p3",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (103,)
            },
            104: {
                "PDGs": (22,),
                "momenta": ("p4",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (104,)
            },
            1: {
                "PDGs": (6, 22, 6),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (1, 101, 4)
            },
            2: {
                "PDGs": (6, 22, 6),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (2, 102, 1)
            },
            3: {
                "PDGs": (6, 22, 6),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (3, 103, 2)
            },
            4: {
                "PDGs": (6, 22, 6),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (4, 104, 3)
            },
        },
        "overall_factor": "1"
    }
)
