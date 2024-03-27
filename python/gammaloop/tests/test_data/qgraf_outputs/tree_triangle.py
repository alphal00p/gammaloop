graphs = []

graphs.append(  # type: ignore
    {
        "edges": {
            101: {
                "name": "p1",
                "PDG": 1000,
                "type": "in",
                "momentum": "p1",
                "indices": (),
                "vertices": (101, 1)
            },
            102: {
                "name": "p2",
                "PDG": 1000,
                "type": "in",
                "momentum": "p2",
                "indices": (),
                "vertices": (102, 1)
            },
            103: {
                "name": "p3",
                "PDG": 1000,
                "type": "out",
                "momentum": "p3",
                "indices": (),
                "vertices": (3, 103)
            },
            104: {
                "name": "p4",
                "PDG": 1000,
                "type": "out",
                "momentum": "p4",
                "indices": (),
                "vertices": (4, 104)
            },
            1: {
                "name": "q1",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (1, 2)
            },
            2: {
                "name": "q2",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (2, 3)
            },
            3: {
                "name": "q3",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (3, 4)
            },
            4: {
                "name": "q4",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (4, 2)
            },
        },
        "nodes": {
            101: {
                "PDGs": (1000,),
                "momenta": ("p1",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (101,)
            },
            102: {
                "PDGs": (1000,),
                "momenta": ("p2",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (102,)
            },
            103: {
                "PDGs": (1000,),
                "momenta": ("p3",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (103,)
            },
            104: {
                "PDGs": (1000,),
                "momenta": ("p4",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (104,)
            },
            1: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (1, 101, 102),
            },
            2: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (2, 1, 4),
            },
            3: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (3, 2, 103),
            },
            4: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (4, 3, 104),
            }
        },
        "overall_factor": "1",
    }
)
