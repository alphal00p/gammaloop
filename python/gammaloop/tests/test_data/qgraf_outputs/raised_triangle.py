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
                "type": "out",
                "momentum": "p2",
                "indices": (),
                "vertices": (2, 102)
            },
            103: {
                "name": "p3",
                "PDG": 1000,
                "type": "out",
                "momentum": "p3",
                "indices": (),
                "vertices": (3, 103)
            },
            1: {
                "name": "q1",
                "PDG": 1000,
                "type": "virtual",
                "indices": (),
                "vertices": (1, 4)
            },
            2: {
                "name": "q2",
                "PDG": 1000,
                "type": "virtual",
                "indices": (),
                "vertices": (4, 5),
            },
            3: {
                "name": "q3",
                "PDG": 1000,
                "type": "virtual",
                "indices": (),
                "vertices": (4, 5)
            },
            4: {
                "name": "q4",
                "PDG": 1000,
                "type": "virtual",
                "indices": (),
                "vertices": (5, 2)
            },
            5: {
                "name": "q5",
                "PDG": 1000,
                "type": "virtual",
                "indices": (),
                "vertices": (2, 3)
            },
            6: {
                "name": "q6",
                "PDG": 1000,
                "type": "virtual",
                "indices": (),
                "vertices": (3, 1)
            }
        },
        "nodes": {
            101: {
                "PDGs": (1000,),
                "momenta": ("p1",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (101,),
            },
            102: {
                "PDGs": (1000,),
                "momenta": ("p2",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (102,),
            },
            103: {
                "PDGs": (1000,),
                "momenta": ("p3",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (103,),
            },
            1: {
                "PDGs": (1000, 1000, 1000),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (1, 6, 101),
            },
            2: {
                "PDGs": (1000, 1000, 1000),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (4, 5, 102),
            },
            3: {
                "PDGs": (1000, 1000, 1000),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (5, 6, 103),
            },
            4: {
                "PDGs": (1000, 1000, 1000),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (1, 2, 3),
            },
            5: {
                "PDGs": (1000, 1000, 1000),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (2, 3, 4),
            },
        },
        "overall_factor": "1"
    }
)
