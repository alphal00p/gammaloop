graphs = []

graphs.append(  # type: ignore
    {
        "edges": {
            101: {
                "name": "p1",
                "PDG": 1000,
                "type": "in",
                "momentum": "",
                "indices": (),
                "vertices": (101, 1),
            },
            102: {
                "name": "p2",
                "PDG": 1000,
                "momentum": "",
                "type": "out",
                "indices": (),
                "vertices": (3, 102),
            },
            1: {
                "name": "q1",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (1, 2)
            },
            2: {
                "name": "q2",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (2, 3)
            },
            3: {
                "name": "q3",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (3, 4)
            },
            4: {
                "name": "q4",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (4, 1)
            },
            5: {
                "name": "q5",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (5, 1)
            },
            6: {
                "name": "q6",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (4, 6)
            },
            7: {
                "name": "q7",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (5, 6)
            },
            8: {
                "name": "q8",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (6, 2)
            },
            9: {
                "name": "q9",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (5, 3)
            },
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
            1: {
                "PDGs": (1000, 1000, 1000, 1000),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (101, 1, 4, 5),
            },
            2: {
                "PDGs": (1000, 1000, 1000),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (1, 2, 8),
            },
            3: {
                "PDGs": (1000, 1000, 1000, 1000),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (2, 3, 9, 102),
            },
            4: {
                "PDGs": (1000, 1000, 1000),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (3, 4, 6),
            },
            5: {
                "PDGs": (1000, 1000, 1000),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (5, 7, 9),
            },
            6: {
                "PDGs": (1000, 1000, 1000),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (6, 7, 8),
            },
        },
        "overall_factor": "1",
    }
)
