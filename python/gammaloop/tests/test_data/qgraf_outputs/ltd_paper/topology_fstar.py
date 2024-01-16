graphs = []

graphs.append(  # type: ignore
    {
        "edges": {
            101: {
                "name": "p1",
                "PDG": 1001,
                "type": "in",
                "momentum": "p1",
                "indices": (),
                "vertices": (101, 1),
            },
            102: {
                "name": "p2",
                "PDG": 1001,
                "type": "out",
                "momentum": "-p1",
                "indices": (),
                "vertices": (2, 102),
            },
            1: {
                "name": "q1",
                "PDG": 1001,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (1, 3),
            },
            2: {
                "name": "q2",
                "PDG": 1001,
                "type": "virtual",
                "indices": (),
                "vertices": (3, 2),
            },
            3: {
                "name": "q3",
                "PDG": 1001,
                "type": "virtual",
                "indices": (),
                "vertices": (2, 4),
            },
            4: {
                "name": "q4",
                "PDG": 1001,
                "type": "virtual",
                "indices": (),
                "vertices": (4, 1),
            },
            5: {
                "name": "q5",
                "PDG": 1001,
                "type": "virtual",
                "indices": (),
                "vertices":  (3, 5),
            },
            6: {
                "name": "q6",
                "PDG": 1001,
                "type": "virtual",
                "indices": (),
                "vertices":  (2, 5),
            },
            7: {
                "name": "q7",
                "PDG": 1001,
                "type": "virtual",
                "indices": (),
                "vertices":  (4, 5),
            },
        },
        "nodes": {
            101: {
                "PDGs": (1001,),
                "momenta": ("p1",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (101,),
            },
            102: {
                "PDGs": (1001,),
                "momenta": ("p2",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (102,),
            },
            1: {
                "PDGs": (1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (101, 1, 4),
            },
            2: {
                "PDGs": (1001, 1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (102, 2, 6, 3),
            },
            3: {
                "PDGs": (1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (1, 2, 5),
            },
            4: {
                "PDGs": (1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (4, 3, 7),
            },
            5: {
                "PDGs": (1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (5, 6, 7),
            },
        },
        "overall_factor": "1"
    }
)
