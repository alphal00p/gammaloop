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
                "vertices": (101, 1)
            },
            102: {
                "name": "p2",
                "PDG": 1001,
                "type": "in",
                "momentum": "p2",
                "indices": (),
                "vertices": (102, 2)
            },
            103: {
                "name": "p3",
                "PDG": 1001,
                "type": "in",
                "momentum": "p3",
                "indices": (),
                "vertices": (103, 3)
            },
            104: {
                "name": "p4",
                "PDG": 1001,
                "type": "in",
                "momentum": "p4",
                "indices": (),
                "vertices": (104, 4)
            },
            105: {
                "name": "p5",
                "PDG": 1001,
                "type": "in",
                "momentum": "p5",
                "indices": (),
                "vertices": (105, 5)
            },
            106: {
                "name": "p6",
                "PDG": 1001,
                "type": "out",
                "momentum": "p1 + p2 + p3 + p4 + p5",
                "indices": (),
                "vertices": (6, 106)
            },
            1: {
                "name": "q1",
                "PDG": 1001,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (1, 6),
            },
            2: {
                "name": "q2",
                "PDG": 1001,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (6, 9),
            },
            3: {
                "name": "q3",
                "PDG": 1001,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (9, 2),
            },
            4: {
                "name": "q4",
                "PDG": 1001,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (2, 3),
            },
            5: {
                "name": "q5",
                "PDG": 1001,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (3, 8),
            },
            6: {
                "name": "q6",
                "PDG": 1001,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (8, 1),
            },
            7: {
                "name": "q7",
                "PDG": 1001,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (9, 5),
            },
            8: {
                "name": "q8",
                "PDG": 1001,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (5, 7),
            },
            9: {
                "name": "q9",
                "PDG": 1001,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (7, 4),
            },
            10: {
                "name": "q10",
                "PDG": 1001,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (4, 8),
            },
            11: {
                "name": "q11",
                "PDG": 1001,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (2, 7),
            }
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
            103: {
                "PDGs": (1001,),
                "momenta": ("p3",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (103,),
            },
            104: {
                "PDGs": (1001,),
                "momenta": ("p4",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (104,),
            },
            105: {
                "PDGs": (1001,),
                "momenta": ("p5",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (105,),
            },
            106: {
                "PDGs": (1001,),
                "momenta": ("p6",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (106,),
            },
            1: {
                "PDGs": (1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (101, 1, 6)
            },
            2: {
                "PDGs": (1001, 1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (102, 11, 3, 4)
            },
            3: {
                "PDGs": (1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (103, 4, 5)
            },
            4: {
                "PDGs": (1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (104, 10, 9)
            },
            5: {
                "PDGs": (1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (105, 7, 8)
            },
            6: {
                "PDGs": (1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (106, 1, 2)
            },
            7: {
                "PDGs": (1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (8, 9, 11)
            },
            8: {
                "PDGs": (1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (6, 10, 5)
            },
            9: {
                "PDGs": (1001, 1001, 1001),
                "momenta": (),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (2, 7, 3)
            },
        },
        "overall_factor": "1"
    }
)
