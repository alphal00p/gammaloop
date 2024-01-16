graphs = []

graphs.append({  # type: ignore
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
            "vertices": (102, 6)
        },
        103: {
            "name": "p3",
            "PDG": 1000,
            "type": "out",
            "momentum": "p3",
            "indices": (),
            "vertices": (5, 103)
        },
        104: {
            "name": "p4",
            "PDG": 1000,
            "type": "out",
            "momentum": "p4",
            "indices": (),
            "vertices": (10, 104)
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
            "vertices": (4, 5)
        },
        5: {
            "name": "q5",
            "PDG": 1000,
            "type": "virtual",
            "momentum": "",
            "indices": (),
            "vertices": (6, 7)
        },
        6: {
            "name": "q6",
            "PDG": 1000,
            "type": "virtual",
            "momentum": "",
            "indices": (),
            "vertices": (7, 8)
        },
        7: {
            "name": "q7",
            "PDG": 1000,
            "type": "virtual",
            "momentum": "",
            "indices": (),
            "vertices": (8, 9)
        },
        8: {
            "name": "q8",
            "PDG": 1000,
            "type": "virtual",
            "momentum": "",
            "indices": (),
            "vertices": (9, 10)
        },
        9: {
            "name": "q9",
            "PDG": 1000,
            "type": "virtual",
            "momentum": "",
            "indices": (),
            "vertices": (1, 6)
        },
        10: {
            "name": "q10",
            "PDG": 1000,
            "type": "virtual",
            "momentum": "",
            "indices": (),
            "vertices": (2, 7)
        },
        11: {
            "name": "q11",
            "PDG": 1000,
            "type": "virtual",
            "momentum": "",
            "indices": (),
            "vertices": (3, 8)
        },
        12: {
            "name": "q12",
            "PDG": 1000,
            "type": "virtual",
            "momentum": "",
            "indices": (),
            "vertices": (4, 9)
        },
        13: {
            "name": "q13",
            "PDG": 1000,
            "type": "virtual",
            "momentum": "",
            "indices": (),
            "vertices": (5, 10)
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
        103: {
            "PDGs": (1000,),
            "momenta": ("p3",),
            "indices": (),
            "vertex_id": -1,
            "edge_ids": (103,),
        },
        104: {
            "PDGs": (1000,),
            "momenta": ("p4",),
            "indices": (),
            "vertex_id": -1,
            "edge_ids": (104,),
        },
        1: {
            "PDGs": (1000, 1000, 1000),
            "momenta": (),
            "indices": (),
            "vertex_id": 0,
            "edge_ids": (101, 1, 9),
        },
        2: {
            "PDGs": (1000, 1000, 1000),
            "momenta": (),
            "indices": (),
            "vertex_id": 0,
            "edge_ids": (2, 1, 10),
        },
        3: {
            "PDGs": (1000, 1000, 1000),
            "momenta": (),
            "indices": (),
            "vertex_id": 0,
            "edge_ids": (2, 3, 11),
        },
        4: {
            "PDGs": (1000, 1000, 1000),
            "momenta": (),
            "indices": (),
            "vertex_id": 0,
            "edge_ids": (3, 4, 12),
        },
        5: {
            "PDGs": (1000, 1000, 1000),
            "momenta": (),
            "indices": (),
            "vertex_id": 0,
            "edge_ids": (4, 13, 103),
        },
        6: {
            "PDGs": (1000, 1000, 1000),
            "momenta": (),
            "indices": (),
            "vertex_id": 0,
            "edge_ids": (102, 9, 5),
        },
        7: {
            "PDGs": (1000, 1000, 1000),
            "momenta": (),
            "indices": (),
            "vertex_id": 0,
            "edge_ids": (5, 10, 6),
        },
        8: {
            "PDGs": (1000, 1000, 1000),
            "momenta": (),
            "indices": (),
            "vertex_id": 0,
            "edge_ids": (6, 7, 11),
        },
        9: {
            "PDGs": (1000, 1000, 1000),
            "momenta": (),
            "indices": (),
            "vertex_id": 0,
            "edge_ids": (7, 8, 12),
        },
        10: {
            "PDGs": (1000, 1000, 1000),
            "momenta": (),
            "indices": (),
            "vertex_id": 0,
            "edge_ids": (104, 9, 13),
        },
    },
    "overall_factor": 1,
})
