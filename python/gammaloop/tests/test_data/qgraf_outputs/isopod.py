graphs = [] 
# this graph looks like an isopod if you draw it with round edges :)
# also I did not want to write triangle box box everywhere

graphs.append(
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
                "momentum": "N/A",
                "indices": (),
                "vertices": (1, 4)
            },
            2: {
                "name": "q2",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (1, 5)
            },
            3: {
                "name": "q3",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (4, 5)
            },
            4: {
                "name": "q4",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (4, 6)
            },
            5: {
                "name": "q5",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (5, 7)
            },
            6: {
                "name": "q6",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (6, 7)
            },
            7: {
                "name": "q7",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (6, 2)
            },
            8: {
                "name": "q8",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (7, 3)
            },
            9: {
                "name": "q9",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (2, 3)
            }
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
            1: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (1, 101, 2)
            },
            2: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (7, 9, 102)
            },
            3: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (8, 103, 9)
            },
            4: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (4, 1, 3)
            },
            5: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (5, 2, 3)
            },
            6: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"), 
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (6, 4, 7)
            },
            7: {
                "PDGs": (1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A"), 
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (6, 5, 8)
            }
        },
        "overall_factor": "1",
    }
)
