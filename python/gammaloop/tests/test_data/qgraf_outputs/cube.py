graphs = []
# To be imported after having the gammaloop UFO model named "scalars"

# edge name convention XY where XY are the vertices connected by the edge
# massless cube diagram
graphs.append(
    {
        "edges": {
            1001: {
                "name": "p1",
                "PDG": 1000,
                "type": "in",
                "momentum": "p1",
                "indices": (),
                "vertices": (101, 1)
            },
            1002: {
                "name": "p2",
                "PDG": 1000,
                "type": "in",
                "momentum": "p2",
                "indices": (),
                "vertices": (102, 2)
            },
            1003: {
                "name": "p3",
                "PDG": 1000,
                "type": "in",
                "momentum": "p3",
                "indices": (),
                "vertices": (103, 3)
            },
            1004: {
                "name": "p4",
                "PDG": 1000,
                "type": "in",
                "momentum": "p4",
                "indices": (),
                "vertices": (104, 4)
            },
            1005: {
                "name": "p5",
                "PDG": 1000,
                "type": "in",
                "momentum": "p5",
                "indices": (),
                "vertices": (105, 5)
            },
            1006: {
                "name": "p6",
                "PDG": 1000,
                "type": "in",
                "momentum": "p6",
                "indices": (),
                "vertices": (106, 6)
            },
            1007: {
                "name": "p7",
                "PDG": 1000,
                "type": "in",
                "momentum": "p7",
                "indices": (),
                "vertices": (107, 7)
            },
            1008: {
                "name": "p8",
                "PDG": 1000,
                "type": "in",
                "momentum": "p8",
                "indices": (),
                "vertices": (108, 8)
            },
            12: {
                "name": "q12",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (1, 2)
            },
            24: {
                "name": "q24",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (2, 4)
            },
            43: {
                "name": "q21",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (2, 1)
            },
            31: {
                "name": "q31",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (3, 1)
            },
            56: {
                "name": "q56",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (5, 6)
            },
            68: {
                "name": "q68",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (6, 8)
            },
            87: {
                "name": "q78",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (7, 8)
            },
            75: {
                "name": "q75",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (7, 5)
            },
            15: {
                "name": "q15",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (1, 5)
            },
            26: {
                "name": "q26",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (2, 6)
            },
            37: {
                "name": "q37",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (3, 7)
            },
            48: {
                "name": "q48",
                "PDG": 1000,
                "type": "virtual",
                "momentum": "N/A",
                "indices": (),
                "vertices": (4, 8)
            },
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
            103: {
                "PDGs": (1000,),
                "momenta": ("p3",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (1003,)
            },
            104: {
                "PDGs": (1000,),
                "momenta": ("p4",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (1004,)
            },
            105: {
                "PDGs": (1000,),
                "momenta": ("p5",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (1005,)
            },
            106: {
                "PDGs": (1000,),
                "momenta": ("p6",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (1006,)
            },
            107: {
                "PDGs": (1000,),
                "momenta": ("p7",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (1007,)
            },
            108: {
                "PDGs": (1000,),
                "momenta": ("p8",),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (1008,)
            },
            1: {
                "PDGs": (1000, 1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (1001, 12, 15, 31),
            },
            2: {
                "PDGs": (1000, 1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (12, 1002, 24, 26),
            },
            3: {
                "PDGs": (1000, 1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (31, 43, 1003, 37),
            },
            4: {
                "PDGs": (1000, 1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (24, 1004, 48, 43),
            },
            5: {
                "PDGs": (1000, 1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (15, 75, 1005, 56),
            },
            6: {
                "PDGs": (1000, 1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (26, 56, 1006, 68),
            },
            7: {
                "PDGs": (1000, 1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (37, 75, 1007, 87),
            },
            8: {
                "PDGs": (1000, 1000, 1000, 1000),
                "momenta": ("N/A", "N/A", "N/A", "N/A"),
                "indices": (),
                "vertex_id": 0,
                "edge_ids": (48, 68, 87, 1008),
            },
        },
        "overall_factor": "1"
    }
)
