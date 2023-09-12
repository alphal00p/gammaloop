graphs=[]


graphs.append(
{
"edges":{
(-1+4):{
    "name":"p1",
    "PDG": -11,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
(-3+4):{
    "name":"p2",
    "PDG": +11,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,1+4)
 },
(-2+4):{
    "name": "p3",
    "PDG": -11,
    "type": "out",
    "momentum": "p1",
    "indices": (-2+4+1,),
    "vertices":(2+4,-2+4+1)
 },
(-4+4):{
    "name": "p4",
    "PDG": +11,
    "type": "out",
    "momentum": "p2",
    "indices": (-4+4+1,),
    "vertices":(2+4,-4+4+1)
 },
(1+4):{
    "name":"q"+str(1),
    "PDG": 22,
    "type": "virtual",
    "momentum": "-p1-p2",
    "indices": (2+4+1,1+4+1,),
    "vertices":(3+4,1+4)
 },
(2+4):{
    "name":"q"+str(2),
    "PDG": 22,
    "type": "virtual",
    "momentum": "p1+p2",
    "indices": (4+4+1,3+4+1,),
    "vertices":(4+4,2+4)
 },
(3+4):{
    "name":"q"+str(3),
    "PDG": 1,
    "type": "virtual",
    "momentum": "-k1",
    "indices": (6+4+1,5+4+1,),
    "vertices":(5+4,3+4)
 },
(4+4):{
    "name":"q"+str(4),
    "PDG": 1,
    "type": "virtual",
    "momentum": "-k1+p1+p2",
    "indices": (8+4+1,7+4+1,),
    "vertices":(3+4,6+4)
 },
(5+4):{
    "name":"q"+str(5),
    "PDG": 1,
    "type": "virtual",
    "momentum": "k2",
    "indices": (10+4+1,9+4+1,),
    "vertices":(4+4,5+4)
 },
(6+4):{
    "name":"q"+str(6),
    "PDG": 1,
    "type": "virtual",
    "momentum": "k2+p1+p2",
    "indices": (12+4+1,11+4+1,),
    "vertices":(6+4,4+4)
 },
(7+4):{
    "name":"q"+str(7),
    "PDG": 21,
    "type": "virtual",
    "momentum": "-k1-k2",
    "indices": (14+4+1,13+4+1,),
    "vertices":(6+4,5+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (-11,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1,
    "edge_ids": (-1+4,)
 },
 -3+4+1:{
    "PDGs": (+11,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1,
    "edge_ids": (-3+4,)
 },
 -2+4+1:{
    "PDGs": (-11,),
    "momenta": ("p1"),
    "indices": (-2+4+1,),
    "vertex_id": -2,
    "edge_ids": (-2+4,)
 },
 -4+4+1:{
    "PDGs": (+11,),
    "momenta": ("p2"),
    "indices": (-4+4+1,),
    "vertex_id": -2,
    "edge_ids": (-4+4,)
 },
 1+4:{
    "PDGs": (-11,22,+11),
    "momenta": ("p1","-p1-p2","p2"),
    "indices": (-1+4+1,1+4+1,-3+4+1),
    "vertex_id": 0,
    "edge_ids": (-1+4,1+4,-3+4)
    },
 2+4:{
    "PDGs": (-11,22,+11),
    "momenta": ("-p2","p1+p2","-p1"),
    "indices": (-4+4+1,3+4+1,-2+4+1),
    "vertex_id": 0,
    "edge_ids": (-4+4,2+4,-2+4)
    },
 3+4:{
    "PDGs": (-1,22,1),
    "momenta": ("k1-p1-p2","p1+p2","-k1"),
    "indices": (8+4+1,2+4+1,5+4+1),
    "vertex_id": 0,
    "edge_ids": (4+4,1+4,3+4)
    },
 4+4:{
    "PDGs": (-1,22,1),
    "momenta": ("-k2","-p1-p2","k2+p1+p2"),
    "indices": (10+4+1,4+4+1,11+4+1),
    "vertex_id": 0,
    "edge_ids": (5+4,2+4,6+4)
    },
 5+4:{
    "PDGs": (-1,21,1),
    "momenta": ("k1","-k1-k2","k2"),
    "indices": (6+4+1,13+4+1,9+4+1),
    "vertex_id": 0,
    "edge_ids": (3+4,7+4,5+4)
    },
 6+4:{
    "PDGs": (-1,21,1),
    "momenta": ("-k2-p1-p2","k1+k2","-k1+p1+p2"),
    "indices": (12+4+1,14+4+1,7+4+1),
    "vertex_id": 0,
    "edge_ids": (6+4,7+4,4+4)
    }
},
"overall_factor": "-1"
}
)


graphs.append(
{
"edges":{
(-1+4):{
    "name":"p1",
    "PDG": -11,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
(-3+4):{
    "name":"p2",
    "PDG": +11,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,1+4)
 },
(-2+4):{
    "name": "p3",
    "PDG": -11,
    "type": "out",
    "momentum": "p1",
    "indices": (-2+4+1,),
    "vertices":(2+4,-2+4+1)
 },
(-4+4):{
    "name": "p4",
    "PDG": +11,
    "type": "out",
    "momentum": "p2",
    "indices": (-4+4+1,),
    "vertices":(2+4,-4+4+1)
 },
(1+4):{
    "name":"q"+str(1),
    "PDG": 22,
    "type": "virtual",
    "momentum": "-p1-p2",
    "indices": (2+4+1,1+4+1,),
    "vertices":(3+4,1+4)
 },
(2+4):{
    "name":"q"+str(2),
    "PDG": 22,
    "type": "virtual",
    "momentum": "p1+p2",
    "indices": (4+4+1,3+4+1,),
    "vertices":(4+4,2+4)
 },
(3+4):{
    "name":"q"+str(3),
    "PDG": 1,
    "type": "virtual",
    "momentum": "k1-p1-p2",
    "indices": (6+4+1,5+4+1,),
    "vertices":(5+4,3+4)
 },
(4+4):{
    "name":"q"+str(4),
    "PDG": 1,
    "type": "virtual",
    "momentum": "k1",
    "indices": (8+4+1,7+4+1,),
    "vertices":(3+4,5+4)
 },
(5+4):{
    "name":"q"+str(5),
    "PDG": 1,
    "type": "virtual",
    "momentum": "k2+p1+p2",
    "indices": (10+4+1,9+4+1,),
    "vertices":(6+4,4+4)
 },
(6+4):{
    "name":"q"+str(6),
    "PDG": 1,
    "type": "virtual",
    "momentum": "k2",
    "indices": (12+4+1,11+4+1,),
    "vertices":(4+4,6+4)
 },
(7+4):{
    "name":"q"+str(7),
    "PDG": 21,
    "type": "virtual",
    "momentum": "-p1-p2",
    "indices": (14+4+1,13+4+1,),
    "vertices":(6+4,5+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (-11,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1,
    "edge_ids": (-1+4,)
 },
 -3+4+1:{
    "PDGs": (+11,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1,
    "edge_ids": (-3+4,)
 },
 -2+4+1:{
    "PDGs": (-11,),
    "momenta": ("p1"),
    "indices": (-2+4+1,),
    "vertex_id": -2,
    "edge_ids": (-2+4,)
 },
 -4+4+1:{
    "PDGs": (+11,),
    "momenta": ("p2"),
    "indices": (-4+4+1,),
    "vertex_id": -2,
    "edge_ids": (-4+4,)
 },
 1+4:{
    "PDGs": (-11,22,+11),
    "momenta": ("p1","-p1-p2","p2"),
    "indices": (-1+4+1,1+4+1,-3+4+1),
    "vertex_id": 0,
    "edge_ids": (-1+4,1+4,-3+4)
    },
 2+4:{
    "PDGs": (-11,22,+11),
    "momenta": ("-p2","p1+p2","-p1"),
    "indices": (-4+4+1,3+4+1,-2+4+1),
    "vertex_id": 0,
    "edge_ids": (-4+4,2+4,-2+4)
    },
 3+4:{
    "PDGs": (-1,22,1),
    "momenta": ("-k1","p1+p2","k1-p1-p2"),
    "indices": (8+4+1,2+4+1,5+4+1),
    "vertex_id": 0,
    "edge_ids": (4+4,1+4,3+4)
    },
 4+4:{
    "PDGs": (-1,22,1),
    "momenta": ("-k2","-p1-p2","k2+p1+p2"),
    "indices": (12+4+1,4+4+1,9+4+1),
    "vertex_id": 0,
    "edge_ids": (6+4,2+4,5+4)
    },
 5+4:{
    "PDGs": (-1,21,1),
    "momenta": ("-k1+p1+p2","-p1-p2","k1"),
    "indices": (6+4+1,13+4+1,7+4+1),
    "vertex_id": 0,
    "edge_ids": (3+4,7+4,4+4)
    },
 6+4:{
    "PDGs": (-1,21,1),
    "momenta": ("-k2-p1-p2","p1+p2","k2"),
    "indices": (10+4+1,14+4+1,11+4+1),
    "vertex_id": 0,
    "edge_ids": (5+4,7+4,6+4)
    }
},
"overall_factor": "+1"
}
)


graphs.append(
{
"edges":{
(-1+4):{
    "name":"p1",
    "PDG": -11,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
(-3+4):{
    "name":"p2",
    "PDG": +11,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,1+4)
 },
(-2+4):{
    "name": "p3",
    "PDG": -11,
    "type": "out",
    "momentum": "p1",
    "indices": (-2+4+1,),
    "vertices":(2+4,-2+4+1)
 },
(-4+4):{
    "name": "p4",
    "PDG": +11,
    "type": "out",
    "momentum": "p2",
    "indices": (-4+4+1,),
    "vertices":(2+4,-4+4+1)
 },
(1+4):{
    "name":"q"+str(1),
    "PDG": 22,
    "type": "virtual",
    "momentum": "-p1-p2",
    "indices": (2+4+1,1+4+1,),
    "vertices":(3+4,1+4)
 },
(2+4):{
    "name":"q"+str(2),
    "PDG": 22,
    "type": "virtual",
    "momentum": "p1+p2",
    "indices": (4+4+1,3+4+1,),
    "vertices":(4+4,2+4)
 },
(3+4):{
    "name":"q"+str(3),
    "PDG": 1,
    "type": "virtual",
    "momentum": "-k1",
    "indices": (6+4+1,5+4+1,),
    "vertices":(4+4,3+4)
 },
(4+4):{
    "name":"q"+str(4),
    "PDG": 1,
    "type": "virtual",
    "momentum": "-k1+p1+p2",
    "indices": (8+4+1,7+4+1,),
    "vertices":(3+4,5+4)
 },
(5+4):{
    "name":"q"+str(5),
    "PDG": 1,
    "type": "virtual",
    "momentum": "-k1+p1+p2",
    "indices": (10+4+1,9+4+1,),
    "vertices":(6+4,4+4)
 },
(6+4):{
    "name":"q"+str(6),
    "PDG": 21,
    "type": "virtual",
    "momentum": "k1+k2-p1-p2",
    "indices": (12+4+1,11+4+1,),
    "vertices":(6+4,5+4)
 },
(7+4):{
    "name":"q"+str(7),
    "PDG": 1,
    "type": "virtual",
    "momentum": "k2",
    "indices": (14+4+1,13+4+1,),
    "vertices":(5+4,6+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (-11,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1,
    "edge_ids": (-1+4,)
 },
 -3+4+1:{
    "PDGs": (+11,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1,
    "edge_ids": (-3+4,)
 },
 -2+4+1:{
    "PDGs": (-11,),
    "momenta": ("p1"),
    "indices": (-2+4+1,),
    "vertex_id": -2,
    "edge_ids": (-2+4,)
 },
 -4+4+1:{
    "PDGs": (+11,),
    "momenta": ("p2"),
    "indices": (-4+4+1,),
    "vertex_id": -2,
    "edge_ids": (-4+4,)
 },
 1+4:{
    "PDGs": (-11,22,+11),
    "momenta": ("p1","-p1-p2","p2"),
    "indices": (-1+4+1,1+4+1,-3+4+1),
    "vertex_id": 0,
    "edge_ids": (-1+4,1+4,-3+4)
    },
 2+4:{
    "PDGs": (-11,22,+11),
    "momenta": ("-p2","p1+p2","-p1"),
    "indices": (-4+4+1,3+4+1,-2+4+1),
    "vertex_id": 0,
    "edge_ids": (-4+4,2+4,-2+4)
    },
 3+4:{
    "PDGs": (-1,22,1),
    "momenta": ("k1-p1-p2","p1+p2","-k1"),
    "indices": (8+4+1,2+4+1,5+4+1),
    "vertex_id": 0,
    "edge_ids": (4+4,1+4,3+4)
    },
 4+4:{
    "PDGs": (-1,22,1),
    "momenta": ("k1","-p1-p2","-k1+p1+p2"),
    "indices": (6+4+1,4+4+1,9+4+1),
    "vertex_id": 0,
    "edge_ids": (3+4,2+4,5+4)
    },
 5+4:{
    "PDGs": (-1,21,1),
    "momenta": ("-k2","k1+k2-p1-p2","-k1+p1+p2"),
    "indices": (14+4+1,11+4+1,7+4+1),
    "vertex_id": 0,
    "edge_ids": (7+4,6+4,4+4)
    },
 6+4:{
    "PDGs": (-1,21,1),
    "momenta": ("k1-p1-p2","-k1-k2+p1+p2","k2"),
    "indices": (10+4+1,12+4+1,13+4+1),
    "vertex_id": 0,
    "edge_ids": (5+4,6+4,7+4)
    }
},
"overall_factor": "-1"
}
)


graphs.append(
{
"edges":{
(-1+4):{
    "name":"p1",
    "PDG": -11,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
(-3+4):{
    "name":"p2",
    "PDG": +11,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,1+4)
 },
(-2+4):{
    "name": "p3",
    "PDG": -11,
    "type": "out",
    "momentum": "p1",
    "indices": (-2+4+1,),
    "vertices":(2+4,-2+4+1)
 },
(-4+4):{
    "name": "p4",
    "PDG": +11,
    "type": "out",
    "momentum": "p2",
    "indices": (-4+4+1,),
    "vertices":(2+4,-4+4+1)
 },
(1+4):{
    "name":"q"+str(1),
    "PDG": 22,
    "type": "virtual",
    "momentum": "-p1-p2",
    "indices": (2+4+1,1+4+1,),
    "vertices":(3+4,1+4)
 },
(2+4):{
    "name":"q"+str(2),
    "PDG": 22,
    "type": "virtual",
    "momentum": "p1+p2",
    "indices": (4+4+1,3+4+1,),
    "vertices":(4+4,2+4)
 },
(3+4):{
    "name":"q"+str(3),
    "PDG": 1,
    "type": "virtual",
    "momentum": "k1",
    "indices": (6+4+1,5+4+1,),
    "vertices":(3+4,4+4)
 },
(4+4):{
    "name":"q"+str(4),
    "PDG": 1,
    "type": "virtual",
    "momentum": "k1-p1-p2",
    "indices": (8+4+1,7+4+1,),
    "vertices":(5+4,3+4)
 },
(5+4):{
    "name":"q"+str(5),
    "PDG": 1,
    "type": "virtual",
    "momentum": "k1-p1-p2",
    "indices": (10+4+1,9+4+1,),
    "vertices":(4+4,6+4)
 },
(6+4):{
    "name":"q"+str(6),
    "PDG": 21,
    "type": "virtual",
    "momentum": "k1+k2-p1-p2",
    "indices": (12+4+1,11+4+1,),
    "vertices":(6+4,5+4)
 },
(7+4):{
    "name":"q"+str(7),
    "PDG": 1,
    "type": "virtual",
    "momentum": "-k2",
    "indices": (14+4+1,13+4+1,),
    "vertices":(6+4,5+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (-11,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1,
    "edge_ids": (-1+4,)
 },
 -3+4+1:{
    "PDGs": (+11,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1,
    "edge_ids": (-3+4,)
 },
 -2+4+1:{
    "PDGs": (-11,),
    "momenta": ("p1"),
    "indices": (-2+4+1,),
    "vertex_id": -2,
    "edge_ids": (-2+4,)
 },
 -4+4+1:{
    "PDGs": (+11,),
    "momenta": ("p2"),
    "indices": (-4+4+1,),
    "vertex_id": -2,
    "edge_ids": (-4+4,)
 },
 1+4:{
    "PDGs": (-11,22,+11),
    "momenta": ("p1","-p1-p2","p2"),
    "indices": (-1+4+1,1+4+1,-3+4+1),
    "vertex_id": 0,
    "edge_ids": (-1+4,1+4,-3+4)
    },
 2+4:{
    "PDGs": (-11,22,+11),
    "momenta": ("-p2","p1+p2","-p1"),
    "indices": (-4+4+1,3+4+1,-2+4+1),
    "vertex_id": 0,
    "edge_ids": (-4+4,2+4,-2+4)
    },
 3+4:{
    "PDGs": (-1,22,1),
    "momenta": ("-k1","p1+p2","k1-p1-p2"),
    "indices": (6+4+1,2+4+1,7+4+1),
    "vertex_id": 0,
    "edge_ids": (3+4,1+4,4+4)
    },
 4+4:{
    "PDGs": (-1,22,1),
    "momenta": ("-k1+p1+p2","-p1-p2","k1"),
    "indices": (10+4+1,4+4+1,5+4+1),
    "vertex_id": 0,
    "edge_ids": (5+4,2+4,3+4)
    },
 5+4:{
    "PDGs": (-1,21,1),
    "momenta": ("-k1+p1+p2","k1+k2-p1-p2","-k2"),
    "indices": (8+4+1,11+4+1,13+4+1),
    "vertex_id": 0,
    "edge_ids": (4+4,6+4,7+4)
    },
 6+4:{
    "PDGs": (-1,21,1),
    "momenta": ("k2","-k1-k2+p1+p2","k1-p1-p2"),
    "indices": (14+4+1,12+4+1,9+4+1),
    "vertex_id": 0,
    "edge_ids": (7+4,6+4,5+4)
    }
},
"overall_factor": "-1"
}
)


