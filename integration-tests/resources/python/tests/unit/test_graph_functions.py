import gammaloop.misc.utils as utils


class TestGraphFunctions:
    def test_graph_fishnet_2x2(self):
        fishnet_2x2_edge_map = [(0, 4), (1, 10), (6, 2), (12, 3), (4, 5), (5, 6), (7, 8), (
            8, 9), (10, 11), (11, 12), (4, 7), (5, 8), (6, 9), (7, 10), (8, 11), (9, 12)]
        assert len(utils.find_longest_cycle(
            fishnet_2x2_edge_map)) == 8  # type: ignore

    def test_graph_triangle(self):
        triangle_edge_map = [
            (101, 0), (0, 1), (103, 2), (2, 1), (2, 0), (102, 1)]
        assert len(utils.find_longest_cycle(
            triangle_edge_map)) == 3  # type: ignore

    def test_cube(self):
        cube_edge_map = [(0, 8), (1, 9), (2, 10), (3, 11), (4, 12), (5, 13), (6, 14), (7, 15), (8, 9), (
            9, 11), (9, 8), (10, 8), (12, 13), (13, 15), (14, 15), (14, 12), (8, 12), (9, 13), (10, 14), (11, 15)]
        assert len(utils.find_longest_cycle(
            cube_edge_map)) == 8  # type: ignore
