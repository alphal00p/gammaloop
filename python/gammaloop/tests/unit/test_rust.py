import pytest
from pathlib import Path
from gammaloop.tests.common import run_rust_test


@pytest.mark.rust
class TestRust:

    def test_generate_rust_test_inputs(self,
                                       sm_model_yaml_file: Path,
                                       scalars_model_yaml_file: Path,
                                       massless_scalar_triangle_export: Path,
                                       scalar_fishnet_2x2_export: Path,
                                       scalar_sunrise_export: Path,
                                       scalar_fishnet_2x3_export: Path,
                                       scalar_cube_export: Path,
                                       scalar_bubble_export: Path,
                                       scalar_massless_box_export: Path,
                                       scalar_double_triangle_export: Path,
                                       scalar_mercedes_export: Path,
                                       scalar_triangle_box_export: Path,
                                       scalar_isopod_export: Path,
                                       scalar_raised_triangle_export: Path,
                                       scalar_hexagon_export: Path,
                                       lbl_box_export: Path,
                                       scalar_ltd_topology_c_export: Path,
                                       scalar_massless_pentabox_export: Path,
                                       scalar_massless_3l_pentabox_export: Path
                                       ):
        assert True

    def test_rust_massless_scalar_triangle(self, compile_rust_tests, scalars_model_yaml_file: Path, massless_scalar_triangle_export: Path):
        assert run_rust_test(compile_rust_tests, massless_scalar_triangle_export,
                             'massless_scalar_triangle')

    def test_rust_scalar_fishnet_2x2(self, compile_rust_tests: Path, scalar_fishnet_2x2_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_fishnet_2x2_export,
                             'scalar_fishnet_2x2')

    def test_rust_scalar_sunrise(self, compile_rust_tests: Path, scalar_sunrise_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_sunrise_export,
                             'scalar_sunrise')

    @pytest.mark.slow
    def test_rust_scalar_fishnet_2x3(self, compile_rust_tests: Path, scalar_fishnet_2x3_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_fishnet_2x3_export,
                             'scalar_fishnet_2x3')

    def test_rust_scalar_cube(self, compile_rust_tests: Path, scalar_cube_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_cube_export,
                             'scalar_cube')

    def test_rust_scalar_bubble(self, compile_rust_tests: Path, scalar_bubble_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_bubble_export,
                             'scalar_bubble')

    def test_rust_massless_scalar_box(self, compile_rust_tests: Path, scalar_massless_box_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_massless_box_export,
                             'massless_scalar_box')

    def test_rust_scalar_double_triangle(self, compile_rust_tests: Path, scalar_double_triangle_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_double_triangle_export,
                             'scalar_double_triangle')

    def test_rust_scalar_mercedes(self, compile_rust_tests: Path, scalar_mercedes_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_mercedes_export,
                             'scalar_mercedes')

    def test_rust_scalar_triangle_box(self, compile_rust_tests: Path, scalar_triangle_box_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_triangle_box_export,
                             'scalar_triangle_box')

    def test_rust_scalar_isopod(self, compile_rust_tests: Path, scalar_isopod_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_isopod_export,
                             'scalar_isopod')

    def test_rust_raised_triangle(self, compile_rust_tests: Path, scalar_raised_triangle_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_raised_triangle_export,
                             'raised_triangle')

    def test_rust_hexagon(self, compile_rust_tests: Path, scalar_hexagon_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_hexagon_export,
                             'hexagon')

    def test_rust_lbl_box(self, compile_rust_tests: Path, lbl_box_export: Path):
        assert run_rust_test(compile_rust_tests, lbl_box_export,
                             'lbl_box')

    def test_rust_topology_c(self, compile_rust_tests: Path, scalar_ltd_topology_c_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_ltd_topology_c_export,
                             'topology_c')

    def test_rust_massless_pentabox(self, compile_rust_tests: Path, scalar_massless_pentabox_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_massless_pentabox_export,
                             'massless_pentabox')

    def test_rust_massless_3l_pentabox(self, compile_rust_tests: Path, scalar_massless_3l_pentabox_export: Path):
        assert run_rust_test(compile_rust_tests, scalar_massless_3l_pentabox_export,
                             'massless_3l_pentabox')
