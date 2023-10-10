import pytest
from gammaloop.tests.common import run_rust_test


class TestRust:

    def test_rust_massless_scalar_triangle(self, compile_rust_tests, scalar_massless_triangle_export):
        assert run_rust_test(compile_rust_tests, scalar_massless_triangle_export,
                             'massless_scalar_triangle')

    def test_rust_scalar_fishnet_2x2(self, compile_rust_tests, scalar_fishnet_2x2_export):
        assert run_rust_test(compile_rust_tests, scalar_fishnet_2x2_export,
                             'scalar_fishnet_2x2')

    @pytest.mark.slow
    def test_rust_scalar_fishnet_2x3(self, compile_rust_tests, scalar_fishnet_2x3_export):
        assert run_rust_test(compile_rust_tests, scalar_fishnet_2x3_export,
                             'scalar_fishnet_2x3')

    def test_rust_scalar_cube(self, compile_rust_tests, scalar_cube_export):
        assert run_rust_test(compile_rust_tests, scalar_cube_export,
                             'scalar_cube')
