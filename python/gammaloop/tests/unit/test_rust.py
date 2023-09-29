from gammaloop.tests.common import run_rust_test


class TestRust:

    def test_rust_massless_scalar_triangle(self, compile_rust_tests, scalar_massless_triangle_export):
        assert run_rust_test(compile_rust_tests, scalar_massless_triangle_export,
                             'massless_scalar_triangle')
