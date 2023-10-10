import pytest
from gammaloop.tests.common import get_gamma_loop_interpreter, RESOURCES_PATH, pjoin, run_drawing
from gammaloop.interface.gammaloop_interface import CommandList


class TestShellCommand:
    def test_load_sm_from_ufo(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string("!ls"))

    def test_set_drawing_config(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "set drawing.combined_graphs_pdf_grid_shape [1,1]"))
        gloop.run(CommandList.from_string(
            "set drawing.feynmp.show_edge_labels True"))
        assert gloop.config['drawing']['combined_graphs_pdf_grid_shape'] == [
            1, 1]
        assert gloop.config['drawing']['feynmp']['show_edge_labels'] == True


class TestLoadModel:

    def test_load_sm_from_ufo(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string("import_model sm --format ufo"))

        assert len(gloop.model.particles) == 43
        assert len(gloop.model.lorentz_structures) == 22
        assert len(gloop.model.couplings) == 108
        assert len(gloop.model.parameters) == 72
        assert len(gloop.model.vertex_rules) == 153
        assert len(gloop.model.orders) == 2

    # This test uses a session-wide fixture defined in conftest.py
    def test_load_sm_from_yaml(self, sm_model_yaml_file):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"import_model {sm_model_yaml_file} --format yaml"))

        assert len(gloop.model.particles) == 43
        assert len(gloop.model.lorentz_structures) == 22
        assert len(gloop.model.couplings) == 108
        assert len(gloop.model.parameters) == 72
        assert len(gloop.model.vertex_rules) == 153
        assert len(gloop.model.orders) == 2

    def test_load_scalars(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string("import_model scalars"))

        assert len(gloop.model.particles) == 3
        assert len(gloop.model.lorentz_structures) == 1
        assert len(gloop.model.couplings) == 1
        assert len(gloop.model.parameters) == 12
        assert len(gloop.model.vertex_rules) == 10
        assert len(gloop.model.orders) == 2


class TestLoadQGraph:

    def test_epem_a_ddx_nlo(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"import_model sm; import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'epem_a_ddx_NLO.py')} -f qgraph --no_compile"))
        assert len(gloop.cross_sections) == 1
        assert len(gloop.cross_sections[0].supergraphs) == 4

    def test_massless_scalar_triangle(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"import_model scalars; import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'massless_triangle.py')} -f qgraph --no_compile"))
        assert len(gloop.amplitudes) == 1
        assert len(gloop.amplitudes[0].amplitude_graphs) == 1

    def test_fishnet_2x2(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"import_model scalars; import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'fishnet_2x2.py')} -f qgraph --no_compile"))
        assert len(gloop.amplitudes) == 1
        assert len(gloop.amplitudes[0].amplitude_graphs) == 1
        assert len(gloop.amplitudes[0].amplitude_graphs[0].graph.edges) == 16
        assert len(
            gloop.amplitudes[0].amplitude_graphs[0].graph.vertices) == 13
        assert len(
            gloop.amplitudes[0].amplitude_graphs[0].graph.external_connections) == 4

    def test_fishnet_2x3(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"import_model scalars; import_graphs {pjoin(RESOURCES_PATH, 'qgraf_outputs', 'fishnet_2x3.py')} -f qgraph --no_compile"))
        assert len(gloop.amplitudes) == 1
        assert len(gloop.amplitudes[0].amplitude_graphs) == 1
        assert len(gloop.amplitudes[0].amplitude_graphs[0].graph.edges) == 21
        assert len(
            gloop.amplitudes[0].amplitude_graphs[0].graph.vertices) == 16
        assert len(
            gloop.amplitudes[0].amplitude_graphs[0].graph.external_connections) == 4


class TestMasslessScalarTriangleAmplitude:

    # This test uses a session-wide fixture defined in conftest.py
    def test_info(self, scalar_massless_triangle_export: str):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {scalar_massless_triangle_export}"))
        assert gloop.model.name == 'scalars'
        assert gloop.get_model_from_rust_worker().name == 'scalars'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 0
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 1
            assert len(amplitudes[0].amplitude_graphs) == 1
            assert amplitudes[0].name == 'massless_triangle'
        gloop.run(CommandList.from_string("info"))

    def test_drawing(self, scalar_massless_triangle_export):
        assert run_drawing(pjoin(scalar_massless_triangle_export, 'sources',
                           'amplitudes', 'massless_triangle', 'drawings'))


class TestScalarFishnet2x2:

    def test_info(self, scalar_fishnet_2x2_export):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {scalar_fishnet_2x2_export}"))
        assert gloop.model.name == 'scalars'
        assert gloop.get_model_from_rust_worker().name == 'scalars'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 0
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 1
            assert len(amplitudes[0].amplitude_graphs) == 1
            assert amplitudes[0].name == 'fishnet_2x2'
        gloop.run(CommandList.from_string("info"))

    def test_drawing(self, scalar_fishnet_2x2_export):
        assert run_drawing(pjoin(scalar_fishnet_2x2_export, 'sources',
                           'amplitudes', 'fishnet_2x2', 'drawings'))


@pytest.mark.slow
class TestScalarFishnet2x3:

    def test_info(self, scalar_fishnet_2x3_export):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {scalar_fishnet_2x3_export}"))
        assert gloop.model.name == 'scalars'
        assert gloop.get_model_from_rust_worker().name == 'scalars'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 0
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 1
            assert len(amplitudes[0].amplitude_graphs) == 1
            assert amplitudes[0].name == 'fishnet_2x3'
        gloop.run(CommandList.from_string("info"))

    def test_drawing(self, scalar_fishnet_2x3_export):
        assert run_drawing(pjoin(scalar_fishnet_2x3_export, 'sources',
                           'amplitudes', 'fishnet_2x3', 'drawings'))


class TestScalarCube:

    def test_info(self, scalar_cube_export):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {scalar_cube_export}"))
        assert gloop.model.name == 'scalars'
        assert gloop.get_model_from_rust_worker().name == 'scalars'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 0
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 1
            assert len(amplitudes[0].amplitude_graphs) == 1
            assert amplitudes[0].name == 'cube'
        gloop.run(CommandList.from_string("info"))

    def test_drawing(self, scalar_cube_export):
        assert run_drawing(pjoin(scalar_cube_export, 'sources',
                           'amplitudes', 'cube', 'drawings'))


class TestEpEmADdxNLOCrossSection:

    # This test uses a session-wide fixture defined in conftest.py
    def test_info(self, epem_a_ddx_nlo_export):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {epem_a_ddx_nlo_export}"))
        assert gloop.model.name == 'sm'
        assert gloop.get_model_from_rust_worker().name == 'sm'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 1
            assert len(cross_sections[0].supergraphs) == 4
            assert cross_sections[0].name == 'epem_a_ddx_NLO'
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 0
        gloop.run(CommandList.from_string("info"))

    def test_drawing(self, epem_a_ddx_nlo_export):
        assert run_drawing(pjoin(epem_a_ddx_nlo_export, 'sources',
                           'cross_sections', 'epem_a_ddx_NLO', 'drawings'))
