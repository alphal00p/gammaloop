import pytest
import os
from gammaloop.interface.gammaloop_interface import CommandList, GammaLoopConfiguration, GammaLoopError, GammaLoop
from gammaloop.misc.common import load_configuration, GL_PATH
from gammaloop.tests.common import get_gamma_loop_interpreter, get_gamma_loop_interpreter_no_compilation, RESOURCES_PATH, pjoin, run_drawing
from pathlib import Path


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


class TestGammaLoopConfigurationDefaults:

    def test_default_config(self):
        default_configuration = GammaLoopConfiguration(path='')
        rust_run_config = load_configuration(os.path.join(
            GL_PATH, 'data', 'run_cards', 'rust_run_config.yaml'))
        gammaloop_config = load_configuration(os.path.join(
            GL_PATH, 'data', 'config', 'gammaloop_config.yaml'))

        rust_run_default_config = default_configuration._config['run_settings']
        del default_configuration._config['run_settings']
        # from pprint import pprint
        # pprint(rust_run_default_config)
        # pprint(rust_run_config)
        assert rust_run_default_config == rust_run_config
        # from pprint import pprint
        # pprint(default_configuration._config)
        # pprint(gammaloop_config)
        assert default_configuration._config == gammaloop_config


class TestLoadModel:

    def test_load_full_sm_from_ufo(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string("import_model sm-full --format ufo"))

        assert len(gloop.model.particles) == 43
        assert len(gloop.model.lorentz_structures) == 22
        assert len(gloop.model.couplings) == 108
        assert len(gloop.model.parameters) == 72
        assert len(gloop.model.vertex_rules) == 153
        assert len(gloop.model.orders) == 2

    def test_load_sm_from_ufo(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string("import_model sm --format ufo"))

        assert len(gloop.model.particles) == 43
        assert len(gloop.model.lorentz_structures) == 22
        assert len(gloop.model.couplings) == 72
        assert len(gloop.model.parameters) == 72
        assert len(gloop.model.vertex_rules) == 119
        assert len(gloop.model.orders) == 2

    def test_load_full_loop_sm_from_ufo(self):
        loop_sm_model_path = os.path.join(RESOURCES_PATH, 'models', 'loop_sm')
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"import_model {loop_sm_model_path}-full --format ufo"))

        assert len(gloop.model.particles) == 35
        assert len(gloop.model.lorentz_structures) == 40
        assert len(gloop.model.couplings) == 266
        assert len(gloop.model.parameters) == 146
        assert len(gloop.model.vertex_rules) == 130
        assert len(gloop.model.orders) == 2

    def test_load_loop_sm_from_ufo(self):
        loop_sm_model_path = os.path.join(RESOURCES_PATH, 'models', 'loop_sm')
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"import_model {loop_sm_model_path} --format ufo"))

        assert len(gloop.model.particles) == 35
        assert len(gloop.model.lorentz_structures) == 40
        assert len(gloop.model.couplings) == 163
        assert len(gloop.model.parameters) == 146
        assert len(gloop.model.vertex_rules) == 96
        assert len(gloop.model.orders) == 2

    # This test uses a session-wide fixture defined in conftest.py
    def test_load_sm_from_yaml(self, sm_model_yaml_file: str):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"import_model {sm_model_yaml_file} --format yaml"))

        assert len(gloop.model.particles) == 43
        assert len(gloop.model.lorentz_structures) == 22
        assert len(gloop.model.couplings) == 72
        assert len(gloop.model.parameters) == 72
        assert len(gloop.model.vertex_rules) == 119
        assert len(gloop.model.orders) == 2

    def test_load_scalars(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string("import_model scalars"))

        assert len(gloop.model.particles) == 3
        assert len(gloop.model.lorentz_structures) == 2
        assert len(gloop.model.couplings) == 1
        assert len(gloop.model.parameters) == 12
        assert len(gloop.model.vertex_rules) == 25
        assert len(gloop.model.orders) == 2


class TestProcessGeneration:

    @staticmethod
    def run_tests(gloop: GammaLoop, tests: list[tuple[str, int]]):
        for test, expected_graph_number in tests:
            try:
                gloop.run(CommandList.from_string(
                    f"generate {test} --clear_existing_processes"))
            except GammaLoopError as e:
                pass
            n_graphs = (
                0 if len(gloop.amplitudes) == 0 else len(
                    gloop.amplitudes[0].amplitude_graphs)
            ) + (
                0 if len(gloop.cross_sections) == 0 else len(
                    gloop.cross_sections[0].supergraphs)
            )
            assert n_graphs == expected_graph_number, f"For process: '{
                test}' | Expected {expected_graph_number} graphs, got {n_graphs}"
            gloop.amplitudes.clear()
            gloop.cross_sections.clear()

    def test_generate_sm_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            # Only d, g and a as particle contents
            ('a > d d~ [{{1}}] | d g a -num_grouping only_detect_zeroes', 1),
            ('a > d d~ g [{{2}}] | d g a -num_grouping only_detect_zeroes', 3),
            ('a > d d~ g [{{2}}] | d g a --symmetrize_left_right_states', 2),
            ('a > d d~ z [{{2}}] | d g a -num_grouping only_detect_zeroes', 0),
            # Full particle contents
            ('a > d d~ [{{1}}] --symmetrize_left_right_states', 1),
            ('a > d d~ [{{2}}] --symmetrize_left_right_states', 9),
            # Only 1-flavour pure QCD corrections
            ('a > d d~ | d g ghG a QED^2=2 [{{1}}] --symmetrize_left_right_states', 1),
            ('a > d d~ | d g ghG a QED^2=2 [{{2}} QCD=1] --symmetrize_left_right_states', 2),
            ('a > d d~ | d g ghG a QED^2=2 [{{3}} QCD=2] --symmetrize_left_right_states', 17),
            ('a > d d~ | d g ghG a QED^2=2 [{{4}} QCD=3] --symmetrize_left_right_states', 258),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_sm_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            # Full particle contents
            ('a > d d~ [{{3}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 658),
            # Only 1-flavour pure QCD corrections
            ('a > d d~ | d g ghG a QED^2=2 [{{5}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 5016),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_sm_full_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm-full"))
        tests = [
            # Full particle contents
            ('a > d d~ [{{1}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 1),
            ('a > d d~ [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 25),
            ('a > d d~ [{{3}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 4060),
            ('a > d d~ [{{3}}] --symmetrize_left_right_states -num_grouping no_grouping', 4142),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_sm_h_n_j(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            # Full particle contents
            ('h > g g | h g b t ghg [{{3}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 6),
            ('h > g g | h g b t ghg [{{3}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 3),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_sm_h_n_j(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            # Full particle contents
            ('h > g g g | h g b t ghg [{{4}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 56),
            # FIX: Weird, this process keeps randomly yielding 36, 37 or 38
            # ('h > g g g | h g b t ghg [{{4}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 36),
            ('h > g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 2),
            ('h > g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 36),
            ('h > g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 2),
            ('h > g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 40),
            ('h > g g | h g b t ghg [{3}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 1169),
            ('h > g g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 6),
            ('h > g g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 192),
            ('h > g g g | h g b t ghg [{3}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 7142),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_epem_a_qqh(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            ('e+ e- > d d~ | d a e- ghg g QED^2=4 [{{2}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 2),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2=6 [{{2}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 2),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2=6 [{{3}} QCD=1] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 16),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2=6 [{{4}} QCD=2] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 314),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2=6 [{{2}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 12),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2=6 [{{3}} QCD=2] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 100),
            # Too slow sadly
            # ('e+ e- > b b~ h | d b h a e- ghg g z QED^2=6 [{{4}} QCD=4] --symmetrize_left_right_states -num_grouping only_detect_zeroes', ?),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_a_qqh(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            ('a > d d~ | d a ghg g QED^2=2 [{{2}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 2),
            ('a > b b~ h | d b h a ghg g QED^2=4 [{{2}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 2),
            ('a > b b~ h | d b h a ghg g QED^2=4 [{{3}} QCD=1] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 16),
            ('a > b b~ h | d b h a ghg g QED^2=4 [{{4}} QCD=2] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 314),
            # Too slow sadly
            # ('e+ e- > b b~ h | d b h a e- ghg g z QED^2=6 [{{4}} QCD=4] --symmetrize_left_right_states -num_grouping only_detect_zeroes', ?),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_amplitude_1l_sm_jets(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        tests = [
            ('d d~ > d d~ | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 20),
            ('d d~ > d d~ g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 192),
            ('d d~ > g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 384),
            ('g g > g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 1139),
            ('d d~ > d d~ g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 2284),
            ('d d~ > u u~ | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 10),
            ('d d~ > u u~ g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 96),
            ('d d~ > u u~ g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 1142),
            ('d d~ > u u~ d d~ | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 414),
            ('d d~ > d d~ d d~ | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 1242),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    # Currently bugged with spenso network parsing for "-num_grouping group_identical_graphs_up_to_sign"
    def test_generate_amplitude_1l_sm_jets_with_grouping(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        tests = [
            ('d d~ > d d~ | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 18),
            ('d d~ > u u~ | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 9),
            # Grouping diagrams with spenso-expand is way too slow with current approach
            # ('d d~ > u u~ g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 88),
            # ('d d~ > d d~ g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 176),
            # ('d d~ > g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 341),
            # ('g g > g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 905),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_amplitude_1l_sm_jets(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        tests = [
            ('d d~ > u u~ d d~ g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 5424),
            # A bit too slow, let's skip for now
            # ('g g > g g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 14875),
            # ('d d~ > d d~ g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 32074),
            # ('d d~ > u u~ g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 16037),
            # ('d d~ > d d~ d d~ g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 16272),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    # A bit too slow, let's skip for now
    def very_slow_generate_amplitude_1l_sm_jets(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        tests = [
            ('g g > g g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 14875),
            ('d d~ > d d~ g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 32074),
            ('d d~ > u u~ g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 16037),
            ('d d~ > d d~ d d~ g | u d g ghg a QED=0 [QCD=1] -a -num_grouping only_detect_zeroes', 16272),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    # Grouping diagrams with spenso-expand is way too slow with current approach
    @pytest.mark.slow
    def no_test_slow_generate_amplitude_1l_sm_jets_with_grouping(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        tests = [
            ('d d~ > d d~ g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 2090),
            ('d d~ > u u~ g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 1045),
            ('d d~ > u u~ d d~ | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 380),
            ('d d~ > d d~ d d~ | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 1140),
            ('g g > g g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 11850),
            ('d d~ > d d~ g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 29210),
            ('d d~ > u u~ g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 14605),
            ('d d~ > d d~ d d~ g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 15030),
            ('d d~ > u u~ d d~ g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 5010),
        ]
        TestProcessGeneration.run_tests(gloop, tests)


class TestLoadQGraph:

    def test_epem_a_ddx_nlo(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"import_model sm; import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'epem_a_ddx_NLO.dot')} --no_compile"))
        assert len(gloop.cross_sections) == 1
        assert len(gloop.cross_sections[0].supergraphs) == 4

    def test_massless_scalar_triangle(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"import_model scalars; import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'massless_triangle.dot')} --no_compile"))
        assert len(gloop.amplitudes) == 1
        assert len(gloop.amplitudes[0].amplitude_graphs) == 1

    def test_fishnet_2x2(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"import_model scalars; import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'fishnet_2x2.dot')} --no_compile"))
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
            f"import_model scalars; import_graphs {pjoin(RESOURCES_PATH, 'graph_inputs', 'fishnet_2x3.dot')} --no_compile"))
        assert len(gloop.amplitudes) == 1
        assert len(gloop.amplitudes[0].amplitude_graphs) == 1
        assert len(gloop.amplitudes[0].amplitude_graphs[0].graph.edges) == 21
        assert len(
            gloop.amplitudes[0].amplitude_graphs[0].graph.vertices) == 16
        assert len(
            gloop.amplitudes[0].amplitude_graphs[0].graph.external_connections) == 4


class TestMasslessScalarTriangleAmplitude:

    # This test uses a session-wide fixture defined in conftest.py
    def test_info(self, massless_scalar_triangle_export: Path):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {massless_scalar_triangle_export}"))
        assert gloop.model.name == 'scalars'
        assert gloop.get_model_from_rust_worker().name == 'scalars'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 0
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 1
            assert len(amplitudes[0].amplitude_graphs) == 1
            assert amplitudes[0].name == 'massless_triangle'
        gloop.run(CommandList.from_string("info"))

    @pytest.mark.drawing
    def test_drawing(self, massless_scalar_triangle_export: str):
        assert run_drawing(pjoin(massless_scalar_triangle_export, 'sources',
                           'amplitudes', 'massless_triangle', 'drawings'))


class TestMasslessScalarBoxAmplitude:

    # This test uses a session-wide fixture defined in conftest.py
    def test_info(self, scalar_massless_box_export: Path):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {scalar_massless_box_export}"))
        assert gloop.model.name == 'scalars'
        assert gloop.get_model_from_rust_worker().name == 'scalars'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 0
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 1
            assert len(amplitudes[0].amplitude_graphs) == 1
            assert amplitudes[0].name == 'massless_box'
        gloop.run(CommandList.from_string("info"))

    @pytest.mark.drawing
    def test_drawing(self, scalar_massless_box_export: str):
        assert run_drawing(pjoin(scalar_massless_box_export, 'sources',
                           'amplitudes', 'massless_box', 'drawings'))


class TestScalarFishnet2x2:

    def test_info(self, scalar_fishnet_2x2_export: Path):
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

    @pytest.mark.drawing
    def test_drawing(self, scalar_fishnet_2x2_export: Path):
        assert run_drawing(pjoin(scalar_fishnet_2x2_export, 'sources',
                           'amplitudes', 'fishnet_2x2', 'drawings'))


@pytest.mark.slow
class TestScalarFishnet2x3:

    def test_info(self, scalar_fishnet_2x3_export: Path):
        gloop = get_gamma_loop_interpreter_no_compilation()
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


class TestScalarCube:

    def test_info(self, scalar_cube_export: Path):
        gloop = get_gamma_loop_interpreter()

        command_list = CommandList.from_string(
            "set externals.data.momenta [[1.,3.,4.,5.],[1.,6.,7.,8.],[1.,9.,10.,11.],[1.,12.,13.,14.],[1.,15.,16.,17.],[1.,18.,19.,20.],[1.,21.,22.,23.]]")
        command_list.add_command(f"launch {scalar_cube_export}")

        gloop.run(command_list)
        assert gloop.model.name == 'scalars'
        assert gloop.get_model_from_rust_worker().name == 'scalars'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 0
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 1
            assert len(amplitudes[0].amplitude_graphs) == 1
            assert amplitudes[0].name == 'cube'
        gloop.run(CommandList.from_string("info"))

    @pytest.mark.drawing
    def test_drawing(self, scalar_cube_export: Path):
        assert run_drawing(pjoin(scalar_cube_export, 'sources',
                           'amplitudes', 'cube', 'drawings'))


class TestEpEmADdxNLOCrossSection:

    # This test uses a session-wide fixture defined in conftest.py
    def NO_TEST_YET_test_info(self, epem_a_ddx_nlo_export: Path):
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

    @pytest.mark.drawing
    def test_drawing(self, epem_a_ddx_nlo_export: Path):
        assert run_drawing(pjoin(epem_a_ddx_nlo_export, 'sources',
                           'cross_sections', 'epem_a_ddx_NLO', 'drawings'))

    def NO_TEST_YET_test_info_massive(self, massive_epem_a_ddx_nlo_export: Path):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {massive_epem_a_ddx_nlo_export}"))
        assert gloop.model.name == 'sm'
        assert gloop.get_model_from_rust_worker().name == 'sm'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 1
            assert len(cross_sections[0].supergraphs) == 4
            assert cross_sections[0].name == 'epem_a_ddx_NLO'
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 0
        gloop.run(CommandList.from_string("info"))

    @pytest.mark.drawing
    def test_drawing_massive(self, massive_epem_a_ddx_nlo_export: Path):
        assert run_drawing(pjoin(massive_epem_a_ddx_nlo_export, 'sources',
                           'cross_sections', 'epem_a_ddx_NLO', 'drawings'))


class TestScalarBubble:

    def test_info(self, scalar_bubble_export: Path):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {scalar_bubble_export}"))
        assert gloop.model.name == 'scalars'
        assert gloop.get_model_from_rust_worker().name == 'scalars'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 0
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 1
            assert len(amplitudes[0].amplitude_graphs) == 1
            assert amplitudes[0].name == 'bubble'
        gloop.run(CommandList.from_string("info"))

    @pytest.mark.drawing
    def test_drawing(self, scalar_bubble_export: Path):
        assert run_drawing(pjoin(scalar_bubble_export, 'sources',
                           'amplitudes', 'bubble', 'drawings'))


class TestScalarDoubleTriangle:

    def test_info(self, scalar_double_triangle_export: Path):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {scalar_double_triangle_export}"))
        assert gloop.model.name == 'scalars'
        assert gloop.get_model_from_rust_worker().name == 'scalars'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 0
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 1
            assert len(amplitudes[0].amplitude_graphs) == 1
            assert amplitudes[0].name == 'double_triangle'
        gloop.run(CommandList.from_string("info"))

    @pytest.mark.drawing
    def test_drawing(self, scalar_double_triangle_export: Path):
        assert run_drawing(pjoin(scalar_double_triangle_export, 'sources',
                           'amplitudes', 'double_triangle', 'drawings'))


class TestScalarMercedes:

    def test_info(self, scalar_mercedes_export: Path):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {scalar_mercedes_export}"))
        assert gloop.model.name == 'scalars'
        assert gloop.get_model_from_rust_worker().name == 'scalars'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 0
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 1
            assert len(amplitudes[0].amplitude_graphs) == 1
            assert amplitudes[0].name == 'mercedes'
        gloop.run(CommandList.from_string("info"))

    @pytest.mark.drawing
    def test_drawing(self, scalar_mercedes_export: Path):
        assert run_drawing(pjoin(scalar_mercedes_export, 'sources',
                           'amplitudes', 'mercedes', 'drawings'))


class TestScalarTriangleBox:

    def test_info(self, scalar_triangle_box_export: Path):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {scalar_triangle_box_export}"))
        assert gloop.model.name == 'scalars'
        assert gloop.get_model_from_rust_worker().name == 'scalars'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 0
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 1
            assert len(amplitudes[0].amplitude_graphs) == 1
            assert amplitudes[0].name == 'triangle_box'
        gloop.run(CommandList.from_string("info"))

    @pytest.mark.drawing
    def test_drawing(self, scalar_triangle_box_export: Path):
        assert run_drawing(pjoin(scalar_triangle_box_export, 'sources',
                           'amplitudes', 'triangle_box', 'drawings'))


class TestScalarIsopod:

    def test_info(self, scalar_isopod_export: Path):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            f"launch {scalar_isopod_export}"))
        assert gloop.model.name == 'scalars'
        assert gloop.get_model_from_rust_worker().name == 'scalars'
        for cross_sections in [gloop.cross_sections, gloop.get_cross_sections_from_rust_worker()]:
            assert len(cross_sections) == 0
        for amplitudes in [gloop.amplitudes, gloop.get_amplitudes_from_rust_worker()]:
            assert len(amplitudes) == 1
            assert len(amplitudes[0].amplitude_graphs) == 1
            assert amplitudes[0].name == 'isopod'
        gloop.run(CommandList.from_string("info"))

    @pytest.mark.drawing
    def test_drawing(self, scalar_isopod_export: Path):
        assert run_drawing(pjoin(scalar_isopod_export, 'sources',
                           'amplitudes', 'isopod', 'drawings'))
