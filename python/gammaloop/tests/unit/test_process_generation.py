import pytest
import os
from gammaloop.interface.gammaloop_interface import CommandList, GammaLoop
from gammaloop.misc.common import GammaLoopError
from gammaloop.tests.common import get_gamma_loop_interpreter
from pathlib import Path


class TestProcessGeneration:

    @staticmethod
    def run_tests(gloop: GammaLoop, tests: list[tuple[str, int]]):
        n_parallel_threads = 1
        for test, expected_graph_number in tests:
            try:
                gloop.run(CommandList.from_string(
                    f"generate {test} --clear_existing_processes -nt {n_parallel_threads}"))
            except GammaLoopError as e:
                pass
            n_graphs = (
                0 if len(gloop.amplitudes) == 0 else len(
                    gloop.amplitudes[0].amplitude_graphs)
            ) + (
                0 if len(gloop.cross_sections) == 0 else len(
                    gloop.cross_sections[0].supergraphs)
            )
            assert n_graphs == expected_graph_number, f"For process: '{test}' | Expected {expected_graph_number} graphs, got {n_graphs}"  # nopep8
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
            # # Full particle contents
            ('a > d d~ [{{1}}] --symmetrize_left_right_states', 1),
            ('a > d d~ [{{2}}] --symmetrize_left_right_states', 10),
            # # Only 1-flavour pure QCD corrections
            ('a > d d~ | d g ghG a QED^2=2 [{{1}}] --symmetrize_left_right_states', 1),
            ('a > d d~ | d g ghG a QED^2=2 [{{2}} QCD=1] --symmetrize_left_right_states', 2),
            ('a > d d~ | d g ghG a QED^2=2 [{{3}} QCD=2] --symmetrize_left_right_states', 16),
            ('a > d d~ | d g ghG a QED^2=2 [{{4}} QCD=3] --symmetrize_left_right_states', 166),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_sm_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            # Full particle contents
            # Adding --symmetrize_left_right_states below would give:
            # Could not find the CP conjugate of this vertex in the Feynman rules of the model: (ghWm~, G-, ghA). Consider generating without the option '--symmetrize_left_right_states'.
            ('a > d d~ [{{3}}] -num_grouping only_detect_zeroes', 1249),
            # Only 1-flavour pure QCD corrections
            ('a > d d~ | d g ghG a QED^2=2 [{{5}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 9739),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_sm_full_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm-full"))
        tests = [
            # Full particle contents
            ('a > d d~ [{{1}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 1),
            ('a > d d~ [{{2}}] -num_grouping only_detect_zeroes', 47),
            ('a > d d~ [{{2}}] -num_grouping group_identical_graphs_up_to_sign', 45),
            ('a > d d~ [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 40),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_sm_full_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm-full"))
        tests = [
            # Full particle contents
            ('a > d d~ QED^2=4 [{{3}}] -num_grouping group_identical_graphs_up_to_sign', 339),
            ('a > d d~ [{{3}}] -num_grouping no_grouping', 8189),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_sm_h_n_j(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            # Full particle contents
            ('h > g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 2),
            ('h > g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 30),
            ('h > g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 2),
            ('h > g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 40),
            ('h > g g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 6),
            ('h > g g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 188),
            ('h > g g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 4),
            ('h > g g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 132),
            ('h > g g | h g b t ghg [{{3}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 8),
            ('h > g g | h g b t ghg [{{3}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 3),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_sm_h_n_j_cross_section(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            ('h > g g g | h g b t ghg [{{4}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 96),
            ('h > g g g | h g b t ghg [{{4}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 26),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_sm_h_n_j(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            ('h > g g | h g b t ghg [{3}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 1153),
            ('h > g g g | h g b t ghg [{3}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 6876),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_epem_a_qqh(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            ('e+ e- > d d~ | d a e- ghg g QED^2=4 [{{2}}] -num_grouping only_detect_zeroes', 3),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2=6 [{{2}}] -num_grouping only_detect_zeroes', 3),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2=6 [{{3}} QCD=1] -num_grouping only_detect_zeroes', 30),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2=6 [{{4}} QCD=2] -num_grouping only_detect_zeroes', 594),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2=6 [{{2}}] -num_grouping only_detect_zeroes', 21),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2=6 [{{3}} QCD=2] -num_grouping only_detect_zeroes', 171),
            # Enable grouping
            ('e+ e- > d d~ | d a e- ghg g QED^2=4 [{{2}}] -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2=6 [{{2}}] -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2=6 [{{3}} QCD=1] -num_grouping group_identical_graphs_up_to_scalar_rescaling', 16),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2=6 [{{4}} QCD=2] -num_grouping group_identical_graphs_up_to_scalar_rescaling', 266),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2=6 [{{2}}] -num_grouping group_identical_graphs_up_to_scalar_rescaling', 20),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2=6 [{{3}} QCD=1] -num_grouping group_identical_graphs_up_to_scalar_rescaling', 157),
            # Enable left-right-symmetry
            ('e+ e- > d d~ | d a e- ghg g QED^2=4 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2=6 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2=6 [{{3}} QCD=1] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 11),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2=6 [{{4}} QCD=2] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 166),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2=6 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 16),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2=6 [{{3}} QCD=1] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 104),
            # Too slow sadly
            # ('e+ e- > b b~ h | d b h a e- ghg g z QED^2=6 [{{4}} QCD=4] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', ?),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_a_qqh(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            ('a > d d~ | d a ghg g QED^2=2 [{{2}}] -num_grouping only_detect_zeroes', 3),
            ('a > b b~ h | d b h a ghg g QED^2=4 [{{2}}] -num_grouping only_detect_zeroes', 3),
            ('a > b b~ h | d b h a ghg g QED^2=4 [{{3}} QCD=1] -num_grouping only_detect_zeroes', 30),
            ('a > b b~ h | d b h a ghg g QED^2=4 [{{4}} QCD=2] -num_grouping only_detect_zeroes', 594),
            # # Enable left-right-symmetry and grouping
            ('a > d d~ | d a ghg g QED^2=2 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2),
            ('a > b b~ h | d b h a ghg g QED^2=4 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2),
            ('a > b b~ h | d b h a ghg g QED^2=4 [{{3}} QCD=1] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 11),
            ('a > b b~ h | d b h a ghg g QED^2=4 [{{4}} QCD=2] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 166),
            ('a > b b~ h | b h a ghg g QED^2=4 [{{4}} QCD=2] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 151),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_aa_ttx_amplitude(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            ('a a > t t~ | a t g b ghg QED=2 [{0} QCD=0] -a -num_grouping only_detect_zeroes', 2),
            ('a a > t t~ | a t g b ghg QED=2 [{1} QCD=1] -a -num_grouping only_detect_zeroes', 8),
            ('a a > t t~ | a t g b ghg QED=2 [{2} QCD=2] -a -num_grouping only_detect_zeroes', 144),
            ('a a > t t~ | a t g b ghg QED=2 [{3} QCD=3] -a -num_grouping only_detect_zeroes', 3319),
            ('a a > t t~ | a t g b ghg QED=2 [{0} QCD=0] -a --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 1),
            ('a a > t t~ | a t g b ghg QED=2 [{1} QCD=1] -a --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 4),
            ('a a > t t~ | a t g b ghg QED=2 [{2} QCD=2] -a --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 58),
            ('a a > t t~ | a t g b ghg QED=2 [{3} QCD=3] -a --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 1184),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_aa_ttx_cross_section(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            ('a a > t t~ | a t g b ghg QED^2=4 [{{1}} QCD=0] -num_grouping only_detect_zeroes', 4),
            ('a a > t t~ | a t g b ghg QED^2=4 [{{2}} QCD=1] -num_grouping only_detect_zeroes', 36),
            ('a a > t t~ | a t g b ghg QED^2=4 [{{3}} QCD=2] -num_grouping only_detect_zeroes', 750),
            ('a a > t t~ | a t g b ghg QED^2=4 [{{1}} QCD=0] --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2),
            ('a a > t t~ | a t g b ghg QED^2=4 [{{2}} QCD=1] --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 14),
            ('a a > t t~ | a t g b ghg QED^2=4 [{{3}} QCD=2] --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 219),
            ('a a > t t~ | a t g b ghg QED^2=4 [{{1}} QCD=0] --symmetrize_initial_states --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2),
            ('a a > t t~ | a t g b ghg QED^2=4 [{{2}} QCD=1] --symmetrize_initial_states --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 10),
            ('a a > t t~ | a t g b ghg QED^2=4 [{{3}} QCD=2] --symmetrize_initial_states --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 143),
            ('a a > t t~ | a t g ghg QED^2=4 [{{3}} QCD=2] --symmetrize_initial_states --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 128),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_aa_ttx_cross_section_slow_filter(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            ('a a > t t~ | a t g b ghg QED^2=4 [{{1}} QCD=0] -num_grouping only_detect_zeroes --max_multiplicity_for_fast_cut_filter 0', 4),
            ('a a > t t~ | a t g b ghg QED^2=4 [{{2}} QCD=1] -num_grouping only_detect_zeroes --max_multiplicity_for_fast_cut_filter 0', 36),
            ('a a > t t~ | a t g b ghg QED^2=4 [{{3}} QCD=2] -num_grouping only_detect_zeroes --max_multiplicity_for_fast_cut_filter 0', 750),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_vaccuum_amplitude_generation(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            # Including all graphs
            ('{} > {} | g ghg t u d QED=0 [{2}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 6),
            ('{} > {} | g ghg t u d QED=0 [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 5),
            ('{} > {} | g ghg t u d QED=0 [{3}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 36),
            ('{} > {} | g ghg t u d QED=0 [{3}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 27),
            ('{} > {} | g ghg t u d QED=0 [{4}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 264),
            ('{} > {} | g ghg t u d QED=0 [{4}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 179),
            # Only non-factorizable ones
            ('{} > {} | g ghg t u d QED=0 [{2}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 5),
            ('{} > {} | g ghg t u d QED=0 [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 4),
            ('{} > {} | g ghg t u d QED=0 [{3}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 29),
            ('{} > {} | g ghg t u d QED=0 [{3}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 21),
            ('{} > {} | g ghg t u d QED=0 [{4}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 212),
            ('{} > {} | g ghg t u d QED=0 [{4}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 138),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_vaccuum_amplitude_generation(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            ('{} > {} | g ghg t u d QED=0 [{5}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 2560),
            ('{} > {} | g ghg t u d QED=0 [{5}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 1438),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_vaccuum_amplitude_generation_full_sm(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm-full"))
        tests = [
            # Including all graphs
            ('{} > {} [{2}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 448),
            ('{} > {} [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 243),
            ('{} > {} [{3}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 36362),
            ('{} > {} [{3}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 13845),
            # Only non-factorizable ones
            ('{} > {} [{2}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 92),
            ('{} > {} [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 69),
            ('{} > {} [{3}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 3085),
            ('{} > {} [{3}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 2009),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_dis_isr(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            ('e- d > e- d | e- a d g ghg QED=2 [{0}] -a -num_grouping only_detect_zeroes', 1),
            ('e- d > e- d | e- a d g ghg QED=2 [{1}] -a -num_grouping only_detect_zeroes', 1),
            ('e- d > e- d | e- a d g ghg QED=2 [{2}] -a -num_grouping only_detect_zeroes', 13),
            ('e- d > e- d | e- a d g ghg QED=2 [{3}] -a -num_grouping only_detect_zeroes', 236),
            ('e- d > e- d | e- a d g ghg QED=2 [{0}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 1),
            ('e- d > e- d | e- a d g ghg QED=2 [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 1),
            ('e- d > e- d | e- a d g ghg QED=2 [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 11),
            ('e- d > e- d | e- a d g ghg QED=2 [{3}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 186),
            # 3>N no symmetrization
            ('e- d g > e- d | e- a d g ghg QED=2 [{0}] -a -num_grouping only_detect_zeroes', 2),
            ('e- d g > e- d | e- a d g ghg QED=2 [{1}] -a -num_grouping only_detect_zeroes', 13),
            ('e- d g > e- d | e- a d g ghg QED=2 [{2}] -a -num_grouping only_detect_zeroes', 223),
            # 3>N with symmetrization
            ('e- d g > e- d | e- a d g ghg QED=2 [{0}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2),
            ('e- d g > e- d | e- a d g ghg QED=2 [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 11),
            ('e- d g > e- d | e- a d g ghg QED=2 [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 179),
            # 4>N no symmetrization
            ('e- d u u~ > e- d | e- a d u g ghg QED=2 [{0}] -a -num_grouping only_detect_zeroes', 4),
            ('e- d u u~ > e- d | e- a d u g ghg QED=2 [{1}] -a -num_grouping only_detect_zeroes', 62),
            ('e- d u u~ > e- d | e- a d u g ghg QED=2 [{2}] -a -num_grouping only_detect_zeroes', 1428),
            # 4>N with symmetrization
            ('e- d u u~ > e- d | e- a d u g ghg QED=2 [{0}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 4),
            ('e- d u u~ > e- d | e- a d u g ghg QED=2 [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 54),
            ('e- d u u~ > e- d | e- a d u g ghg QED=2 [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 1107)
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_a_qqh(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests = [
            ('a > b b~ h | b h a ghg g QED^2=4 [{{5}} QCD=3] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2763),
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
            ('d d~ > u u~ g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 88),
            ('d d~ > d d~ g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 176),
            ('d d~ > g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 341),
            ('g g > g g g | g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 671),
            # For fun, symmetrize it all
            ('g g > g g g | g ghg a QED=0 [QCD=1] -a --symmetrize_initial_states --symmetrize_final_states -num_grouping group_identical_graphs_up_to_sign', 109),
            ('g g > g g g | g ghg a QED=0 [QCD=1] -a --symmetrize_left_right_states --symmetrize_initial_states --symmetrize_final_states -num_grouping group_identical_graphs_up_to_sign', 25),
            ('g g > g g g | u d g ghg a QED=0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 905),
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
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    # A bit too slow, let's skip for now
    @pytest.mark.skip
    @pytest.mark.slow
    def test_very_slow_generate_amplitude_1l_sm_jets(self):
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

    @pytest.mark.skip
    @pytest.mark.slow
    def test_slow_generate_amplitude_1l_sm_jets_with_grouping(self):
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
