import pytest
import os
from gammaloop.interface.gammaloop_interface import CommandList, GammaLoop
from gammaloop.misc.common import GammaLoopError
from gammaloop.tests.common import get_gamma_loop_interpreter
from pathlib import Path
from gammaloop.misc.utils import parse_python_expression, expression_to_string
from symbolica import Expression, E


class TestProcessGeneration:

    @staticmethod
    def evaluate_overall_factor(overall_factor: str) -> Expression:
        e = parse_python_expression(overall_factor)
        assert e is not None, f"Could not parse overall factor: '{overall_factor}'"  # nopep8
        for header in ["AutG",            "CouplingsMultiplicity",
                       "InternalFermionLoopSign",
                       "ExternalFermionOrderingSign",
                       "NumeratorIndependentSymmetryGrouping"]:
            e = e.replace_all(E(f"{header}(x_)"), E("x_"))
        e = e.replace_all(
            E("NumeratorDependentGrouping(GraphId_,ratio_,GraphSymmetryFactor_)"), E("ratio_*GraphSymmetryFactor_"))
        return e.expand()

    @staticmethod
    def run_tests(gloop: GammaLoop, tests: list[tuple[str, int, str | None]]):
        n_parallel_threads = 1
        for test, expected_graph_number, expected_total_overall_factor in tests:
            try:
                gloop.run(CommandList.from_string(
                    f"generate {test} --clear_existing_processes -nt {n_parallel_threads}"))
            except GammaLoopError as e:
                pass
            all_graphs = (
                ([] if len(gloop.amplitudes) == 0 else [
                    g.graph.overall_factor for g in gloop.amplitudes[0].amplitude_graphs
                ])
                + ([] if len(gloop.cross_sections) == 0 else [
                    g.cuts[0].forward_scattering_graph.graph.overall_factor for g in gloop.cross_sections[0].supergraphs
                ])
            )

            total_overall_factor = E("0")
            all_overall_factors: list[Expression] = [TestProcessGeneration.evaluate_overall_factor(of) for of in all_graphs]  # nopep8 # type: ignore
            for of in all_overall_factors:
                total_overall_factor += of
            total_overall_factor_str = expression_to_string(total_overall_factor.expand())  # nopep8
            assert total_overall_factor_str is not None, f"Could not convert total overall factor to string: '{total_overall_factor}'"  # nopep8

            n_graphs = len(all_graphs)
            assert n_graphs == expected_graph_number, f"For process: '{test}' | Expected {expected_graph_number} graphs, got {n_graphs} (total_overall_factor: {total_overall_factor_str})"  # nopep8
            if expected_total_overall_factor is not None:
                expected_total_overall_factor_expr = parse_python_expression(
                    expected_total_overall_factor)
                assert expected_total_overall_factor_expr is not None, f"Could not parse expected total overall factor: '{expected_total_overall_factor}'"  # nopep8
                difference = (total_overall_factor -
                              expected_total_overall_factor_expr).expand()
                assert (difference == E("0")).eval(), f"For process: '{test}' | Expected total overall factor: '{expected_total_overall_factor}', got '{total_overall_factor_str}'"  # nopep8
            gloop.amplitudes.clear()
            gloop.cross_sections.clear()

    def test_generate_sm_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            # Only d, g and a as particle contents
            # Test the validity of multiple final-state specifications for cross-section generation
            ('a > {d d~, g g} [{{1}}] | d g a -num_grouping only_detect_zeroes', 1, "-1"),
            ('a > d d~ g [{{2}}] | d g a -num_grouping only_detect_zeroes', 3, "-3"),
            ('a > d d~ g [{{2}}] | d g a --symmetrize_left_right_states', 2, "-3"),
            ('a > d d~ z [{{2}}] | d g a -num_grouping only_detect_zeroes', 0, "0"),
            # # Full particle contents
            ('a > d d~ [{{1}}] --symmetrize_left_right_states', 1, "-1"),
            ('a > d d~ [{{2}}] --symmetrize_left_right_states', 10, "-3*(9*ee^-2*G^2*Nc*TR-9*ee^-2*G^2*Nc^-1*TR)-12"),  # nopep8
            # # Only 1-flavour pure QCD corrections
            ('a > d d~ | d g ghG a QED^2==2 [{{1}}] --symmetrize_left_right_states', 1, "-1"),
            ('a > d d~ | d g ghG a QED^2==2 [{{2}} QCD=1] --symmetrize_left_right_states', 2, "-3"),
            ('a > d d~ | d g ghG a QED^2==2 [{{3}} QCD=2] --symmetrize_left_right_states', 16, "-41/2"),
            ('a > d d~ | d g ghG a QED^2==2 [{{4}} QCD=3] --symmetrize_left_right_states', 166, "-3107/14"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_sm_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            # Full particle contents
            # Adding --symmetrize_left_right_states below would give:
            # Could not find the CP conjugate of this vertex in the Feynman rules of the model: (ghWm~, G-, ghA). Consider generating without the option '--symmetrize_left_right_states'.
            ('a > d d~ [{{3}}] -num_grouping only_detect_zeroes', 1321, "-1103/2"),  # nopep8
            # Only 1-flavour pure QCD corrections
            ('a > d d~ | d g ghG a QED^2==2 [{{5}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 6303, "-51683/24"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_sm_full_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm-full"))
        tests: list[tuple[str, int, str | None]] = [
            # Full particle contents
            ('a > d d~ [{{1}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 1, "-1"),
            ('a > d d~ [{{2}}] -num_grouping only_detect_zeroes', 47, "-47"),
            ('a > d d~ [{{2}}] -num_grouping group_identical_graphs_up_to_sign', 45, "-47"),
            ('a > d d~ [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 40, "-47"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_sm_full_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm-full"))
        tests: list[tuple[str, int, str | None]] = [
            # Full particle contents
            ('a > d d~ QED^2==4 [{{3}}] -num_grouping group_identical_graphs_up_to_sign', 339, "-368"),
            ('a > d d~ [{{3}}] -num_grouping no_grouping', 8549, "-9111/2"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_sm_h_n_j(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            # Full particle contents
            ('h > g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 2, "-4"),
            ('h > g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 30, "-74"),
            ('h > g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 2, "-4"),
            ('h > g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 40, "-74"),
            ('h > g g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 6, "-24"),
            ('h > g g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 188, "-618"),
            ('h > g g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 4, "-24"),
            ('h > g g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 132, "-618"),
            ('h > g g | h g b t ghg [{{3}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 8, "8"),
            ('h > g g | h g b t ghg [{{3}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 3, "8"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_sm_h_n_j_cross_section(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            ('h > g g g | h g b t ghg [{{4}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', 72, "88"),
            ('h > g g g | h g b t ghg [{{4}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', 26, "88"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_sm_h_n_j(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            ('h > g g | h g b t ghg [{3}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 1194, "-1260"),
            ('h > g g g | h g b t ghg [{3}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', 6946, "-12429"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_epem_a_qqh(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            ('e+ e- > d d~ | d a e- ghg g QED^2==4 [{{2}}] -num_grouping only_detect_zeroes', 3, "-3"),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{2}}] -num_grouping only_detect_zeroes', 3, "-3"),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{3}} QCD=1] -num_grouping only_detect_zeroes', 30, "-30"),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{4}} QCD=2] -num_grouping only_detect_zeroes', 594, "-261"),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{2}}] -num_grouping only_detect_zeroes', 21, "-21"),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{3}} QCD=2] -num_grouping only_detect_zeroes', 171, "-171"),
            # Enable grouping
            ('e+ e- > d d~ | d a e- ghg g QED^2==4 [{{2}}] -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2, "-3"),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{2}}] -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2, "-3"),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{3}} QCD=1] -num_grouping group_identical_graphs_up_to_scalar_rescaling', 16, "-30"),
            # This grouping involves 21 new cancellations between graphs, so that modifies the total overall factor
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{4}} QCD=2] -num_grouping group_identical_graphs_up_to_scalar_rescaling', 266, "-327"),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{2}}] -num_grouping group_identical_graphs_up_to_scalar_rescaling', 20, "-21"),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{3}} QCD=1] -num_grouping group_identical_graphs_up_to_scalar_rescaling', 157, "-171"),
            # Enable left-right-symmetry
            ('e+ e- > d d~ | d a e- ghg g QED^2==4 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2, "-3"),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2, "-3"),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{3}} QCD=1] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 11, "-30"),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{4}} QCD=2] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 166, "-327"),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 16, "-21"),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{3}} QCD=1] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 104, "-171"),
            # Too slow sadly
            # ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{4}} QCD=4] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', ?),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_a_qqh(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            ('a > d d~ | d a ghg g QED^2==2 [{{2}}] -num_grouping only_detect_zeroes', 3, "-3"),
            ('a > b b~ h | d b h a ghg g QED^2==4 [{{2}}] -num_grouping only_detect_zeroes', 3, "-3"),
            ('a > b b~ h | d b h a ghg g QED^2==4 [{{3}} QCD=1] -num_grouping only_detect_zeroes', 30, "-30"),
            ('a > b b~ h | d b h a ghg g QED^2==4 [{{4}} QCD=2] -num_grouping only_detect_zeroes', 594, "-261"),
            # Enable left-right-symmetry and grouping
            ('a > d d~ | d a ghg g QED^2==2 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2, "-3"),
            ('a > b b~ h | d b h a ghg g QED^2==4 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2, "-3"),
            ('a > b b~ h | d b h a ghg g QED^2==4 [{{3}} QCD=1] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 11, "-30"),
            # Additional cancellations, so that modifies the total overall factor
            ('a > b b~ h | d b h a ghg g QED^2==4 [{{4}} QCD=2] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 166, "-327"),
            ('a > b b~ h | b h a ghg g QED^2==4 [{{4}} QCD=2] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 151, "-366"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_aa_ttx_amplitude(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            ('a a > t t~ | a t g b ghg QED==2 [{0} QCD=0] -a -num_grouping only_detect_zeroes', 2, "2"),
            ('a a > t t~ | a t g b ghg QED==2 [{1} QCD=1] -a -num_grouping only_detect_zeroes', 8, "8"),
            ('a a > t t~ | a t g b ghg QED==2 [{2} QCD=2] -a -num_grouping only_detect_zeroes', 144, "36"),
            ('a a > t t~ | a t g b ghg QED==2 [{3} QCD=3] -a -num_grouping only_detect_zeroes', 3424, "28/3"),
            ('a a > t t~ | a t g b ghg QED==2 [{0} QCD=0] -a --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 1, "2"),
            ('a a > t t~ | a t g b ghg QED==2 [{1} QCD=1] -a --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 4, "8"),
            # Six additional cancellations imply that the total overall factor is different
            ('a a > t t~ | a t g b ghg QED==2 [{2} QCD=2] -a --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 58, "60"),
            ('a a > t t~ | a t g b ghg QED==2 [{3} QCD=3] -a --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 1221, "1600/3"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_aa_ttx_cross_section(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            ('a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] -num_grouping only_detect_zeroes', 4, "-4"),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] -num_grouping only_detect_zeroes', 40, "-40"),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] -num_grouping only_detect_zeroes', 874, "-266"),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2, "-4"),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 14, "-40"),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 219, "-426"),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] --symmetrize_initial_states --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2, "-4"),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] --symmetrize_initial_states --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 10, "-40"),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] --symmetrize_initial_states --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 143, "-426"),
            ('a a > t t~ | a t g ghg QED^2==4 [{{3}} QCD=2] --symmetrize_initial_states --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 128, "-502"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_aa_ttx_cross_section_slow_filter(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            ('a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] -num_grouping only_detect_zeroes --max_multiplicity_for_fast_cut_filter 0', 4, "-4"),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] -num_grouping only_detect_zeroes --max_multiplicity_for_fast_cut_filter 0', 40, "-40"),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] -num_grouping only_detect_zeroes --max_multiplicity_for_fast_cut_filter 0', 874, "-266"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_vaccuum_amplitude_generation(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            # Including all graphs
            ('{} > {} | g ghg t u d QED==0 [{2}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 6, "-19/24"),
            ('{} > {} | g ghg t u d QED==0 [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 5, "-19/24"),
            ('{} > {} | g ghg t u d QED==0 [{3}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 36, "-71/48"),
            ('{} > {} | g ghg t u d QED==0 [{3}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 27, "-71/48"),
            ('{} > {} | g ghg t u d QED==0 [{4}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 264, "-787/120"),
            ('{} > {} | g ghg t u d QED==0 [{4}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 179, "-4159/840"),
            # Only non-factorizable ones
            ('{} > {} | g ghg t u d QED==0 [{2}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 5, "-11/12"),
            ('{} > {} | g ghg t u d QED==0 [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 4, "-11/12"),
            ('{} > {} | g ghg t u d QED==0 [{3}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 29, "-19/16"),
            ('{} > {} | g ghg t u d QED==0 [{3}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 21, "-19/16"),
            ('{} > {} | g ghg t u d QED==0 [{4}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 212, "-811/160"),
            ('{} > {} | g ghg t u d QED==0 [{4}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 138, "-3877/1120"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_vaccuum_amplitude_generation(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            ('{} > {} | g ghg t u d QED==0 [{5}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 2560, "-45055/2688"),
            ('{} > {} | g ghg t u d QED==0 [{5}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 1438, "226763/51072"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_vaccuum_amplitude_generation_full_sm(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm-full"))
        tests: list[tuple[str, int, str | None]] = [
            # Including all graphs
            ('{} > {} [{2}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 448, "-19/12"),
            ('{} > {} [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 243,
             "-81/8*ee^-2*G^2*Nc*TR+81/8*ee^-2*G^2*Nc^-1*TR-CKM2x1^-1*CKM2x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-CKM3x1^-1*CKM3x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-I2x12^-1*I2x22*I3x21^-1*I3x22-I2x13^-1*I2x23*I3x31^-1*I3x32-CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+313/24"),
            ('{} > {} [{3}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 36362, "-257/4"),
            ('{} > {} [{3}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', 13845,
             "-27/2*ee^-2*G^2+7533/32*ee^-4*G^4*TR^2-5/2*CKM1x3*I1x31*I4x13^-1*complexconjugate(CKM1x3)^-1-5/2*CKM2x1*I2x12*I3x21^-1*complexconjugate(CKM2x1)^-1-5/2*CKM3x1*I2x13*I3x31^-1*complexconjugate(CKM3x1)^-1-5/2*CKM2x2*I2x22*I3x21^-1*complexconjugate(CKM2x1)^-1-5/2*CKM3x2*I2x23*I3x31^-1*complexconjugate(CKM3x1)^-1-729/4*ee^-4*G^4*Nc^-2*TR^2-1701/32*ee^-4*G^4*Nc^2*TR^2+2133/32*ee^-2*G^2*Nc*TR-2133/32*ee^-2*G^2*Nc^-1*TR-31/2*ee^2*cw^-1*vev*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1+ee^4*cw^-2*vev^2*(ee^4*cw^-2*vev^2-2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1+2*ee^4*sw^-2*vev^2*(ee^4*cw^-2*vev^2-2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1-1/2*CKM2x1^-2*CKM2x2^2*complexconjugate(CKM2x1)^-2*complexconjugate(CKM2x2)^2+31/4*CKM2x1^-1*CKM2x2*I2x12^-1*I2x22-2*CKM2x1^-1*CKM2x2*I2x13^-1*I2x23+5/2*CKM2x1^-1*CKM2x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-6*CKM2x1^-1*CKM2x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-5*CKM2x1^-1*CKM2x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-1/2*CKM3x1^-2*CKM3x2^2*complexconjugate(CKM3x1)^-2*complexconjugate(CKM3x2)^2-2*CKM3x1^-1*CKM3x2*I2x12^-1*I2x22+31/4*CKM3x1^-1*CKM3x2*I2x13^-1*I2x23-6*CKM3x1^-1*CKM3x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+5/2*CKM3x1^-1*CKM3x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-5*CKM3x1^-1*CKM3x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-1/2*I2x12^-2*I2x22^2*I3x21^-2*I3x22^2+7/2*I2x12^-1*I2x22*I3x21^-1*I3x22-6*I2x12^-1*I2x22*I3x31^-1*I3x32-2*I2x12^-1*I2x22*CKM1x1^-1*CKM1x2-1/2*I2x13^-2*I2x23^2*I3x31^-2*I3x32^2-6*I2x13^-1*I2x23*I3x21^-1*I3x22+7/2*I2x13^-1*I2x23*I3x31^-1*I3x32-2*I2x13^-1*I2x23*CKM1x1^-1*CKM1x2+11/2*I3x21^-1*I3x22*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-2*I3x21^-1*I3x22*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-2*I3x21^-1*I3x22*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-2*I3x31^-1*I3x32*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+11/2*I3x31^-1*I3x32*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-2*I3x31^-1*I3x32*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-1/2*CKM1x1^-2*CKM1x2^2*complexconjugate(CKM1x1)^-2*complexconjugate(CKM1x2)^2-5*CKM1x1^-1*CKM1x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-5*CKM1x1^-1*CKM1x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-1/2*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-31/2*ee^2*cw*sw^-2*vev*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1+ee^4*cw^2*sw^-4*vev^2*(ee^4*cw^-2*vev^2-2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1+2*ee^8*cw^-4*vev^4*(ee^4*cw^-2*vev^2-2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1*(ee^4*cw^-2*vev^2+2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1-4*ee^8*sw^-4*vev^4*(ee^4*cw^-2*vev^2-2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1*(ee^4*cw^-2*vev^2+2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1-3*ee^-2*G^2*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+2*ee^8*cw^4*sw^-8*vev^4*(ee^4*cw^-2*vev^2-2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1*(ee^4*cw^-2*vev^2+2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1-27/4*ee^-2*CKM2x1^-1*CKM2x2*G^2*Nc*TR*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+27/4*ee^-2*CKM2x1^-1*CKM2x2*G^2*Nc^-1*TR*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-27/4*ee^-2*CKM3x1^-1*CKM3x2*G^2*Nc*TR*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+27/4*ee^-2*CKM3x1^-1*CKM3x2*G^2*Nc^-1*TR*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-27/4*ee^-2*G^2*I2x12^-1*I2x22*I3x21^-1*I3x22*Nc*TR+27/4*ee^-2*G^2*I2x12^-1*I2x22*I3x21^-1*I3x22*Nc^-1*TR-27/4*ee^-2*G^2*I2x13^-1*I2x23*I3x31^-1*I3x32*Nc*TR+27/4*ee^-2*G^2*I2x13^-1*I2x23*I3x31^-1*I3x32*Nc^-1*TR+9/4*ee^-2*G^2*CKM1x1^-1*CKM1x2*Nc*TR*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-9/4*ee^-2*G^2*CKM1x1^-1*CKM1x2*Nc^-1*TR*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+1/4*ee^2*cw^-1*vev*CKM2x1^-1*CKM2x2*I2x12^-1*I2x22*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1+1/4*ee^2*cw^-1*vev*CKM3x1^-1*CKM3x2*I2x13^-1*I2x23*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1+1/4*ee^2*cw^-1*vev*I3x21^-1*I3x22*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+1/4*ee^2*cw^-1*vev*I3x31^-1*I3x32*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-CKM2x1^-1*CKM2x2*I2x12^-1*I2x22*I3x21^-1*I3x22*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-CKM2x1^-1*CKM2x2*I2x13^-1*I2x23*I3x31^-1*I3x32*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-CKM2x1^-1*CKM2x2*CKM1x1^-1*CKM1x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-CKM2x1^-1*CKM3x1^-1*CKM2x2*CKM3x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM3x1)^-1*complexconjugate(CKM2x2)*complexconjugate(CKM3x2)-CKM3x1^-1*CKM3x2*I2x12^-1*I2x22*I3x21^-1*I3x22*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-CKM3x1^-1*CKM3x2*I2x13^-1*I2x23*I3x31^-1*I3x32*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-CKM3x1^-1*CKM3x2*CKM1x1^-1*CKM1x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-I2x12^-1*I2x22*I3x21^-1*I3x22*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-I2x12^-1*I2x13^-1*I2x22*I2x23*I3x21^-1*I3x22*I3x31^-1*I3x32-I2x13^-1*I2x23*I3x31^-1*I3x32*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+1/4*ee^2*cw*sw^-2*vev*CKM2x1^-1*CKM2x2*I2x12^-1*I2x22*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1+1/4*ee^2*cw*sw^-2*vev*CKM3x1^-1*CKM3x2*I2x13^-1*I2x23*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1+1/4*ee^2*cw*sw^-2*vev*I3x21^-1*I3x22*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+1/4*ee^2*cw*sw^-2*vev*I3x31^-1*I3x32*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+1169/4"),
            # Only non-factorizable ones
            ('{} > {} [{2}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 92, "-67/3"),
            ('{} > {} [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 69,
             "-81/8*ee^-2*G^2*Nc*TR+81/8*ee^-2*G^2*Nc^-1*TR-CKM2x1^-1*CKM2x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-CKM3x1^-1*CKM3x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-I2x12^-1*I2x22*I3x21^-1*I3x22-I2x13^-1*I2x23*I3x31^-1*I3x32-CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-163/12"),
            ('{} > {} [{3}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 3085, "-6703/12"),
            ('{} > {} [{3}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', 2009, "-3*ee^-2*G^2+7533/32*ee^-4*G^4*TR^2-5/2*CKM1x3*I1x31*I4x13^-1*complexconjugate(CKM1x3)^-1-5/2*CKM2x1*I2x12*I3x21^-1*complexconjugate(CKM2x1)^-1-5/2*CKM3x1*I2x13*I3x31^-1*complexconjugate(CKM3x1)^-1-5/2*CKM2x2*I2x22*I3x21^-1*complexconjugate(CKM2x1)^-1-5/2*CKM3x2*I2x23*I3x31^-1*complexconjugate(CKM3x1)^-1-729/4*ee^-4*G^4*Nc^-2*TR^2-1701/32*ee^-4*G^4*Nc^2*TR^2-2907/16*ee^-2*G^2*Nc*TR+2907/16*ee^-2*G^2*Nc^-1*TR-9*ee^2*cw^-1*vev*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1+ee^4*cw^-2*vev^2*(ee^4*cw^-2*vev^2-2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1+2*ee^4*sw^-2*vev^2*(ee^4*cw^-2*vev^2-2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1-1/2*CKM2x1^-2*CKM2x2^2*complexconjugate(CKM2x1)^-2*complexconjugate(CKM2x2)^2+7/2*CKM2x1^-1*CKM2x2*I2x12^-1*I2x22-2*CKM2x1^-1*CKM2x2*I2x13^-1*I2x23-81/4*CKM2x1^-1*CKM2x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-6*CKM2x1^-1*CKM2x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-5*CKM2x1^-1*CKM2x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-1/2*CKM3x1^-2*CKM3x2^2*complexconjugate(CKM3x1)^-2*complexconjugate(CKM3x2)^2-2*CKM3x1^-1*CKM3x2*I2x12^-1*I2x22+7/2*CKM3x1^-1*CKM3x2*I2x13^-1*I2x23-6*CKM3x1^-1*CKM3x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-81/4*CKM3x1^-1*CKM3x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-5*CKM3x1^-1*CKM3x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-1/2*I2x12^-2*I2x22^2*I3x21^-2*I3x22^2-77/4*I2x12^-1*I2x22*I3x21^-1*I3x22-6*I2x12^-1*I2x22*I3x31^-1*I3x32-2*I2x12^-1*I2x22*CKM1x1^-1*CKM1x2-1/2*I2x13^-2*I2x23^2*I3x31^-2*I3x32^2-6*I2x13^-1*I2x23*I3x21^-1*I3x22-77/4*I2x13^-1*I2x23*I3x31^-1*I3x32-2*I2x13^-1*I2x23*CKM1x1^-1*CKM1x2+5/4*I3x21^-1*I3x22*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-2*I3x21^-1*I3x22*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-2*I3x21^-1*I3x22*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-2*I3x31^-1*I3x32*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+5/4*I3x31^-1*I3x32*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-2*I3x31^-1*I3x32*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-1/2*CKM1x1^-2*CKM1x2^2*complexconjugate(CKM1x1)^-2*complexconjugate(CKM1x2)^2-5*CKM1x1^-1*CKM1x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-5*CKM1x1^-1*CKM1x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-63/4*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-9*ee^2*cw*sw^-2*vev*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1+ee^4*cw^2*sw^-4*vev^2*(ee^4*cw^-2*vev^2-2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1+2*ee^8*cw^-4*vev^4*(ee^4*cw^-2*vev^2-2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1*(ee^4*cw^-2*vev^2+2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1-4*ee^8*sw^-4*vev^4*(ee^4*cw^-2*vev^2-2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1*(ee^4*cw^-2*vev^2+2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1-3*ee^-2*G^2*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+2*ee^8*cw^4*sw^-8*vev^4*(ee^4*cw^-2*vev^2-2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1*(ee^4*cw^-2*vev^2+2*ee^4*sw^-2*vev^2+ee^4*cw^2*sw^-4*vev^2)^-1-27/4*ee^-2*CKM2x1^-1*CKM2x2*G^2*Nc*TR*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+27/4*ee^-2*CKM2x1^-1*CKM2x2*G^2*Nc^-1*TR*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-27/4*ee^-2*CKM3x1^-1*CKM3x2*G^2*Nc*TR*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+27/4*ee^-2*CKM3x1^-1*CKM3x2*G^2*Nc^-1*TR*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-27/4*ee^-2*G^2*I2x12^-1*I2x22*I3x21^-1*I3x22*Nc*TR+27/4*ee^-2*G^2*I2x12^-1*I2x22*I3x21^-1*I3x22*Nc^-1*TR-27/4*ee^-2*G^2*I2x13^-1*I2x23*I3x31^-1*I3x32*Nc*TR+27/4*ee^-2*G^2*I2x13^-1*I2x23*I3x31^-1*I3x32*Nc^-1*TR+9/4*ee^-2*G^2*CKM1x1^-1*CKM1x2*Nc*TR*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-9/4*ee^-2*G^2*CKM1x1^-1*CKM1x2*Nc^-1*TR*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+1/4*ee^2*cw^-1*vev*CKM2x1^-1*CKM2x2*I2x12^-1*I2x22*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1+1/4*ee^2*cw^-1*vev*CKM3x1^-1*CKM3x2*I2x13^-1*I2x23*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1+1/4*ee^2*cw^-1*vev*I3x21^-1*I3x22*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+1/4*ee^2*cw^-1*vev*I3x31^-1*I3x32*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-CKM2x1^-1*CKM2x2*I2x12^-1*I2x22*I3x21^-1*I3x22*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-CKM2x1^-1*CKM2x2*I2x13^-1*I2x23*I3x31^-1*I3x32*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)-CKM2x1^-1*CKM2x2*CKM1x1^-1*CKM1x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-CKM2x1^-1*CKM3x1^-1*CKM2x2*CKM3x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM3x1)^-1*complexconjugate(CKM2x2)*complexconjugate(CKM3x2)-CKM3x1^-1*CKM3x2*I2x12^-1*I2x22*I3x21^-1*I3x22*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-CKM3x1^-1*CKM3x2*I2x13^-1*I2x23*I3x31^-1*I3x32*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-CKM3x1^-1*CKM3x2*CKM1x1^-1*CKM1x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-I2x12^-1*I2x22*I3x21^-1*I3x22*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)-I2x12^-1*I2x13^-1*I2x22*I2x23*I3x21^-1*I3x22*I3x31^-1*I3x32-I2x13^-1*I2x23*I3x31^-1*I3x32*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+1/4*ee^2*cw*sw^-2*vev*CKM2x1^-1*CKM2x2*I2x12^-1*I2x22*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1+1/4*ee^2*cw*sw^-2*vev*CKM3x1^-1*CKM3x2*I2x13^-1*I2x23*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1+1/4*ee^2*cw*sw^-2*vev*I3x21^-1*I3x22*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+1/4*ee^2*cw*sw^-2*vev*I3x31^-1*I3x32*(ee^2*cw^-1*vev-ee^2*cw*sw^-2*vev)^-1*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)-25619/96"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_dis_isr(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            ('e- d > e- d | e- a d g ghg QED==2 [{0}] -a -num_grouping only_detect_zeroes', 1, "-1"),
            ('e- d > e- d | e- a d g ghg QED==2 [{1}] -a -num_grouping only_detect_zeroes', 1, "-1"),
            ('e- d > e- d | e- a d g ghg QED==2 [{2}] -a -num_grouping only_detect_zeroes', 13, "-13/2"),
            ('e- d > e- d | e- a d g ghg QED==2 [{3}] -a -num_grouping only_detect_zeroes', 236, "-194/3"),
            ('e- d > e- d | e- a d g ghg QED==2 [{0}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 1, "-1"),
            ('e- d > e- d | e- a d g ghg QED==2 [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 1, "-1"),
            ('e- d > e- d | e- a d g ghg QED==2 [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 11, "-17/2"),
            ('e- d > e- d | e- a d g ghg QED==2 [{3}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 186, "-2183/21"),
            # 3>N no symmetrization
            ('e- d g > e- d | e- a d g ghg QED==2 [{0}] -a -num_grouping only_detect_zeroes', 2, "-2"),
            ('e- d g > e- d | e- a d g ghg QED==2 [{1}] -a -num_grouping only_detect_zeroes', 13, "-9"),
            ('e- d g > e- d | e- a d g ghg QED==2 [{2}] -a -num_grouping only_detect_zeroes', 223, "-169/2"),
            # 3>N with symmetrization
            ('e- d g > e- d | e- a d g ghg QED==2 [{0}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2, "-2"),
            ('e- d g > e- d | e- a d g ghg QED==2 [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 11, "-11"),
            ('e- d g > e- d | e- a d g ghg QED==2 [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 179, "-842/7"),
            # 4>N no symmetrization
            ('e- d u u~ > e- d | e- a d u g ghg QED==2 [{0}] -a -num_grouping only_detect_zeroes', 4, "4"),
            ('e- d u u~ > e- d | e- a d u g ghg QED==2 [{1}] -a -num_grouping only_detect_zeroes', 62, "36"),
            ('e- d u u~ > e- d | e- a d u g ghg QED==2 [{2}] -a -num_grouping only_detect_zeroes', 1428, "1298/3"),
            # 4>N with symmetrization
            ('e- d u u~ > e- d | e- a d u g ghg QED==2 [{0}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 4, "4"),
            ('e- d u u~ > e- d | e- a d u g ghg QED==2 [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 54, "40"),
            ('e- d u u~ > e- d | e- a d u g ghg QED==2 [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 1107, "3343/6")
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_a_qqh(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            ('a > b b~ h | b h a ghg g QED^2==4 [{{5}} QCD=3] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', 2763, "-38712/7"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_amplitude_1l_sm_jets(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        tests: list[tuple[str, int, str | None]] = [
            ('d d~ > d d~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 20, "0"),
            ('d d~ > d d~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 192, "0"),
            ('d d~ > g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 384, "-186"),
            ('g g > g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 1139, "171/2"),
            ('d d~ > d d~ g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 2284, "0"),
            ('d d~ > u u~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 10, "-11/2"),
            ('d d~ > u u~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 96, "-119/2"),
            ('d d~ > u u~ g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 1142, "-698"),
            ('d d~ > u u~ d d~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 414, "0"),
            ('d d~ > d d~ d d~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 1242, "0"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    # Currently bugged with spenso network parsing for "-num_grouping group_identical_graphs_up_to_sign"
    def test_generate_amplitude_1l_sm_jets_with_grouping(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        tests: list[tuple[str, int, str | None]] = [
            ('d d~ > d d~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 18, "0"),
            ('d d~ > u u~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 9, "-11/2"),
            ('d d~ > u u~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 88, "-119/2"),
            ('d d~ > d d~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 176, "0"),
            ('d d~ > g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 341, "-186"),
            ('g g > g g g | g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 671, "1107/2"),
            # For fun, symmetrize it all
            ('g g > g g g | g ghg a QED==0 [QCD=1] -a --symmetrize_initial_states --symmetrize_final_states -num_grouping group_identical_graphs_up_to_sign', 109, "1107/2"),
            ('g g > g g g | g ghg a QED==0 [QCD=1] -a --symmetrize_left_right_states --symmetrize_initial_states --symmetrize_final_states -num_grouping group_identical_graphs_up_to_sign', 25, "1107/2"),
            ('g g > g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 905, "171/2"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.slow
    def test_slow_generate_amplitude_1l_sm_jets(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        tests: list[tuple[str, int, str | None]] = [
            ('d d~ > u u~ d d~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 5424, "0"),
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
        tests: list[tuple[str, int, str | None]] = [
            ('g g > g g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 14875, "1375"),
            ('d d~ > d d~ g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 32074, "0"),
            ('d d~ > u u~ g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 16037, "-18993/2"),
            ('d d~ > d d~ d d~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', 16272, "0"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @pytest.mark.skip
    @pytest.mark.slow
    def test_slow_generate_amplitude_1l_sm_jets_with_grouping(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        tests: list[tuple[str, int, str | None]] = [
            ('d d~ > d d~ g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 2090, "0"),
            ('d d~ > u u~ g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 1045, "-698"),
            ('d d~ > u u~ d d~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 380, "0"),
            ('d d~ > d d~ d d~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 1140, "0"),
            ('g g > g g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 11850, "1375"),
            ('d d~ > d d~ g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 29210, "0"),
            ('d d~ > u u~ g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 14605, "-18993/2"),
            ('d d~ > d d~ d d~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 15030, "0"),
            ('d d~ > u u~ d d~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 5010, "0"),
        ]
        TestProcessGeneration.run_tests(gloop, tests)
