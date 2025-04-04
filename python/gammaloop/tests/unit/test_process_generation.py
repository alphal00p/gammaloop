# nopep8
import pytest
from gammaloop.interface.gammaloop_interface import CommandList, GammaLoop
from gammaloop.misc.common import GammaLoopError
from gammaloop.tests.common import get_gamma_loop_interpreter
from pathlib import Path
from gammaloop.misc.utils import parse_python_expression, expression_to_string, evaluate_graph_overall_factor
from symbolica import Expression, E

USE_SNAPSHOT_TESTS = True
if USE_SNAPSHOT_TESTS:
    try:
        from inline_snapshot import snapshot  # type: ignore
    except ImportError:
        USE_SNAPSHOT_TESTS = False
        def snapshot(target): return target  # type: ignore
else:
    def snapshot(target): return target  # type: ignore


class TestProcessGeneration:

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
            all_overall_factors: list[Expression] = [evaluate_graph_overall_factor(of) for of in all_graphs]  # nopep8 # type: ignore
            for of in all_overall_factors:
                total_overall_factor += of
            total_overall_factor_str = expression_to_string(total_overall_factor.expand(), canonical=True)  # nopep8
            assert total_overall_factor_str is not None, f"Could not convert total overall factor to string: '{total_overall_factor}'"  # nopep8

            n_graphs = len(all_graphs)
            assert n_graphs == expected_graph_number, f"For process: '{test}' | Expected {expected_graph_number} graphs, got {n_graphs} (total_overall_factor: {total_overall_factor_str})"  # nopep8
            if expected_total_overall_factor is not None:
                # expected_total_overall_factor_expr = parse_python_expression(
                #     expected_total_overall_factor)
                # assert expected_total_overall_factor_expr is not None, f"Could not parse expected total overall factor: '{expected_total_overall_factor}'"  # nopep8
                # difference = (total_overall_factor -
                #               expected_total_overall_factor_expr).expand()
                # assert (difference == E("0")).eval(), f"For process: '{test}' | Expected total overall factor: '{expected_total_overall_factor}', got '{total_overall_factor_str}'"  # nopep8
                assert expected_total_overall_factor == total_overall_factor_str, f"For process: '{test}' | Expected total overall factor: '{expected_total_overall_factor}', got '{total_overall_factor_str}'"  # nopep8

            gloop.amplitudes.clear()
            gloop.cross_sections.clear()

    def test_generate_sm_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        tests: list[tuple[str, int, str | None]] = [
            # Only d, g and a as particle contents
            # Test the validity of multiple final-state specifications for cross-section generation
            ('a > {d d~, g g} [{{1}}] | d g a -num_grouping only_detect_zeroes', snapshot(1), snapshot("-1")),  # nopep8
            ('a > d d~ g [{{2}}] | d g a -num_grouping only_detect_zeroes', snapshot(3), snapshot("-3")),  # nopep8
            ('a > d d~ g [{{2}}] | d g a --symmetrize_left_right_states', snapshot(2), snapshot("-3")),  # nopep8
            ('a > d d~ z [{{2}}] | d g a -num_grouping only_detect_zeroes', snapshot(0), snapshot("0")),  # nopep8
            # Full particle contents
            ('a > d d~ [{{1}}] --symmetrize_left_right_states', snapshot(1), snapshot("-1")),  # nopep8
            ('a > d d~ [{{2}}] --symmetrize_left_right_states --compare_canonized_numerator --number_of_samples_for_numerator_comparisons 0 --no-fully_numerical_substitution_when_comparing_numerators', snapshot(10), snapshot('-12+-27*G^2*Nc*TR*ee^-2+27*G^2*Nc^-1*TR*ee^-2')),  # nopep8
            ('a > d d~ [{{2}}] --symmetrize_left_right_states --compare_canonized_numerator --number_of_samples_for_numerator_comparisons 0 --fully_numerical_substitution_when_comparing_numerators', snapshot(10), snapshot('-12+-36*G^2*ee^-2')),  # nopep8
            ('a > d d~ [{{2}}] --symmetrize_left_right_states --no-compare_canonized_numerator --number_of_samples_for_numerator_comparisons 3 --fully_numerical_substitution_when_comparing_numerators', snapshot(10), snapshot('-12+-36*G^2*ee^-2')),  # nopep8
            ('a > d d~ [{{2}}] --symmetrize_left_right_states --no-compare_canonized_numerator --number_of_samples_for_numerator_comparisons 3 --no-fully_numerical_substitution_when_comparing_numerators', snapshot(10), snapshot('-12+-3*Nc*TR+3*Nc^-1*TR')),  # nopep8
            # # Only 1-flavour pure QCD corrections
            ('a > d d~ | d g ghG a QED^2==2 [{{1}}] --symmetrize_left_right_states', snapshot(1), snapshot("-1")),  # nopep8
            ('a > d d~ | d g ghG a QED^2==2 [{{2}} QCD=1] --symmetrize_left_right_states', snapshot(2), snapshot("-3")),  # nopep8
            ('a > d d~ | d g ghG a QED^2==2 [{{3}} QCD=2] --symmetrize_left_right_states', snapshot(16), snapshot("-41/2")),  # nopep8
            ('a > d d~ | d g ghG a QED^2==2 [{{4}} QCD=3] --symmetrize_left_right_states', snapshot(166), snapshot("-3107/14")),  # nopep8
        ]
        TestProcessGeneration.run_tests(gloop, tests)

    @ pytest.mark.slow
    def test_slow_generate_sm_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            # Full particle contents
            # Adding --symmetrize_left_right_states below would give:
            # Could not find the CP conjugate of this vertex in the Feynman rules of the model: (ghWm~, G-, ghA). Consider generating without the option '--symmetrize_left_right_states'.
            ('a > d d~ [{{3}}] -num_grouping only_detect_zeroes', snapshot(1321), snapshot("-1103/2")),
            # Only 1-flavour pure QCD corrections
            ('a > d d~ | d g ghG a QED^2==2 [{{5}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', snapshot(6303), snapshot("-51683/24")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_sm_full_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm-full"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            # Full particle contents
            ('a > d d~ [{{1}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', snapshot(1), snapshot("-1")),
            ('a > d d~ [{{2}}] -num_grouping only_detect_zeroes', snapshot(47), snapshot("-47")),   # nopep8
            ('a > d d~ [{{2}}] -num_grouping group_identical_graphs_up_to_sign', snapshot(45), snapshot("-47")),   # nopep8
            ('a > d d~ [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', snapshot(40), snapshot("-47")),  # nopep8
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    @ pytest.mark.slow
    def test_slow_generate_sm_full_a_ddx(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm-full"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            # Full particle contents
            ('a > d d~ QED^2==4 [{{3}}] -num_grouping group_identical_graphs_up_to_sign', snapshot(339), snapshot("-368")),
            ('a > d d~ [{{3}}] -num_grouping no_grouping', snapshot(8549), snapshot("-9111/2")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_sm_h_n_j(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            # Full particle contents
            ('h > g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', snapshot(2), snapshot("-4")),
            ('h > g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', snapshot(30), snapshot("-74")),
            ('h > g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', snapshot(2), snapshot("-4")),
            ('h > g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', snapshot(40), snapshot("-74")),
            ('h > g g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', snapshot(6), snapshot("-24")),
            ('h > g g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', snapshot(188), snapshot("-618")),
            ('h > g g g | h g b t ghg [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', snapshot(4), snapshot("-24")),
            ('h > g g g | h g b t ghg [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', snapshot(132), snapshot("-618")),
            ('h > g g | h g b t ghg [{{3}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', snapshot(8), snapshot("8")),
            ('h > g g | h g b t ghg [{{3}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', snapshot(3), snapshot("8")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    @ pytest.mark.slow
    def test_slow_generate_sm_h_n_j_cross_section(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('h > g g g | h g b t ghg [{{4}}] --symmetrize_left_right_states -num_grouping only_detect_zeroes', snapshot(72), snapshot("88")),
            ('h > g g g | h g b t ghg [{{4}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign', snapshot(26), snapshot("88")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    @ pytest.mark.slow
    def test_slow_generate_sm_h_n_j(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('h > g g | h g b t ghg [{3}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', snapshot(1194), snapshot("-1260")),
            ('h > g g g | h g b t ghg [{3}] -a --symmetrize_left_right_states -num_grouping only_detect_zeroes', snapshot(6946), snapshot("-12429")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_epem_a_qqh(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('e+ e- > d d~ | d a e- ghg g QED^2==4 [{{2}}] -num_grouping only_detect_zeroes', snapshot(3), snapshot("-3")),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{2}}] -num_grouping only_detect_zeroes', snapshot(3), snapshot("-3")),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{3}} QCD=1] -num_grouping only_detect_zeroes', snapshot(30), snapshot("-30")),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{4}} QCD=2] -num_grouping only_detect_zeroes', snapshot(594), snapshot("-261")),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{2}}] -num_grouping only_detect_zeroes', snapshot(21), snapshot("-21")),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{3}} QCD=2] -num_grouping only_detect_zeroes', snapshot(171), snapshot("-171")),
            # Enable grouping
            ('e+ e- > d d~ | d a e- ghg g QED^2==4 [{{2}}] -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(2), snapshot("-3")),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{2}}] -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(2), snapshot("-3")),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{3}} QCD=1] -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(16), snapshot("-30")),
            # This grouping involves 21 new cancellations between graphs, so that modifies the total overall factor
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{4}} QCD=2] -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(266), snapshot("-327")),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{2}}] -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(20), snapshot("-21")),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{3}} QCD=1] -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(157), snapshot("-171")),
            # Enable left-right-symmetry
            ('e+ e- > d d~ | d a e- ghg g QED^2==4 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(2), snapshot("-3")),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(2), snapshot("-3")),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{3}} QCD=1] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(11), snapshot("-30")),
            ('e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{4}} QCD=2] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(166), snapshot("-327")),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(16), snapshot("-21")),
            ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{3}} QCD=1] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(104), snapshot("-171")),
            # Too slow and RAM intensive
            # ('e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{4}} QCD=2] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', ?),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_a_qqh(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('a > d d~ | d a ghg g QED^2==2 [{{2}}] -num_grouping only_detect_zeroes', snapshot(3), snapshot("-3")),
            ('a > b b~ h | d b h a ghg g QED^2==4 [{{2}}] -num_grouping only_detect_zeroes', snapshot(3), snapshot("-3")),
            ('a > b b~ h | d b h a ghg g QED^2==4 [{{3}} QCD=1] -num_grouping only_detect_zeroes', snapshot(30), snapshot("-30")),
            ('a > b b~ h | d b h a ghg g QED^2==4 [{{4}} QCD=2] -num_grouping only_detect_zeroes', snapshot(594), snapshot("-261")),
            # Enable left-right-symmetry and grouping
            ('a > d d~ | d a ghg g QED^2==2 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(2), snapshot("-3")),
            ('a > b b~ h | d b h a ghg g QED^2==4 [{{2}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(2), snapshot("-3")),
            ('a > b b~ h | d b h a ghg g QED^2==4 [{{3}} QCD=1] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(11), snapshot("-30")),
            # Additional cancellations, so that modifies the total overall factor
            ('a > b b~ h | d b h a ghg g QED^2==4 [{{4}} QCD=2] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(166), snapshot("-327")),
            ('a > b b~ h | b h a ghg g QED^2==4 [{{4}} QCD=2] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(151), snapshot("-366")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_aa_ttx_amplitude(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('a a > t t~ | a t g b ghg QED==2 [{0} QCD=0] -a -num_grouping only_detect_zeroes', snapshot(2), snapshot("2")),
            ('a a > t t~ | a t g b ghg QED==2 [{1} QCD=1] -a -num_grouping only_detect_zeroes', snapshot(8), snapshot("8")),
            ('a a > t t~ | a t g b ghg QED==2 [{2} QCD=2] -a -num_grouping only_detect_zeroes', snapshot(144), snapshot("36")),
            ('a a > t t~ | a t g b ghg QED==2 [{3} QCD=3] -a -num_grouping only_detect_zeroes', snapshot(3424), snapshot("28/3")),
            ('a a > t t~ | a t g b ghg QED==2 [{0} QCD=0] -a --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(1), snapshot("2")),
            ('a a > t t~ | a t g b ghg QED==2 [{1} QCD=1] -a --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(4), snapshot("8")),
            # Six additional cancellations imply that the total overall factor is different
            ('a a > t t~ | a t g b ghg QED==2 [{2} QCD=2] -a --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(58), snapshot("60")),
            ('a a > t t~ | a t g b ghg QED==2 [{3} QCD=3] -a --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(1221), snapshot("1600/3")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_aa_ttx_cross_section(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] -num_grouping only_detect_zeroes', snapshot(4), snapshot("-4")),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] -num_grouping only_detect_zeroes', snapshot(40), snapshot("-40")),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] -num_grouping only_detect_zeroes', snapshot(874), snapshot("-266")),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(2), snapshot("-4")),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(14), snapshot("-40")),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] --symmetrize_initial_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(223), snapshot("-426")),  # was 219 # nopep8
            ('a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] --symmetrize_initial_states --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(2), snapshot("-4")),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] --symmetrize_initial_states --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(10), snapshot("-40")),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] --symmetrize_initial_states --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(143), snapshot("-426")),
            ('a a > t t~ | a t g ghg QED^2==4 [{{3}} QCD=2] --symmetrize_initial_states --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(128), snapshot("-502")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_aa_ttx_cross_section_slow_filter(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] -num_grouping only_detect_zeroes --max_multiplicity_for_fast_cut_filter 0', snapshot(4), snapshot("-4")),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] -num_grouping only_detect_zeroes --max_multiplicity_for_fast_cut_filter 0', snapshot(40), snapshot("-40")),
            ('a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] -num_grouping only_detect_zeroes --max_multiplicity_for_fast_cut_filter 0', snapshot(874), snapshot("-266")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    def test_vaccuum_amplitude_generation(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            # Including all graphs
            ('{} > {} | g ghg t u d QED==0 [{2}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', snapshot(6), snapshot("-19/24")),
            ('{} > {} | g ghg t u d QED==0 [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', snapshot(5), snapshot("-19/24")),
            ('{} > {} | g ghg t u d QED==0 [{3}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', snapshot(36), snapshot('-83/48')),
            ('{} > {} | g ghg t u d QED==0 [{3}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', snapshot(27), snapshot('-83/48')),
            ('{} > {} | g ghg t u d QED==0 [{4}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', snapshot(264), snapshot('-143/24')),
            ('{} > {} | g ghg t u d QED==0 [{4}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', snapshot(179), snapshot('-731/168')),
            # Only non-factorizable ones
            ('{} > {} | g ghg t u d QED==0 [{2}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', snapshot(5), snapshot("-11/12")),
            ('{} > {} | g ghg t u d QED==0 [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', snapshot(4), snapshot("-11/12")),
            ('{} > {} | g ghg t u d QED==0 [{3}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', snapshot(29), snapshot('-23/16')),
            ('{} > {} | g ghg t u d QED==0 [{3}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', snapshot(21), snapshot('-23/16')),
            ('{} > {} | g ghg t u d QED==0 [{4}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', snapshot(212), snapshot('-135/32')),
            ('{} > {} | g ghg t u d QED==0 [{4}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', snapshot(138), snapshot('-585/224')),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    @ pytest.mark.slow
    def test_slow_vaccuum_amplitude_generation(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            # CURRENTLY BUGGED BECAUSE OF SYMBOLICA INCORRECT SYMMETRY FACTORS. POSSIBLY OTHER SCENARIOS BUGGED TOO.
            ('{} > {} | g ghg t u d QED==0 [{5}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', snapshot(2560), snapshot("-5785/384")),
            ('{} > {} | g ghg t u d QED==0 [{5}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', snapshot(1440), snapshot("233015/51072")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    def test_vaccuum_amplitude_generation_full_sm(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm-full"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            # Including all graphs
            ('{} > {} [{2}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', snapshot(448), snapshot("-19/12")),

            # Compare various numerator comparison options below:
            # a) Using only symbolic tensor canonization of the numerator and without substituting color group factors, some groupings are missed (i.e. 293 graphs instead of 243 when grouping as much as possible)
            ('{} > {} [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1 --compare_canonized_numerator --number_of_samples_for_numerator_comparisons 0 --fully_numerical_substitution_when_comparing_numerators', snapshot(281),
             snapshot('-1*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-1*CKM2x1^-1*CKM2x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-1*CKM3x1^-1*CKM3x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*I2x12^-1*I2x22*I3x21^-1*I3x22+-1*I2x13^-1*I2x23*I3x31^-1*I3x32+-27/2*G^2*ee^-2+319/24')),
            # b) Using only symbolic tensor canonization of the numerator but this time substituting color group factors allows for more groupings, but still not the maximum (i.e. 281 graphs instead of 243 when grouping as much as possible)
            ('{} > {} [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1 --compare_canonized_numerator --number_of_samples_for_numerator_comparisons 0 --no-fully_numerical_substitution_when_comparing_numerators', snapshot(281),
             snapshot('-1*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-1*CKM2x1^-1*CKM2x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-1*CKM3x1^-1*CKM3x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*I2x12^-1*I2x22*I3x21^-1*I3x22+-1*I2x13^-1*I2x23*I3x31^-1*I3x32+-81/8*G^2*Nc*TR*ee^-2+319/24+81/8*G^2*Nc^-1*TR*ee^-2')),
            # c) Now instead rely on numerical sample evalutions but still with symbolic color group factors
            # Generation below yields "Could not numerically evaluate numerator: more than one node in the graph", which is a current issue in EW feynman rules when kept with symbolic parameters 
            # ('{} > {} [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1 --no-compare_canonized_numerator --number_of_samples_for_numerator_comparisons 3 --no-fully_numerical_substitution_when_comparing_numerators', snapshot(?), snapshot('?')),
            # d) Finally using fully numerical evaluations of the numerator, yielding the largest amount of groupings (i.e. 243 graphs)
            ('{} > {} [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1 --no-compare_canonized_numerator --number_of_samples_for_numerator_comparisons 3 --fully_numerical_substitution_when_comparing_numerators', snapshot(243),
             snapshot('-1*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-1*CKM2x1^-1*CKM2x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-1*CKM3x1^-1*CKM3x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*I2x12^-1*I2x22*I3x21^-1*I3x22+-1*I2x13^-1*I2x23*I3x31^-1*I3x32+-27/2*G^2*ee^-2+313/24')),

            ('{} > {} [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', snapshot(243),
             snapshot('-1*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-1*CKM2x1^-1*CKM2x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-1*CKM3x1^-1*CKM3x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*I2x12^-1*I2x22*I3x21^-1*I3x22+-1*I2x13^-1*I2x23*I3x31^-1*I3x32+-27/2*G^2*ee^-2+313/24')),
            ('{} > {} [{3}] -a -num_grouping only_detect_zeroes --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', snapshot(36362), snapshot('-68')),
            ('{} > {} [{3}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges -1 --number_of_factorized_loop_subtopologies -1', snapshot(13869),
             snapshot('(-1*CKM2x1*I2x12*cw^2+CKM2x1*I2x12*sw^2)^-1*1/4*CKM2x2*I2x22*cw^2+(-1*CKM2x1*I2x12*cw^2+CKM2x1*I2x12*sw^2)^-1*1/4*CKM2x2*I2x22*sw^2+(-1*CKM3x1*I2x13*cw^2+CKM3x1*I2x13*sw^2)^-1*1/4*CKM3x2*I2x23*cw^2+(-1*CKM3x1*I2x13*cw^2+CKM3x1*I2x13*sw^2)^-1*1/4*CKM3x2*I2x23*sw^2+(-1*I3x21*complexconjugate(CKM2x1)*cw^2+I3x21*complexconjugate(CKM2x1)*sw^2)^-1*1/4*I3x22*complexconjugate(CKM2x2)*cw^2+(-1*I3x21*complexconjugate(CKM2x1)*cw^2+I3x21*complexconjugate(CKM2x1)*sw^2)^-1*1/4*I3x22*complexconjugate(CKM2x2)*sw^2+(-1*I3x31*complexconjugate(CKM3x1)*cw^2+I3x31*complexconjugate(CKM3x1)*sw^2)^-1*1/4*I3x32*complexconjugate(CKM3x2)*cw^2+(-1*I3x31*complexconjugate(CKM3x1)*cw^2+I3x31*complexconjugate(CKM3x1)*sw^2)^-1*1/4*I3x32*complexconjugate(CKM3x2)*sw^2+(-1*cw^2+sw^2)^-1*-31/2*cw^2+(-1*cw^2+sw^2)^-1*-31/2*sw^2+(-1/2*cw^4+-1/2*sw^4+cw^2*sw^2)^-1*-1*cw^2*sw^2+(-1/2*cw^4+-1/2*sw^4+cw^2*sw^2)^-1*-1/2*cw^4+(-1/2*cw^4+-1/2*sw^4+cw^2*sw^2)^-1*-1/2*sw^4+-1*CKM1x1^-1*CKM1x2*CKM2x1^-1*CKM2x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-1*CKM1x1^-1*CKM1x2*CKM3x1^-1*CKM3x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*CKM1x1^-1*CKM1x2*I2x12^-1*I2x22*I3x21^-1*I3x22*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-1*CKM1x1^-1*CKM1x2*I2x13^-1*I2x23*I3x31^-1*I3x32*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-1*CKM2x1^-1*CKM2x2*CKM3x1^-1*CKM3x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*CKM2x1^-1*CKM2x2*I2x12^-1*I2x22*I3x21^-1*I3x22*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-1*CKM2x1^-1*CKM2x2*I2x13^-1*I2x23*I3x31^-1*I3x32*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-1*CKM3x1^-1*CKM3x2*I2x12^-1*I2x22*I3x21^-1*I3x22*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*CKM3x1^-1*CKM3x2*I2x13^-1*I2x23*I3x31^-1*I3x32*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*I2x12^-1*I2x13^-1*I2x22*I2x23*I3x21^-1*I3x22*I3x31^-1*I3x32+-1/2*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-1/2*CKM1x1^-2*CKM1x2^2*complexconjugate(CKM1x1)^-2*complexconjugate(CKM1x2)^2+-1/2*CKM2x1^-2*CKM2x2^2*complexconjugate(CKM2x1)^-2*complexconjugate(CKM2x2)^2+-1/2*CKM3x1^-2*CKM3x2^2*complexconjugate(CKM3x1)^-2*complexconjugate(CKM3x2)^2+-1/2*I2x12^-2*I2x22^2*I3x21^-2*I3x22^2+-1/2*I2x13^-2*I2x23^2*I3x31^-2*I3x32^2+-1053/16*G^4*ee^-4+-2*CKM1x1^-1*CKM1x2*I2x12^-1*I2x22+-2*CKM1x1^-1*CKM1x2*I2x13^-1*I2x23+-2*CKM2x1^-1*CKM2x2*I2x13^-1*I2x23+-2*CKM3x1^-1*CKM3x2*I2x12^-1*I2x22+-2*I3x21^-1*I3x22*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-2*I3x21^-1*I3x22*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-2*I3x31^-1*I3x32*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-2*I3x31^-1*I3x32*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-5*CKM1x1^-1*CKM1x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-5*CKM1x1^-1*CKM1x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-5*CKM2x1^-1*CKM2x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-5*CKM3x1^-1*CKM3x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-5/2*CKM1x3*I1x31*I4x13^-1*complexconjugate(CKM1x3)^-1+-5/2*CKM2x1*I2x12*I3x21^-1*complexconjugate(CKM2x1)^-1+-5/2*CKM2x2*I2x22*I3x21^-1*complexconjugate(CKM2x1)^-1+-5/2*CKM3x1*I2x13*I3x31^-1*complexconjugate(CKM3x1)^-1+-5/2*CKM3x2*I2x23*I3x31^-1*complexconjugate(CKM3x1)^-1+-6*CKM2x1^-1*CKM2x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-6*CKM3x1^-1*CKM3x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-6*I2x12^-1*I2x22*I3x31^-1*I3x32+-6*I2x13^-1*I2x23*I3x21^-1*I3x22+-9*CKM2x1^-1*CKM2x2*G^2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)*ee^-2+-9*CKM3x1^-1*CKM3x2*G^2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)*ee^-2+-9*G^2*I2x12^-1*I2x22*I3x21^-1*I3x22*ee^-2+-9*G^2*I2x13^-1*I2x23*I3x31^-1*I3x32*ee^-2+11/2*I3x21^-1*I3x22*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+11/2*I3x31^-1*I3x32*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+31/4*CKM2x1^-1*CKM2x2*I2x12^-1*I2x22+31/4*CKM3x1^-1*CKM3x2*I2x13^-1*I2x23+4525/16+5/2*CKM2x1^-1*CKM2x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+5/2*CKM3x1^-1*CKM3x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+7/2*I2x12^-1*I2x22*I3x21^-1*I3x22+7/2*I2x13^-1*I2x23*I3x31^-1*I3x32+81*G^2*ee^-2')),
            # Only non-factorizable ones
            ('{} > {} [{2}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', snapshot(92), snapshot("-67/3")),
            ('{} > {} [{2}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', snapshot(69),
             snapshot( '-1*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-1*CKM2x1^-1*CKM2x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-1*CKM3x1^-1*CKM3x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*I2x12^-1*I2x22*I3x21^-1*I3x22+-1*I2x13^-1*I2x23*I3x31^-1*I3x32+-163/12+-27/2*G^2*ee^-2')),
            ('{} > {} [{3}] -a -num_grouping only_detect_zeroes --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', snapshot(3085), snapshot('-1687/3')),
            ('{} > {} [{3}] -a -num_grouping group_identical_graphs_up_to_scalar_rescaling --max_n_bridges 0 --number_of_factorized_loop_subtopologies 1', snapshot(2009), snapshot('(-1*CKM2x1*I2x12*cw^2+CKM2x1*I2x12*sw^2)^-1*1/4*CKM2x2*I2x22*cw^2+(-1*CKM2x1*I2x12*cw^2+CKM2x1*I2x12*sw^2)^-1*1/4*CKM2x2*I2x22*sw^2+(-1*CKM3x1*I2x13*cw^2+CKM3x1*I2x13*sw^2)^-1*1/4*CKM3x2*I2x23*cw^2+(-1*CKM3x1*I2x13*cw^2+CKM3x1*I2x13*sw^2)^-1*1/4*CKM3x2*I2x23*sw^2+(-1*I3x21*complexconjugate(CKM2x1)*cw^2+I3x21*complexconjugate(CKM2x1)*sw^2)^-1*1/4*I3x22*complexconjugate(CKM2x2)*cw^2+(-1*I3x21*complexconjugate(CKM2x1)*cw^2+I3x21*complexconjugate(CKM2x1)*sw^2)^-1*1/4*I3x22*complexconjugate(CKM2x2)*sw^2+(-1*I3x31*complexconjugate(CKM3x1)*cw^2+I3x31*complexconjugate(CKM3x1)*sw^2)^-1*1/4*I3x32*complexconjugate(CKM3x2)*cw^2+(-1*I3x31*complexconjugate(CKM3x1)*cw^2+I3x31*complexconjugate(CKM3x1)*sw^2)^-1*1/4*I3x32*complexconjugate(CKM3x2)*sw^2+(-1*cw^2+sw^2)^-1*-9*cw^2+(-1*cw^2+sw^2)^-1*-9*sw^2+(-1/2*cw^4+-1/2*sw^4+cw^2*sw^2)^-1*-1*cw^2*sw^2+(-1/2*cw^4+-1/2*sw^4+cw^2*sw^2)^-1*-1/2*cw^4+(-1/2*cw^4+-1/2*sw^4+cw^2*sw^2)^-1*-1/2*sw^4+-1*CKM1x1^-1*CKM1x2*CKM2x1^-1*CKM2x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-1*CKM1x1^-1*CKM1x2*CKM3x1^-1*CKM3x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*CKM1x1^-1*CKM1x2*I2x12^-1*I2x22*I3x21^-1*I3x22*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-1*CKM1x1^-1*CKM1x2*I2x13^-1*I2x23*I3x31^-1*I3x32*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-1*CKM2x1^-1*CKM2x2*CKM3x1^-1*CKM3x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*CKM2x1^-1*CKM2x2*I2x12^-1*I2x22*I3x21^-1*I3x22*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-1*CKM2x1^-1*CKM2x2*I2x13^-1*I2x23*I3x31^-1*I3x32*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-1*CKM3x1^-1*CKM3x2*I2x12^-1*I2x22*I3x21^-1*I3x22*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*CKM3x1^-1*CKM3x2*I2x13^-1*I2x23*I3x31^-1*I3x32*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-1*I2x12^-1*I2x13^-1*I2x22*I2x23*I3x21^-1*I3x22*I3x31^-1*I3x32+-1/2*CKM1x1^-2*CKM1x2^2*complexconjugate(CKM1x1)^-2*complexconjugate(CKM1x2)^2+-1/2*CKM2x1^-2*CKM2x2^2*complexconjugate(CKM2x1)^-2*complexconjugate(CKM2x2)^2+-1/2*CKM3x1^-2*CKM3x2^2*complexconjugate(CKM3x1)^-2*complexconjugate(CKM3x2)^2+-1/2*I2x12^-2*I2x22^2*I3x21^-2*I3x22^2+-1/2*I2x13^-2*I2x23^2*I3x31^-2*I3x32^2+-1053/16*G^4*ee^-4+-2*CKM1x1^-1*CKM1x2*I2x12^-1*I2x22+-2*CKM1x1^-1*CKM1x2*I2x13^-1*I2x23+-2*CKM2x1^-1*CKM2x2*I2x13^-1*I2x23+-2*CKM3x1^-1*CKM3x2*I2x12^-1*I2x22+-2*I3x21^-1*I3x22*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-2*I3x21^-1*I3x22*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-2*I3x31^-1*I3x32*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-2*I3x31^-1*I3x32*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-25787/96+-5*CKM1x1^-1*CKM1x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-5*CKM1x1^-1*CKM1x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-5*CKM2x1^-1*CKM2x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-5*CKM3x1^-1*CKM3x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-5/2*CKM1x3*I1x31*I4x13^-1*complexconjugate(CKM1x3)^-1+-5/2*CKM2x1*I2x12*I3x21^-1*complexconjugate(CKM2x1)^-1+-5/2*CKM2x2*I2x22*I3x21^-1*complexconjugate(CKM2x1)^-1+-5/2*CKM3x1*I2x13*I3x31^-1*complexconjugate(CKM3x1)^-1+-5/2*CKM3x2*I2x23*I3x31^-1*complexconjugate(CKM3x1)^-1+-6*CKM2x1^-1*CKM2x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-6*CKM3x1^-1*CKM3x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-6*I2x12^-1*I2x22*I3x31^-1*I3x32+-6*I2x13^-1*I2x23*I3x21^-1*I3x22+-63/4*CKM1x1^-1*CKM1x2*complexconjugate(CKM1x1)^-1*complexconjugate(CKM1x2)+-77/4*I2x12^-1*I2x22*I3x21^-1*I3x22+-77/4*I2x13^-1*I2x23*I3x31^-1*I3x32+-81/4*CKM2x1^-1*CKM2x2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+-81/4*CKM3x1^-1*CKM3x2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+-9*CKM2x1^-1*CKM2x2*G^2*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)*ee^-2+-9*CKM3x1^-1*CKM3x2*G^2*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)*ee^-2+-9*G^2*I2x12^-1*I2x22*I3x21^-1*I3x22*ee^-2+-9*G^2*I2x13^-1*I2x23*I3x31^-1*I3x32*ee^-2+-981/4*G^2*ee^-2+5/4*I3x21^-1*I3x22*complexconjugate(CKM2x1)^-1*complexconjugate(CKM2x2)+5/4*I3x31^-1*I3x32*complexconjugate(CKM3x1)^-1*complexconjugate(CKM3x2)+7/2*CKM2x1^-1*CKM2x2*I2x12^-1*I2x22+7/2*CKM3x1^-1*CKM3x2*I2x13^-1*I2x23')),

        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    def test_dis_isr(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('e- d > e- d | e- a d g ghg QED==2 [{0}] -a -num_grouping only_detect_zeroes', snapshot(1), snapshot("-1")),
            ('e- d > e- d | e- a d g ghg QED==2 [{1}] -a -num_grouping only_detect_zeroes', snapshot(1), snapshot("-1")),
            ('e- d > e- d | e- a d g ghg QED==2 [{2}] -a -num_grouping only_detect_zeroes', snapshot(13), snapshot("-13/2")),
            ('e- d > e- d | e- a d g ghg QED==2 [{3}] -a -num_grouping only_detect_zeroes', snapshot(241), snapshot("-194/3")),
            ('e- d > e- d | e- a d g ghg QED==2 [{0}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(1), snapshot("-1")),
            ('e- d > e- d | e- a d g ghg QED==2 [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(1), snapshot("-1")),
            ('e- d > e- d | e- a d g ghg QED==2 [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(11), snapshot("-17/2")),
            ('e- d > e- d | e- a d g ghg QED==2 [{3}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(191), snapshot("-2183/21")),
            # 3>N no symmetrization
            ('e- d g > e- d | e- a d g ghg QED==2 [{0}] -a -num_grouping only_detect_zeroes', snapshot(2), snapshot("-2")),
            ('e- d g > e- d | e- a d g ghg QED==2 [{1}] -a -num_grouping only_detect_zeroes', snapshot(13), snapshot("-9")),
            ('e- d g > e- d | e- a d g ghg QED==2 [{2}] -a -num_grouping only_detect_zeroes', snapshot(225), snapshot("-169/2")),
            # 3>N with symmetrization
            ('e- d g > e- d | e- a d g ghg QED==2 [{0}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(2), snapshot("-2")),
            ('e- d g > e- d | e- a d g ghg QED==2 [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(11), snapshot("-11")),
            ('e- d g > e- d | e- a d g ghg QED==2 [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(181), snapshot("-842/7")),
            # 4>N no symmetrization
            ('e- d u u~ > e- d | e- a d u g ghg QED==2 [{0}] -a -num_grouping only_detect_zeroes', snapshot(4), snapshot("4")),
            ('e- d u u~ > e- d | e- a d u g ghg QED==2 [{1}] -a -num_grouping only_detect_zeroes', snapshot(62), snapshot("36")),
            ('e- d u u~ > e- d | e- a d u g ghg QED==2 [{2}] -a -num_grouping only_detect_zeroes', snapshot(1428), snapshot("1298/3")),
            # 4>N with symmetrization
            ('e- d u u~ > e- d | e- a d u g ghg QED==2 [{0}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(4), snapshot("4")),
            ('e- d u u~ > e- d | e- a d u g ghg QED==2 [{1}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(54), snapshot("40")),
            ('e- d u u~ > e- d | e- a d u g ghg QED==2 [{2}] -a --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(1087), snapshot("22357/42"))
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    @ pytest.mark.slow
    def test_slow_generate_a_qqh(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('a > b b~ h | b h a ghg g QED^2==4 [{{5}} QCD=3] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_scalar_rescaling', snapshot(2763), snapshot("-38712/7")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    def test_generate_amplitude_1l_sm_jets(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('d d~ > d d~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(20), snapshot("0")),
            ('d d~ > d d~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(192), snapshot("0")),
            ('d d~ > g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(384), snapshot("-186")),
            ('g g > g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(1139), snapshot("171/2")),
            ('d d~ > d d~ g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(2284), snapshot("0")),
            ('d d~ > u u~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(10), snapshot("-11/2")),
            ('d d~ > u u~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(96), snapshot("-119/2")),
            ('d d~ > u u~ g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(1142), snapshot("-698")),
            ('d d~ > u u~ d d~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(414), snapshot("0")),
            ('d d~ > d d~ d d~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(1242), snapshot("0")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    # Currently bugged with spenso network parsing for "-num_grouping group_identical_graphs_up_to_sign"
    def test_generate_amplitude_1l_sm_jets_with_grouping(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('d d~ > d d~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(18), snapshot("0")),
            ('d d~ > u u~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(9), snapshot("-11/2")),
            ('d d~ > u u~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(88), snapshot("-119/2")),
            ('d d~ > d d~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(176), snapshot("0")),
            ('d d~ > g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(341), snapshot("-186")),
            ('g g > g g g | g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', 671, snapshot("1107/2")),
            # For fun, symmetrize it all
            ('g g > g g g | g ghg a QED==0 [QCD=1] -a --symmetrize_initial_states --symmetrize_final_states -num_grouping group_identical_graphs_up_to_sign', snapshot(109), snapshot("1107/2")),
            ('g g > g g g | g ghg a QED==0 [QCD=1] -a --symmetrize_left_right_states --symmetrize_initial_states --symmetrize_final_states -num_grouping group_identical_graphs_up_to_sign', snapshot(25), snapshot("1107/2")),
            ('g g > g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(905), snapshot("171/2")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    @ pytest.mark.slow
    def test_slow_generate_amplitude_1l_sm_jets(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('d d~ > u u~ d d~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(5424), snapshot("0")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    # A bit too slow, let's skip for now
    @ pytest.mark.skip
    @ pytest.mark.slow
    def test_very_slow_generate_amplitude_1l_sm_jets(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('g g > g g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(14875), snapshot("1375")),
            ('d d~ > d d~ g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(32074), snapshot("0")),
            ('d d~ > u u~ g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(16037), snapshot("-18993/2")),
            ('d d~ > d d~ d d~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping only_detect_zeroes', snapshot(16272), snapshot("0")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)

    @ pytest.mark.skip
    @ pytest.mark.slow
    def test_slow_generate_amplitude_1l_sm_jets_with_grouping(self):
        gloop = get_gamma_loop_interpreter()
        gloop.run(CommandList.from_string(
            "import_model sm"))
        # Targets confirmed by MadGraph
        # autopep8: off
        tests: list[tuple[str, int, str | None]] = [
            ('d d~ > d d~ g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(2090), snapshot("0")),
            ('d d~ > u u~ g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(1045), snapshot("-698")),
            ('d d~ > u u~ d d~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(380), snapshot("0")),
            ('d d~ > d d~ d d~ | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(1140), snapshot("0")),
            ('g g > g g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(11850), snapshot("1375")),
            ('d d~ > d d~ g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(29210), snapshot("0")),
            ('d d~ > u u~ g g g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(14605), snapshot("-18993/2")),
            ('d d~ > d d~ d d~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(15030), snapshot("0")),
            ('d d~ > u u~ d d~ g | u d g ghg a QED==0 [QCD=1] -a -num_grouping group_identical_graphs_up_to_sign', snapshot(5010), snapshot("0")),
        ]
        # autopep8: on
        TestProcessGeneration.run_tests(gloop, tests)
