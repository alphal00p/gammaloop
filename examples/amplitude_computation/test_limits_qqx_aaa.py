#!/usr/bin/env python3

import argparse
import polars as pl
import subprocess
import os

import re

from math import log10

_NUM = r"(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

# Complex: ( re  ±  [optional sign] im i )   OR   ( re ,  [optional sign] im i )
RX_COMPLEX = re.compile(
    rf"""
    \(\s*
      (?P<re>[+-]?{_NUM})
      \s*
      (?:
        (?P<pm>[+-])\s*(?P<imsign>[+-]?)\s*(?P<im_a>{_NUM})\s*[ij]   # a±±b i
        |
        ,\s*(?P<im_b_sign>[+-]?)\s*(?P<im_b>{_NUM})\s*[ij]           # a, ±b i
      )
    \s*\)?                                                           # optional ')'
    """,
    re.VERBOSE | re.IGNORECASE,
)

# Real number after '=' at line end, e.g. "... = 1.23e4"
RX_REAL_EQ_END = re.compile(rf"=\s*(?P<real>[+-]?{_NUM})\s*$")

# Any standalone real number (used as last-resort fallback)
RX_ANY_NUM = re.compile(rf"[+-]?{_NUM}")


def extract_complex(s: str) -> complex:
    m = RX_COMPLEX.search(s)
    if m:
        re_part = float(m.group("re"))
        if m.group("im_a") is not None:
            im_val = float(m.group("im_a"))
            neg = (m.group("pm") == "-") ^ (m.group("imsign") == "-")
            if neg:
                im_val = -im_val
        else:
            im_val = float(m.group("im_b"))
            if m.group("im_b_sign") == "-":
                im_val = -im_val
        return complex(re_part, im_val)

    m = RX_REAL_EQ_END.search(s)
    if m:
        return complex(float(m.group("real")), 0.0)

    nums = RX_ANY_NUM.findall(s)
    if nums:
        return complex(float(nums[-1]), 0.0)

    raise ValueError("No number found")


def get_integrand_results_old(stdout):
    all_evals = []
    total = None
    total_read = None
    jac = None
    for l in stdout.split('\n'):
        if 'evaluated integrand:' in l:
            all_evals.append(extract_complex(l))
        if "f128 sampling jacobian for this point" in l:
            if jac is not None:
                raise ValueError("Jacobian already set, but found another one")
            jac = extract_complex(l)
        if l.startswith("The evaluation of integrand '"):
            total_read = 0
        if total_read is True:
            total = extract_complex(l)
            total_read = None
        elif total_read is not None:
            total_read += 1
            if total_read == 2:
                total_read = True

    return all_evals, total, jac


graph_names = [
    "qqx_aaa_pentagon",
    "qqx_aaa_box_A",
    "qqx_aaa_box_B",
    "qqx_aaa_tri_A",
    "uv_triangle_A",
    "qqx_aaa_tri_B",
    "uv_triangle_B",
    "qqx_aaa_tri_C",
    "uv_triangle_C",
    "qqx_aaa_bub_A",
    "uv_bubble_A",
    "qqx_aaa_bub_B",
    "uv_bubble_B",
    "ir_ct",
    "uv_ir_ct",
    "int_ir_bare_ct",
    "int_ir_dressed_ct",
    "int_uv_ir_ct",
    "int_uv_ct",
]


def get_integrand_results(log_file_path):
    df = pl.read_ndjson(log_file_path)

    has_span = "span" in df.columns
    has_spans = "spans" in df.columns

    all_evals = {}
    for g_name in graph_names:
        all_evals[g_name] = {'I': complex(
            0.0, 0.0), 'threshold_CT': complex(0.0, 0.0)}
        pred_span = (
            pl.col("span").struct.field("term.name") == g_name
        ) if has_span else pl.lit(False)

        pred_spans = (
            pl.col("spans")
            .list.eval(pl.element().struct.field("term.name") == g_name)
            .list.any()
        ) if has_spans else pl.lit(False)

        g_df = df.filter(pred_span | pred_spans)

        g_original_integrand = g_df.filter(
            pl.col("message") == "Original integrand value")
        if g_original_integrand.is_empty():
            raise ValueError(
                f"No original integrand entries found for graph '{g_name}'")
        if g_original_integrand.height > 1:
            raise ValueError(
                f"Multiple original integrand entries found for graph '{g_name}', expected only one")
        eval_str = g_original_integrand[0, "bare_eval"]
        all_evals[g_name]['I'] = extract_complex(eval_str)

        g_ct_eval = g_df.filter(pl.col("message") ==
                                "evaluated sum of threshold counterterms")

        if g_ct_eval.is_empty():
            raise ValueError(
                f"No evaluated threshold counterterm entries found for graph '{g_name}'")
        if g_ct_eval.height > 1:
            raise ValueError(
                f"Multiple evaluated threshold counterterm entries found for graph '{g_name}', expected only one")
        eval_str = g_ct_eval[0, "value"]
        all_evals[g_name]['threshold_CT'] = -extract_complex(eval_str)

    # inspect_df = df.filter(pl.col("target") == "gammalooprs::inspect")
    inspect_df = df.filter(pl.col("target") == "status")

    total_df = inspect_df.filter(
        pl.col("message").str.contains("The evaluation of integrand"))
    if total_df.is_empty():
        raise ValueError("No entries found for total integrand")
    if total_df.height > 1:
        raise ValueError(
            "Multiple entries found for total integrand, expected only one")
    total_str = total_df[0, "message"]
    total = extract_complex(total_str)

    # inspect_df = df.filter(pl.col("target") == "gammalooprs::inspect")
    inspect_df = df.filter(pl.col("target") == "status")

    jac_df = inspect_df.filter(pl.col("message").str.contains(
        "f128 sampling jacobian for this point"))
    if jac_df.is_empty():
        raise ValueError("No entries found for Jacobian")
    if jac_df.height > 1:
        raise ValueError(
            "Multiple entries found for Jacobian, expected only one")
    jac_str = jac_df[0, "message"]
    jac = extract_complex(jac_str)
    # print("jac ", jac)
    # print("total ", total)

    return all_evals, total, jac


def eval_integrand(k_input, debug=False, gammaloop_state='./GL_QQX_AAA_euclidean', gammaloop_executable='../../target/dev-optim/gammaloop'):
    if os.path.exists(f'{gammaloop_state}/logs/gammalog-inspect.jsonl'):
        os.remove(f'{gammaloop_state}/logs/gammalog-inspect.jsonl')
    cmd = f'{gammaloop_executable} -s {gammaloop_state} -n -t inspect inspect --process-id 0 --name default --discrete-dim 0 -m -p {
        k_input[0]:.16f} {k_input[1]:.16f} {k_input[2]:.16f}'
    if debug:
        print(f"Running command:\n{cmd}")
    tmp_env = os.environ.copy()
    tmp_env['GL_LOGFILE_FILTER'] = 'debug'
    if debug:
        tmp_env['GL_DISPLAY_FILTER'] = 'debug'
    else:
        tmp_env['GL_DISPLAY_FILTER'] = 'info'
    r = subprocess.run(cmd, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, shell=True, env=tmp_env)
    if debug:
        print(f"Output:\n{r.stdout}")
    all_res, total, jac = get_integrand_results(
        f'{gammaloop_state}/logs/gammalog-inspect.jsonl')

    reconstructed_total = sum(v['I']+v['threshold_CT']
                              for v in all_res.values())*jac
    if debug:
        for g_name in graph_names:
            res = all_res[g_name]['I'] + all_res[g_name]['threshold_CT']
            print(f"{g_name+' [I]':30s}: {all_res[g_name]
                  ['I'].real:+.16e} {all_res[g_name]['I'].imag:+.16e}")
            print(f"{g_name+' [th.CT]':30s}: {all_res[g_name]['threshold_CT'].real:+.16e} {
                  all_res[g_name]['threshold_CT'].imag:+.16e}")
            print(f"{g_name+' [tot]':30s}: {res.real:+.16e} {res.imag:+.16e}")
        print(f"{'Jacobian':30s}: {jac.real:+.16e} {jac.imag:+.16e}")

        print(f"{'Reconstructed total':30s}: {reconstructed_total.real:+.16e} {reconstructed_total.imag:+.16e}")  # nopep8
        print(f"{'Total':30s}: {total.real:+.16e} {total.imag:+.16e}")
    return all_res, jac, total, reconstructed_total


def run_limit_test(elements_to_display=None, tests_to_run=None, gammaloop_state='./GL_QQX_AAA_euclidean', gammaloop_executable='../../target/dev-optim/gammaloop', debug=False, threshold_to_study='A'):

    p1: list[float] = [1.0, 0.0, 0.0, 1.0]
    p2: list[float] = [1.0, 0.0, 0.0, -1.0]
    p3: list[float] = [
        0.9903067841547669,
        0.3734458315384239,
        0.6393671561700484,
        -0.657613394982612,
    ]
    p4: list[float] = [
        0.2053558277500182,
        -0.1412408546154445,
        -0.0869119059139323,
        0.121112995127698,
    ]
    p5: list[float] = [
        p1[i]+p2[i]-p3[i]-p4[i] for i in range(4)
    ]
    if elements_to_display is None:
        elements_to_display = ['I', 'th.CT', 'tot']
    if tests_to_run is None:
        tests_to_run = ['IR', 'UV']

    scaling_factor = 0.1

    match threshold_to_study:
        case 'A':
            on_threshold_point = [0.0, 1.0, 0.0]
            # Add the -p1 shift since the LMB locates the defining is between p1 and p2
            on_threshold_point = [
                on_threshold_point[i] + p1[i+1] for i in range(3)]
            # on_threshold_point = [on_threshold_point[i] * 2. for i in range(3)]
        case _:
            raise ValueError(
                f"Unknown threshold to study '{threshold_to_study}'")

    x = 0.1
    for test_name, k_base, scaling, pref_power in [
        ('IR', [0.0001, 0.0001, 1.0*x], 0.1, 0.0),
        ('UV', [1.0e3, 1.0e3, 1.0e3], 10.0, 1./2.),
        ('THRES', on_threshold_point, 0.1, 0.0),
    ]:
        if test_name not in tests_to_run:
            continue
        running_sum_non_ct = complex(0.0, 0.0)
        running_sum_ct = complex(0.0, 0.0)
        running_sum_non_ct_scaled = complex(0.0, 0.0)
        running_sum_ct_scaled = complex(0.0, 0.0)
        print('')
        print(f'Investigating {test_name}...')
        scaling_factor = scaling
        if test_name == 'IR':
            k = k_base
            scaled_k = [k[0]*scaling_factor,
                        k[1]*scaling_factor, k[2]]
        elif test_name == 'UV':
            k = k_base
            scaled_k = [k[0]+scaling_factor,
                        k[1]*scaling_factor, k[2]*scaling_factor]
        elif test_name == 'THRES':
            eps = 0.001
            shift = [eps, eps, eps]
            k = [k_base[0] + shift[0],
                 k_base[1] + shift[1], k_base[2] + shift[2]]
            scaled_k = [k_base[0] + shift[0]*scaling_factor,
                        k_base[1] + shift[1]*scaling_factor, k_base[2] + shift[2]*scaling_factor]
            # print(k)
            # print(scaled_k)
        else:
            raise ValueError(f"Unknown limit type to approach '{test_name}'")
        all_res, jac, total, reconstructed_total = eval_integrand(
            k, debug, gammaloop_state, gammaloop_executable)
        scaled_all_res, scaled_jac, scaled_total, scaled_reconstructed_total = eval_integrand(
            scaled_k, debug, gammaloop_state, gammaloop_executable)
        for i_g, g_name in enumerate(graph_names):
            eval_res = all_res[g_name]
            eval_res_scaled = scaled_all_res[g_name]

            eval_res_scaled = {k: v*(scaling_factor**pref_power)
                               for k, v in eval_res_scaled.items()}

            if test_name == 'THRES':
                running_sum_ct += eval_res['I']
                running_sum_ct_scaled += eval_res_scaled['I']
                running_sum_non_ct += eval_res['threshold_CT']
                running_sum_non_ct_scaled += eval_res_scaled['threshold_CT']
            else:
                if g_name in ['ir_ct', 'uv_triangle_A', 'uv_triangle_B', 'uv_triangle_C', 'uv_bubble_A', 'uv_bubble_B', 'uv_ir_ct']:
                    running_sum_ct += eval_res['I'] + eval_res['threshold_CT']
                    running_sum_ct_scaled += eval_res_scaled['I'] + \
                        eval_res_scaled['threshold_CT']
                else:
                    running_sum_non_ct += eval_res['I'] + \
                        eval_res['threshold_CT']
                    running_sum_non_ct_scaled += eval_res_scaled['I'] + \
                        eval_res_scaled['threshold_CT']
            res = eval_res['I']
            res_scaled = eval_res_scaled['I']
            if 'I' in elements_to_display:
                print(f"{'Evaluation result for graph "'+g_name+'" [I]':60s}: {res:50.16e} | scaled: {
                    res_scaled:50.16e} | log10(abs(scaled / orig)) = {log10(abs(res_scaled / res)) if res != 0. else float('nan'):.3f} ")
            res = eval_res['threshold_CT']
            res_scaled = eval_res_scaled['threshold_CT']
            if 'th.CT' in elements_to_display:
                print(f"{'Evaluation result for graph "'+g_name+'" [th.CT]':60s}: {res:50.16e} | scaled: {
                    res_scaled:50.16e} | log10(abs(scaled / orig)) = {log10(abs(res_scaled / res)) if res != 0. else float('nan'):.3f} ")
            res = eval_res['I'] + eval_res['threshold_CT']
            res_scaled = eval_res_scaled['I'] + eval_res_scaled['threshold_CT']
            if 'tot' in elements_to_display:
                print(f"{'Evaluation result for graph "'+g_name+'" [tot]':60s}: {res:50.16e} | scaled: {
                    res_scaled:50.16e} | log10(abs(scaled / orig)) = {log10(abs(res_scaled / res)) if res != 0. else float('nan'):.3f} ")

        print('-'*100)
        print(f"{'Total evaluation result from non-CT':60s}: {running_sum_non_ct:50.16e} | scaled: {running_sum_non_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {
              log10(abs(running_sum_non_ct_scaled / running_sum_non_ct)) if running_sum_non_ct != 0. else float('nan'):.3f} ")
        print(f"{'Total evaluation result from CT':60s}: {running_sum_ct:50.16e} | scaled: {running_sum_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {
              log10(abs(running_sum_ct_scaled / running_sum_ct)) if running_sum_non_ct != 0. else float('nan'):.3f} ")
        print(f"{'Total evaluation result':60s}: {running_sum_ct+running_sum_non_ct:50.16e} | scaled: {running_sum_ct_scaled+running_sum_non_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {
              log10(abs((running_sum_ct_scaled+running_sum_non_ct_scaled) / (running_sum_ct+running_sum_non_ct))) if running_sum_ct+running_sum_non_ct != 0. else float('nan'):.3f} ")
        # print(f"{'Ratio non-CT / CT ':60s}: {running_sum_non_ct / running_sum_ct:50.16e} | scaled: {running_sum_non_ct_scaled / running_sum_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {log10(abs((running_sum_non_ct_scaled / running_sum_ct_scaled) / (running_sum_non_ct / running_sum_ct))) if running_sum_ct_scaled!=0. and running_sum_ct!=0. and running_sum_non_ct!=0. else float('nan'):.3f} ")

        print('')


parser = argparse.ArgumentParser("Run limit tests for qqx_aaa process")
parser.add_argument('--display_elements', '-d', nargs='*', choices=['I', 'th.CT', 'tot'],
                    help="Which elements to display per graph", default=['tot'])
parser.add_argument('--tests', '-t', nargs='+', choices=['IR', 'UV', 'THRES'],
                    help="Which limits to test", default=['IR', 'UV', 'THRES'])
parser.add_argument('--threshold_to_study', '-thres', type=str, choices=['A'],
                    help="Which threshold to study", default='A')
parser.add_argument('--gammaloop_state', '-s', type=str,
                    help="Path to the GammaLoop state directory. Default=%(default)s", default='./GL_QQX_AAA_euclidean')
parser.add_argument('--gammaloop_executable', '-e', type=str,
                    help="Path to the GammaLoop executable. Default=%(default)s", default='../../target/dev-optim/gammaloop')
parser.add_argument('--debug', action='store_true',
                    help="Enable debug output", default=False)
if __name__ == "__main__":
    args = parser.parse_args()
    run_limit_test(args.display_elements, args.tests,
                   args.gammaloop_state, args.gammaloop_executable, args.debug, args.threshold_to_study)
