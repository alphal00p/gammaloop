#!/usr/bin/env python3

import subprocess

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


def get_integrand_results(stdout):
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


def eval_integrand(k_input, debug=False):
    cmd = f'../../target/release/cli -n inspect --process-id 0 --name qqx_aaa_subtracted -m -p {
        k_input[0]:.16f} {k_input[1]:.16f} {k_input[2]:.16f}'
    if debug:
        print(f"Running command:\n{cmd}")
    r = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, shell=True)
    if debug:
        print(f"Output:\n{r.stdout}")
    all_res, total, jac = get_integrand_results(r.stdout)
    reconstructed_total = sum(all_res)*jac
    if debug:
        for g_name, res in zip(graph_names, all_res):
            print(f"{g_name:30s}: {res.real:+.16e} {res.imag:+.16e}")
            print(f"{'Jacobian':30s}: {jac.real:+.16e} {jac.imag:+.16e}")

        print(f"{'Reconstructed total':30s}: {reconstructed_total.real:+.16e} {reconstructed_total.imag:+.16e}")  # nopep8
        print(f"{'Total':30s}: {total.real:+.16e} {total.imag:+.16e}")
    return all_res, jac, total, reconstructed_total


def run_limit_test():
    scaling_factor = 0.1

    x = 0.1
    for test_name, k_base, scaling, pref_power in [
        ('IR limit', [0.0001, 0.0001, 1.0*x], 0.1, 0.0),
        ('UV limit', [100.0, 100.0, 100.0], 10.0, 0.0),
    ]:
        running_sum_non_ct = complex(0.0, 0.0)
        running_sum_ct = complex(0.0, 0.0)
        running_sum_non_ct_scaled = complex(0.0, 0.0)
        running_sum_ct_scaled = complex(0.0, 0.0)
        print('')
        print(f'Investigating {test_name}...')
        k = k_base
        scaling_factor = scaling
        all_res, jac, total, reconstructed_total = eval_integrand(k)
        if test_name == 'IR limit':
            scaled_k = [k[0]*scaling_factor,
                        k[1]*scaling_factor, k[2]]
        elif test_name == 'UV limit':
            scaled_k = [k[0]*scaling_factor,
                        k[1]*scaling_factor, k[2]*scaling_factor]
        scaled_all_res, scaled_jac, scaled_total, scaled_reconstructed_total = eval_integrand(
            scaled_k)
        for i_g, g_name in enumerate(graph_names):
            eval_res = all_res[i_g]
            eval_res_scaled = scaled_all_res[i_g]

            eval_res_scaled *= scaling_factor**pref_power
            if g_name in ['ir_ct', 'uv_triangle_A', 'uv_triangle_B', 'uv_triangle_C', 'uv_bubble_A', 'uv_bubble_B', 'uv_ir_ct']:
                running_sum_ct += eval_res
                running_sum_ct_scaled += eval_res_scaled
            else:
                running_sum_non_ct += eval_res
                running_sum_non_ct_scaled += eval_res_scaled
            print(f"{'Evaluation result for graph "'+g_name+'"':50s}: {eval_res:50.16e} | scaled: {
                  eval_res_scaled:50.16e} | log10(abs(scaled / orig)) = {log10(abs(eval_res_scaled / eval_res)) if eval_res != 0. else float('nan'):.3f} ")

        print('-'*100)
        print(f"{'Total evaluation result from non-CT':50s}: {running_sum_non_ct:50.16e} | scaled: {running_sum_non_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {
              log10(abs(running_sum_non_ct_scaled / running_sum_non_ct)) if running_sum_non_ct != 0. else float('nan'):.3f} ")
        print(f"{'Total evaluation result from CT':50s}: {running_sum_ct:50.16e} | scaled: {running_sum_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {
              log10(abs(running_sum_ct_scaled / running_sum_ct)) if running_sum_non_ct != 0. else float('nan'):.3f} ")
        print(f"{'Total evaluation result':50s}: {running_sum_ct+running_sum_non_ct:50.16e} | scaled: {running_sum_ct_scaled+running_sum_non_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {
              log10(abs((running_sum_ct_scaled+running_sum_non_ct_scaled) / (running_sum_ct+running_sum_non_ct))) if running_sum_ct+running_sum_non_ct != 0. else float('nan'):.3f} ")
        # print(f"{'Ratio non-CT / CT ':50s}: {running_sum_non_ct / running_sum_ct:50.16e} | scaled: {running_sum_non_ct_scaled / running_sum_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {log10(abs((running_sum_non_ct_scaled / running_sum_ct_scaled) / (running_sum_non_ct / running_sum_ct))) if running_sum_ct_scaled!=0. and running_sum_ct!=0. and running_sum_non_ct!=0. else float('nan'):.3f} ")

        print('')


run_limit_test()
