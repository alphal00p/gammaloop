# BNL Gamma Simplification And Chain Cleanup Captures

Source logs:

```text
/tmp/bnl-integrated-state-chainlog-fields/logs/gammalog-2026-05-27T12-53-40.076-02-00.jsonl
/tmp/bnl-integrated-state-postmetricgamma/logs/gammalog-2026-05-27T13-07-10.451-02-00.jsonl
```

This file only shows `log_print(Some(120))` fields. Exact `to_plain_string()`
payloads are intentionally left out because they are too noisy for inspection.

## `integrated_uv_start_after_simplify_gamma.before_gamma`

```text
1𝑖·GC_2·GC_11²·q₁₃(αᵉ¹³.¹)·q₁₄(αᵉ¹⁴.¹)
    ·γ(α₆,b₁₃ b₂₂)
    ·γ(α₂₄,b₂₁ b₁₂)
    ·γ(α₂₄,b₂₅ b₁₄)
    ·γ(αᵉ¹⁴.¹,b₁₄ b₁₃)
    ·γ(αᵉ¹³.¹,b₂₂ b₂₁)
    ·t(c₂₄ г₁₂)(г₁₄)
    ·t(c₂₄ г₁₄)(г₂₅)
+1𝑖·MT·GC_2·GC_11²·q₁₃(αᵉ¹³.¹)·γ(α₆,b₁₄ b₂₂)
        ·γ(α₂₄,b₂₁ b₁₂)
        ·γ(α₂₄,b₂₅ b₁₄)
        ·γ(αᵉ¹³.¹,b₂₂ b₂₁)
        ·t(c₂₄ г₁₂)(г₁₄)
        ·t(c₂₄ г₁₄)(г₂₅)
+1𝑖·MT·GC_2·GC_11²·q₁₄(αᵉ¹⁴.¹)·γ(α₆,b₁₃ b₂₂)
        ·γ(α₂₄,b₂₂ b₁₂)
        ·γ(α₂₄,b₂₅ b₁₄)
        ·γ(αᵉ¹⁴.¹,b₁₄ b₁₃)
        ·t(c₂₄ г₁₂)(г₁₄)
        ·t(c₂₄ г₁₄)(г₂₅)
+1𝑖·MT²·GC_2·GC_11²·γ(α₆,b₁₄ b₂₂)
        ·γ(α₂₄,b₂₂ b₁₂)
        ·γ(α₂₄,b₂₅ b₁₄)
        ·t(c₂₄ г₁₂)(г₁₄)
        ·t(c₂₄ г₁₄)(г₂₅)
```

## `integrated_uv_start_after_simplify_gamma.after_gamma`

```text
2𝑖·GC_2·GC_11²
    ·b₂₅[γ(q₁₃(mink(dim)))γ(q₁₄(mink(dim)))γ(α₆)]b₁₂·t(c₂₄ г₁₂)(г₁₄)
    ·t(c₂₄ г₁₄)(г₂₅)
-2𝑖·GC_2·GC_11²
        ·b₂₅[γ(α₆)γ(q₁₄(mink(dim)))γ(q₁₃(mink(dim)))]b₁₂·t(c₂₄ г₁₂)(г₁₄)
        ·t(c₂₄ г₁₄)(г₂₅)
+2𝑖·MT·GC_2·GC_11²·b₂₅[γ(q₁₃(mink(dim)))γ(α₆)]b₁₂
        ·t(c₂₄ г₁₂)(г₁₄)
        ·t(c₂₄ г₁₄)(г₂₅)
+2𝑖·MT·GC_2·GC_11²·b₂₅[γ(α₆)γ(q₁₄(mink(dim)))]b₁₂
        ·t(c₂₄ г₁₂)(г₁₄)
        ·t(c₂₄ г₁₄)(г₂₅)
+(-1𝑖·dim·GC_2·GC_11²+2𝑖·GC_2·GC_11²)
        ·b₂₅[γ(q₁₄(mink(dim)))γ(α₆)γ(q₁₃(mink(dim)))]b₁₂·t(c₂₄ г₁₂)(г₁₄)
        ·t(c₂₄ г₁₄)(г₂₅)
+(1𝑖·dim·MT·GC_2·GC_11²-2𝑖·MT·GC_2·GC_11²)
        ·b₂₅[γ(q₁₄(mink(dim)))γ(α₆)]b₁₂
        ·t(c₂₄ г₁₂)(г₁₄)
        ·t(c₂₄ г₁₄)(г₂₅)
+(1𝑖·dim·MT·GC_2·GC_11²-2𝑖·MT·GC_2·GC_11²)
        ·b₂₅[γ(α₆)γ(q₁₃(mink(dim)))]b₁₂
        ·t(c₂₄ г₁₂)(г₁₄)
        ·t(c₂₄ г₁₄)(г₂₅)
+(-1𝑖·dim·MT²·GC_2·GC_11²+2𝑖·MT²·GC_2·GC_11²)
        ·b₂₅[γ(α₆)]b₁₂
        ·t(c₂₄ г₁₂)(г₁₄)
        ·t(c₂₄ г₁₄)(г₂₅)
```

## `integrate_and_truncate_after_collect_chains.before_log_print`

```text
2·ε·(-1/32·GC_2·GC_11²/𝜋²-1/16·GC_2·GC_11²·log(mUVI)/𝜋²+1/32·GC_2·GC_11²·log(μᵣ²)/𝜋²)
    ·b₂₅[γ(α₆)]b₁₂
    ·t(c₂₄ г₁₂)(г₁₄)
    ·t(c₂₄ г₁₄)(г₂₅)
+(-2·(-1/32·GC_2·GC_11²/𝜋²-1/16·GC_2·GC_11²·log(mUVI)/𝜋²+1/32·GC_2·GC_11²·log(μᵣ²)/𝜋²)
        +1/16·GC_2·GC_11²/𝜋²
    )·b₂₅[γ(α₆)]b₁₂
        ·t(c₂₄ г₁₂)(г₁₄)
        ·t(c₂₄ г₁₄)(г₂₅)
```

## `integrate_and_truncate_after_collect_chains.after_collect_log_print`

```text
(2·ε·(-1/32·GC_2·GC_11²/𝜋²-1/16·GC_2·GC_11²·log(mUVI)/𝜋²+1/32·GC_2·GC_11²·log(μᵣ²)/𝜋²)
    ·t(c₂₄ г₁₂)(г₁₄)
    ·t(c₂₄ г₁₄)(г₂₅)
+(-2·(-1/32·GC_2·GC_11²/𝜋²-1/16·GC_2·GC_11²·log(mUVI)/𝜋²+1/32·GC_2·GC_11²·log(μᵣ²)/𝜋²)
    +1/16·GC_2·GC_11²/𝜋²
    )·t(c₂₄ г₁₂)(г₁₄)
        ·t(c₂₄ г₁₄)(г₂₅)
)·b₂₅[γ(α₆)]b₁₂
```

## `integrate_and_truncate_after_collect_chains.after_undo_single_length_log_print`

`after_undo_single_length_log_print` is identical to `after_collect_log_print`
for this capture.
