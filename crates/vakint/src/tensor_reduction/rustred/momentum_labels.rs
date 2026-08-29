use std::collections::BTreeSet;

use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use symbolica::function;

use super::contains_function_head;
use crate::get_integer_from_atom;
use crate::symbols::S;
use crate::tensor_reduction::TensorReductionError;
use crate::utils::vakint_macros::vk_symbol;

/// Adapter-private momentum identities passed through RustRed's tensor API.
///
/// Vakint numbers loop and external momenta independently, so `k(1)` and
/// `p(1)` are different vectors. RustRed deliberately keeps both kinds in one
/// label namespace. Distinct inert wrapper heads preserve Vakint's type as
/// part of the label without nesting any of RustRed's reserved tensor heads.
#[derive(Clone, Copy)]
pub(super) struct MomentumLabels {
    loop_head: Symbol,
    external_head: Symbol,
}

impl MomentumLabels {
    pub(super) fn new() -> Self {
        Self {
            loop_head: vk_symbol!("rustred_loop_label"),
            external_head: vk_symbol!("rustred_external_label"),
        }
    }

    pub(super) fn loop_labels(
        &self,
        term_index: usize,
        loop_count: usize,
    ) -> Result<Vec<Atom>, TensorReductionError> {
        (1..=loop_count)
            .map(|id| {
                i64::try_from(id)
                    .map(|id| self.loop_label(Atom::num(id)))
                    .map_err(|_| TensorReductionError::RustRedUnsupportedOutput {
                        term: term_index,
                        detail: format!("canonical loop label {id} does not fit in i64"),
                    })
            })
            .collect()
    }

    /// Reject caller-spelled occurrences of either adapter-private head.
    ///
    /// Symbolica symbols are globally interned, so privacy must be enforced at
    /// the value boundary in addition to Rust module visibility. This scans
    /// recursively and therefore also catches a tag hidden inside a `dot`
    /// argument or nested as the first argument of a vector.
    pub(super) fn reject_reserved_input(
        &self,
        term_index: usize,
        numerator: AtomView<'_>,
    ) -> Result<(), TensorReductionError> {
        let head = if contains_function_head(numerator, self.loop_head) {
            Some("rustred_loop_label")
        } else if contains_function_head(numerator, self.external_head) {
            Some("rustred_external_label")
        } else {
            None
        };
        match head {
            Some(head) => Err(TensorReductionError::RustRedReservedMomentumLabel {
                term: term_index,
                head: head.to_owned(),
            }),
            None => Ok(()),
        }
    }

    /// Convert Vakint vectors to RustRed's globally disjoint label namespace.
    ///
    /// `replace_map` is top-down and does not visit a replacement, so this is
    /// one simultaneous pass: a newly inserted wrapper cannot be wrapped
    /// again. The caller first converts indexed contractions to `dot`, leaving
    /// every remaining vector in one- or two-argument Vakint form.
    pub(super) fn encode(&self, numerator: AtomView<'_>) -> (Atom, Vec<Atom>) {
        let mut external_ids = BTreeSet::new();
        let encoded = numerator.replace_map(|atom, _context, output| {
            let AtomView::Fun(vector) = atom else {
                return;
            };
            let head = vector.get_symbol();
            if head == S.dot && vector.get_nargs() == 2 {
                let mut arguments = vector.iter();
                let Some(left) = arguments
                    .next()
                    .and_then(|argument| self.encode_dot_argument(argument, &mut external_ids))
                else {
                    return;
                };
                let Some(right) = arguments
                    .next()
                    .and_then(|argument| self.encode_dot_argument(argument, &mut external_ids))
                else {
                    return;
                };
                **output = function!(S.dot, left, right);
                return;
            }
            if (head != S.k && head != S.p) || !(vector.get_nargs() == 1 || vector.get_nargs() == 2)
            {
                return;
            }
            let mut arguments = vector.iter();
            let id = arguments.next().and_then(get_integer_from_atom);
            let Some(id) = id else {
                return;
            };
            let label = if head == S.k {
                self.loop_label(Atom::num(id))
            } else {
                external_ids.insert(id);
                self.external_label(Atom::num(id))
            };
            let mut replacement = FunctionBuilder::new(head).add_arg(label);
            if let Some(index) = arguments.next() {
                replacement = replacement.add_arg(index);
            }
            **output = replacement.finish();
        });
        let external_labels = external_ids
            .into_iter()
            .map(|id| self.external_label(Atom::num(id)))
            .collect();
        (encoded, external_labels)
    }

    /// Restore Vakint's typed vector syntax and reject any leaked private tag.
    pub(super) fn decode(
        &self,
        term_index: usize,
        expression: AtomView<'_>,
    ) -> Result<Atom, TensorReductionError> {
        let decoded = expression.replace_map(|atom, _context, output| {
            let AtomView::Fun(vector) = atom else {
                return;
            };
            let head = vector.get_symbol();
            if head == S.dot && vector.get_nargs() == 2 {
                let mut arguments = vector.iter();
                let Some(left) = arguments
                    .next()
                    .and_then(|argument| self.decode_dot_argument(argument))
                else {
                    return;
                };
                let Some(right) = arguments
                    .next()
                    .and_then(|argument| self.decode_dot_argument(argument))
                else {
                    return;
                };
                **output = function!(S.dot, left, right);
                return;
            }
            if (head != S.k && head != S.p) || !(vector.get_nargs() == 1 || vector.get_nargs() == 2)
            {
                return;
            }
            let mut arguments = vector.iter();
            let Some(label) = arguments.next() else {
                return;
            };
            let expected_label_head = if head == S.k {
                self.loop_head
            } else {
                self.external_head
            };
            let Some(id) = tagged_id(label, expected_label_head) else {
                return;
            };
            let mut replacement = FunctionBuilder::new(head).add_arg(Atom::num(id));
            if let Some(index) = arguments.next() {
                replacement = replacement.add_arg(index);
            }
            **output = replacement.finish();
        });
        if contains_function_head(decoded.as_view(), self.loop_head)
            || contains_function_head(decoded.as_view(), self.external_head)
        {
            return Err(TensorReductionError::RustRedUnsupportedOutput {
                term: term_index,
                detail: "adapter-private momentum label escaped RustRed output".to_owned(),
            });
        }
        Ok(decoded)
    }

    fn loop_label(&self, id: Atom) -> Atom {
        function!(self.loop_head, id)
    }

    fn external_label(&self, id: Atom) -> Atom {
        function!(self.external_head, id)
    }

    fn encode_dot_argument(
        &self,
        argument: AtomView<'_>,
        external_ids: &mut BTreeSet<i64>,
    ) -> Option<Atom> {
        let AtomView::Fun(vector) = argument else {
            return None;
        };
        if vector.get_nargs() != 1 || (vector.get_symbol() != S.k && vector.get_symbol() != S.p) {
            return None;
        }
        let id = vector.iter().next().and_then(get_integer_from_atom)?;
        if vector.get_symbol() == S.k {
            Some(self.loop_label(Atom::num(id)))
        } else {
            external_ids.insert(id);
            Some(self.external_label(Atom::num(id)))
        }
    }

    fn decode_dot_argument(&self, argument: AtomView<'_>) -> Option<Atom> {
        if let Some(id) = tagged_id(argument, self.loop_head) {
            Some(function!(S.k, Atom::num(id)))
        } else {
            tagged_id(argument, self.external_head).map(|id| function!(S.p, Atom::num(id)))
        }
    }
}

fn tagged_id(label: AtomView<'_>, expected_head: Symbol) -> Option<i64> {
    let AtomView::Fun(label) = label else {
        return None;
    };
    if label.get_symbol() != expected_head || label.get_nargs() != 1 {
        return None;
    }
    label.iter().next().and_then(get_integer_from_atom)
}
