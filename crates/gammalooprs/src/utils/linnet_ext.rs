use linnet::permutation::Permutation;

pub trait FromMappings: Sized {
    fn from_mappings<I>(mappings: I, size: usize) -> Result<Self, String>
    where
        I: IntoIterator<Item = (usize, usize)>;
}

impl FromMappings for Permutation {
    /// Creates a permutation from a set of mapping pairs, minimizing changes to other elements.
    ///
    /// Takes a collection of (from, to) pairs and constructs a valid permutation that satisfies
    /// these mappings while keeping other elements as close to their identity positions as possible.
    ///
    /// # Arguments
    ///
    /// * `mappings` - An iterator of (from, to) pairs where each pair represents `from -> to`
    /// * `size` - The size of the permutation to create
    ///
    /// # Returns
    ///
    /// A `Result` containing the permutation if successful, or an error string if the mappings
    /// are inconsistent or invalid.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// // Create a permutation where 5->0 and 0->1, others stay in place where possible
    /// let mappings = vec![(5, 0), (0, 1)];
    /// let p = Permutation::from_mappings(mappings.iter().copied(), 6).unwrap();
    /// assert_eq!(p.map(), &[1, 5, 2, 3, 4, 0]);
    ///
    /// // Verify the specific mappings work
    /// assert_eq!(p[5], 0);  // 5 -> 0
    /// assert_eq!(p[0], 1);  // 0 -> 1
    /// ```
    fn from_mappings<I>(mappings: I, size: usize) -> Result<Self, String>
    where
        I: IntoIterator<Item = (usize, usize)>,
    {
        if size == 0 {
            return Ok(Self::id(0));
        }

        // Start with identity
        let mut map: Vec<usize> = (0..size).collect();

        // owner_of_target[t] = Some(s) iff target t is already taken by source s (forced or assigned)
        let mut owner_of_target: Vec<Option<usize>> = vec![None; size];

        // forced_sources[s] = true iff s has a forced mapping in the input
        let mut forced_sources = vec![false; size];

        // 1) Apply and validate forced mappings in O(1) each
        for (from, to) in mappings {
            if from >= size || to >= size {
                return Err(format!(
                    "Index out of bounds: mapping ({}, {}) for size {}",
                    from, to, size
                ));
            }

            if let Some(prev_from) = owner_of_target[to]
                && prev_from != from
            {
                return Err(format!(
                    "Target {} is mapped from multiple sources ({} and {})",
                    to, prev_from, from
                ));
            }
            if forced_sources[from] && map[from] != to {
                return Err(format!(
                    "Source {} has conflicting mappings ({} vs {})",
                    from, map[from], to
                ));
            }

            // Record the forced mapping
            map[from] = to;
            owner_of_target[to] = Some(from);
            forced_sources[from] = true;
        }

        // 2) Keep identities whenever possible (minimal changes)
        //    If a source s is not forced and its identity target s is free, keep s->s.
        //    Otherwise, it becomes "displaced" and will be assigned later.
        let mut displaced_sources: Vec<usize> = Vec::new();
        for s in 0..size {
            if forced_sources[s] {
                continue;
            }
            if owner_of_target[s].is_none() {
                // take the identity slot
                map[s] = s;
                owner_of_target[s] = Some(s);
            } else {
                // identity target already used by some forced mapping -> displace s
                displaced_sources.push(s);
            }
        }

        // 3) Gather remaining free targets (those without an owner)
        let mut free_targets: Vec<usize> = Vec::new();
        for (t, owner) in owner_of_target.iter().enumerate().take(size) {
            if owner.is_none() {
                free_targets.push(t);
            }
        }

        // Sanity: counts must match, otherwise constraints were inconsistent
        if displaced_sources.len() != free_targets.len() {
            return Err("Failed to create valid permutation (inconsistent counts)".to_string());
        }

        // 4) Assign displaced sources to remaining free targets
        for (s, t) in displaced_sources.into_iter().zip(free_targets.into_iter()) {
            // Note: t != s by construction (if t==s, owner_of_target[s] would have been None earlier,
            // and we would have kept the identity)
            map[s] = t;
            owner_of_target[t] = Some(s);
        }

        // 5) Final uniqueness check (linear)
        let mut seen = vec![false; size];
        for &t in &map {
            if t >= size || std::mem::replace(&mut seen[t], true) {
                return Err("Failed to create valid permutation".to_string());
            }
        }

        Ok(Self::from_map(map))
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_snapshot;
    use linnet::permutation::Permutation;

    #[test]
    fn test_from_mappings_complex() {
        // Test more complex displacement
        let mappings = [(0, 3), (1, 0)];
        let p = Permutation::from_mappings(mappings.iter().copied(), 4).unwrap();
        assert_eq!(p[0], 3);
        assert_eq!(p[1], 0);

        assert_snapshot!( format!("{:?}",p.apply_slice([0,1,2,3])),@"[1, 3, 2, 0]");
    }
}
