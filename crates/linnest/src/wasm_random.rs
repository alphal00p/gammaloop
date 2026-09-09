use std::sync::{Mutex, OnceLock};

const DEFAULT_SEED: u32 = 0x5eed_1234;

static RNG_STATE: OnceLock<Mutex<[u32; 3]>> = OnceLock::new();

#[inline(always)]
fn tausworthe(s: u32, a: u8, b: u8, c: u32, d: u8) -> u32 {
    let s1 = (s & c) << d;
    let s2 = ((s << a) ^ s) >> b;
    s1 ^ s2
}

#[inline(always)]
fn taus_get(state: [u32; 3]) -> ([u32; 3], u32) {
    let s1 = tausworthe(state[0], 13, 19, 0xFFFF_FFFE, 12);
    let s2 = tausworthe(state[1], 2, 25, 0xFFFF_FFF8, 4);
    let s3 = tausworthe(state[2], 3, 11, 0xFFFF_FFF0, 17);
    let value = s1 ^ s2 ^ s3;
    ([s1, s2, s3], value)
}

fn taus_set(seed: u32) -> [u32; 3] {
    let mut s1 = 69069 * seed;
    if s1 < 2 {
        s1 += 2;
    }

    let mut s2 = 69069 * s1;
    if s2 < 8 {
        s2 += 8;
    }

    let mut s3 = 69069 * s2;
    if s3 < 16 {
        s3 += 16;
    }

    let mut state = [s1, s2, s3];
    for _ in 0..6 {
        (state, _) = taus_get(state);
    }
    state
}

#[inline(always)]
fn next_u32(state: &mut [u32; 3]) -> u32 {
    let (next_state, value) = taus_get(*state);
    *state = next_state;
    value
}

pub(crate) fn custom_getrandom(buf: &mut [u8]) -> Result<(), getrandom::Error> {
    let state = RNG_STATE.get_or_init(|| Mutex::new(taus_set(DEFAULT_SEED)));
    let mut state = state
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());

    for chunk in buf.chunks_mut(std::mem::size_of::<u32>()) {
        let random_bytes = next_u32(&mut state).to_le_bytes();
        let len = chunk.len();
        chunk.copy_from_slice(&random_bytes[..len]);
    }

    Ok(())
}
