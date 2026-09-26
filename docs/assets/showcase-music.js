/* Procedural slow-blues piano score for the GammaLoop showcase.
   A seeded composer writes a twelve-bar jazz blues in F on a triplet (shuffle) grid: a
   two-feel upright bass that walks in the second chorus, rootless left-hand voicings, and a
   sparse right-hand melody with blue notes, grace-note slides, and one tremolo fill per
   chorus. A small additive piano/bass synthesizer with a Schroeder reverb renders it. The
   module runs unchanged in Node (WAV for the video) and in the browser (Web Audio chunks),
   so both hear the same notes. */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.GammaLoopShowcaseMusic = factory();
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';

  const BPM = 66;
  const BEAT = 60 / BPM;
  const BAR = 4 * BEAT;
  const SLOT = BEAT / 3; // twelve triplet slots per bar
  const INTRO_BARS = 2;

  function rng(seed) {
    let s = seed >>> 0;
    return () => {
      s = (s + 0x6d2b79f5) | 0;
      let t = Math.imul(s ^ (s >>> 15), 1 | s);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  const pick = (next, list) => list[Math.floor(next() * list.length)];
  const midiHz = (m) => 440 * Math.pow(2, (m - 69) / 12);
  const pc = (midi) => ((midi % 12) + 12) % 12;

  // Chord vocabulary in F: bass root, rootless voicing, chord tones, and extra colour tones.
  const CHORDS = {
    F7: { root: 41, voicing: [57, 63, 67, 74], tones: [5, 9, 0, 3], colour: [7, 2] },
    Bb7: { root: 46, voicing: [62, 68, 72], tones: [10, 2, 5, 8], colour: [0, 7] },
    Bdim7: { root: 47, voicing: [62, 65, 68], tones: [11, 2, 5, 8], colour: [] },
    Cm7: { root: 48, voicing: [63, 67, 70], tones: [0, 3, 7, 10], colour: [2, 5] },
    D7b9: { root: 50, voicing: [66, 72, 75], tones: [2, 6, 9, 0], colour: [3] },
    Gm7: { root: 43, voicing: [58, 62, 65, 69], tones: [7, 10, 2, 5], colour: [9, 0] },
    C7: { root: 48, voicing: [64, 70, 74], tones: [0, 4, 7, 10], colour: [2, 9, 3] },
  };
  const BLUES = [5, 8, 10, 11, 0, 3]; // F Ab Bb B C Eb
  const CHORUS = [
    ['F7'], ['Bb7'], ['F7'], ['Cm7', 'F7'],
    ['Bb7'], ['Bdim7'], ['F7'], ['D7b9'],
    ['Gm7'], ['C7'], ['F7', 'D7b9'], ['Gm7', 'C7'],
  ];
  // Comping patterns as [slot, length in slots] on the triplet grid.
  const COMP = [
    [[0, 5], [6, 4]],
    [[0, 4], [5, 3], [9, 2]],
    [[2, 4], [6, 5]],
    [[0, 6], [8, 3]],
    [[0, 3], [3, 2], [6, 5]],
    [[0, 5], [8, 2], [11, 1]],
  ];
  const COMP_SHORT = [[[0, 4]], [[0, 3], [3, 2]], [[2, 4]]];
  // Melody rhythm cells in triplet slots; the last value is usually a held note.
  const CELLS = [
    [2, 1, 3, 6],
    [1, 1, 1, 3, 6],
    [3, 3, 6],
    [2, 1, 2, 1, 6],
    [4, 2, 6],
    [1, 2, 3, 6],
    [2, 1, 2, 1, 2, 4],
    [3, 1, 2, 6],
  ];

  const slotTime = (bar, slot) => bar * BAR + slot * SLOT;

  // ---------------------------------------------------------------- composer
  function compose(duration, seed) {
    const next = rng(seed);
    const notes = [];
    const push = (inst, midi, start, dur, vel) => {
      if (start >= duration - 0.05 || start < 0) return;
      notes.push({ inst, midi, start, dur: Math.max(0.05, Math.min(dur, duration - start)), vel });
    };
    const bars = Math.floor((duration - 4.5) / BAR);

    // Chord timeline: an intro on the tonic, then choruses of the blues form.
    const segments = [];
    for (let b = 0; b < bars; b++) {
      const bar = b < INTRO_BARS ? ['F7'] : CHORUS[(b - INTRO_BARS) % CHORUS.length];
      bar.forEach((name, i) => {
        segments.push({
          name,
          chord: CHORDS[name],
          bar: b,
          slot: (i * 12) / bar.length,
          slots: 12 / bar.length,
          chorus: b < INTRO_BARS ? -1 : Math.floor((b - INTRO_BARS) / CHORUS.length),
        });
      });
    }
    const segmentAt = (bar, slot) =>
      segments.find((s) => s.bar === bar && slot >= s.slot && slot < s.slot + s.slots) || segments[segments.length - 1];
    const bassRange = (m) => {
      while (m > 55) m -= 12;
      while (m < 36) m += 12;
      return m;
    };

    // Bass: a two-feel in the first chorus, walking quarters afterwards.
    segments.forEach((seg, i) => {
      const nextRoot = segments[i + 1] ? segments[i + 1].chord.root : 41;
      const r = seg.chord.root;
      const third = r + (seg.name.includes('m7') || seg.name.includes('dim') ? 3 : 4);
      const fifth = r + (seg.name.includes('dim') ? 6 : 7);
      const seventh = r + 10;
      const start = slotTime(seg.bar, seg.slot);
      if (seg.bar === 0) return; // solo piano intro bar
      if (seg.chorus < 1) {
        push('bass', bassRange(r), start, BEAT * 1.9, 0.5 + next() * 0.06);
        if (seg.slots === 12) {
          const second = next() < 0.6 ? fifth : pick(next, [third, seventh, nextRoot - 1, nextRoot + 1]);
          push('bass', bassRange(second), start + 2 * BEAT, BEAT * 1.85, 0.44 + next() * 0.06);
        }
      } else {
        const line = [r];
        if (seg.slots === 12) {
          line.push(pick(next, [third, fifth, seventh]));
          line.push(pick(next, [fifth, third, nextRoot + 2, nextRoot - 2]));
        }
        line.push(pick(next, [nextRoot - 1, nextRoot + 1, nextRoot - 5, nextRoot + 7]));
        line.forEach((m, k) => push('bass', bassRange(m), start + k * BEAT, BEAT * 0.92, 0.46 + next() * 0.08 + (k === 0 ? 0.05 : 0)));
      }
      if (seg.slots === 12 && next() < 0.3) {
        // triplet pickup into the next chord
        push('bass', bassRange(nextRoot + (next() < 0.5 ? -1 : 2)), start + 11 * SLOT, SLOT * 0.9, 0.3);
      }
    });

    // Left hand: rolled voicings on the shuffle grid.
    const roll = (voicing, start, dur, vel, spread) =>
      voicing.forEach((m, v) => push('piano', m, start + v * spread + next() * 0.004, dur, vel * (v === 0 ? 0.92 : 1)));
    segments.forEach((seg) => {
      const start = slotTime(seg.bar, seg.slot);
      if (seg.bar === 0) {
        roll([...seg.chord.voicing, 79], start + 0.4, BEAT * 2.4, 0.3, 0.07);
        roll(seg.chord.voicing, start + 2.5 * BEAT, BEAT * 1.4, 0.24, 0.05);
        return;
      }
      const pattern = pick(next, seg.slots === 12 ? COMP : COMP_SHORT);
      const dynamics = seg.chorus < 1 ? 0.88 : 1;
      pattern.forEach(([slot, len]) => {
        if (slot >= seg.slots) return;
        const spread = next() < 0.25 ? 0.045 : 0.014;
        roll(seg.chord.voicing, start + slot * SLOT, len * SLOT * 0.95, (0.27 + next() * 0.1) * dynamics, spread);
      });
    });

    // Right hand: call-and-response phrases with blue notes and grace-note slides.
    const grace = (midi, time, vel) => push('piano', midi - 1, time - 0.15, 0.17, vel * 0.7);
    const melodyNote = (midi, time, dur, vel) => {
      const p = pc(midi);
      if ((p === 9 || p === 0 || p === 2) && next() < 0.45) grace(midi, time, vel);
      push('piano', midi, time, dur, vel);
    };
    const choosePitch = (current, seg, strong) => {
      const pool = [];
      for (let m = 65; m <= 89; m++) {
        const p = pc(m);
        let w = 0;
        if (seg.chord.tones.includes(p)) w += 3;
        if (BLUES.includes(p)) w += 2;
        if (seg.chord.colour.includes(p)) w += 1;
        if (w === 0) continue;
        if (strong && !seg.chord.tones.includes(p) && p !== 8) continue;
        const dist = current === null ? Math.abs(m - 79) : Math.abs(m - current);
        if (current !== null && dist === 0) continue;
        w *= dist <= 2 ? 3 : dist <= 5 ? 1.6 : dist <= 7 ? 0.6 : 0.15;
        pool.push([m, w]);
      }
      const total = pool.reduce((s, [, w]) => s + w, 0);
      let r = next() * total;
      for (const [m, w] of pool) {
        r -= w;
        if (r <= 0) return m;
      }
      return pool[pool.length - 1][0];
    };
    let previous = null;
    let current = null;
    for (let b = 1; b < bars; b++) {
      const inChorus = b - INTRO_BARS;
      if (b === 1) {
        // pickup lick into the first chorus: C Eb F rising into the downbeat
        [[8, 72], [9, 75], [10, 77]].forEach(([slot, m], k) => push('piano', m, slotTime(1, slot), SLOT * 0.9, 0.3 + k * 0.03));
        current = 77;
        continue;
      }
      if (inChorus % 12 === 4 && next() < 0.8) {
        // slow-blues tremolo between a third, resolving down to the tonic
        const [lo, hi] = pick(next, [[81, 84], [79, 82], [77, 81]]);
        const t0 = slotTime(b, 6);
        const period = 1 / 7.5;
        for (let k = 0; k * period < 5 * SLOT; k++) {
          push('piano', k % 2 ? hi : lo, t0 + k * period, period * 1.4, 0.3 - k * 0.012);
        }
        push('piano', 77, t0 + 5.2 * SLOT, SLOT * 3, 0.3);
        current = 77;
        previous = null;
        continue;
      }
      if (next() < 0.38) continue; // rests are part of the phrasing
      const response = previous && next() < 0.5;
      const cell = response ? previous.cell : pick(next, CELLS);
      const offset = response ? previous.offset : pick(next, [0, 2, 3, 6]);
      previous = { cell, offset };
      let slot = offset;
      cell.forEach((len, i) => {
        if (slot >= 12) return;
        const seg = segmentAt(b, slot);
        const strong = slot % 6 === 0;
        const midi = choosePitch(current, seg, strong);
        current = midi;
        const held = i === cell.length - 1;
        melodyNote(midi, slotTime(b, slot), len * SLOT * (held ? 0.85 : 0.9), 0.28 + next() * 0.12 + (held ? 0.03 : 0));
        slot += len;
      });
    }

    // Ending: a rolled F6/9 with the bass on the tonic, then the master fade takes over.
    const t = bars * BAR;
    push('bass', 41, t, duration - t, 0.5);
    [57, 60, 62, 67, 74, 79].forEach((m, v) => push('piano', m, t + v * 0.045, duration - t, 0.34 - v * 0.01));

    notes.sort((a, b) => a.start - b.start);
    return notes;
  }

  // --------------------------------------------------------------- synthesis
  function makeVoice(note, sampleRate) {
    const f0 = midiHz(note.midi);
    const piano = note.inst === 'piano';
    const partials = [];
    const count = piano ? Math.min(14, Math.floor(7500 / f0)) : Math.min(7, Math.floor(2600 / f0));
    const brightness = piano ? 0.95 + Math.max(0, note.midi - 60) / 90 - note.vel * 0.25 : 1.7;
    const d1 = piano ? 0.3 + Math.max(0, note.midi - 36) * 0.03 : 1.1 + Math.max(0, note.midi - 36) * 0.04;
    const B = piano ? 0.00016 : 0.00004;
    const add = (f, amp, decay, phase) => {
      if (f > sampleRate * 0.45) return;
      const w = (2 * Math.PI * f) / sampleRate;
      partials.push({ re: Math.cos(phase), im: Math.sin(phase), cos: Math.cos(w), sin: Math.sin(w), amp, g: Math.exp(-decay / sampleRate) });
    };
    for (let n = 1; n <= count; n++) {
      const f = n * f0 * Math.sqrt(1 + B * n * n);
      const amp = Math.pow(n, -brightness) * Math.exp(-(n - 1) * 0.1) * (!piano && n % 2 === 0 ? 0.75 : 1);
      add(f, amp, d1 * (1 + 0.3 * (n - 1)), 0);
      if (piano && n <= 2) {
        // Unison strings beat slightly against each other and sustain longer.
        add(f * (1 + 0.0011 * n), amp * 0.55, d1 * 0.65 * (1 + 0.3 * (n - 1)), 1.3);
        add(f * (1 - 0.0008 * n), amp * 0.4, d1 * 0.75 * (1 + 0.3 * (n - 1)), 2.6);
      }
    }
    const pan = piano ? Math.max(-0.5, Math.min(0.5, (note.midi - 63) / 44)) : -0.12;
    return {
      note,
      partials,
      startSample: Math.round(note.start * sampleRate),
      endSample: Math.round((note.start + note.dur) * sampleRate),
      attack: Math.max(1, Math.round((piano ? 0.005 : 0.012) * sampleRate)),
      release: Math.exp(-1 / ((piano ? 0.09 : 0.06) * sampleRate)),
      gainL: Math.cos(((pan + 1) * Math.PI) / 4),
      gainR: Math.sin(((pan + 1) * Math.PI) / 4),
      level: note.vel * (piano ? 0.13 : 0.3),
      noise: (piano ? 0.01 : 0.018) * note.vel,
      env: 0,
      releasing: false,
      done: false,
    };
  }

  class Reverb {
    constructor(sampleRate) {
      const scale = sampleRate / 44100;
      const mk = (len) => ({ buf: new Float32Array(Math.round(len * scale)), idx: 0, filt: 0 });
      this.combsL = [1116, 1188, 1277, 1356].map(mk);
      this.combsR = [1139, 1211, 1300, 1379].map(mk);
      this.allL = [556, 441].map(mk);
      this.allR = [579, 464].map(mk);
      this.feedback = 0.84;
      this.damp = 0.4;
    }
    comb(c, x) {
      const y = c.buf[c.idx];
      c.filt = y * (1 - this.damp) + c.filt * this.damp;
      c.buf[c.idx] = x + c.filt * this.feedback;
      c.idx = (c.idx + 1) % c.buf.length;
      return y;
    }
    allpass(a, x) {
      const y = a.buf[a.idx];
      const out = y - x;
      a.buf[a.idx] = x + y * 0.5;
      a.idx = (a.idx + 1) % a.buf.length;
      return out;
    }
    process(input, combs, alls) {
      let y = 0;
      for (const c of combs) y += this.comb(c, input);
      for (const a of alls) y = this.allpass(a, y);
      return y;
    }
  }

  function createSynth({ sampleRate = 44100, duration = 100, seed = 20260926 } = {}) {
    const notes = compose(duration, seed);
    const totalSamples = Math.round(duration * sampleRate);
    let cursor = 0;
    let nextNote = 0;
    let active = [];
    let reverb = new Reverb(sampleRate);
    let lpL = 0;
    let lpR = 0;
    const lpCoefficient = 1 - Math.exp((-2 * Math.PI * 4200) / sampleRate);
    const noise = rng(seed ^ 0x9e3779b9);
    const fadeIn = 1.2 * sampleRate;
    const fadeOut = 2.6 * sampleRate;

    function seek(time) {
      cursor = Math.max(0, Math.min(totalSamples, Math.round(time * sampleRate)));
      active = [];
      reverb = new Reverb(sampleRate);
      lpL = lpR = 0;
      nextNote = 0;
      // Notes that started before the seek point and still sound are re-entered mid-flight.
      while (nextNote < notes.length && Math.round(notes[nextNote].start * sampleRate) < cursor) {
        const v = makeVoice(notes[nextNote], sampleRate);
        if (v.endSample > cursor - 0.5 * sampleRate) {
          const elapsed = cursor - v.startSample;
          for (const p of v.partials) p.amp *= Math.pow(p.g, elapsed);
          v.env = 1;
          active.push(v);
        }
        nextNote++;
      }
    }

    function render(left, right) {
      const n = left.length;
      const end = Math.min(totalSamples, cursor + n);
      left.fill(0);
      right.fill(0);
      if (cursor >= totalSamples) return false;
      while (nextNote < notes.length && Math.round(notes[nextNote].start * sampleRate) < end) {
        active.push(makeVoice(notes[nextNote], sampleRate));
        nextNote++;
      }
      for (const v of active) {
        const from = Math.max(0, v.startSample - cursor);
        for (let i = from; i < end - cursor; i++) {
          const s = cursor + i;
          if (s >= v.endSample) v.releasing = true;
          if (v.releasing) v.env *= v.release;
          else if (v.env < 1) v.env = Math.min(1, v.env + 1 / v.attack);
          let x = 0;
          let energy = 0;
          for (const p of v.partials) {
            const re = p.re * p.cos - p.im * p.sin;
            p.im = p.re * p.sin + p.im * p.cos;
            p.re = re;
            p.amp *= p.g;
            x += p.im * p.amp;
            energy += p.amp;
          }
          const age = s - v.startSample;
          if (age < v.attack * 2) x += (noise() * 2 - 1) * v.noise * (1 - age / (v.attack * 2));
          const y = x * v.env * v.level;
          left[i] += y * v.gainL;
          right[i] += y * v.gainR;
          if (energy * v.env < 1e-4) {
            v.done = true;
            break;
          }
        }
      }
      active = active.filter((v) => !v.done);
      for (let i = 0; i < end - cursor; i++) {
        const s = cursor + i;
        // A gentle soundboard low-pass takes the synthetic edge off the partials.
        lpL += lpCoefficient * (left[i] - lpL);
        lpR += lpCoefficient * (right[i] - lpR);
        const mono = (lpL + lpR) * 0.5;
        const wetL = reverb.process(mono, reverb.combsL, reverb.allL) * 0.19;
        const wetR = reverb.process(mono, reverb.combsR, reverb.allR) * 0.19;
        let gain = 1;
        if (s < fadeIn) gain = s / fadeIn;
        if (totalSamples - s < fadeOut) gain *= (totalSamples - s) / fadeOut;
        left[i] = Math.tanh((lpL + wetL) * gain * 1.15);
        right[i] = Math.tanh((lpR + wetR) * gain * 1.15);
      }
      cursor = end;
      return cursor < totalSamples;
    }

    seek(0);
    return { seek, render, notes, sampleRate, duration, get time() { return cursor / sampleRate; } };
  }

  return { BPM, BAR, compose, createSynth };
});
