/**
 * Single-sequence k-mer coverage calculator.
 */
public class MaskingSequence {

    /** Minimum k-mer occurrence count for a region to be treated as repetitive.
     *  Configurable via 'repeatMinCount='; 2 reproduces the original behaviour. */
    public static int MIN_COUNT = 2;

    private final int kmer;

    public MaskingSequence() {
        this(18);
    }

    public MaskingSequence(int kmer) {
        this.kmer = kmer;
    }

    // ── Public API ──────────────────────────────────────────────────────────
    /**
     * Compute per-base repeat coverage for a single sequence.
     *
     * @param seq DNA sequence
     * @return int array of per-base coverage counts from repeated k-mers
     */
    public int[] mask(String seq) {
        final byte[] fwd = normalise(seq);
        final byte[] rev = normalise(dna.ComplementDNA(seq));
        final int len = fwd.length;

        // Phase 1: count k-mers on both strands
        LongIntHashMap kmerMap = new LongIntHashMap();
        countStrand(kmerMap, fwd, len);
        countStrand(kmerMap, rev, len);

        // Phase 2: compute per-base coverage — kmerMap  
        // filtering through get(code) >= MIN_COUNT instead of Set.contains()
        int[] coverage = new int[len];
        markStrandForward(kmerMap, fwd, coverage, len);
        markStrandReverse(kmerMap, rev, coverage, len);

        return coverage;
    }

    // ── Phase 1: k-mer counting ─────────────────────────────────────────────
    /**
     * Slide a k-mer window, counting gap-free k-mers.
     */
    private void countStrand(LongIntHashMap map, byte[] bases, int len) {
        int gapCount = 0;

        for (int i = 0; i < kmer - 1; i++) {
            if (bases[i] == 4) {
                gapCount++;
            }
        }

        for (int i = kmer - 1; i < len; i++) {
            if (bases[i] == 4) {
                gapCount++;
            }

            if (gapCount == 0) {
                long code = encodeKmer(bases, i - kmer + 1);
                int current = map.get(code);
                if (current == LongIntHashMap.NO_VALUE) {
                    map.put(code, 1);
                } else if (current < MIN_COUNT) {
                    map.put(code, current + 1);
                }
            }

            if (bases[i + 1 - kmer] == 4) {
                gapCount--;
            }
        }
    }

    // ── Phase 2: coverage marking ───────────────────────────────────────────
    /**
     * get(code) >= MIN_COUNT
     */
    private void markStrandForward(LongIntHashMap kmerMap, byte[] bases, int[] coverage, int len) {
        int gapCount = 0;

        for (int i = 0; i < kmer - 1; i++) {
            if (bases[i] == 4) {
                gapCount++;
            }
        }

        for (int i = kmer - 1; i < len; i++) {
            if (bases[i] == 4) {
                gapCount++;
            }
            int start = i - kmer + 1;

            if (gapCount == 0) {
                long code = encodeKmer(bases, start);
                if (kmerMap.get(code) >= MIN_COUNT) {
                    for (int j = start; j <= i; j++) {
                        coverage[j]++;
                    }
                }
            }

            if (bases[start] == 4) {
                gapCount--;
            }
        }
    }

    private void markStrandReverse(LongIntHashMap kmerMap, byte[] bases, int[] coverage, int len) {
        int gapCount = 0;

        for (int i = 0; i < kmer - 1; i++) {
            if (bases[i] == 4) {
                gapCount++;
            }
        }

        for (int i = kmer - 1; i < len; i++) {
            if (bases[i] == 4) {
                gapCount++;
            }
            int start = i - kmer + 1;

            if (gapCount == 0) {
                long code = encodeKmer(bases, start);
                if (kmerMap.get(code) >= MIN_COUNT) {
                    // Map reverse-complement window → forward-strand coordinates
                    int fwdStart = len - start - kmer;
                    int fwdEnd = len - start;
                    for (int j = fwdStart; j < fwdEnd; j++) {
                        coverage[j]++;
                    }
                }
            }

            if (bases[start] == 4) {
                gapCount--;
            }
        }
    }

 
    private static byte[] normalise(String seq) {
        byte[] raw = seq.getBytes();
        for (int i = 0; i < raw.length; i++) {
            raw[i] = tables.dx2[raw[i]];
        }
        return raw;
    }

    /**
     * Base-10 k-mer encoding. Collision-free for the DNA+gap alphabet; stays  positive for k ≤ 19.
     */
    private long encodeKmer(byte[] bases, int start) {
        long code = 0L;
        for (int i = start; i < start + kmer; i++) {
            code = code * 5 + bases[i];
        }
        return code;
    }
}
